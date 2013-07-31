/********************************************************************************
    CalcJKTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Fri Aug 15 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Custom_Tasks/CalcJKTask.hpp"
#include "System/ArchiveData.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Elements/ElementBase.hpp"
#include "Patches/GridPatch.hpp"
#include "NairnMPM_Class/MPMTask.hpp"

// Global
CalcJKTask *theJKTask=NULL;

#pragma mark INITIALIZE

// Constructors
CalcJKTask::CalcJKTask()
{	theJKTask=this;
	// allocate J integral data on each particle
    int p;
	for(p=0;p<nmpmsNR;p++)
		mpm[p]->AllocateJStructures();
}

// Return name of this task
const char *CalcJKTask::TaskName(void) { return "Calculate J and K Task"; }

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *CalcJKTask::Initialize(void)
{
    cout << "J Integral and Stress Intensity calculation activated." << endl;
    cout << "   Rectangular contour 2*(" << JGridSize << "X" << JGridSize << ") elements" << endl;
	if(JTerms<0) JTerms = fmobj->IsAxisymmetric() ? 2 : 1 ;
	switch(JTerms)
	{	case 1:
			cout << "   Contour integral only";
			break;
		case 2:
			cout << "   Contour and volume integrals";
			break;
		default:
			cout << "   Unknown option, will use contour integral only";
			JTerms=1;
			break;
	}
	cout << endl;
    return nextTask;
}

// called when MPM step is getting ready to do custom tasks
// has its own extrapolations for speed
CustomTask *CalcJKTask::PrepareForStep(bool &needExtraps)
{
	getJKThisStep=archiver->WillArchiveJK(TRUE);
    return nextTask;
}

// Called when custom tasks are all done on a step
CustomTask *CalcJKTask::FinishForStep(void)
{
    // skip if not done
    if(!getJKThisStep) return nextTask;
    
    // finished with strain fields
    int i;
    for(i=1;i<=nnodes;i++)
        nd[i]->DeleteDisp();
        
    int totalPatches = fmobj->GetTotalNumberOfPatches();
    if(totalPatches>1)
    {	for(i=0;i<totalPatches;i++)
            patches[i]->DeleteDisp();
    }
    
    return nextTask;
}

// Calculate J and K at crack tips
CustomTask *CalcJKTask::StepCalculation(void)
{
    // skip if not needed
    if(!getJKThisStep) return nextTask;
    
    int nds[maxShapeNodes];
    double fn[maxShapeNodes];
	//double xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
    
    // set up strain fields for crack extrapolations
//#pragma omp parallel private(nds,fn,xDeriv,yDeriv,zDeriv)
#pragma omp parallel private(nds,fn)
    {	// in case 2D planar
        //for(int i=0;i<maxShapeNodes;i++) zDeriv[i] = 0.;
        
#pragma omp for
        for(int i=1;i<=nnodes;i++)
            nd[i]->ZeroDisp();
	
        // zero displacement fields on ghost nodes
        int pn = MPMTask::GetPatchNumber();
        patches[pn]->ZeroDisp();
        
        // loop over only non-rigid particles in patch
        MPMBase *mpnt = patches[pn]->GetFirstBlockPointer(FIRST_NONRIGID);
        while(mpnt!=NULL)
        {   // material reference
            const MaterialBase *matref = theMaterials[mpnt->MatID()];
		
            // find shape functions and derviatives
            int numnds;
            const ElementBase *elref = theElements[mpnt->ElemID()];
            //elref->GetShapeGradients(&numnds,fn,nds,xDeriv,yDeriv,zDeriv,mpnt);
			elref->GetShapeFunctions(&numnds,fn,nds,mpnt);
		
            // Add particle property to each node in the element
            NodalPoint *ndmi;
            short vfld;
            double fnmp;
            for(int i=1;i<=numnds;i++)
            {   // global mass matrix
                vfld=(short)mpnt->vfld[i];				// velocity field to use
                fnmp=fn[i]*mpnt->mp;
                
                // get node pointer
                ndmi = MPMTask::GetNodePointer(pn,nds[i]);
			
                // get 2D gradient terms (dimensionless) and track material (if needed)
                int activeMatField = matref->GetActiveField();
                ndmi->AddUGradient(vfld,fnmp,mpnt->GetDuDx(),mpnt->GetDuDy(),mpnt->GetDvDx(),mpnt->GetDvDy(),activeMatField,mpnt->mp);
			
                // get a nodal stress (rho*stress has units N/m^2)
                fnmp *= matref->rho;
                Tensor sp = mpnt->ReadStressTensor();
                ndmi->AddStress(vfld,fnmp,&sp);
			
                // get energy and rho*energy has units J/m^3 = N/m^2
                // In axisymmetric, energy density is 2 pi m U/(2 pi rp Ap), but since m = rho rp Ap
                //		energy density it still rho*energy
                ndmi->AddEnergy(vfld,fnmp,mpnt->vel.x,mpnt->vel.y,mpnt->GetStrainEnergy());
			
            }
            
            // next non-rigid material point
            mpnt = (MPMBase *)mpnt->GetNextObject();
        }
    }
        
    // copy ghost to real nodes
    int totalPatches = fmobj->GetTotalNumberOfPatches();
    if(totalPatches>1)
    {	for(int j=0;j<totalPatches;j++)
            patches[j]->JKTaskReduction();
    }
 	
    // finish strain fields
#pragma omp parallel for
    for(int i=1;i<=nnodes;i++)
        nd[i]->CalcStrainField();
    
    // No Do the J Integral calculations
    
	int inMat;
    Vector d,C;
    CrackSegment *crkTip;

    CrackHeader *nextCrack=firstCrack;
    while(nextCrack!=NULL)
    {   nextCrack->JIntegral();         // crack-axis components of J-integral
        
        // if material known, find KI and KII for crack tips
        if(getJKThisStep & NEED_K)
        {   for(int i=START_OF_CRACK;i<=END_OF_CRACK;i++)
            {   crkTip=nextCrack->GetCrackTip(i);
                inMat=crkTip->tipMatnum;
                if(inMat>0)
                {   //C=crkTip->C;		// will be crack segment property
                    C.x=0.0; C.y=0.0; ; C.z=0.;
                    
                    // find normal and shear COD
					nextCrack->GetCOD(crkTip,d,true);
					d.z=0.;

					// convert to K
                    crkTip->sif=theMaterials[inMat-1]->ConvertJToK(d,C,crkTip->Jint,fmobj->np);
                }
                else
                {   crkTip->sif.x=0;
                    crkTip->sif.y=0;
                }
            }
        }

        // next crack
        nextCrack=(CrackHeader *)nextCrack->GetNextObject();
    }
    
    return nextTask;
}

#pragma mark CalcJKTask METHODS

// Another task can invoke J and K if needed
//   use FALSE, NEED_J, NEED_JANDK, or NEED_J+NEED_K (NEED_K alone no good)
void CalcJKTask::ScheduleJK(int newNeed) { getJKThisStep|=newNeed; }

