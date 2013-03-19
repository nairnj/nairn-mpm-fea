/********************************************************************************
    CalcJKTask.cpp
    NairnMPM
    
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

// Global
CalcJKTask *theJKTask=NULL;

#pragma mark INITIALIZE

// Constructors
CalcJKTask::CalcJKTask()
{	theJKTask=this;
	// allocate J integral data on each particle
    int p;
	for(p=0;p<nmpms;p++)
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
CustomTask *CalcJKTask::PrepareForStep(bool &doCrackExtraps)
{
	if((getJKThisStep=archiver->WillArchiveJK(TRUE)))
    {	// if archiving J or K, need this task's extrapolations
    	doCrackExtraps=TRUE;
    }
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
        
    return nextTask;
}

// Calculate J and K at crack tips
CustomTask *CalcJKTask::StepCalculation(void)
{
    // skip if not needed
    if(!getJKThisStep) return nextTask;
    
    int i,inMat;
    Vector d,C;
    CrackSegment *crkTip;

    CrackHeader *nextCrack=firstCrack;
    while(nextCrack!=NULL)
    {   nextCrack->JIntegral();         // crack-axis components of J-integral
        
        // if material known, find KI and KII for crack tips
        if(getJKThisStep & NEED_K)
        {   for(i=START_OF_CRACK;i<=END_OF_CRACK;i++)
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

#pragma mark TASK EXTRAPOLATION METHODS

// initialize for crack extrapolations
CustomTask *CalcJKTask::BeginExtrapolations(void)
{
    // skip if already set up
    if(!getJKThisStep) return nextTask;
	
    // set up strain fields for crack extrapolations
    int i;
    for(i=1;i<=nnodes;i++)
        nd[i]->ZeroDisp();
    
    return nextTask;
}

// add particle data to a node
CustomTask *CalcJKTask::NodalExtrapolation(NodalPoint *ndmi,MPMBase *mpnt,short vfld,int matfld,double wt,short isRigid)
{
    // skip if already set up
    if(!getJKThisStep || isRigid) return nextTask;
    
	// get 2D gradient terms (dimensionless) and track material (if needed)
	int matid = mpnt->MatID();
	ndmi->AddUGradient(vfld,wt,mpnt->GetDuDx(),mpnt->GetDuDy(),mpnt->GetDvDx(),mpnt->GetDvDy(),matid,mpnt->mp);
	
	// get a nodal stress (rho*stress has units N/m^2)
    wt *= theMaterials[matid]->rho;
    Tensor sp = mpnt->ReadStressTensor();
    ndmi->AddStress(vfld,wt,&sp);
    
	// get energy and rho*energy has units J/m^3 = N/m^2
	// In axisymmetric, energy density is 2 pi m U/(2 pi rp Ap), but since m = rho rp Ap
	//		energy density it still rho*energy
    ndmi->AddEnergy(vfld,wt,mpnt->vel.x,mpnt->vel.y,mpnt->GetStrainEnergy());
	
    return nextTask;
}

// initialize for crack extrapolations
CustomTask *CalcJKTask::EndExtrapolations(void)
{
    // skip if already set up
    if(!getJKThisStep) return nextTask;
    
    // finish strain fields
	int i;
    for(i=1;i<=nnodes;i++)
        nd[i]->CalcStrainField();
    
    return nextTask;
}
