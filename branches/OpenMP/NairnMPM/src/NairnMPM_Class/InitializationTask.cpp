/********************************************************************************
	InitializationTask.cpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/InitializationTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/MaterialInterfaceNode.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "Cracks/CrackNode.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Patches/GridPatch.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Materials/MaterialBase.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark CONSTRUCTORS

InitializationTask::InitializationTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
void InitializationTask::Execute(void)
{
	CommonException *initErr = NULL;
	
	// Zero Mass Matrix and vectors
	warnings.BeginStep();

	int tp = fmobj->GetTotalNumberOfPatches();
#pragma omp parallel
	{
        // zero all nodal variables on real nodes
#pragma omp for
		for(int i=1;i<=nnodes;i++)
			nd[i]->InitializeForTimeStep();
		
        // zero ghost nodes in patch for this thread
        int pn = GetPatchNumber();
		patches[pn]->InitializeForTimeStep();

		// particle calculations
#pragma omp for
		for(int p=0;p<nmpmsRC;p++)
		{	MPMBase *mpmptr = mpm[p];										// pointer
			const ElementBase *elref = theElements[mpmptr->ElemID()];		// element containing this particle
			try
			{	elref->GetShapeFunctionData(mpmptr);
			}
			catch(CommonException err)
			{	if(initErr==NULL)
				{
#pragma omp critical
					initErr = new CommonException(err);
				}
			}			
		}
	}
	
	// was there an error?
	if(initErr!=NULL) throw *initErr;
	
	// allocate crack and material velocity fields needed for time step on real nodes
    // tried critical sections when nodes changes, but it was slower
    // can't use ghost nodes, because need to test all on real nodes
	if(firstCrack!=NULL || maxMaterialFields>1)
	{   int nds[maxShapeNodes];

		for(int pn=0;pn<tp;pn++)
		{
           // do non-rigid and rigid contact materials in patch pn
            for(int block=FIRST_NONRIGID;block<=FIRST_RIGID_CONTACT;block++)
            {   MPMBase *mpmptr = patches[pn]->GetFirstBlockPointer(block);
                while(mpmptr!=NULL)
                {	const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
                    int matfld = matID->GetField();									// material velocity field
                    
                    // get nodes and shape function for material point p
                    int i,numnds;
                    const ElementBase *elref = theElements[mpmptr->ElemID()];		// element containing this particle
                    elref->GetShapeFunctionNodes(&numnds,nds,mpmptr->GetNcpos(),mpmptr);
                    
                    // Add particle property to each node in the element
                    short vfld;
                    NodalPoint *ndptr;
                    for(i=1;i<=numnds;i++)
                    {	// use real node in this loop
                        ndptr = nd[nds[i]];
                        
                        if(firstCrack==NULL)
                        {	vfld=0;
                        }
                        else
                        {	// in CRAMP, find crack crossing and appropriate velocity field
                            CrackField cfld[2];
                            cfld[0].loc = NO_CRACK;			// NO_CRACK, ABOVE_CRACK, or BELOW_CRACK
                            cfld[1].loc = NO_CRACK;
                            int cfound=0;
                            Vector norm;
                            
                            CrackHeader *nextCrack = firstCrack;
                            while(nextCrack!=NULL)
                            {	vfld = nextCrack->CrackCross(mpmptr->pos.x,mpmptr->pos.y,ndptr->x,ndptr->y,&norm);
                                if(vfld!=NO_CRACK)
                                {	cfld[cfound].loc=vfld;
                                    cfld[cfound].norm=norm;
#ifdef IGNORE_CRACK_INTERACTIONS
                                    cfld[cfound].crackNum=1;	// appears to always be same crack, and stop when found one
                                    break;
#else
                                    cfld[cfound].crackNum=nextCrack->GetNumber();
                                    cfound++;
                                    if(cfound>1) break;			// stop if found two, if there are more then two, physics will be off
#endif
                                }
                                nextCrack=(CrackHeader *)nextCrack->GetNextObject();
                            }
                                
                            
                            // momentum vector (and allocate velocity field if needed)
                            vfld = ndptr->AddCrackVelocityField(matfld,cfld);
                            mpmptr->vfld[i] = vfld;
                        }
                        
                        // make sure material velocity field is created too
                        if(maxMaterialFields>1)
                            ndptr->AddMatVelocityField(vfld,matfld);
                    }
                    
                    // next material point
                    mpmptr = (MPMBase *)mpmptr->GetNextObject();
                }
            }
		}
    
		// copy fields on real nodes to ghost nodes
		if(tp>1)
        {   for(int pn=0;pn<tp;pn++)
				patches[pn]->InitializationReduction();
		}
	}
	
    // Update forces applied to particles
	MatPtLoadBC::SetParticleFext(mtime);
	
	// remove contact conditions
	CrackNode::RemoveCrackNodes();
	MaterialInterfaceNode::RemoveInterfaceNodes();
	
    // turn off isothermal ramp when done and ramp step initialization
	thermal.CheckDone(mtime);
	
}	
