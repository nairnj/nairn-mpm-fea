/********************************************************************************
	InitializationTask.cpp
	nairn-mpm-fea

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

		// particle calculations get xipos for particles and if doing CPDI
        // precalculate CPDI info needed for subsequent shape functions
#pragma omp for nowait
		for(int p=0;p<nmpmsRC;p++)
        {   MPMBase *mpmptr = mpm[p];                                       // pointer
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
    // tried critical sections when nodes changed, but it was slower
    // can't use ghost nodes, because need to test all on real nodes
	if(firstCrack!=NULL || maxMaterialFields>1)
	{
#pragma omp parallel
        {
            int nds[maxShapeNodes];
            double fn[maxShapeNodes];
            
            //for(int pn=0;pn<tp;pn++)
            int pn = GetPatchNumber();
		
            // do non-rigid and rigid contact materials in patch pn
            for(int block=FIRST_NONRIGID;block<=FIRST_RIGID_CONTACT;block++)
            {   // get material point (only in this patch)
                MPMBase *mpmptr = patches[pn]->GetFirstBlockPointer(block);
                
                while(mpmptr!=NULL)
                {	const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
                    const int matfld = matID->GetField();                           // material velocity field
                    
                    // get nodes and shape function for material point p
                    const ElementBase *elref = theElements[mpmptr->ElemID()];		// element containing this particle
                    
                    // don't actually need shape functions, but need to screen out zero shape function
                    // like done in subsequent tasks, otherwise node numbers will not align correctly
                    // only think used from return are numnds and nds
                    int numnds;
                    elref->GetShapeFunctions(&numnds,fn,nds,mpmptr);
                    
                    // Add particle property to each node in the element
                    for(int i=1;i<=numnds;i++)
                    {	// use real node in this loop
                        NodalPoint *ndptr = nd[nds[i]];
                        
                        // always zero when no cracks
                        short vfld = 0;
#ifdef COMBINE_RIGID_MATERIALS
                        // when combining rigid particles, extrapolate all to field 0 and later
                        // copy to other active fields
                        if(firstCrack!=NULL && block!=FIRST_RIGID_CONTACT)
#else
                        if(firstCrack!=NULL)
#endif
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
                                
                            
                            // find (and allocate if needed) the velocity field
                            // Use vfld=0 if no cracks found
                            if(cfound>0)
                            {   // In parallel, this is critical code
#pragma omp critical
                                {   try
                                    {   vfld = ndptr->AddCrackVelocityField(matfld,cfld);
                                    }
                                    catch(CommonException err)
                                    {   if(initErr==NULL)
                                            initErr = new CommonException(err);
                                    }
                                }
                            }
                            
                            // set material point velocity field for this node
                            mpmptr->vfld[i] = vfld;
                        }
                        
                        // make sure material velocity field is created too
                        if(maxMaterialFields>1 && ndptr->NeedsMatVelocityField(vfld,matfld))
                        {   // If parallel, this is critical code
#pragma omp critical
                            {   try
                                {   ndptr->AddMatVelocityField(vfld,matfld);
                                }
                                catch(CommonException err)
                                {   if(initErr==NULL)
                                        initErr = new CommonException(err);
                                }
                            }
                            
                        }
                    }
                    
                    // next material point
                    mpmptr = (MPMBase *)mpmptr->GetNextObject();
                }
            }
		}
    
        // was there an error?
        if(initErr!=NULL) throw *initErr;

		// copy crack and material fields on real nodes to ghost nodes
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
