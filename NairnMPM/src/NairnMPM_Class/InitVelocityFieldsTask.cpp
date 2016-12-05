/********************************************************************************
	InitVelocityFieldsTask.cpp
	nairn-mpm-fea

	Created by John Nairn on March 6, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/InitVelocityFieldsTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Patches/GridPatch.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Exceptions/CommonException.hpp"

// to ignore crack interactions (only valid if 1 crack or non-interacting cracks)
//#define IGNORE_CRACK_INTERACTIONS

#pragma mark CONSTRUCTORS

InitVelocityFieldsTask::InitVelocityFieldsTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// allocate crack and material velocity fields needed for time step on real nodes
// tried critical sections when nodes changed, but it was slower
// can't use ghost nodes, because need to test all on real nodes
//
// This task only used if have cracks or in multimaterial mode
// throws CommonException()
void InitVelocityFieldsTask::Execute(void)
{
	CommonException *initErr = NULL;
	
	int tp = fmobj->GetTotalNumberOfPatches();
	
#pragma omp parallel
	{
#ifndef LOAD_GIMP_INFO
#ifdef CONST_ARRAYS
		int ndsArray[MAX_SHAPE_NODES];
		double fn[MAX_SHAPE_NODES];
#else
		int ndsArray[maxShapeNodes];
		double fn[maxShapeNodes];
#endif
#endif
		
		int pn = GetPatchNumber();
		
		// do non-rigid and rigid contact materials in patch pn
		for(int block=FIRST_NONRIGID;block<=FIRST_RIGID_CONTACT;block++)
		{   // get material point (only in this patch)
			MPMBase *mpmptr = patches[pn]->GetFirstBlockPointer(block);

			while(mpmptr!=NULL)
			{	const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
				const int matfld = matID->GetField();                           // material velocity field
				
#ifdef LOAD_GIMP_INFO
				GIMPNodes *gimp = mpmptr->GetGIMPInfo();
				int *nds = gimp->nds;
				int numnds = nds[0];
#else
				// get nodes and shape function for material point p
				const ElementBase *elref = theElements[mpmptr->ElemID()];		// element containing this particle
				
				// don't actually need shape functions, but need to screen out zero shape function
				// like done in subsequent tasks, otherwise node numbers will not align correctly
				// only thing used from return are numnds and nds
				int *nds = ndsArray;
				elref->GetShapeFunctions(fn,&nds,mpmptr);
				int numnds = nds[0];
#endif
				
				// Only need to decipher crack velocity field if has cracks (firstCrack!=NULL)
				//      and if this material allows cracks.
				bool decipherCVF = firstCrack!=NULL && block!=FIRST_RIGID_CONTACT;
				
				// Check each node
				for(int i=1;i<=numnds;i++)
				{	// use real node in this loop
					NodalPoint *ndptr = nd[nds[i]];
					
					// always zero when no cracks (or when ignoring cracks)
					short vfld = 0;
					
					// If need, find vlocity field and for each field set location
					// (above or below crack) and crack number (1 based) or 0 for NO_CRACK
					if(decipherCVF)
					{	// in CRAMP, find crack crossing and appropriate velocity field
						CrackField cfld[2];
						cfld[0].loc = NO_CRACK;			// NO_CRACK=0, ABOVE_CRACK=1, or BELOW_CRACK=2
						cfld[1].loc = NO_CRACK;
						int cfound=0;
						Vector norm;					// track normal vector for crack plane
						
						CrackHeader *nextCrack = firstCrack;
						while(nextCrack!=NULL)
						{	vfld = nextCrack->CrackCross(mpmptr->pos.x,mpmptr->pos.y,ndptr->x,ndptr->y,&norm);
							if(vfld!=NO_CRACK)
							{	cfld[cfound].loc=vfld;
								cfld[cfound].norm=norm;
								
#ifdef IGNORE_CRACK_INTERACTIONS
								// appears to always be same crack, and stop when found one
								cfld[cfound].crackNum=1;
								break;
#endif
								
								// Get crack number (default code does not ignore interactions)
								cfld[cfound].crackNum=nextCrack->GetNumber();
								cfound++;
								
								// stop if found two because code can only handle two interacting cracks
								// It exits loop now to go ahead with the first two found, by physics may be off
								if(cfound>1) break;
							}
							nextCrack=(CrackHeader *)nextCrack->GetNextObject();
						}
						
						
						// find (and allocate if needed) the velocity field
						// Use vfld=0 if no cracks found
						if(cfound>0)
						{   // In parallel, this is critical code
#pragma omp critical (addcvf)
							{   try
								{   vfld = ndptr->AddCrackVelocityField(matfld,cfld);
								}
								catch(std::bad_alloc&)
								{   if(initErr==NULL)
										initErr = new CommonException("Memory error","InitVelocityFieldsTask::Execute");
								}
								catch(...)
								{	if(initErr==NULL)
										initErr = new CommonException("Unexpected error","InitVelocityFieldsTask::Execute");
								}
							}
						}
						
						// set material point velocity field for this node
						mpmptr->vfld[i] = (char)vfld;
					}
					
					// make sure material velocity field is created too
					// (Note: when maxMaterialFields==1 (Singe Mat Mode), mvf[0] is always there
					//        so no need to create it here)
					// When some materials ignore cracks, those materials always use [0]
					if(maxMaterialFields>1 && ndptr->NeedsMatVelocityField(vfld,matfld))
					{   // If parallel, this is critical code
#pragma omp critical (addcvf)
						{   try
							{   ndptr->AddMatVelocityField(vfld,matfld);
							}
							catch(std::bad_alloc&)
							{   if(initErr==NULL)
									initErr = new CommonException("Memory error","InitVelocityFieldsTask::Execute");
							}
							catch(...)
						 	{	if(initErr==NULL)
									initErr = new CommonException("Unexpected error","InitVelocityFieldsTask::Execute");
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
