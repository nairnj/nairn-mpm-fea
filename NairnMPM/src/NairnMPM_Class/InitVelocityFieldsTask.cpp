/********************************************************************************
	InitVelocityFieldsTask.cpp
	nairn-mpm-fea

	Created by John Nairn on March 6, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.
 
	The tasks are:
	--------------
	* If no crack and single material, set each MP to use vfld=0
    * If cracks or multimaterial mode
		- Find vfld for CVF (create CVF if will be needed)
		- In multimaterial mode, and MVF if will be needed
	* Reduction to copy field info to ghost nodes
********************************************************************************/
#if defined ( _MSC_VER) || defined (__APPLE__) 
#include "stdafx.h"
#endif
#include "System/MPMPrefix.hpp"
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

#ifdef PREHASH_CRACKS
bool InitVelocityFieldsTask::prehashed = true;
#endif

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
bool InitVelocityFieldsTask::Execute(int taskOption)
{
#ifdef PREHASH_CRACKS
	// If cracks, prehashing maintains a list of cracks seen
    // by each element (including nearest neighbors for GIMP or CPDI).
	// This section clears lists from previous step
    // .. and then rebuilds them for this step.
	if(prehashed && numberOfCracks>0)
    {   // clear elements first (not parallel)
		for (int i = 0; i < nelems; i++)
			theElements[i]->DeleteCrackList();

		// loop over all cracks (not parallel)
		for (int cn = 0; cn < numberOfCracks; cn++)
			crackList[cn]->UpdateElementCrackList(cn);
	}
#endif
	
	CommonException *initErr = NULL;
	
	int tp = fmobj->GetTotalNumberOfPatches();

#pragma omp parallel
	{
#ifdef CONST_ARRAYS
		int ndsArray[MAX_SHAPE_NODES];
		double fn[MAX_SHAPE_NODES];
#else
		int ndsArray[maxShapeNodes];
		double fn[maxShapeNodes];
#endif
		
		int pn = GetPatchNumber();
		
        // This is the "main particle loop" to initialize crack velocity fields.
        // It is in the appendix of our paper on interacting cracks. The steps
        // in this loop are labeled with "paper" steps in that paper's appendix.
		// Loop over non-rigid, rigid block, and rigid contact particles in patch pn
		for(int block=FIRST_NONRIGID;block<=FIRST_RIGID_CONTACT;block++)
		{   // get material point (only in this patch)
			MPMBase *mpmptr = patches[pn]->GetFirstBlockPointer(block);

			while(mpmptr!=NULL)
			{	const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
				const int matfld = matID->GetField();                           // material velocity field
				
				// Get element containing this particle (paper step 1)
				const ElementBase *elref = theElements[mpmptr->ElemID()];
				
                // Get nodes and shape function for material point p (paper step 2)
				// don't actually need shape functions, but need to screen out zero shape function
				// like done in subsequent tasks, otherwise node numbers will not align correctly
				// only thing used from return are numnds and nds
				int *nds = ndsArray;
				try
				{	elref->GetShapeFunctions(fn, &nds, mpmptr);
				}
				catch(CommonException& err)
                {   if(initErr==NULL)
                    {
    #pragma omp critical (error)
                        initErr = new CommonException(err);
                    }
                    break;
                }
				catch(...)
				{	if(initErr==NULL)
						initErr = new CommonException("Unexpected error","InitVelocityFieldsTask::Execute");
                    break;
				}
				int numnds = nds[0];

				// Only need to decipher crack velocity field if has cracks (firstCrack!=NULL)
				//      and if this material allows cracks.
                bool decipherCVF = (firstCrack!=NULL) && matID->AllowsCracks();
#ifdef PREHASH_CRACKS
                // paper step 3 - skip looking for fields if element sees no cracks
                if(prehashed && elref->SeesCrack==nullptr) decipherCVF = false;
#endif
				// Check each node seen by this material point
				for(int i=1;i<=numnds;i++)
				{	// use real node in this loop
					NodalPoint *ndptr = nd[nds[i]];
					Vector ndpt = MakeVector(ndptr->x,ndptr->y,ndptr->z);
					
					// always zero when no cracks (or when ignoring cracks)
					short vfld = 0;
					
					// If needed, find velocity field and for each field set location
					// (above or below crack) and crack number (1 based) or 0 for NO_CRACK
					if(decipherCVF)
                    {   // in CRAMP, find crack crossing and appropriate velocity field
						CrackField cfld[2];
						cfld[0].loc = NO_CRACK;			// NO_CRACK=0, ABOVE_CRACK=1, or BELOW_CRACK=2
						cfld[1].loc = NO_CRACK;
                        
                        // paper step 3.1 for alpha=cfound to count intersectin cracks
						int cfound = 0;
						Vector norm;					// track normal vector for crack plane

						// Loop over potential crack crossings
#ifdef PREHASH_CRACKS
                        // When prehashed, we only need to check cracks
                        // seen be element elref
                        vector<int> *checkCracks;
                        int numCheck;
                        if(prehashed)
                        {   checkCracks = elref->SeesCrack;
                            numCheck = (int)checkCracks->size();
                        }
						else
						{	checkCracks = NULL;
							numCheck = numberOfCracks;
						}
#else
                        // We need to check all cracks
                        int numCheck = numberOfCracks;
#endif
                        CrackHeader *nextCrack;

                        // Paper step 3.b to loop of each possible cracks
                        for(int cn=0;cn<numCheck;cn++)
                        {
#ifdef PREHASH_CRACKS
                            nextCrack =  prehashed ? crackList[(*checkCracks)[cn]] : crackList[cn];
#else
                            nextCrack = crackList[cn];
#endif
                            // Paper step 3.b.i : trace line from particle to node and return
                            // 0 (no cross), 1 (cross from above), or 2 (cross from below)
							vfld = nextCrack->CrackCross(&(mpmptr->pos), &ndpt, &norm, nds[i]);
                            
                            // Paper step 3.b.ii : when vfld is not zero, save vfld (in loc),
                            // crack number, and normal. When done increment cfound.
							if(vfld!=NO_CRACK)
							{   cfld[cfound].loc = vfld;
								cfld[cfound].norm = norm;
#ifdef IGNORE_CRACK_INTERACTIONS
								// appears to always be same crack, and stop when found one
								cfld[cfound].crackNum = 1;
								break;
#endif
								// Get crack number (default code does not ignore interactions)
								cfld[cfound].crackNum = nextCrack->GetNumber();
								cfound++;
                                
                                // Paper 3.b.iii : exit loop if now have two crack because this code
                                // can only handle two interacting cracks. We exit to below and handle
                                // the first two cracks that were found.
								if(cfound>1) break;
							}
						}
                        
                        // Paper 3.c : if cfound=0 skip all crack velocity code and set v(p,i)=0
                        // But if cfound>0, we enter the code to assign crack velocity fields
						if(cfound>0)
						{   // Some stuff in below needs critical. Two options are to make it all critical
							// (use here comment out all pragma's inside the method) or comment out here and keep
							// all in the method
//#pragma omp critical (addcvf)
							{   try
								{   // This call code that implement the crack velocity field block
                                    // in the interacting cracks paper.
                                    vfld = ndptr->AddCrackVelocityField(matfld,cfld);
								}
								catch(std::bad_alloc&)
								{   if(initErr==NULL)
										initErr = new CommonException("Memory error","InitVelocityFieldsTask::Execute");
                                    break;
								}
								catch(...)
								{	if(initErr==NULL)
										initErr = new CommonException("Unexpected error","InitVelocityFieldsTask::Execute");
                                    break;
								}
							}
						}
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
						
					// set material point velocity field for this node
					mpmptr->vfld[i] = (char)vfld;
				}
				
				// next material point
				mpmptr = (MPMBase *)mpmptr->GetNextObject();
			}
		}
	}


	// was there an error?
	if(initErr!=NULL) throw initErr;
	
	// copy crack and material fields on real nodes to ghost nodes
	if(tp>1)
	{   for(int pn=0;pn<tp;pn++)
			patches[pn]->InitializationReduction();
	}
 
    return true;
}
