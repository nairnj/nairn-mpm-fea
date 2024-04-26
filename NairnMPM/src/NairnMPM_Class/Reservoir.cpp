/*********************************************************************
	Reservoir.hpp
	Nairn Research Group MPM Code
	
	Created by John Nairn on 6/14/2021.
	Copyright (c) 2021, All rights reserved.
*********************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/Reservoir.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Patches/GridPatch.hpp"
#include "Exceptions/CommonException.hpp"
#include "Elements/ElementBase.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Materials/MaterialBase.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "Boundary_Conditions/MatPtHeatFluxBC.hpp"

#include  <iostream>

Reservoir *mpmReservoir = NULL;

#pragma mark Reservoir:Constructors and Destructors

Reservoir::Reservoir(Vector *centroid,Vector *elemSize)
{	store = *centroid;
	storeSize = *elemSize;
	count = 0;
}

// print details on the reservoir
void Reservoir::output(void)
{	cout << "Material points in the reservoir: " << count << endl;
	if(count==0) return;
	for(int i=0;i<matIDs.size();i++)
	{	cout << "   Material " << matIDs[i] << ", size=(";
		Vector lp = resSizes[i];
		lp.x *= storeSize.x;
		lp.y *= storeSize.y;
		cout << lp.x << "," << lp.y;
		if(fmobj->IsThreeD())
		{	lp.z *= storeSize.z;
			cout << "," << lp.z;
		}
		cout << ")" << endl;
	}
}

#pragma mark Reservoir:Methods

// Add material point to the reservoir
// Not thread safe because class data is changed
void Reservoir::AddParticle(MPMBase *mpmptr)
{
	// get the dimensionless size in reservoir
	Vector lp;
	if(mpmptr->InReservoir())
		mpmptr->GetDimensionlessSize(lp);
	else
	{	// adjust to reservoir element size
		lp = mpmptr->GetParticleSize();			// actual particle radius (z=1 in 2D)
		lp.x *= 2./storeSize.x;
		lp.y *= 2./storeSize.y;
		if(fmobj->IsThreeD()) lp.z *= 2./storeSize.z;
	}

	// check for material and size
	int matid = mpmptr->MatID()+1;		// convert to 1 based ID
	int i = IndexForClass(&lp,matid);
	
	// if no match and this call is deleting particle, look for material match and adjust size if found
	if(i<0 && !mpmptr->InReservoir())
	{	// check for material only
		i = IndexForClass(NULL,matid);
		if(i>=0)
		{	// found material, but need to change size to match class
			
			// change mass using sizes in the reservoir
			double rescale = (resSizes[i].x/lp.x)*(resSizes[i].y/lp.y);
			if(fmobj->IsThreeD()) rescale *= (resSizes[i].z/lp.z);
			mpmptr->mp *= rescale;

			// set size in current element to match resSize[i] in reservoir
			// (only matters for Tartan grid)
			ElementBase *elref = theElements[mpmptr->ElemID()];
			lp = resSizes[i];
			lp.x *= storeSize.x/elref->GetDeltaX();
			lp.y *= storeSize.y/elref->GetDeltaY();
			if(fmobj->IsThreeD()) lp.z *= storeSize.z/elref->GetDeltaZ();
			mpmptr->SetDimensionlessSize(&lp);
			
			// save the resizings - number and new absolute size
			resizings.push_back((double)mpmptr->GetNum());
			resizings.push_back(lp.x*elref->GetDeltaX());
			resizings.push_back(lp.y*elref->GetDeltaY());
			if(fmobj->IsThreeD()) resizings.push_back(lp.z*elref->GetDeltaZ());
		}
	}
	
	// Finally, create list (if i<0), start list (if lastPt[i]==NULL), or put at end of list i
	if(i<0)
	{	// create new list for size and material type
		resSizes.push_back(lp);
		matIDs.push_back(matid);
		firstPt.push_back(mpmptr);
		lastPt.push_back(mpmptr);
		mpmptr->SetNextObject(NULL);
	}
	else if(lastPt[i]==NULL)
	{	// insert into existing empty list
		firstPt[i] = mpmptr;
		lastPt[i] = mpmptr;
		mpmptr->SetNextObject(NULL);
	}
	else
	{	// put at end of existing list
		lastPt[i]->SetNextObject(mpmptr);
		mpmptr->SetNextObject(NULL);
		lastPt[i] = mpmptr;
	}

	// increment count of particles in the reservoir
	count++;
}

// Inject particle to calculations
// if lpart>=NULL or matid>0, requires material point of that class
//		lpart is actual size (can be NULL to ignore size and 2D only checks x and y size)
//		matid>0 to match material type (1 based) can be zero to accept any material)
// Warning: when called in any other parallel loop (such as loop of material points),
//		the call must in in an omp critical (delparticle) to avoid race condition
//		in code to add the particle to the reservoir
// Warning: do not call in patch loop
// throws CommonException() is location not valid location in the grid for a particle
MPMBase *Reservoir::InjectParticle(Vector *location,Vector *lpart,int matid)
{
	MPMBase *mptr = NULL;
	
	// is a specfic type wanted?
	if(lpart!=NULL || matid>0)
	{	// convert to dimensionless size in the reservoir
		Vector lp,*lpptr=NULL;
		if(lpart!=NULL)
		{	lp.x = lpart->x/storeSize.x;
			lp.y = lpart->y/storeSize.y;
			lp.z = fmobj->IsThreeD() ? lpart->z/storeSize.z : 1. ;
			lpptr = &lp;
		}
		int i = IndexForClass(lpptr,matid);
		
		// grab first particle in list with that size if that list has particles
		if(i>=0)
		{	if(firstPt[i]!=NULL)
			{	mptr = firstPt[i];
				firstPt[i] = (MPMBase *)mptr->GetNextObject();
				if(firstPt[i]==NULL) lastPt[i]=NULL;
			}
		}
		
		// if not found try other lists
		if(mptr==NULL)
		{	// get first particle of requested material type or any type if don't care
			for(i=0;i<firstPt.size();i++)
			{	if(firstPt[i]!=NULL && ((matid==matIDs[i]) || matid<=0))
				{	// found particle, extract and start this list with next particle
					mptr = firstPt[i];
					firstPt[i] = (MPMBase *)mptr->GetNextObject();
					if(firstPt[i]==NULL) lastPt[i]=NULL;
					break;
				}
			}

			// found one, but if wrong size, change the size now
			if(mptr!=NULL && lpart!=NULL)
			{	// get original size in the reservoir
				Vector lp0;
				mptr->GetDimensionlessSize(lp0);
				
				// set new size in reservoir and scale mass to accommodate it
				mptr->SetDimensionlessSize(&lp);
				double rescale = (lp.x/lp0.x)*(lp.y/lp0.y);
				if(fmobj->IsThreeD()) rescale *= (lp.z/lp0.z);
				mptr->mp *= rescale;
				
				// save the resizings - number and new absolute size
				resizings.push_back((double)mptr->GetNum());
				resizings.push_back(lpart->x);
				resizings.push_back(lpart->y);
				if(fmobj->IsThreeD()) resizings.push_back(lpart->z);
			}
		}
	}
	
	else
	{	// if don't care, first particle of any size and any material
		for(int i=0;i<=firstPt.size();i++)
		{	if(firstPt[i]!=NULL)
			{	// found particle, extract and start this list with next particle
				mptr = firstPt[i];
				firstPt[i] = (MPMBase *)mptr->GetNextObject();
				if(firstPt[i]==NULL) lastPt[i]=NULL;
				break;
			}
		}
	}

	// No particle found matching input conditions
	if(mptr==NULL) return NULL;
	
	// Set position
	mptr->SetPosition(location);
	mptr->SetOrigin(location);

	// set initial particle state (some materials need it)
	// Called instead of when deleted in case initial state depends on position
	theMaterials[mptr->MatID()]->SetInitialParticleState(mptr,fmobj->np,0);
	
	int j;
	try
	{   // get element (note this method throw exception if has left the grid)
		j = mpmgrid.FindElementFromPoint(location,mptr)-1;		// elem ID (0 based)
	}
	catch(...)
	{   throw CommonException("Injection location is off the grid","Reservoir::InjectParticle");
	}
	
	// move to an edge element, GIMP has effective left the grid
	if(theElements[j]->OnTheEdge())
		throw CommonException("Injection location on the edge of the grid","Reservoir::InjectParticle");
	
	// IF axsymmetric, cannot got to element with r<0 (even if not on the edge)
	if(fmobj->IsAxisymmetric() && location->x<=0.)
		throw CommonException("Injection location r<0 in axisymmetric grid","Reservoir::InjectParticle");

	// change particle to new element
	mptr->ChangeElemID(j,!mpmgrid.IsStructuredEqualElementsGrid());
	
	// add to a patch
	int pn = mpmgrid.GetPatchForElement(mptr->ElemID());
	patches[pn]->AddParticle(mptr);
	
	// return particle
	count--;
	return mptr;
}

// Delete particle and add it to the reservoir
// Set enough to zero and add to start of list with matching size
// Note: do not call this in loops over patches because particle deletion changes
//		pointers needed by patch looping. Instead, build list of particle to delete
//		and then delete then in the reduction phase for that loop
//		(see ResetElementsTask() for an example of this process)
// Warning: when called in any other parallel loop (such as loop of material points),
//		the call must in in an omp critical (delparticle) to avoid race condition
//		in code to add the particle to the reservoir
void Reservoir::DeleteParticle(MPMBase *mpmptr)
{
	// remove it from current patch, but only if not called from patches loop
	int pn = mpmgrid.GetPatchForElement(mpmptr->ElemID());
	patches[pn]->RemoveParticle(mpmptr);
	AddParticle(mpmptr);
	
	// reset particle
	// clear all, but not particle spin or history data
	mpmptr->ResetMaterialPoint();
	
	// clear history data (always call SetConcentration in case other diffusions are active)
    double resetConc = diffusion!=NULL ? diffusion->reference : 0.;
	mpmptr->SetConcentration(resetConc,true);
	mpmptr->SetTemperature(thermal.reference,thermal.reference);
    
    // Default is to zero NumberOfHistoryDoubles(). If needed a material can override
    // and do other history reset tasks (such as those that should not be zero or
    // when history is not a list of doubles.
	theMaterials[mpmptr->MatID()]->ResetHistoryData(mpmptr->GetHistoryPtr(0),mpmptr);
	
	// Remove particle BC on the deleted particle
	MatPtLoadBC::DeleteParticleBCs(mpmptr,&firstLoadedPt);
	MatPtLoadBC::DeleteParticleBCs(mpmptr,(MatPtLoadBC **)&firstTractionPt);
    for(int i=0;i<NUM_DUFFUSION_OPTIONS;i++)
        MatPtLoadBC::DeleteParticleBCs(mpmptr,(MatPtLoadBC **)&firstDiffFluxBC[i]);
	MatPtLoadBC::DeleteParticleBCs(mpmptr,(MatPtLoadBC **)&firstHeatFluxPt);

	// move it to element 1 and update elements list (if active)
	mpmptr->SetOrigin(&store);
	mpmptr->SetPosition(&store);
	mpmptr->ChangeElemID(0,!mpmgrid.IsStructuredEqualElementsGrid());
	
	// pTemp is for thermal gradients and zeroed each time step when needed
	
}

// If code wants to change particle size, it should call this method
// with the new actual size. This method scales mass by mass change
// and registers the sie change to be written to PtDims.txt file at
// the next archive
void Reservoir::ChangeParticleSize(MPMBase *mptr,Vector *lpart)
{	// get original size in the reservoir
	Vector lp0;
	mptr->GetDimensionlessSize(lp0);
	
	// get new size in current element
	ElementBase *elref = theElements[mptr->ElemID()];
	Vector lp;
	lp.x = lpart->x/elref->GetDeltaX();
	lp.y = lpart->y/elref->GetDeltaY();
	lp.z = fmobj->IsThreeD() ? lpart->z/elref->GetDeltaZ() : 1. ;
	mptr->SetDimensionlessSize(&lp);
	
	// scale the mass
	double rescale = (lp.x/lp0.x)*(lp.y/lp0.y);
	if(fmobj->IsThreeD()) rescale *= (lp.z/lp0.z);
	mptr->mp *= rescale;
	
	// save the resizings - number and new absolute size
	resizings.push_back((double)mptr->GetNum());
	resizings.push_back(lpart->x);
	resizings.push_back(lpart->y);
	if(fmobj->IsThreeD()) resizings.push_back(lpart->z);
}

#pragma mark Reservoir:Accessors

// Return index particle by class or -1 if none found
// If lp==NULL, match first that matches material ID
// If matid<=0, match first that matches size
// If both provided, try to match both
// If none found return -1
int Reservoir::IndexForClass(Vector *lp,int matid) const
{	// give up if nothing to match
	if(lp==NULL && matid<=0) return -1;
	
	for(int j=0;j<resSizes.size();j++)
	{	// does it match size
		if(lp!=NULL)
		{	Vector lstore = resSizes[j];
			if(!DbleEqual(lp->x,lstore.x) || !DbleEqual(lp->y,lstore.y)) continue;
			if(fmobj->IsThreeD() && !DbleEqual(lp->z,lstore.z)) continue;
		}
		
		// here whenever size matched or was not checked
		
		// if material not checked, it is a match
		if(matid<=0) return j;
		
		// does material match?
		if(matid==matIDs[j]) return j;
	}
	
	// no matches
	return -1;
}

// running count of particles
long Reservoir::currentCount(void) const { return count; }

// get list of resizings
vector <double> Reservoir::GetResizings(void) { return resizings; }
void Reservoir::ClearResizings(void) { resizings.clear(); }
