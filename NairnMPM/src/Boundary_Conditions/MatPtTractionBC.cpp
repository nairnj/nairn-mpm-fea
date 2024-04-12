/********************************************************************************
	MatPtTractionBC.cpp
	nairn-mpm-fea

	Created by John Nairn on 9/13/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/CommonException.hpp"

// Local globals for BCs
// xi1 starts with x12i[1] and xi2 starts with x12i[0]]
static double x12i[5]={-1.,-1.,1.,1.,-1.};

// global
MatPtTractionBC *firstTractionPt=NULL;

#pragma mark MatPtTractionBC::Constructors and Destructors

// Constructors
MatPtTractionBC::MatPtTractionBC(int num,int dof,int edge,int sty)
							: MatPtLoadBC(num,dof,sty)
{
	face = edge;
	// direction set in super class can be 11, 12, 13 for normal, tangential 1, or tangential 2 loading
}

#pragma mark MatPtTractionBC: Methods

// print to output
BoundaryCondition *MatPtTractionBC::PrintBC(ostream &os)
{
    char nline[200];
	size_t nlsize=200;
	
	snprintf(nline,nlsize,"%7d %2d   %2d  %2d %15.7e %15.7e",ptNum,direction,face,style,
			UnitsController::Scaling(1.e-6)*GetBCValueOut(),GetBCFirstTimeOut());
    os << nline;
	PrintFunction(os);
	
	// not allowed on rigid particles
	MaterialBase *matref = theMaterials[mpm[ptNum-1]->MatID()];
	if(matref->IsRigid())
	{	throw CommonException("Cannot set traction force on rigid particles",
							  "MatPtTractionBC::PrintBC");
	}
	
 	// scale function output
	if(style==FUNCTION_VALUE) scale = UnitsController::Scaling(1.e6);
	
	return (BoundaryCondition *)GetNextObject();
}

// increment external load on a particle
// input is analysis time in seconds
MatPtLoadBC *MatPtTractionBC::AddMPFluxBC(double bctime)
{
    // condition value
	MPMBase *mpmptr = mpm[ptNum-1];
	double tmag = BCValue(bctime);

	// Particle information about field
	const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
	int matfld = matID->GetField();									// material velocity field

    // only in 2D and when activated
	if(fmobj->np!=THREED_MPM && fmobj->exactTractions)
	{	int cElem[4],numElements=4;
		Vector ends[8],nscale[4];
		mpmptr->GetExactTractionInfo(face,direction,cElem,ends,nscale,&numElements);
		
		int j=0;
		double rmid=1.,dr=1.;
		bool isAxisymmetric = fmobj->np==AXISYMMETRIC_MPM;
		double thick = isAxisymmetric ? 1. : mpmptr->thickness();
		Vector xi1,xi2;
		int ndsT[5];
		double fnT[4];

		for(int e=0;e<numElements;e++)
		{	// element containing endpoint j and j+1
			const ElementBase *elref = theElements[cElem[e]];		// element containing this particle
			
			// Get dimensionless points
			elref->GetXiPos(&ends[j],&xi1);
			elref->GetXiPos(&ends[j+1],&xi2);
			
			// get radial terms if needed
			if(isAxisymmetric)
			{	rmid = 0.5*(ends[j].x+ends[j+1].x);
				dr = 0.5*(ends[j+1].x-ends[j].x);
			}
		
			// Get shape functions
			elref->GridTractionFunction(&xi1,&xi2,isAxisymmetric,fnT,ndsT,rmid,dr);
			Vector theFrc;
			for(int i=1;i<=ndsT[0];i++)
			{	// external force vector - nscale has direction and length, scale by function and thickness
				CopyScaleVector(&theFrc,&nscale[e],thick*tmag*fnT[i-1]);
				
				// may need changes to handle cracks
				
				short vfld = 0;
				nd[ndsT[i]]->AddTractionTask3(mpmptr,vfld,matfld,&theFrc);
			}
			
			// next corner pair
			j += 2;
		}
		
		// next boundary condition
		return (MatPtTractionBC *)GetNextObject();
	}
	
	// otherwise try to continue with regular tractions

	// get corners and radii from deformed material point (2 in 2D and 4 in 3D)
	int cElem[4],numDnds;
	Vector corners[4],radii[3];
	double redge;
	Vector wtNorm = mpmptr->GetSurfaceInfo(face,direction,cElem,corners,radii,&numDnds,&redge);
	
#ifdef CONST_ARRAYS
	int nds[MAX_SHAPE_NODES],sndsArray[MAX_SHAPE_NODES];
	double fn[MAX_SHAPE_NODES];
#else
	int nds[maxShapeNodes],sndsArray[maxShapeNodes];
	double fn[maxShapeNodes];
#endif
	int numnds;
	
	// get crack velocity fields, if they are needed
	int cnumnds,*snds=NULL;
	if(firstCrack!=NULL)
	{	const ElementBase *elref = theElements[mpmptr->ElemID()];		// element containing this particle

		// The list of GIMP nodes will let code below associate a node with a velocity field
		// Actually shape functions are not needed (but function needs them to find non-zero nodes)
		snds = sndsArray;
		elref->GetShapeFunctions(fn,&snds,mpmptr);
		numnds = snds[0];
	}

	double efffn;
	for(int c=0;c<numDnds;c++)
	{	// get straight grid shape functions
		theElements[cElem[c]]->GetShapeFunctionsForTractions(fn,nds,&corners[c]);
		cnumnds = nds[0];
		
		// track total force due to each corner and rescale in a second pass if needed
		double netshape = 0.;
		double shapeScale = 1.;
		for(int pass=0;pass<1;pass++)
		{	// add force to each node
			short vfld = 0;
			for(int i=1;i<=cnumnds;i++)
			{   // skip empty nodes
				if(nd[nds[i]]->NodeHasNonrigidParticles())
				{   // external force vector
					Vector theFrc;
					ZeroVector(&theFrc);		// otherwise VS optimization fails
					if(fmobj->IsAxisymmetric())
					{	// wtNorm has direction and |r1|. Also scale by axisymmetric term
						efffn = fn[i]*shapeScale;
						CopyScaleVector(&theFrc,&wtNorm,(redge-x12i[c+1]*radii[0].x/3.)*tmag*efffn);
					}
					else
					{	// wtNorm has direction and Area/2 (2D) or Area/4 (3D) to average the nodes
						efffn = fn[i]*shapeScale;
						CopyScaleVector(&theFrc,&wtNorm,tmag*efffn);
					}
				
					// Find the matching velocity field
					bool addedForce = false;
					if(firstCrack!=NULL)
					{	vfld = -1;
						for(int ii=1;ii<=numnds;ii++)
						{	if(nds[i] == snds[ii])
							{	vfld = mpmptr->vfld[ii];
								if(nd[nds[i]]->AddTractionTask3(mpmptr,vfld,matfld,&theFrc))
									addedForce = true;
								break;
							}
						}
						
						// no field possible if highly deformed particle with uGIMP shape functions
						//  add force to field 0 in this case
						if(vfld<0)
						{	vfld = 0;
							if(nd[nds[i]]->AddTractionTask3(mpmptr,vfld,matfld,&theFrc))
								addedForce = true;
						}
					}
					else
					{	if(nd[nds[i]]->AddTractionTask3(mpmptr,vfld,matfld,&theFrc))
							addedForce = true;
					}
					
					// if added force track net shape functions
					if(addedForce) netshape += efffn;
				}
			}
			
			// skip second pass if caught them on first pass
			if(fabs(netshape-1.0)<1.e-6) break;
			
			// normalize to get correct total force from this corner
			// redo each one with effn/netshape so new sum will be 1
			// but we already added effn, so now add effn/netshape - effn = (1/netshape-1)*effn
			shapeScale = 1./netshape - 1.;
			
			// should write this as a warning
			//cout << "# rescale corner " << c << " of particle at (" << mpmptr->pos.x << "," << mpmptr->pos.y << ")"
			//			" by " << shapeScale << endl;
		}
	}

    // next boundary condition
    return (MatPtTractionBC *)GetNextObject();
}

#pragma mark MatPtTractionBC: Class Methods

// Calculate traction forces applied to particles and add to nodal force
void MatPtTractionBC::SetParticleSurfaceTractions(double stepTime)
{
    MatPtLoadBC *nextLoad = firstTractionPt;
    while(nextLoad!=NULL)
    	nextLoad = nextLoad->AddMPFluxBC(stepTime);
}

