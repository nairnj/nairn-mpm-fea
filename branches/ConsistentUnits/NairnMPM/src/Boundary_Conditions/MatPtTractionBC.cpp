/********************************************************************************
	MatPtTractionBC.cpp
	nairn-mpm-fea

	Created by John Nairn on 9/13/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Nodes/NodalPoint.hpp"
#ifdef LOG_PROGRESS
#include "System/ArchiveData.hpp"
#endif

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
    
    sprintf(nline,"%7d %2d   %2d  %2d %15.7e %15.7e",ptNum,direction,face,style,value,ftime);
    os << nline;
	PrintFunction(os);
	
	// rescale
	value*=1.e6;		// Multiply by 1e6 to get N/mm^2 (kg-m/sec^2/mm^2) to g-mm/sec^2/mm^2
    scale*=1.e6;		// ... same in case using a function
	
    return (BoundaryCondition *)GetNextObject();
}

// increment external load on a particle
// input is analysis time in seconds
MatPtTractionBC *MatPtTractionBC::AddMPTraction(double bctime)
{
    // condition value
	double mstime=1000.*bctime;
	MPMBase *mpmptr = mpm[ptNum-1];
	double tmag = BCValue(mstime);
	
	// get corners and direction from material point
	// note 3D has four corners on face and 2D has 2
	// direction is x, y, z, N, or T1 (T1 not allowed in 3D)
	int cElem[4],numDnds;
	Vector corners[4],tscaled;
	double ratio = mpmptr->GetTractionInfo(face,direction,cElem,corners,&tscaled,&numDnds);
    
    // compact CPDI nodes into list of nodes (nds) and final shape function term (fn)
    // May need up to 8 (in 3D) for each of the numDnds (2 in 2D or 4 in 3D)
    int nds[8*numDnds+1];
    double fn[8*numDnds+1];
    int numCnds = CompactCornerNodes(numDnds,corners,cElem,ratio,nds,fn);
    
    // Particle information about field
	const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
	int matfld = matID->GetField();									// material velocity field
		
	// get crack velocity fields, if they are needed
	int numnds;
    int snds[maxShapeNodes];
	if(firstCrack!=NULL)
	{	const ElementBase *elref = theElements[mpmptr->ElemID()];		// element containing this particle
		double fn[maxShapeNodes];
		elref->GetShapeFunctions(&numnds,fn,snds,mpmptr);
	}
		
    // add force to each node
    Vector theFrc;
	short vfld = 0;
    for(int i=1;i<=numCnds;i++)
    {   // skip empty nodes
        if(nd[nds[i]]->NodeHasNonrigidParticles())
        {   // external force vector - tscaled has direction, surface area, and factor 1/2 (2D) or 1/4 (3D) to average the nodes
            CopyScaleVector(&theFrc,&tscaled,tmag*fn[i]);
			
			// Find the matching field
			if(firstCrack!=NULL)
			{	vfld = -1;
				for(int ii=1;ii<=numnds;ii++)
				{	if(nds[i] == snds[ii])
					{	vfld = mpmptr->vfld[ii];
						break;
					}
				}
			}
			
			// add to velocity field
            nd[nds[i]]->AddTractionTask3(mpmptr,vfld,matfld,&theFrc);
        }
    }
   
    // next boundary condition
    return (MatPtTractionBC *)GetNextObject();
}

#pragma mark MatPtTractionBC: Class Methods

// Calculate traction forces applied to particles and add to nodal force
void MatPtTractionBC::SetParticleSurfaceTractions(double stepTime)
{
    MatPtTractionBC *nextLoad=firstTractionPt;
    while(nextLoad!=NULL)
    	nextLoad=nextLoad->AddMPTraction(stepTime);
}

