/********************************************************************************
	MatPtTractionBC.cpp
	NairnMPM

	Created by John Nairn on 9/13/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Nodes/NodalPoint.hpp"

// global
MatPtTractionBC *firstTractionPt=NULL;

#pragma mark MatPtTractionBC::Constructors and Destructors

// Constructors
MatPtTractionBC::MatPtTractionBC(int num,int dof,int edge,int sty)
							: MatPtLoadBC(num,dof,sty)
{
	face = edge;
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
    int i,j;
    
	double mstime=1000.*bctime;
	MPMBase *mpmptr = mpm[ptNum-1];
	double tmag = BCValue(mstime);
	
	// get corners and direction from material point
	int cElem[4],numDnds;
	Vector corners[4],tscaled;
	mpmptr->GetTractionInfo(face,direction,cElem,corners,&tscaled,&numDnds);
	
	// loop over corner finding all nodes and add to fext
    int numnds,nds[8*numDnds+1],ncnds=0;
    double fn[8*numDnds+1],cnodes[8*numDnds],twt[8*numDnds];              // allows 3D which can have 8 nodes each
    for(i=0;i<numDnds;i++)
	{	// get straight grid shape functions
		theElements[cElem[i]]->GridShapeFunctions(&numnds,nds,&corners[i],fn);
		
		// loop over shape grid shape functions and collect in arrays
		for(j=1;j<=numnds;j++)
		{   cnodes[ncnds] = nds[j];
			twt[ncnds] = fn[j];
			ncnds++;
		}
	}
    
    /*
    cout << "# Initial:" << endl;
    for(i=0;i<ncnds;i++)
    {   cout << "# node = " << cnodes[i] << ", Si = " << twt[i] << endl;
    }
    */
    
 	// shell sort by node numbers in cnodes[] (always 16 for linear CPDI)
	int lognb2=(int)(log((double)ncnds)*1.442695022+1.0e-5);	// log base 2
	int k=ncnds,l,cmpNode;
	double cmpSi;
	for(l=1;l<=lognb2;l++)
	{	k>>=1;		// divide by 2
		for(j=k;j<ncnds;j++)
		{	i=j-k;
			cmpNode = cnodes[j];
			cmpSi = twt[j];
			
			// Back up until find insertion point
			while(i>=0 && cnodes[i]>cmpNode)
			{	cnodes[i+k] = cnodes[i];
				twt[i+k] = twt[i];
				i-=k;
			}
			
			// Insert point
			cnodes[i+k]=cmpNode;
			twt[i+k]=cmpSi;
		}
	}
    
    /*
    cout << "# Sorted:" << endl;
    for(i=0;i<ncnds;i++)
    {   cout << "# node = " << cnodes[i] << ", Si = " << twt[i] << endl;
    }
    */

 	// compact same node number
	int count = 0;
	nds[0] = -1;
	fn[0] = 1.;
	for(i=0;i<ncnds;i++)
	{   if(cnodes[i] == nds[count])
        {   fn[count] += twt[i];
        }
        else
        {	if(fn[count]>1.e-10) count++;       // keep only if shape is nonzero
            nds[count] = cnodes[i];
            fn[count] = twt[i];
        }
	}
	if(fn[count]<1.e-10) count--;
	numnds = count;
    
    /*
    cout << "# Compacted: tmag = " << tmag;
    PrintVector(", t vec =",&tscaled);
    cout << endl;
    for(i=1;i<=numnds;i++)
    {   cout << "# node = " << nds[i] << ", Total Si = " << fn[i] << endl;
    }
    */
    
    // Particle information about field
    int vfld=0;                                             // To support traction near cracks need to calculate for each node
    MaterialBase *matID=theMaterials[mpmptr->MatID()];		// material class object
    int matfld=matID->GetField();                           // material field
    Vector theFrc;
    
    // add force to each node
    for(i=1;i<=numnds;i++)
    {   // external force vector
        CopyScaleVector(&theFrc,&tscaled,tmag*fn[i]);
        cout << "# ... node = " << nds[i];
        PrintVector(", F = ",&theFrc);
        cout << endl;
        nd[nds[i]]->AddFextTask3(vfld,matfld,&theFrc);
    }
   
    // next boundary condition
    return (MatPtTractionBC *)GetNextObject();
}

#pragma mark MatPtTractionBC: Class Methods

// Calculate traction forces applied to particles and add to nodal fext
void MatPtTractionBC::SetParticleSurfaceTractions(double stepTime)
{
    MatPtTractionBC *nextLoad=firstTractionPt;
    while(nextLoad!=NULL)
    	nextLoad=nextLoad->AddMPTraction(stepTime);
}

