/********************************************************************************
    MatPtHeatFluxBC.cpp
    nairn-mpm-fea

    Created by John Nairn on May 28, 2013.
    Copyright (c) 2013 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/MatPtHeatFluxBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Read_XML/mathexpr.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "System/UnitsController.hpp"

// global
MatPtHeatFluxBC *firstHeatFluxPt=NULL;

#pragma mark MatPtHeatFluxBC: Constructors and Destructors

// Constructors
MatPtHeatFluxBC::MatPtHeatFluxBC(int num,int dof,int sty,int edge) : MatPtLoadBC(num,dof,sty)
{
	face = edge;
}

#pragma mark MatPtHeatFluxBC: Methods

// print to output
BoundaryCondition *MatPtHeatFluxBC::PrintBC(ostream &os)
{
    char nline[200];
    
    sprintf(nline,"%7d %2d   %2d  %2d %15.7e %15.7e",ptNum,direction,face,style,
			UnitsController::Scaling(1.e-3)*GetBCValueOut(),GetBCFirstTimeOut());
    os << nline;
	PrintFunction(os);
	
    return (BoundaryCondition *)GetNextObject();
}

// increment external heat flux on a particle
// input is analysis time in seconds
// (only called when conduction is active)
MatPtHeatFluxBC *MatPtHeatFluxBC::AddMPHeatFlux(double bctime)
{	
    // condition value
	// flux BC in W/m^2 = N/(m-sec), but need mJ/(mm^2-sec) = N/(mm-sec) = (1/1000) N/(m-sec)
	MPMBase *mpmptr = mpm[ptNum-1];
    MaterialBase *matptr = theMaterials[mpmptr->MatID()];
	
	// Flux is a scalar and we need int_(face) F Ni(x) dA
	// Since F is constant, only need integral which is done by CPDI methods
	//		which has be generalized to work for GIMP too
	// We use X_DIRECTION for bcDIR for efficiency. For Silent BC, change to
	//      Normal direction to all caculation of n
	Vector fluxMag;
	ZeroVector(&fluxMag);
    int bcDir=X_DIRECTION;
	
	if(style==SILENT)
	{	TransportProperties t;
		matptr->GetTransportProps(mpm[ptNum-1],fmobj->np,&t);
		Tensor *k = &(t.kCondTensor);
        
        // k is k/rho0 (nJ mm^2/(sec-K-g)), Dt in K/mm, k Dt (nJ mm/(sec-g))
		if(fmobj->IsThreeD())
		{	fluxMag.x = k->xx*mpmptr->pTemp->DT.x + k->xy*mpmptr->pTemp->DT.y + k->xz*mpmptr->pTemp->DT.z;
			fluxMag.y = k->xy*mpmptr->pTemp->DT.x + k->yy*mpmptr->pTemp->DT.y + k->yz*mpmptr->pTemp->DT.z;
			fluxMag.x = k->xz*mpmptr->pTemp->DT.x + k->yz*mpmptr->pTemp->DT.y + k->zz*mpmptr->pTemp->DT.z;
		}
		else
		{	fluxMag.x = k->xx*mpmptr->pTemp->DT.x + k->xy*mpmptr->pTemp->DT.y;
			fluxMag.y = k->xy*mpmptr->pTemp->DT.x + k->yy*mpmptr->pTemp->DT.y;
		}
		
		// remove 1/rho0 scaling on k to get N/(mm-sec)
		ScaleVector(&fluxMag,matptr->rho);
		
		// need to get normal vector from cpdi functions below
        bcDir = N_DIRECTION;
	}
	else if(direction==EXTERNAL_FLUX)
	{	// user-supplied constant
		fluxMag.x = BCValue(bctime);
	}
	else
    {   // coupled surface flux
		// time variable (t) is replaced by particle temperature
		varTime = mpmptr->pPreviousTemperature;
		GetPosition(&varXValue,&varYValue,&varZValue,&varRotValue);
		
		// Legacy scaling of W/m^2 to nW/mm^2
		fluxMag.x = UnitsController::Scaling(1.e3)*function->Val();
	}
	
	// get corners and direction from material point
	int cElem[4],numDnds;
	Vector corners[4],tscaled;
	double ratio = mpmptr->GetTractionInfo(face,bcDir,cElem,corners,&tscaled,&numDnds);
	
    // compact CPDI nodes into list of nodes (nds) and final shape function term (fn)
    // May need up to 8 (in 3D) for each of the numDnds (2 in 2D or 4 in 3D)
    int nds[8*numDnds+1];
    double fn[8*numDnds+1];
    int numnds = CompactCornerNodes(numDnds,corners,cElem,ratio,nds,fn);
	
    // add force to each node
	// tscaled has units mm^s for final flux is (N-mm)/sec
	int i;
    for(i=1;i<=numnds;i++)
    {   // skip empty nodes
        if(nd[nds[i]]->NodeHasNonrigidParticles())
			conduction->AddFluxCondition(nd[nds[i]],DotVectors(&fluxMag,&tscaled)*fn[i],false);
    }
	
    return (MatPtHeatFluxBC *)GetNextObject();
}

#pragma mark MatPtHeatFluxBC:Accessors

// set value (and scale legacy W/m^2 to nW/mm^2)
void MatPtHeatFluxBC::SetBCValue(double bcvalue)
{	BoundaryCondition::SetBCValue(UnitsController::Scaling(1.e3)*bcvalue);
}



