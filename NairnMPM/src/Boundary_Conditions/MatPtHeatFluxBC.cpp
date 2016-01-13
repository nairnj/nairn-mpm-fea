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
	
	// for function input scale for Legacy units
	if(style==FUNCTION_VALUE)
		scale = UnitsController::Scaling(1.e3);
	
    return (BoundaryCondition *)GetNextObject();
}

// increment external heat flux on a particle
// input is analysis time in seconds
// (only called when conduction is active)
MatPtHeatFluxBC *MatPtHeatFluxBC::AddMPHeatFlux(double bctime)
{	
    // condition value
	// flux BC in nW/mm^2
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
		matptr->GetTransportProps(mpmptr,fmobj->np,&t);
		Tensor *k = &(t.kCondTensor);
        
        // k is k/rho0 (nJ mm^2/(sec-K-g)), Dt in K/mm, k Dt (nJ mm/(sec-g))
		if(fmobj->IsThreeD())
		{	fluxMag.x = k->xx*mpmptr->pTemp[gGRADx] + k->xy*mpmptr->pTemp[gGRADy] + k->xz*mpmptr->pTemp[gGRADz];
			fluxMag.y = k->xy*mpmptr->pTemp[gGRADx] + k->yy*mpmptr->pTemp[gGRADy] + k->yz*mpmptr->pTemp[gGRADz];
			fluxMag.z = k->xz*mpmptr->pTemp[gGRADx] + k->yz*mpmptr->pTemp[gGRADy] + k->zz*mpmptr->pTemp[gGRADz];
		}
		else
		{	fluxMag.x = k->xx*mpmptr->pTemp[gGRADx] + k->xy*mpmptr->pTemp[gGRADy];
			fluxMag.y = k->xy*mpmptr->pTemp[gGRADx] + k->yy*mpmptr->pTemp[gGRADy];
		}
		
		// remove 1/rho0 scaling on k to get nW/mm^2
		ScaleVector(&fluxMag,matptr->GetRho(mpmptr));
		
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
		fluxMag.x = scale*function->Val();
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
		conduction->AddFluxCondition(nd[nds[i]],DotVectors(&fluxMag,&tscaled)*fn[i],false);
	
    return (MatPtHeatFluxBC *)GetNextObject();
}

#pragma mark MatPtHeatFluxBC:Accessors

// set value (and scale legacy W/m^2 to nW/mm^2)
void MatPtHeatFluxBC::SetBCValue(double bcvalue)
{	BoundaryCondition::SetBCValue(UnitsController::Scaling(1.e3)*bcvalue);
}



