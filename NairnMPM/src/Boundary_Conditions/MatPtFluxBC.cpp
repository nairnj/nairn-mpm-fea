/********************************************************************************
    MatPtFluxBC.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Mar 17 2004.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Read_XML/mathexpr.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "System/UnitsController.hpp"

// global
MatPtFluxBC *firstFluxPt=NULL;

#pragma mark MatPtFluxBC: Constructors and Destructors

// Constructors
MatPtFluxBC::MatPtFluxBC(int num,int dof,int sty,int edge) : MatPtLoadBC(num,dof,sty)
{
	face = edge;
}

#pragma mark MatPtFluxBC: Methods

// print to output
BoundaryCondition *MatPtFluxBC::PrintBC(ostream &os)
{
    char nline[200];
    
    sprintf(nline,"%7d %2d   %2d  %2d %15.7e %15.7e",ptNum,direction,face,style,
			UnitsController::Scaling(1.e3)*GetBCValueOut(),GetBCFirstTimeOut());
    os << nline;
	PrintFunction(os);
	
	// for function input scale for Legacy units
	if(style==FUNCTION_VALUE)
		scale = UnitsController::Scaling(1.e-3);
	
    return (BoundaryCondition *)GetNextObject();
}

// increment external load on a particle
// input is analysis time in seconds
// (only called when diffusion is active)
MatPtLoadBC *MatPtFluxBC::AddMPFluxBC(double bctime)
{
    // condition value is g/(mm^2-sec), Divide by rho*csat to get potential flux in mm/sec
	// find this flux and then add (times area) to get mm^3-potential/sec
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
		Tensor *D = &(t.diffusionTensor);
        
        // D in mm^2/sec, Dc in 1/mm
		if(fmobj->IsThreeD())
		{	fluxMag.x = D->xx*mpmptr->pDiffusion[gGRADx] + D->xy*mpmptr->pDiffusion[gGRADy] + D->xz*mpmptr->pDiffusion[gGRADz];
			fluxMag.y = D->xy*mpmptr->pDiffusion[gGRADx] + D->yy*mpmptr->pDiffusion[gGRADy] + D->yz*mpmptr->pDiffusion[gGRADz];
			fluxMag.z = D->xz*mpmptr->pDiffusion[gGRADx] + D->yz*mpmptr->pDiffusion[gGRADy] + D->zz*mpmptr->pDiffusion[gGRADz];
		}
		else
		{	fluxMag.x = D->xx*mpmptr->pDiffusion[gGRADx] + D->xy*mpmptr->pDiffusion[gGRADy];
			fluxMag.y = D->xy*mpmptr->pDiffusion[gGRADx] + D->yy*mpmptr->pDiffusion[gGRADy];
		}
        bcDir = N_DIRECTION;
	}
	else if(direction==EXTERNAL_FLUX)
	{	// csatrho = rho0 V0 csat/V (units g/mm^3)
		double csatrho = mpmptr->GetRho()*mpmptr->GetConcSaturation()/mpmptr->GetRelativeVolume();
		fluxMag.x = BCValue(bctime)/csatrho;
	}
	else
    {   // coupled surface flux and ftime is bath concentration
		// time variable (t) is replaced by c-cbath, where c is the particle potential and cbath is bath potential
		varTime = mpmptr->pPreviousConcentration-GetBCFirstTime();
		GetPosition(&varXValue,&varYValue,&varZValue,&varRotValue);
		double currentValue = fabs(scale*function->Val());
		if(varTime>0.) currentValue=-currentValue;
		
		// csatrho = rho0 V0 csat/V (units g/mm^3)
		double csatrho = mpmptr->GetRho()*mpmptr->GetConcSaturation()/mpmptr->GetRelativeVolume();
		fluxMag.x = currentValue/csatrho;
	}
	
	// get corners and direction from material point
	int cElem[4],numDnds;
	Vector corners[4],tscaled;
	double ratio = mpmptr->GetTractionInfo(face,bcDir,cElem,corners,&tscaled,&numDnds);
	
    // compact CPDI nodes into list of nodes (nds) and final shape function term (fn)
    // May need up to 8 (in 3D) for each of the numDnds (2 in 2D or 4 in 3D)
#ifdef CONST_ARRAYS
	int nds[8*4+1];
	double fn[8*4+1];
#else
	int nds[8*numDnds+1];
    double fn[8*numDnds+1];
#endif
    int numnds = CompactCornerNodes(numDnds,corners,cElem,ratio,nds,fn);
	
    // add force to each node
	int i;
    for(i=1;i<=numnds;i++)
		diffusion->AddFluxCondition(nd[nds[i]],DotVectors(&fluxMag,&tscaled)*fn[i],false);
	
    return (MatPtFluxBC *)GetNextObject();
}

#pragma mark MatPtFluxBC:Accessors

// set value (and scale legacy kg/(m^2-sec) to g/(mm^2-sec))
void MatPtFluxBC::SetBCValue(double bcvalue)
{	BoundaryCondition::SetBCValue(UnitsController::Scaling(1.e-3)*bcvalue);
}

// check coupled flux which uses first time as unscaled concentration potential
void MatPtFluxBC::SetBCFirstTime(double bcftime)
{	if(direction==COUPLED_FLUX)
		BoundaryCondition::SetBCFirstTimeCU(bcftime);
	else
		BoundaryCondition::SetBCFirstTime(bcftime);
}

// check coupled flux which uses first time as unscaled concentration potential
double MatPtFluxBC::GetBCFirstTimeOut(void)
{	if(direction==COUPLED_FLUX)
		return BoundaryCondition::GetBCFirstTime();
	else
		return BoundaryCondition::GetBCFirstTimeOut();
}


