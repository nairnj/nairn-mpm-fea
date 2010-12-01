/********************************************************************************
    MatPtFluxBC.cpp
    NairnMPM
    
    Created by John Nairn on Wed Mar 17 2004.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Read_XML/mathexpr.hpp"

// global
MatPtFluxBC *firstFluxPt=NULL;

/*******************************************************************
	MatPtFluxBC: Constructors and Destructors
*******************************************************************/

// Constructors
MatPtFluxBC::MatPtFluxBC(int num,int dof,int sty) : MatPtLoadBC(num,dof,sty)
{
}

/*******************************************************************
	MatPtFluxBC: Methods
*******************************************************************/

// print to output
BoundaryCondition *MatPtFluxBC::PrintBC(ostream &os)
{
    char nline[200];
    
    sprintf(nline,"%5d %2d %2d %15.7e %15.7e",ptNum,direction,style,value,ftime);
    os << nline;
	PrintFunction(os);
	
    return (BoundaryCondition *)GetNextObject();
}

// set component of surface flux to zero (only called when diffusion is active)
MatPtFluxBC *MatPtFluxBC::ZeroMPFlux(void)
{
	mpm[ptNum-1]->pDiffusion->flux=0.;
    return (MatPtFluxBC *)GetNextObject();
}

// increment external load on a particle
// input is analysis time in seconds
// (only called when diffusion is active)
MatPtFluxBC *MatPtFluxBC::AddMPFlux(double bctime)
{
	double volume=mpm[ptNum-1]->volume;						// in mm^3
	
	if(style==SILENT)
	{	// silent assumes isotropic material
		int matnum=mpm[ptNum-1]->MatID();
		theMaterials[matnum]->LoadTransportProps(mpm[ptNum-1],fmobj->np);
		Tensor *D=theMaterials[matnum]->GetDiffusionTensor();			// in mm^2/sec
		
		// get flux force in pot mm^3/s
		double sign = value>=0. ? 1. : -1.;
		if(direction==X_DIRECTION)
			mpm[ptNum-1]->pDiffusion->flux+=sign*(D->xx*mpm[ptNum-1]->pDiffusion->Dc.x+D->xy*mpm[ptNum-1]->pDiffusion->Dc.y)*volume/mpmgrid.partx;
		else
			mpm[ptNum-1]->pDiffusion->flux+=sign*(D->xy*mpm[ptNum-1]->pDiffusion->Dc.x+D->yy*mpm[ptNum-1]->pDiffusion->Dc.y)*volume/mpmgrid.party;
	}
	else if(direction==EXTERNAL_FLUX)
	{	double mstime=1000.*bctime;
		mpm[ptNum-1]->pDiffusion->flux+=BCValue(mstime)*volume;
	}
	else
	{	varTime=mpm[ptNum-1]->pConcentration-ftime;
		GetPosition(&varXValue,&varYValue,&varZValue,&varRotValue);
		double currentValue=fabs(function->Val());
		if(varTime>0.) currentValue=-currentValue;
		mpm[ptNum-1]->pDiffusion->flux+=currentValue*volume;
	}
		
    return (MatPtFluxBC *)GetNextObject();
}
