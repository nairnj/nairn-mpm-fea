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
#include "Nodes/NodalPoint.hpp"				// +AS

// global
MatPtFluxBC *firstFluxPt=NULL;

/*******************************************************************
	MatPtFluxBC: Constructors and Destructors
*******************************************************************/

// Constructors
MatPtFluxBC::MatPtFluxBC(int num,int dof,int sty,int edge) : MatPtLoadBC(num,dof,sty)
{
	face = edge;
}

/*******************************************************************
	MatPtFluxBC: Methods
*******************************************************************/

// print to output
BoundaryCondition *MatPtFluxBC::PrintBC(ostream &os)
{
    char nline[200];
    
    sprintf(nline,"%7d %2d   %2d  %2d %15.7e %15.7e",ptNum,face,direction,style,value,ftime);
    os << nline;
	PrintFunction(os);
	
    return (BoundaryCondition *)GetNextObject();
}

// increment external load on a particle
// input is analysis time in seconds
// (only called when diffusion is active)
MatPtFluxBC *MatPtFluxBC::AddMPFlux(double bctime)
{
	/*
	double volume = mpm[ptNum-1]->GetVolume(DEFORMED_VOLUME);						// in mm^3
	
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
		mpm[ptNum-1]->pDiffusion->flux += BCValue(mstime)*volume;
	}
	else
    {   // coupled surface flux and ftime is bath concentration
		varTime = mpm[ptNum-1]->pConcentration-ftime;
		GetPosition(&varXValue,&varYValue,&varZValue,&varRotValue);
		double currentValue = fabs(function->Val());
		if(varTime>0.) currentValue=-currentValue;
		mpm[ptNum-1]->pDiffusion->flux += currentValue*volume;
	}
	*/
		
    // condition value
	MPMBase *mpmptr = mpm[ptNum-1];
	double fluxMag;
	
	if(direction==EXTERNAL_FLUX)
	{	double mstime=1000.*bctime;
		fluxMag = BCValue(mstime);
	}
	else
    {   // coupled surface flux and ftime is bath concentration
		varTime = mpmptr->pConcentration-ftime;
		GetPosition(&varXValue,&varYValue,&varZValue,&varRotValue);
		double currentValue = fabs(function->Val());
		if(varTime>0.) currentValue=-currentValue;
		fluxMag = currentValue;
	}
	
	// fluxMag in g/(mm^2-sec) - find J/rho in units of mm/sec
	MaterialBase *matptr = theMaterials[mpmptr->MatID()];
	double rho = 0.001*matptr->rho*matptr->concSaturation/mpmptr->GetRelativeVolume();
	fluxMag /= rho;
	
	// get corners and direction from material point
	int cElem[4],numDnds;
	Vector corners[4],tscaled;
	double ratio = mpmptr->GetTractionInfo(face,X_DIRECTION,cElem,corners,&tscaled,&numDnds);
	
    // compact CPDI nodes into list of nodes (nds) and final shape function term (fn)
    // May need up to 8 (in 3D) for each of the numDnds (2 in 2D or 4 in 3D)
    int nds[8*numDnds+1];
    double fn[8*numDnds+1];
    int numnds = CompactCornerNodes(numDnds,corners,cElem,ratio,nds,fn);
	
    // add force to each node
	int i;
    for(i=1;i<=numnds;i++)
    {   // skip empty nodes
        if(nd[nds[i]]->NumberNonrigidParticles())
		{	nd[nds[i]]->fdiff += fluxMag*tscaled.x*fn[i];
        }
    }
	
    return (MatPtFluxBC *)GetNextObject();
}
