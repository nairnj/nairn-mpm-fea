/********************************************************************************
    MatPtFluxBC.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Mar 17 2004.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Read_XML/mathexpr.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"

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
    
    sprintf(nline,"%7d %2d   %2d  %2d %15.7e %15.7e",ptNum,direction,face,style,GetBCValueOut(),GetBCFirstTimeOut());
    os << nline;
	PrintFunction(os);
	
    return (BoundaryCondition *)GetNextObject();
}

// increment external load on a particle
// input is analysis time in seconds
// (only called when diffusion is active)
MatPtFluxBC *MatPtFluxBC::AddMPFlux(double bctime)
{
    // condition value
	// flux BC in kg/(m^2-sec) - find Flux/rho in units of mm/sec
	MPMBase *mpmptr = mpm[ptNum-1];
    MaterialBase *matptr = theMaterials[mpmptr->MatID()];
	double csatrho;
	
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
		Tensor *D = &(t.diffusionTensor);
        
        // D in mm^2/sec, Dc in 1/mm
		if(fmobj->IsThreeD())
		{	fluxMag.x = D->xx*mpmptr->pDiffusion->Dc.x + D->xy*mpmptr->pDiffusion->Dc.y + D->xz*mpmptr->pDiffusion->Dc.z;
			fluxMag.y = D->xy*mpmptr->pDiffusion->Dc.x + D->yy*mpmptr->pDiffusion->Dc.y + D->yz*mpmptr->pDiffusion->Dc.z;
			fluxMag.x = D->xz*mpmptr->pDiffusion->Dc.x + D->yz*mpmptr->pDiffusion->Dc.y + D->zz*mpmptr->pDiffusion->Dc.z;
		}
		else
		{	fluxMag.x = D->xx*mpmptr->pDiffusion->Dc.x + D->xy*mpmptr->pDiffusion->Dc.y;
			fluxMag.x = D->xy*mpmptr->pDiffusion->Dc.x + D->yy*mpmptr->pDiffusion->Dc.y;
		}
        bcDir = N_DIRECTION;
	}
	else if(direction==EXTERNAL_FLUX)
	{	// csatrho = rho V csat/V0 = solvent mass per reference volume
		// units are 1000 kg mm^3/(m^2-g-sec) = mm/sec
		csatrho = 1000.*matptr->rho*matptr->concSaturation/mpmptr->GetRelativeVolume();
		fluxMag.x = BCValue(bctime)/csatrho;
	}
	else
    {   // coupled surface flux and ftime is bath concentration
		// time variable (t) is replaced by c-cbath, where c is the particle potention and cbath and bath potential
		varTime = mpmptr->pPreviousConcentration-GetBCFirstTime();
		GetPosition(&varXValue,&varYValue,&varZValue,&varRotValue);
		double currentValue = fabs(function->Val());
		if(varTime>0.) currentValue=-currentValue;
		// csatrho = rho V csat/V0 = solvent mass per reference volume
		// units are 1000 kg mm^3/(m^2-g-sec) = mm/sec
 		csatrho = 1000.*matptr->rho*matptr->concSaturation/mpmptr->GetRelativeVolume();
		fluxMag.x = currentValue/csatrho;
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
	int i;
    for(i=1;i<=numnds;i++)
    {   // skip empty nodes
        if(nd[nds[i]]->NodeHasNonrigidParticles())
			diffusion->AddFluxCondition(nd[nds[i]],DotVectors(&fluxMag,&tscaled)*fn[i],false);
    }
	
    return (MatPtFluxBC *)GetNextObject();
}
