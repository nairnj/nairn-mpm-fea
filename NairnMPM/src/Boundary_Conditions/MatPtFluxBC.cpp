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
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/CommonException.hpp"
#include "Read_XML/Expression.hpp"
#include "Elements/ElementBase.hpp"

// Local globals for BCs
// c1 starts with c12i[1] and c2 starts with c12i[o]
static double c12i[5]={1.,1.,-1.,-1.,1.};

// global
MatPtFluxBC *firstFluxPt=NULL;
MatPtFluxBC *firstDiffFluxBC[NUM_DUFFUSION_OPTIONS];

#pragma mark MatPtFluxBC: Constructors and Destructors

// Constructors
MatPtFluxBC::MatPtFluxBC(int num,int dof,int sty,int edge,int phase) : MatPtLoadBC(num,dof,sty)
{
	face = edge;
    phaseStyle = phase;
    if(phaseStyle<FRACTURE_PHASE_FIELD)
    {   phaseStyle = MOISTURE_DIFFUSION;
#ifdef POROELASTICITY
        if(DiffusionTask::HasPoroelasticity())
            phaseStyle = POROELASTICITY_DIFFUSION;
#endif
    }
}

#pragma mark MatPtFluxBC: Methods

// print to output
BoundaryCondition *MatPtFluxBC::PrintBC(ostream &os)
{
    char nline[200];
	size_t nlsize=200;
	
    // =1 for all except moisture in Legacy is 1e+3
    double rescale = 1./DiffusionTask::RescaleFlux(phaseStyle);
    snprintf(nline,nlsize,"%7d %2d   %2d  %2d %15.7e %15.7e %3d",ptNum,direction,face,style,
			rescale*GetBCValueOut(),GetBCFirstTimeOut(),phaseStyle);
    os << nline;
	PrintFunction(os);
	
	// not allowed on rigid particles
	MaterialBase *matref = theMaterials[mpm[ptNum-1]->MatID()];
	if(matref->IsRigid())
	{	throw CommonException("Cannot set solvent/transport/pore pressure flux on rigid particles",
							  "MatPtFluxBC::PrintBC");
	}
	
	// for function input scale for Legacy units if needed
	if(style==FUNCTION_VALUE)
	{	scale = 1./rescale;
	}
	
    return (BoundaryCondition *)GetNextObject();
}

// Add transport flux on a particle surface
// input is analysis time in seconds
// (called for diffusion, poroelasticity, ot other diffusion methods)
MatPtLoadBC *MatPtFluxBC::AddMPFluxBC(double bctime)
{
    // Moisture: condition value is g/(mm^2-sec), Divide by rho*csat to get potential flux in mm/sec
	// Poroelasticity: condition value is Pa/sec = µN/(mm^2-sec) - does not account for area (could if needed)
    // Other transport:
	//    find this flux and then add (times area) to get mm^3-potential/sec
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
    
    // for extra diffusion
    DiffusionTask *diffusionTask = DiffusionTask::FindDiffusionTask(phaseStyle);
    if(diffusionTask==NULL) return (MatPtFluxBC *)GetNextObject();
    int diffusionNumber = diffusionTask->GetNumber();
	
	if(style==SILENT)
	{   // for moisture and poroelasticity only
        TransportProperties t;
		matptr->GetTransportProps(mpmptr,fmobj->np,&t);
		Tensor *D = &(t.diffusionTensor);
		
		// D in L^2/T (moisture) or L^2/(P-T) (poroelasticity), grad C in 1/L (moisture), P/L (poroelasticity)
		// Flux is L/T (both)
		if(fmobj->IsThreeD())
		{	fluxMag.x = D->xx*mpmptr->pDiff[0]->grad.x + D->xy*mpmptr->pDiff[0]->grad.y + D->xz*mpmptr->pDiff[0]->grad.z;
			fluxMag.y = D->xy*mpmptr->pDiff[0]->grad.x + D->yy*mpmptr->pDiff[0]->grad.y + D->yz*mpmptr->pDiff[0]->grad.z;
			fluxMag.z = D->xz*mpmptr->pDiff[0]->grad.x + D->yz*mpmptr->pDiff[0]->grad.y + D->zz*mpmptr->pDiff[0]->grad.z;
		}
		else
		{	fluxMag.x = D->xx*mpmptr->pDiff[0]->grad.x + D->xy*mpmptr->pDiff[0]->grad.y;
			fluxMag.y = D->xy*mpmptr->pDiff[0]->grad.x + D->yy*mpmptr->pDiff[0]->grad.y;
		}
		
        bcDir = N_DIRECTION;
	}
	else if(direction==EXTERNAL_FLUX)
    {   if(phaseStyle>POROELASTICITY_DIFFUSION)
        {   fluxMag.x = BCValue(bctime)*mpmptr->GetVolume(DEFORMED_VOLUME);
        }
		else if(fmobj->HasDiffusion())
		{	// provided value in M/(L^2-T), scale by current density to get L/T
			fluxMag.x = BCValue(bctime)*mpmptr->GetVolume(DEFORMED_VOLUME)/mpmptr->mp;
		}
#ifdef POROELASTICITY
		else
		{	// provided value in (dV/V)/(L^2-T), scale by V to get L/T
			fluxMag.x = BCValue(bctime)*mpmptr->GetVolume(DEFORMED_VOLUME);
		}
#endif
	}
	else
    {   double cmcres = 0.;
        if(phaseStyle>POROELASTICITY_DIFFUSION)
        {   // function of the current transport quantity
            cmcres = mpmptr->pDiff[diffusionNumber]->prevConc;
        }
		else if(fmobj->HasDiffusion())
		{	// moisture f(c-cres) (units potential) and function should give flux in M/(L^2-T)
			// time variable (t) is replaced by c-cres, where c is the particle value and cres is reservoir
			cmcres = mpmptr->pDiff[0]->prevConc-GetBCFirstTime();
		}
#ifdef POROELASTICITY
		else if(fmobj->HasPoroelasticity())
		{	// poroelasticity f(p-pres) (units P) and function should give flux in (dV/V)/(L^2-T)
			// Poroelasticity uses particle value to support changed flux when void space
			cmcres = mpmptr->pDiff[0]->conc-GetBCFirstTime();
		}
#endif
        // (see Expression vmap)
		double vars[7];
		vars[0] = 6.5;
		vars[1] = cmcres;		//t
		GetPositionVars(vars);
		
		// scaling only used when in Legacy units
		double currentValue = fabs(scale*function->EvaluateFunction(vars));
		
		// change direction to match sign of the difference
		if(cmcres>0.) currentValue=-currentValue;
		
        if(phaseStyle>POROELASTICITY_DIFFUSION)
        {   // value should be flux for that quantity
            fluxMag.x = currentValue;
        }
		else if(fmobj->HasDiffusion())
		{	// provided value is M/(L^2-T), scale by current density to get L/T
			fluxMag.x = currentValue*mpmptr->GetVolume(DEFORMED_VOLUME)/mpmptr->mp;
		}
#ifdef POROELASTICITY
		else
		{	// provided value is (dV/V)/(L^2-T), scale by V to get L/T
			fluxMag.x = currentValue*mpmptr->GetVolume(DEFORMED_VOLUME);
		}
#endif
	}
	
	// get corners and radii from deformed material point (2 in 2D and 4 in 3D)
	// Note that bcDir will be X_DIRECTION to scale fluxMag, but if silent, it will be N_DIRECTION
	int cElem[4],numDnds;
	Vector corners[4],radii[3];
	double redge;
	Vector wtNorm = mpmptr->GetSurfaceInfo(face,bcDir,cElem,corners,radii,&numDnds,&redge);
	
#ifdef CONST_ARRAYS
	int nds[MAX_SHAPE_NODES];
	double fn[MAX_SHAPE_NODES];
#else
	int nds[maxShapeNodes];
	double fn[maxShapeNodes];
#endif
		
	for(int c=0;c<numDnds;c++)
	{	// get straight grid shape functions
		theElements[cElem[c]]->GetShapeFunctionsForTractions(fn,nds,&corners[c]);
		int numnds = nds[0];
			
		// add flux to each node
		double flux;
		for(int i=1;i<=numnds;i++)
		{	if(fmobj->IsAxisymmetric())
			{	// wtNorm has direction and |r1|. Also scale by axisymmetric term
				flux = DotVectors(&fluxMag,&wtNorm)*(redge-c12i[c+1]*radii[0].x/3.)*fn[i];
			}
			else
			{	// wtNorm has direction and Area/2 (2D) or Area/4 (3D) to average the nodes
				flux = DotVectors(&fluxMag,&wtNorm)*fn[i];
			}
			diffusionTask->AddFluxCondition(nd[nds[i]],flux,false);
		}
	}
	
    return (MatPtFluxBC *)GetNextObject();
}

#pragma mark MatPtFluxBC:Accessors

// set value (and scale legacy kg/(m^2-sec) to g/(mm^2-sec) for concentration
// or (dV/V)/sec to (dV/V)/sec for pore pressure or 1 for all others
void MatPtFluxBC::SetBCValue(double bcvalue)
{   double rescale = DiffusionTask::RescaleFlux(phaseStyle);
	BoundaryCondition::SetBCValue(rescale*bcvalue);
}

// check coupled flux which uses first time as unscaled concentration potential
//		or reservoir pressure (Legacy MPa or stress units)
void MatPtFluxBC::SetBCFirstTime(double bcftime)
{	if(direction==COUPLED_FLUX)
    {   // change parameter when update flux to multiple types to style of this BC
        double rescale = DiffusionTask::RescalePotential(MOISTURE_DIFFUSION);
		BoundaryCondition::SetBCFirstTimeCU(rescale*bcftime);
	}
	else
		BoundaryCondition::SetBCFirstTime(bcftime);
}

// check coupled flux which uses first time as concentration potential or reservoir pressure
double MatPtFluxBC::GetBCFirstTimeOut(void)
{	if(direction==COUPLED_FLUX)
		return BoundaryCondition::GetBCFirstTime();
	else
		return BoundaryCondition::GetBCFirstTimeOut();
}


