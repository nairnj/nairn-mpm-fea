/********************************************************************************
    MatPtHeatFluxBC.cpp
    nairn-mpm-fea

    Created by John Nairn on May 28, 2013.
    Copyright (c) 2013 John A. Nairn, All rights reserved.
********************************************************************************/

#ifdef _MSC_VER
#include "stdafx.h"
#endif

#include "Boundary_Conditions/MatPtHeatFluxBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/CommonException.hpp"
#include "Read_XML/Expression.hpp"
#include "Elements/ElementBase.hpp"
#include "System/MPMPrefix.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"

// Local globals for BCs
// c1 starts with c12i[1] and c2 starts with c12i[o]
static double c12i[5]={1.,1.,-1.,-1.,1.};

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
	size_t nlsize=200;
    
    snprintf(nline,nlsize,"%7d %2d   %2d  %2d %15.7e %15.7e",ptNum,direction,face,style,
			UnitsController::Scaling(1.e-3)*GetBCValueOut(),GetBCFirstTimeOut());
    os << nline;
	PrintFunction(os);
	
	// not allowed on rigid particles
	MaterialBase *matref = theMaterials[mpm[ptNum-1]->MatID()];
	if(matref->IsRigid())
	{	throw CommonException("Cannot set heat flux on rigid particles",
							  "MatPtHeatFluxBC::PrintBC");
	}

	// for function input scale for Legacy units
	if(style==FUNCTION_VALUE)
		scale = UnitsController::Scaling(1.e3);
	
    return (BoundaryCondition *)GetNextObject();
}

// increment external heat flux on a particle
// input is analysis time in seconds
// (only called when conduction is active)
MatPtLoadBC *MatPtHeatFluxBC::AddMPFluxBC(double bctime)
{	
    // condition value
	// flux BC in nW/mm^2
	MPMBase *mpmptr = mpm[ptNum-1];
    MaterialBase *matptr = theMaterials[mpmptr->MatID()];
	
	// Flux is a scalar and we need int_(face) F Ni(x) dA
	// Since F is constant, only need integral which is done by CPDI methods
	//		which has be generalized to work for GIMP too
	// We use X_DIRECTION for bcDIR for efficiency. For Silent BC, change to
	//      Normal direction to all calculation of n
	Vector fluxMag;
	ZeroVector(&fluxMag);
    int bcDir=X_DIRECTION;
	
	if(style==SILENT)
	{	TransportProperties t;
		matptr->GetTransportProps(mpmptr,fmobj->np,&t);
		Tensor *k = &(t.kCondTensor);
		
		// k is k/rho (E-L^2/(K-T-M)), grat T in K/L, k grad T in E-L/(T-M)
		if(fmobj->IsThreeD())
		{	fluxMag.x = k->xx*mpmptr->pTemp[gGRADx] + k->xy*mpmptr->pTemp[gGRADy] + k->xz*mpmptr->pTemp[gGRADz];
			fluxMag.y = k->xy*mpmptr->pTemp[gGRADx] + k->yy*mpmptr->pTemp[gGRADy] + k->yz*mpmptr->pTemp[gGRADz];
			fluxMag.z = k->xz*mpmptr->pTemp[gGRADx] + k->yz*mpmptr->pTemp[gGRADy] + k->zz*mpmptr->pTemp[gGRADz];
		}
		else
		{	fluxMag.x = k->xx*mpmptr->pTemp[gGRADx] + k->xy*mpmptr->pTemp[gGRADy];
			fluxMag.y = k->xy*mpmptr->pTemp[gGRADx] + k->yy*mpmptr->pTemp[gGRADy];
		}
		
		// remove 1/rho0 scaling on k to get E/(T-L^2)
		ScaleVector(&fluxMag,matptr->GetRho(mpmptr));
		
		// need to get normal vector from cpdi functions below
        bcDir = N_DIRECTION;
	}
	else if(direction==EXTERNAL_FLUX)
	{	// user-supplied constant (E/(T-L^2))
		fluxMag.x = BCValue(bctime);
	}
	else
    {	// coupled surface flux
		if(bctime>=GetBCFirstTime())
		{	// time variable (t) is replaced by particle temperature, result should be E/(T-L^2)
            // (see Expression vmap)
			double vars[7];
			vars[0] = 6.5;
			vars[1] = mpmptr->pPreviousTemperature;		//t
			GetPositionVars(vars);
			
			// Legacy scaling of W/m^2 to nW/mm^2
			fluxMag.x = scale*function->EvaluateFunction(vars);
		}
	}
	
	// get corners and radii from deformed material point (2 in 2D and 4 in 3D)
	// Not that bcDir will be X_DIRECTION to scale fluxMag, but if silent, it will be N_DIRECTION
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
		
		// add flux to each node. fluxMag has units Energy/(sec-L^2)
		// wtNorm has units L^2 for final flux is Energy/sec
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
			conduction->AddFluxCondition(nd[nds[i]],flux,false);
		}
	}
	
    return (MatPtHeatFluxBC *)GetNextObject();
}

#pragma mark MatPtHeatFluxBC:Accessors

// set value (and scale legacy W/m^2 to nW/mm^2)
void MatPtHeatFluxBC::SetBCValue(double bcvalue)
{	BoundaryCondition::SetBCValue(UnitsController::Scaling(1.e3)*bcvalue);
}



