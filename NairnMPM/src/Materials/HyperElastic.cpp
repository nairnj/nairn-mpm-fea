/********************************************************************************
    HyperElastic.cpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/HyperElastic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Global_Quantities/ThermalRamp.hpp"

#pragma mark HyperElastic::Constructors and Destructors

// Constructors
HyperElastic::HyperElastic() {}

HyperElastic::HyperElastic(char *matName) : MaterialBase(matName)
{
	aI=0.;
}

#pragma mark HyperElastic::Methods

/*  Get new deformation gradient from current one using dF.F where dF = I + gradV * dt and F is current
        deformation gradient (i.e., two steps with two gradients F1 = F and F2 = dF, the total
        gradient is product of F1 and F2 in reverse order)
    dvij are elements of gradV * time step
    if storeInParticle is true, transfer new gradient to particle strain and rotation tensors
    Note: This assumes plane strain to find F[2][2]=1. If the 2D calculation is plane stress, the
        caller must replace F[2][2] with 1 + dw/dz
*/
void HyperElastic::GetDeformationGrad(double F[][3],MPMBase *mptr,double dvxx,double dvyy,
												double dvxy,double dvyx,bool storeInParticle)
{
	// deformation gradient (found from current strains and rotations)
	Tensor *ep=mptr->GetStrainTensor();
	TensorAntisym *wrot = mptr->GetRotationStrainTensor();
	
	// current deformation gradient in 2D
	double Fxx = 1. + ep->xx;;
	double Fxy = (ep->xy - wrot->xy)/2.;
	double Fyx = (ep->xy + wrot->xy)/2.;
	double Fyy = 1. + ep->yy;
	
	// get new 2D deformation gradient
	F[0][0] = (1. + dvxx)*Fxx + dvxy*Fyx;		// 1 + du/dx
	F[0][1] = (1. + dvxx)*Fxy + dvxy*Fyy;		// du/dy
	//F[0][2] = 0.;								// du/dz
	F[1][0] = dvyx*Fxx + (1. + dvyy)*Fyx;		// dv/dx
	F[1][1] = dvyx*Fxy + (1. + dvyy)*Fyy;		// 1 + dv/dy
	//F[1][2] = 0.;								// dv/dz
	//F[2][0] = 0.;								// dw/dx
	//F[2][1] = 0.;								// dw/dy
	F[2][2] = 1. ;								// 1 + dw/dz (assumes plane strain)
	
	// store in total strain and rotation tensors
	// (assumes plane strain so ep->zz = 0)
	if(storeInParticle)
	{	// strain increments
    	ep->xx = F[0][0] - 1.;
		ep->yy = F[1][1] - 1.;
		ep->xy = F[1][0] + F[0][1];
	
		// rotational strain increments
		wrot->xy = F[1][0] - F[0][1];
	}
}

// Get new deformation gradient from current one using dF.F where dF = I + gradV * dt and F is current
// dvij are elements of gradV * time step
// if storeInParticle is true, transfer new gradient to particle strain and rotation tensors
void HyperElastic::GetDeformationGrad(double F[][3],MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
						 double dvxz,double dvzx,double dvyz,double dvzy,bool storeInParticle)
{
	// get new doformation gradient from current one using dF.F where dF = I + gradV * dt and F is current
	// deformation gradient (stored int current strains and rotations)
	Tensor *ep=mptr->GetStrainTensor();
	TensorAntisym *wrot = mptr->GetRotationStrainTensor();
	
	// current deformation gradient in 3D
	double Fxx = 1. + ep->xx;;
	double Fxy = (ep->xy - wrot->xy)/2.;
	double Fxz = (ep->xz - wrot->xz)/2.;
	double Fyx = (ep->xy + wrot->xy)/2.;
	double Fyy = 1. + ep->yy;
	double Fyz = (ep->yz - wrot->yz)/2.;
	double Fzx = (ep->xz + wrot->xz)/2.;
	double Fzy = (ep->yz + wrot->yz)/2.;
	double Fzz = 1. + ep->zz;
	
	// get new deformation gradient
	F[0][0] = (1. + dvxx)*Fxx + dvxy*Fyx + dvxz*Fzx;		// 1 + du/dx
	F[0][1] = (1. + dvxx)*Fxy + dvxy*Fyy + dvxz*Fzy;		// du/dy
	F[0][2] = (1. + dvxx)*Fxz + dvxy*Fyz + dvxz*Fzz;		// du/dz
	F[1][0] = dvyx*Fxx + (1. + dvyy)*Fyx + dvyz*Fzx;		// dv/dx
	F[1][1] = dvyx*Fxy + (1. + dvyy)*Fyy + dvyz*Fzy;		// 1 + dv/dy
	F[1][2] = dvyx*Fxz + (1. + dvyy)*Fyz + dvyz*Fzz;		// dv/dz
	F[2][0] = dvzx*Fxx + dvzy*Fyx + (1. + dvzz)*Fzx;		// dw/dx
	F[2][1] = dvzx*Fxy + dvzy*Fyy + (1. + dvzz)*Fzy;		// dw/dy
	F[2][2] = dvzx*Fxz + dvzy*Fyz + (1. + dvzz)*Fzz;		// 1 + dw/dz
	
	// store in total strain and rotation tensors
	if(storeInParticle)
	{	ep->xx = F[0][0] - 1.;
		ep->yy = F[1][1] - 1.;
		ep->zz = F[2][2] - 1.;
		ep->xy = F[1][0] + F[0][1];
		ep->xz = F[2][0] + F[0][2];
		ep->yz = F[2][1] + F[1][2];
		
		// rotational strain increments
		wrot->xy = F[1][0] - F[0][1];
		wrot->xz = F[2][0] - F[0][2];
		wrot->yz = F[2][1] - F[1][2];
	}
}

// Find Left-Cauchy Green Tensor B = F.F^T for 2D calculations
// Note: This assumes plane strain to set B.zz=1. If the 2D calculation is plane stress, the
//	caller must replace B.zz with the plane stress result
Tensor HyperElastic::GetLeftCauchyTensor2D(double F[][3])
{
	Tensor B;
	ZeroTensor(&B);
	B.xx = F[0][0]*F[0][0] + F[0][1]*F[0][1];
	B.yy = F[1][0]*F[1][0] + F[1][1]*F[1][1];
	B.xy = F[0][0]*F[1][0] + F[0][1]*F[1][1];
	B.zz = 1.;
	return B;
}

// Find Left-Cauchy Green Tensor B = F.F^T for 3D calculations
Tensor HyperElastic::GetLeftCauchyTensor3D(double F[][3])
{
	// left Cauchy deformation tensor B = F F^T
	Tensor B;
	B.xx = F[0][0]*F[0][0] + F[0][1]*F[0][1] + F[0][2]*F[0][2];
	B.yy = F[1][0]*F[1][0] + F[1][1]*F[1][1] + F[1][2]*F[1][2];
	B.zz = F[2][0]*F[2][0] + F[2][1]*F[2][1] + F[2][2]*F[2][2];
	B.xy = F[0][0]*F[1][0] + F[0][1]*F[1][1] + F[0][2]*F[1][2];
	B.xz = F[0][0]*F[2][0] + F[0][1]*F[2][1] + F[0][2]*F[2][2];
	B.yz = F[1][0]*F[2][0] + F[1][1]*F[2][1] + F[1][2]*F[2][2];
	return B;
}

// Find isotropic stretch for thermal and moisture expansion
// total residual stretch (1 + alpha dT + beta csat dConcentration)
double HyperElastic::GetResidualStretch(MPMBase *mptr)
{
	// total residual stretch (1 + alpha dT + beta csat dConcentration)
	double resStretch = 1.0;
	double dTemp=mptr->pPreviousTemperature-thermal.reference;
	resStretch += CTE1*dTemp;
	if(DiffusionTask::active)
	{	double dConc=mptr->pPreviousConcentration-DiffusionTask::reference;
		resStretch += CME1*dConc;
	}
	return resStretch;
}

#ifndef CONSTANT_RHO
// Get current relative volume change = J = det F = lam1 lam2 lam3
// If using CONSTANT_RHO, then skip, which will pass to base class result of 1
double HyperElastic::GetCurrentRelativeVolume(MPMBase *mptr,bool threeD)
{
	// deformation gradient (found from current strains and rotations)
	Tensor *ep=mptr->GetStrainTensor();
	TensorAntisym *wrot = mptr->GetRotationStrainTensor();
	
	// current deformation gradient in 2D
	double Fxx = 1. + ep->xx;;
	double Fxy = (ep->xy - wrot->xy)/2.;
	double Fyx = (ep->xy + wrot->xy)/2.;
	double Fyy = 1. + ep->yy;
    double Fzz = 1. + ep->zz;
    
    if(threeD)
    {   // add 3D terms
        double Fxz = (ep->xz - wrot->xz)/2.;
        double Fyz = (ep->yz - wrot->yz)/2.;
        double Fzx = (ep->xz + wrot->xz)/2.;
        double Fzy = (ep->yz + wrot->yz)/2.;
        
        return Fxx*(Fyy*Fxx-Fzy*Fyz) - Fxy*(Fyx*Fzz-Fzx*Fyz) + Fxz*(Fyx*Fzy-Fzx*Fyy);
    }
    else
        return Fzz*(Fxx*Fyy - Fxy*Fyx);
}
#endif

// Convert to nominal stress in 2D by premultiply with J F^(-1)
// JFi is transpose of matrix of cofactors for F (since J = det F)
void HyperElastic::ConvertToNominalStress2D(MPMBase *mptr,double F[][3])
{
    Tensor *sp=mptr->GetStressTensor();
	double JFi[3][3];
	JFi[0][0] = F[1][1]*F[2][2];
	JFi[0][1] = -F[0][1]*F[2][2];
	//JFi[0][2] = 0.;
	JFi[1][0] = -F[1][0]*F[2][2];
	JFi[1][1] = F[0][0]*F[2][2];
	//JFi[1][2] = 0.;
	//JFi[2][0] = 0.;
	//JFi[2][1] = 0.;
	JFi[2][2] = F[0][0]*F[1][1] - F[1][0]*F[0][1];
	Tensor sp0=*sp;
	sp->xx = JFi[0][0]*sp0.xx + JFi[0][1]*sp0.xy;
	sp->xy = JFi[0][0]*sp0.xy + JFi[0][1]*sp0.yy;
	sp->yy = JFi[1][0]*sp0.xy + JFi[1][1]*sp0.yy;
	sp->zz = JFi[2][2]*sp0.zz;
}

// In future, may need to convert to Kirchoff or Nominal Stress in 2D or 3D 

// Convert to Kirchoff Stress in 2D (times J)
/*
	sp->xx *= J;
	sp->yy *= J;
	sp->zz *= J;
	sp->xy *= J;
*/

// Convert to nominal stress in 2D by premultiply with J F^(-1)
// JFi is transpose of matrix of cofactors for F
/*
	double JFi[3][3];
	JFi[0][0] = F[1][1]*F[2][2];
	JFi[0][1] = -F[0][1]*F[2][2];
	JFi[0][2] = 0.;
	JFi[1][0] = -F[1][0]*F[2][2];
	JFi[1][1] = F[0][0]*F[2][2];
	JFi[1][2] = 0.;
	JFi[2][0] = 0.;
	JFi[2][1] = 0.;
	JFi[2][2] = F[0][0]*F[1][1] - F[1][0]*F[0][1];
	Tensor sp0=*sp;
	sp->xx = JFi[0][0]*sp0.xx + JFi[0][1]*sp0.xy;
	sp->xy = JFi[0][0]*sp0.xy + JFi[0][1]*sp0.yy;
	sp->yy = JFi[1][0]*sp0.xy + JFi[1][1]*sp0.yy;
	sp->zz = JFi[2][2]*sp0.zz;
*/

// Convert to Kirchoff Stress 3D (times J)
/*
	sp->xx *= J;
	sp->yy *= J;
	sp->zz *= J;
	sp->xy *= J;
	sp->xz *= J;
	sp->yz *= J;
*/

// Convert to nominal stress in 3D by premultiply with J F^(-1)
// JFi is transpose of matrix of cofactors for F
/*
	double JFi[3][3];
	JFi[0][0] = F[1][1]*F[2][2] - F[2][1]*F[1][2];
	JFi[0][1] = -(F[0][1]*F[2][2] - F[2][1]*F[0][2]);
	JFi[0][2] = F[0][1]*F[1][2] - F[1][1]*F[0][2];
	JFi[1][0] = -(F[1][0]*F[2][2] - F[2][0]*F[1][2]);
	JFi[1][1] = F[0][0]*F[2][2] - F[2][0]*F[0][2];
	JFi[1][2] = -(F[0][0]*F[1][2] - F[1][0]*F[0][2]);
	JFi[2][0] = F[1][0]*F[2][1] - F[2][0]*F[1][1];
	JFi[2][1] = -(F[0][0]*F[2][1] - F[2][0]*F[0][1]);
	JFi[2][2] = F[0][0]*F[1][1] - F[1][0]*F[0][1];
	Tensor sp0=*sp;
	sp->xx = JFi[0][0]*sp0.xx + JFi[0][1]*sp0.xy + JFi[0][2]*sp0.xz;
	sp->xy = JFi[0][0]*sp0.xy + JFi[0][1]*sp0.yy + JFi[0][2]*sp0.yz;
	sp->xz = JFi[0][0]*sp0.xz + JFi[0][1]*sp0.yz + JFi[0][2]*sp0.zz;
	sp->yy = JFi[1][0]*sp0.xy + JFi[1][1]*sp0.yy + JFi[1][2]*sp0.yz;
	sp->yz = JFi[1][0]*sp0.xz + JFi[1][1]*sp0.yz + JFi[1][2]*sp0.zz;
	sp->zz = JFi[2][0]*sp0.xz + JFi[2][1]*sp0.yz + JFi[2][2]*sp0.zz;
*/



