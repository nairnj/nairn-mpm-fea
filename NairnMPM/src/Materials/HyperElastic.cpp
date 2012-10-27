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
    if detIncrement is true, return det(dF), otherwise return 0.0 (some materials mighe make use of this result)
    Note: This assumes plane strain to find F[2][2]=1. If the 2D calculation is plane stress, the
        caller must replace F[2][2] with 1 + dw/dz. Likewise, caller must multiply det(dF) by 1 + dw/dz
*/
double HyperElastic::GetDeformationGrad(double F[][3],MPMBase *mptr,double dvxx,double dvyy,
			double dvxy,double dvyx,bool storeInParticle,bool detIncrement)
{
	// current deformation gradient in 2D
    double pF[3][3];
    mptr->GetDeformationGradient(pF);
	
	// get new 2D deformation gradient
	F[0][0] = (1. + dvxx)*pF[0][0] + dvxy*pF[1][0];		// 1 + du/dx
	F[0][1] = (1. + dvxx)*pF[0][1] + dvxy*pF[1][1];		// du/dy
	//F[0][2] = 0.;								// du/dz
	F[1][0] = dvyx*pF[0][0] + (1. + dvyy)*pF[1][0];		// dv/dx
	F[1][1] = dvyx*pF[0][1] + (1. + dvyy)*pF[1][1];		// 1 + dv/dy
	//F[1][2] = 0.;								// dv/dz
	//F[2][0] = 0.;								// dw/dx
	//F[2][1] = 0.;								// dw/dy
	F[2][2] = 1. ;								// 1 + dw/dz (assumes plane strain)
	
	// store in total strain and rotation tensors
	// (assumes plane strain so ep->zz = 0)
	if(storeInParticle)
	{	// strain increments
        Tensor *ep=mptr->GetStrainTensor();
        TensorAntisym *wrot = mptr->GetRotationStrainTensor();
        
    	ep->xx = F[0][0] - 1.;
		ep->yy = F[1][1] - 1.;
		ep->xy = F[1][0] + F[0][1];         // du/dy + dv/dx
	
		// rotational strain increments
		wrot->xy = F[1][0] - F[0][1];		// dv/dx - du/dy
	}
    
    // calculate incremental determinant of (I+ gradV*dt) if desired
    if(!detIncrement) return 0.0;
    
    // assumes plain strain; multiply by 1+dwdz when known to get det(dF) in plane stress
    return (1. + dvxx)*(1. + dvyy)- dvyx*dvxy;
}

/*  Get new deformation gradient from current one using dF.F where dF = I + gradV * dt and F is current
        deformation gradient (i.e., two steps with two gradients F1 = F and F2 = dF, the total
        gradient is product of F1 and F2 in reverse order)
    dvij are elements of gradV * time step
    if Bptr!=NULL, find new B from (I + gradV*dt)pB(I + gradV*dt)^T where pB is current particle B
        This is not stored on the particle unless Bptr is pointer to particle B
    if storeInParticle is true, transfer new gradient to particle strain and rotation tensors
    if detIncrement is true, return det(dF), otherwise return 0.0 (some materials might make use of this result)
    Note: This assumes plane strain to find F[2][2]=1. If the 2D calculation is plane stress, the
        caller must replace F[2][2] with 1 + dw/dz. Likewise, caller must multiply det(dF) by 1 + dw/dz.
        Finally, caller must get B.zz as needed and store on particle
 */
double HyperElastic::IncrementDeformation(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx)
{
	// current deformation gradient in 2D
    double pF[3][3];
    mptr->GetDeformationGradient(pF);
	
	// get new 2D deformation gradient to increment particle strain
    // plaine stress will need to add ep->zz when known
    Tensor *ep=mptr->GetStrainTensor();
    TensorAntisym *wrot = mptr->GetRotationStrainTensor();
	ep->xx += dvxx*(1.+ep->xx) + dvxy*pF[1][0];             // 1 + du/dx - 1
	ep->yy += dvyx*pF[0][1] + dvyy*(1.+ep->yy);             // 1 + dv/dy - 1
	double F01 = (1. + dvxx)*pF[0][1] + dvxy*pF[1][1];		// du/dy
	double F10 = dvyx*pF[0][0] + (1. + dvyy)*pF[1][0];		// dv/dx
    ep->xy = F10 + F01;                                     // du/dy + dv/dx
    wrot->xy = F10 - F01;                                   // dv/dx - du/dy
	
    // increment Left Cauchy tensor B = F.F^T
    // plain stress will need to update B.zz when known
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    double b1 = (1.+dvxx)*pB->xx +      dvxy*pB->xy;
    double b2 = (1.+dvxx)*pB->xy +      dvxy*pB->yy;
    double b3 =      dvyx*pB->xx + (1.+dvyy)*pB->xy;
    double b4 =      dvyx*pB->xy + (1.+dvyy)*pB->yy;
    pB->xx = (1. + dvxx)*b1 +      dvxy*b2;
    pB->xy =        dvyx*b1 + (1.+dvyy)*b2;
    pB->yy =        dvyx*b3 + (1.+dvyy)*b4;
    
    // det(I + grad v * dt), which assumes plain strain
    // multiply by 1+dvzz when known to get det(dF) in plane stress
    return (1. + dvxx)*(1. + dvyy)- dvyx*dvxy;
}

// Get new deformation gradient from current one using dF.F where dF = I + gradV * dt and F is current
// dvij are elements of gradV * time step
// if storeInParticle is true, transfer new gradient to particle strain and rotation tensors
// if detIncrement is true, return det(dF), otherwise return 0.0 (some materials mighe make use of this result)
double HyperElastic::GetDeformationGrad(double F[][3],MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
						 double dvxz,double dvzx,double dvyz,double dvzy,bool storeInParticle,bool detIncrement)
{
	// current deformation gradient in 3D
    double pF[3][3];
    mptr->GetDeformationGradient(pF);
	
	// get new deformation gradient
	F[0][0] = (1. + dvxx)*pF[0][0] + dvxy*pF[1][0] + dvxz*pF[2][0];		// 1 + du/dx
	F[0][1] = (1. + dvxx)*pF[0][1] + dvxy*pF[1][1] + dvxz*pF[2][1];		// du/dy
	F[0][2] = (1. + dvxx)*pF[0][2] + dvxy*pF[1][2] + dvxz*pF[2][2];		// du/dz
	F[1][0] = dvyx*pF[0][0] + (1. + dvyy)*pF[1][0] + dvyz*pF[2][0];		// dv/dx
	F[1][1] = dvyx*pF[0][1] + (1. + dvyy)*pF[1][1] + dvyz*pF[2][1];		// 1 + dv/dy
	F[1][2] = dvyx*pF[0][2] + (1. + dvyy)*pF[1][2] + dvyz*pF[2][2];		// dv/dz
	F[2][0] = dvzx*pF[0][0] + dvzy*pF[1][0] + (1. + dvzz)*pF[2][0];		// dw/dx
	F[2][1] = dvzx*pF[0][1] + dvzy*pF[1][1] + (1. + dvzz)*pF[2][1];		// dw/dy
	F[2][2] = dvzx*pF[0][2] + dvzy*pF[1][2] + (1. + dvzz)*pF[2][2];		// 1 + dw/dz
	
	// store in total strain and rotation tensors
	if(storeInParticle)
    {   Tensor *ep=mptr->GetStrainTensor();
        TensorAntisym *wrot = mptr->GetRotationStrainTensor();
        
		ep->xx = F[0][0] - 1.;
		ep->yy = F[1][1] - 1.;
		ep->zz = F[2][2] - 1.;
		ep->xy = F[1][0] + F[0][1];
		ep->xz = F[2][0] + F[0][2];
		ep->yz = F[2][1] + F[1][2];
		
		// rotational strain increments
		wrot->xy = F[1][0] - F[0][1];			// dv/dx - du/dy
		wrot->xz = F[2][0] - F[0][2];			// dw/dx - du/dz
		wrot->yz = F[2][1] - F[1][2];			// dw/dy - dv/dz
	}
    
    // calculate incremental determinant of (I+ gradV*dt) if desired
    if(!detIncrement) return 0.0;
    
    return (1. + dvxx)*((1. + dvyy)*(1. + dvzz)-dvzy*dvyz)
                    - dvyx*(dvxy*(1. + dvzz)-dvzy*dvxz)
                    + dvzx*(dvxy*dvyz-(1. + dvyy)*dvxz);
}

// Find Left-Cauchy Green Tensor B = F.F^T for 2D calculations from a provided F[][]
// Note: This assumes plane strain and sets B.zz=1. If the 2D calculation is plane stress, the
//	caller must replace B.zz with the plane stress result
Tensor HyperElastic::GetLeftCauchyTensor2D(double F[][3])
{
	// left Cauchy deformation tensor B = F F^T
	Tensor B;
	ZeroTensor(&B);
	B.xx = F[0][0]*F[0][0] + F[0][1]*F[0][1];
	B.yy = F[1][0]*F[1][0] + F[1][1]*F[1][1];
	B.xy = F[0][0]*F[1][0] + F[0][1]*F[1][1];
	B.zz = 1.;
	return B;
}

// Find Left-Cauchy Green Tensor B = F.F^T for 3D calculations from a provided F[][]
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

// Get current relative volume change = J = det F = lam1 lam2 lam3
// Need to have this call in material classes to allow small and large deformation material laws
//  to handle it differently. It is used on archiving to convert Kirchoff Stress/rho0 to Cauchy stress
double HyperElastic::GetCurrentRelativeVolume(MPMBase *mptr)
{   return mptr->GetRelativeVolume();
}

