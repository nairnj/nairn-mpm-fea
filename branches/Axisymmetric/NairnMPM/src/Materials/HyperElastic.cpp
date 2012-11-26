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

/*  Given components of incremental deformation dF = (I + gradV *dt), increment particle strain,
        rotation, and LeftCauchy Green strain (latter is assumed to be stored in the particle's
        plastic strain tensor (which is accessed also with GetElasticLeftCauchyTensor().
    New new F is dF.F, which is used to find new strain
    New B = dF.(Old B).dF^T
    dvij are elements of gradV * time step (except for dvzz)
    Returns |dF|
    Note: This assumes plane strain or F[2][2]=1, B.zz and e.zz unchanged. If the 2D calculation
        is plane stress, the caller will need to update B.zz and e.zz on particle when known,
        where New B.zz = (1+dvzz)^2 * (Old B.zz) = (1+e.zz)^2
        Likewise, caller must multiply det(dF) by 1 + dvzz = sqrt(New B.zz/Old B.zz).
 */
double HyperElastic::IncrementDeformation(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,double dvzz)
{
	// get new 2D deformation gradient to increment particle strain
    // plaine stress will need to add ep->zz when known
    Tensor *ep=mptr->GetStrainTensor();
    TensorAntisym *wrot = mptr->GetRotationStrainTensor();
    
    // exx  = (1+dvxx)*pF[0][0] + dvxy*pF[1][0] - 1
    //      = (1+dvxx)*(1+exx) + dvxy*0.5*(exy + wxy) - 1
    //      = exx + dvxx*(1+exx) + dvxy*pF10
    double pF10 = 0.5*(ep->xy + wrot->xy);                  // previous particle gradient
	double F10 = dvyx*(1. + ep->xx) + (1. + dvyy)*pF10;		// dv/dx in new deformation gradient
	ep->xx += dvxx*(1.+ep->xx) + dvxy*pF10;
    
    // eyy = dvyx*pF[0][1] + (1+dvyy)*(1+eyy) - 1
    //      = eyy + dvyy*(1+eyy) + dvyx*pF01
    double pF01 = 0.5*(ep->xy - wrot->xy);                  // previous particle gradient
	double F01 = (1. + dvxx)*pF01 + dvxy*(1. + ep->yy);		// du/dy in new deformation gradient
	ep->yy += dvyx*pF01 + dvyy*(1.+ep->yy);
    
    // shear and rotation for new deformation gradient
    ep->xy = F10 + F01;                                     // du/dy + dv/dx
    wrot->xy = F10 - F01;                                   // dv/dx - du/dy
	
    // ezz  = (1+dvzz)*pF[2][2] - 1							// axisymmetric only, otherwise dvzz=0
    //      = (1+dvzz)*(1+ezz) - 1
    //      = ezz + dvzz*(1+exx)
	ep->zz += dvzz*(1.+ep->zz);

    // increment Left Cauchy tensor B = F.F^T = dF.old B.dF^T
    // plain stress will need to update B.zz when known
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    double b1 = (1.+dvxx)*pB->xx +      dvxy*pB->xy;
    double b2 = (1.+dvxx)*pB->xy +      dvxy*pB->yy;
    double b3 =      dvyx*pB->xx + (1.+dvyy)*pB->xy;
    double b4 =      dvyx*pB->xy + (1.+dvyy)*pB->yy;
    double b5 = (1.+dvzz)*pB->zz;
    pB->xx = b1*(1. + dvxx) + b2*dvxy;
    pB->xy = b1*dvyx        + b2*(1.+dvyy);
    pB->yy = b3*dvyx        + b4*(1.+dvyy);
	pB->zz = b5*(1. + dvzz);
    
    // det(I + grad v * dt), which assumes plain strain or axisymmetric
    // multiply by 1+dvzz when known to get det(dF) in plane stress
    return (1. + dvzz)*((1. + dvxx)*(1. + dvyy)- dvyx*dvxy);
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

/*  Given components of incremental deformation dF = (I + gradV *dt), increment particle strain,
        rotation, and LeftCauchy Green strain (latter is assumed to be stored in the particle's
        plastic strain tensor (which is accessed also with GetElasticLeftCauchyTensor().
    New new F is dF.F, which is used to find new strain
    New B = dF.(Old B).dF^T
    dvij are elements of gradV * time step
    Returns |dF|
*/
double HyperElastic::IncrementDeformation(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
                                        double dvxz,double dvzx,double dvyz,double dvzy)
{
	// current deformation gradient in 3D
    double pF[3][3];
    mptr->GetDeformationGradient(pF);
	
	// get new deformation gradient for off axis components
	double F01 = (1. + dvxx)*pF[0][1] + dvxy*pF[1][1] + dvxz*pF[2][1];		// du/dy
	double F02 = (1. + dvxx)*pF[0][2] + dvxy*pF[1][2] + dvxz*pF[2][2];		// du/dz
	double F10 = dvyx*pF[0][0] + (1. + dvyy)*pF[1][0] + dvyz*pF[2][0];		// dv/dx
	double F12 = dvyx*pF[0][2] + (1. + dvyy)*pF[1][2] + dvyz*pF[2][2];		// dv/dz
	double F20 = dvzx*pF[0][0] + dvzy*pF[1][0] + (1. + dvzz)*pF[2][0];		// dw/dx
	double F21 = dvzx*pF[0][1] + dvzy*pF[1][1] + (1. + dvzz)*pF[2][1];		// dw/dy
	
	// store in total strain and rotation tensors
    Tensor *ep=mptr->GetStrainTensor();
    TensorAntisym *wrot = mptr->GetRotationStrainTensor();
        
    ep->xx = (1. + dvxx)*pF[0][0] + dvxy*pF[1][0] + dvxz*pF[2][0] - 1.;     // 1 + du/dx - 1
    ep->yy = dvyx*pF[0][1] + (1. + dvyy)*pF[1][1] + dvyz*pF[2][1] - 1.;		// 1 + dv/dy - 1
    ep->zz = dvzx*pF[0][2] + dvzy*pF[1][2] + (1. + dvzz)*pF[2][2] - 1.;		// 1 + dw/dz - 1
    ep->xy = F10 + F01;
    ep->xz = F20 + F02;
    ep->yz = F21 + F12;
    
    // rotational strain increments
    wrot->xy = F10 - F01;			// dv/dx - du/dy
    wrot->xz = F20 - F02;			// dw/dx - du/dz
    wrot->yz = F21 - F12;			// dw/dy - dv/dz
    
    // increment Left Cauchy tensor B = F.F^T = dF.old B.dF^T
    // plain stress will need to update B.zz when known
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    double b1 = (1.+dvxx)*pB->xx +      dvxy*pB->xy +      dvxz*pB->xz;
    double b2 = (1.+dvxx)*pB->xy +      dvxy*pB->yy +      dvxz*pB->yz;
    double b3 = (1.+dvxx)*pB->xz +      dvxy*pB->yz +      dvxz*pB->zz;
    double b4 =      dvyx*pB->xx + (1.+dvyy)*pB->xy +      dvyz*pB->xz;
    double b5 =      dvyx*pB->xy + (1.+dvyy)*pB->yy +      dvyz*pB->yz;
    double b6 =      dvyx*pB->xz + (1.+dvyy)*pB->yz +      dvyz*pB->zz;
    double b7 =      dvzx*pB->xx +      dvzy*pB->xy + (1.+dvzz)*pB->xz;
    double b8 =      dvzx*pB->xy +      dvzy*pB->yy + (1.+dvzz)*pB->yz;
    double b9 =      dvzx*pB->xz +      dvzy*pB->yz + (1.+dvzz)*pB->zz;
    pB->xx = (1+dvxx)*b1 +      dvxy*b2 +      dvxz*b3;
    pB->xy =     dvyx*b1 + (1.+dvyy)*b2 +      dvyz*b3;
    pB->xz =     dvzx*b1 +      dvzy*b2 + (1.+dvzz)*b3;
    pB->yy =     dvyx*b4 + (1.+dvyy)*b5 +      dvyz*b6;
    pB->yz =     dvzx*b4 +      dvzy*b5 + (1.+dvzz)*b6;
    pB->zz =     dvzx*b7 +      dvzy*b8 + (1.+dvzz)*b9;
    
    // return |dF|
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

