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
	Kbulk = -1.;                                        // required (check >0 before starting)
    UofJOption = HALF_J_SQUARED_MINUS_1_MINUS_LN_J;     // default U(J) functino
	aI=0.;
}

// Read material properties
char *HyperElastic::InputMat(char *xName,int &input)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"K")==0)
        return((char *)&Kbulk);
    
    else if(strcmp(xName,"alpha")==0)
        return((char *)&aI);
    
    return(MaterialBase::InputMat(xName,input));
}

#pragma mark HyperElastic::Initialize

// Set intial particle Left Cauchy strain tensor to identity
void HyperElastic::SetInitialParticleState(MPMBase *mptr,int np)
{
    // get previous particle B
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    
    ZeroTensor(pB);
    pB->xx = pB->yy = pB->zz = 1.;
}

// Constant properties used in constitutive law
void HyperElastic::InitialLoadMechProps(int makeSpecific,int np)
{
	// Kbulk in Specific units using initial rho
	// for MPM (units N/m^2 cm^3/g)
	Ksp=Kbulk*1.0e+06/rho;
	
	// expansion coefficients
	CTE1 = 1.e-6*aI;
	CME1 = betaI*concSaturation;
	
    // call superclass
    MaterialBase::InitialLoadMechProps(makeSpecific,np);
}
#pragma mark HyperElastic::Methods

//#define OLD_METHOD
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
double HyperElastic::IncrementDeformation(MPMBase *mptr,double dvxx,double dvyy,double dvxy,
                                            double dvyx,double dvzz,Tensor *Btrial)
{
#ifndef OLD_METHOD
	Matrix3 du(dvxx,dvxy,dvyx,dvyy,dvzz);
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
	double dFxx = dF(0,0);
	double dFxy = dF(0,1);
	double dFyx = dF(1,0);
	double dFyy = dF(1,1);
	double dFzz = dF(2,2);
#else
	// get elements of dF = exp(dt*grad v)
	double c0 = dvxy*dvyx - dvxx*dvyy;					// -det(dt*grad V)
	double c1 = dvxx + dvyy;							// Tr(dt * grad V)
	double alpha0 = 1. + 0.5*c0, alpha1 = 1. + 0.5*c1;
	double alphaz = dvzz*dvzz, dFzz = 1. + dvzz + 0.5*alphaz;
	
	// evalate k terms in the expansion
	int k,kmax=2;
	double factorial = 2.,temp,beta0 = c0,beta1 = c1;
	for(k=3;k<=kmax;k++)
	{	// update beta using beta(k,0) = c0 beta(k-1,1)
		//			     and beta(k,1) = c1 beta(k-1,1) + beta(k-1,0)
		temp = c0*beta1;
		beta1 = c1*beta1 + beta0;
		beta0 = temp;
		factorial *= (double)k;
		alpha0 += beta0/factorial;
		alpha1 += beta1/factorial;
		alphaz *= dvzz;
		dFzz += alphaz/factorial;
	}
	
	// to try first order
	//alpha0 = 1.;
	//alpha1 = 1.;
	//dFzz = 1. + dvzz;
	
	// Now dF = alpha0 I + alpha1 dt grad V
	double dFxx = alpha0 + alpha1*dvxx;
	double dFxy = alpha1*dvxy;
	double dFyx = alpha1*dvyx;
	double dFyy = alpha0 + alpha1*dvyy;
#endif
	
	// get new 2D deformation gradient to increment particle strain
    // plaine stress will need to add ep->zz when known
    Tensor *ep=mptr->GetStrainTensor();
    TensorAntisym *wrot = mptr->GetRotationStrainTensor();
    
    // exx  = dFxx*pFxx + dFxy*pFyx - 1
    double pFyx = 0.5*(ep->xy + wrot->xy);                  // previous particle gradient
	double Fyx = dFyx*(1. + ep->xx) + dFyy*pFyx;		    // dv/dx in new deformation gradient
	ep->xx = dFxx*(1. + ep->xx) + dFxy*pFyx - 1.;
    
    // eyy = dFyx*pFxy + dFyy*(1+eyy) - 1
    double pFxy = 0.5*(ep->xy - wrot->xy);                  // previous particle gradient
	double Fxy = dFxx*pFxy + dFxy*(1. + ep->yy);			// du/dy in new deformation gradient
	ep->yy = dFyx*pFxy + dFyy*(1. + ep->yy) - 1.;
	
    // shear and rotation for new deformation gradient
    ep->xy = Fyx + Fxy;                                     // du/dy + dv/dx
    wrot->xy = Fyx - Fxy;                                   // dv/dx - du/dy
	
    // ezz  = dFzz*pFzz - 1                                 // axisymmetric only, otherwise dvzz=0
	ep->zz = dFzz*(1. + ep->zz) - 1.;

    // increment Left Cauchy tensor B = F.F^T = dF.old B.dF^T
    // plain stress will need to update B.zz when known
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    
    // elements of dF.B
    double dFBxx = dFxx*pB->xx + dFxy*pB->xy;
    double dFBxy = dFxx*pB->xy + dFxy*pB->yy;
    double dFByx = dFyx*pB->xx + dFyy*pB->xy;
    double dFByy = dFyx*pB->xy + dFyy*pB->yy;
    double dFBzz = dFzz*pB->zz;
    
    // return trial B (if provided) or store new B on the particle
    if(Btrial != NULL) pB = Btrial;
    
    // find (dF.B).dF^T
    pB->xx = dFBxx*dFxx + dFBxy*dFxy;
    pB->xy = dFBxx*dFyx + dFBxy*dFyy;
    pB->yy = dFByx*dFyx + dFByy*dFyy;
    pB->zz = dFBzz*dFzz;
    
    // det(I + grad v * dt), which assumes plain strain or axisymmetric
    // multiply by dFzz when known to get det(dF) in plane stress
    return dFzz*(dFxx*dFyy - dFyx*dFxy);
}

/*  Given components of incremental deformation dF = exp(dt*grad v), increment particle strain,
        rotation, and LeftCauchy Green strain (latter is assumed to be stored in the particle's
        plastic strain tensor (which is accessed also with GetElasticLeftCauchyTensor().
    New new F is dF.F, which is used to find new strain
    New B = dF.(Old B).dF^T
    dvij are elements of gradV * time step
    Returns |dF|
*/
double HyperElastic::IncrementDeformation(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
                                        double dvxz,double dvzx,double dvyz,double dvzy,Tensor *Btrial)
{
#ifndef OLD_METHOD
	Matrix3 du(dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz);
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
	
	// current deformation gradient in 3D
	const Matrix3 pF = mptr->GetDeformationGradientMatrix();
	
	// new deformation matrix
	const Matrix3 F = dF*pF;
	
	// store in total strain and rotation tensors
    Tensor *ep=mptr->GetStrainTensor();
    TensorAntisym *wrot = mptr->GetRotationStrainTensor();
	
	ep->xx = F(0,0) - 1.;		// 1 + du/dx - 1
    ep->yy = F(1,1) - 1.;		// 1 + dv/dy - 1
    ep->zz = F(2,2) - 1.;		// 1 + dw/dz - 1
    ep->xy = F(1,0) + F(0,1);
    ep->xz = F(2,0) + F(0,2);
    ep->yz = F(2,1) + F(1,2);
    
    // rotational strain increments
    wrot->xy = F(1,0) - F(0,1);			// dv/dx - du/dy
    wrot->xz = F(2,0) - F(0,2);			// dw/dx - du/dz
    wrot->yz = F(2,1) - F(1,2);			// dw/dy - dv/dz

    // increment Left Cauchy tensor B = F.F^T = dF.old B.dF^T
    // plain stress will need to update B.zz when known
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    
    // elements of dF.B
    double dFBxx = dF(0,0)*pB->xx + dF(0,1)*pB->xy + dF(0,2)*pB->xz;
    double dFBxy = dF(0,0)*pB->xy + dF(0,1)*pB->yy + dF(0,2)*pB->yz;
    double dFBxz = dF(0,0)*pB->xz + dF(0,1)*pB->yz + dF(0,2)*pB->zz;
    double dFByx = dF(1,0)*pB->xx + dF(1,1)*pB->xy + dF(1,2)*pB->xz;
    double dFByy = dF(1,0)*pB->xy + dF(1,1)*pB->yy + dF(1,2)*pB->yz;
    double dFByz = dF(1,0)*pB->xz + dF(1,1)*pB->yz + dF(1,2)*pB->zz;
    double dFBzx = dF(2,0)*pB->xx + dF(2,1)*pB->xy + dF(2,2)*pB->xz;
    double dFBzy = dF(2,0)*pB->xy + dF(2,1)*pB->yy + dF(2,2)*pB->yz;
    double dFBzz = dF(2,0)*pB->xz + dF(2,1)*pB->yz + dF(2,2)*pB->zz;
    
    // return trial B (if provided) or store new B on the particle
    if(Btrial != NULL) pB = Btrial;
    pB->xx = dFBxx*dF(0,0) + dFBxy*dF(0,1) + dFBxz*dF(0,2);
    pB->xy = dFBxx*dF(1,0) + dFBxy*dF(1,1) + dFBxz*dF(1,2);
    pB->xz = dFBxx*dF(2,0) + dFBxy*dF(2,1) + dFBxz*dF(2,2);
    pB->yy = dFByx*dF(1,0) + dFByy*dF(1,1) + dFByz*dF(1,2);
    pB->yz = dFByx*dF(2,0) + dFByy*dF(2,1) + dFByz*dF(2,2);
    pB->zz = dFBzx*dF(2,0) + dFBzy*dF(2,1) + dFBzz*dF(2,2);
    
    // return |dF|
    return dF.determinant();
#else
	// current deformation gradient in 3D
    double pF[3][3];
    mptr->GetDeformationGradient(pF);
	
	// get new deformation gradient for off axis components
	double Fxy = (1. + dvxx)*pF[0][1] + dvxy*pF[1][1] + dvxz*pF[2][1];		// du/dy
	double Fxz = (1. + dvxx)*pF[0][2] + dvxy*pF[1][2] + dvxz*pF[2][2];		// du/dz
	double Fyx = dvyx*pF[0][0] + (1. + dvyy)*pF[1][0] + dvyz*pF[2][0];		// dv/dx
	double Fyz = dvyx*pF[0][2] + (1. + dvyy)*pF[1][2] + dvyz*pF[2][2];		// dv/dz
	double Fzx = dvzx*pF[0][0] + dvzy*pF[1][0] + (1. + dvzz)*pF[2][0];		// dw/dx
	double Fzy = dvzx*pF[0][1] + dvzy*pF[1][1] + (1. + dvzz)*pF[2][1];		// dw/dy
	
	// store in total strain and rotation tensors
    Tensor *ep=mptr->GetStrainTensor();
    TensorAntisym *wrot = mptr->GetRotationStrainTensor();
        
    ep->xx = (1. + dvxx)*pF[0][0] + dvxy*pF[1][0] + dvxz*pF[2][0] - 1.;     // 1 + du/dx - 1
    ep->yy = dvyx*pF[0][1] + (1. + dvyy)*pF[1][1] + dvyz*pF[2][1] - 1.;		// 1 + dv/dy - 1
    ep->zz = dvzx*pF[0][2] + dvzy*pF[1][2] + (1. + dvzz)*pF[2][2] - 1.;		// 1 + dw/dz - 1
    ep->xy = Fyx + Fxy;
    ep->xz = Fzx + Fxz;
    ep->yz = Fzy + Fyz;
    
    // rotational strain increments
    wrot->xy = Fyx - Fxy;			// dv/dx - du/dy
    wrot->xz = Fzx - Fxz;			// dw/dx - du/dz
    wrot->yz = Fzy - Fyz;			// dw/dy - dv/dz
    
    // increment Left Cauchy tensor B = F.F^T = dF.old B.dF^T
    // plain stress will need to update B.zz when known
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    
    // elements of dF.B
    double dFBxx = (1.+dvxx)*pB->xx +      dvxy*pB->xy +      dvxz*pB->xz;
    double dFBxy = (1.+dvxx)*pB->xy +      dvxy*pB->yy +      dvxz*pB->yz;
    double dFBxz = (1.+dvxx)*pB->xz +      dvxy*pB->yz +      dvxz*pB->zz;
    double dFByx =      dvyx*pB->xx + (1.+dvyy)*pB->xy +      dvyz*pB->xz;
    double dFByy =      dvyx*pB->xy + (1.+dvyy)*pB->yy +      dvyz*pB->yz;
    double dFByz =      dvyx*pB->xz + (1.+dvyy)*pB->yz +      dvyz*pB->zz;
    double dFBzx =      dvzx*pB->xx +      dvzy*pB->xy + (1.+dvzz)*pB->xz;
    double dFBzy =      dvzx*pB->xy +      dvzy*pB->yy + (1.+dvzz)*pB->yz;
    double dFBzz =      dvzx*pB->xz +      dvzy*pB->yz + (1.+dvzz)*pB->zz;
    
    // return trial B (if provided) or store new B on the particle
    if(Btrial != NULL) pB = Btrial;
    pB->xx = dFBxx*(1+dvxx) + dFBxy*dvxy      + dFBxz*dvxz;
    pB->xy = dFBxx*dvyx     + dFBxy*(1.+dvyy) + dFBxz*dvyz;
    pB->xz = dFBxx*dvzx     + dFBxy*dvzy      + dFBxz*(1.+dvzz);
    pB->yy = dFByx*dvyx     + dFByy*(1.+dvyy) + dFByz*dvyz;
    pB->yz = dFByx*dvzx     + dFByy*dvzy      + dFByz*(1.+dvzz);
    pB->zz = dFBzx*dvzx     + dFBzy*dvzy      + dFBzz*(1.+dvzz);
    
    // return |dF|
    return (1. + dvxx)*((1. + dvyy)*(1. + dvzz)-dvzy*dvyz)
                - dvyx*(dvxy*(1. + dvzz)-dvzy*dvxz)
                + dvzx*(dvxy*dvyz-(1. + dvyy)*dvxz);
#endif
}

// Find isotropic stretch for thermal and moisture expansion
// total residual stretch (1 + alpha dT + beta csat dConcentration)
// Current assumes isotropic with CTE1 and CME1 expansion coefficients
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

// Return normal stress term (due to bulk modulus) and twice the pressure term (i.e. 2U(J)) for strain energy.
// Each block of lines is for a different U(J).
// Any change here must also be made in 2D MPMConstLaw for the numerical solution to find B.zz in plane stress
double HyperElastic::GetVolumetricTerms(double J,double *Kse)
{
    double Kterm;
    
    switch(UofJOption)
    {   case J_MINUS_1_SQUARED:
            // This is for *Kse/2 = U(J) = (K/2)(J-1)^2
            Kterm = Ksp*(J-1.);
            *Kse = Kterm*(J-1);
            break;
        
        case LN_J_SQUARED:
        {   // This is for for *Kse/2 = U(J) = (K/2)(ln J)^2
            // Zienkiewicz & Taylor recommend not using this one
            double lj = log(J);
            Kterm =Ksp*lj;
            *Kse = Kterm*lj;
            Kterm /= J;           // = Ksp*(ln J)/J
            break;
        }
        
        case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
        default:
            // This is for *Kse/2 = U(J) = (K/2)((1/2)(J^2-1) - ln J)
            // Zienkiewicz & Taylor note that stress is infinite as J->0 and J->infinity for this function, while others are not
            // Simo and Hughes also use this form (see Eq. 9.2.3)
            *Kse = Ksp*(0.5*(J*J-1.)-log(J));
            Kterm = 0.5*Ksp*(J - 1./J);      // = (Ksp/2)*(J - 1/J)
            break;
    }
    
    // return result
    return Kterm;
}  


