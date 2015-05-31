/********************************************************************************
    HyperElastic.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/HyperElastic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "System/UnitsController.hpp"

#pragma mark HyperElastic::Constructors and Destructors

// Constructors
HyperElastic::HyperElastic() {}

HyperElastic::HyperElastic(char *matName) : MaterialBase(matName)
{
	Kbulk = -1.;                                        // required (check >0 before starting)
    UofJOption = HALF_J_SQUARED_MINUS_1_MINUS_LN_J;     // default U(J) function
	aI=0.;
}

// Read material properties
char *HyperElastic::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"K")==0)
        return UnitsController::ScaledPtr((char *)&Kbulk,gScaling,1.e6);
    
    else if(strcmp(xName,"alpha")==0)
        return((char *)&aI);
    
    else if(strcmp(xName,"UJOption")==0)
    {   input = INT_NUM;
        return((char *)&UofJOption);
    }
    
    return(MaterialBase::InputMaterialProperty(xName,input,gScaling));
}

#pragma mark HyperElastic::Initialize

// Set intial particle Left Cauchy strain tensor to identity
void HyperElastic::SetInitialParticleState(MPMBase *mptr,int np) const
{
    // get previous particle B
    Tensor *pB = mptr->GetAltStrainTensor();
    
    ZeroTensor(pB);
    pB->xx = pB->yy = pB->zz = 1.;
	
	MaterialBase::SetInitialParticleState(mptr,np);
}

// Constant properties used in constitutive law
const char *HyperElastic::VerifyAndLoadProperties(int np)
{
	// Kbulk in Specific units using initial rho (F-L/mass)
	Ksp = Kbulk/rho;
	
	// expansion coefficients
	CTE1 = 1.e-6*aI;
	CME1 = betaI*concSaturation;
    
    // for Cp-Cv (units nJ/(g-K^2))
    Ka2sp = Ksp*CTE1*CTE1;
	
    // call superclass
    return MaterialBase::VerifyAndLoadProperties(np);
}

#pragma mark HyperElastic::Methods

/*  Given matrix of incremental deformation dF = exp(dt*grad v), increment particle strain,
        rotation, and LeftCauchy Green strain (latter is assumed to be stored in the particle's
        plastic strain tensor (which is accessed also with GetAltStrainTensor().
    New new F is dF.F, which is used to find new strain
    New B = dF.(Old B).dF^T
    Returns |dF|
*/
double HyperElastic::IncrementDeformation(MPMBase *mptr,Matrix3 du,Tensor *Btrial,int np) const
{
    // get incremental deformation gradient
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
	
	// current deformation gradient
	Matrix3 pF = mptr->GetDeformationGradientMatrix();
	
	// new deformation matrix
	const Matrix3 F = dF*pF;
    mptr->SetDeformationGradientMatrix(F);
	
    // increment Left Cauchy tensor B = F.F^T = dF.old B.dF^T
    // plain stress will need to update B.zz when known
    Matrix3 pBold = mptr->GetElasticLeftCauchyMatrix();
    
    // elements of dF.B
    Matrix3 dFoldB = dF*pBold;
    
    // return trial B (if provided) or store new B on the particle
    Tensor *pB = Btrial!=NULL ? Btrial : mptr->GetAltStrainTensor() ;
    pB->xx = dFoldB(0,0)*dF(0,0) + dFoldB(0,1)*dF(0,1) + dFoldB(0,2)*dF(0,2);
    pB->xy = dFoldB(0,0)*dF(1,0) + dFoldB(0,1)*dF(1,1) + dFoldB(0,2)*dF(1,2);
    
    pB->yy = dFoldB(1,0)*dF(1,0) + dFoldB(1,1)*dF(1,1) + dFoldB(1,2)*dF(1,2);
    
    pB->zz = dFoldB(2,0)*dF(2,0) + dFoldB(2,1)*dF(2,1) + dFoldB(2,2)*dF(2,2);
    
    if(np == THREED_MPM)
    {   pB->xz = dFoldB(0,0)*dF(2,0) + dFoldB(0,1)*dF(2,1) + dFoldB(0,2)*dF(2,2);
        pB->yz = dFoldB(1,0)*dF(2,0) + dFoldB(1,1)*dF(2,1) + dFoldB(1,2)*dF(2,2);
    }
    
    // return |dF|
    return dF.determinant();
}

// Find isotropic stretch for thermal and moisture expansion
// total residual stretch exp(alpha dT + beta csat dC)
// Current assumes isotropic with constant CTE1 and CME1 expansion coefficients
double HyperElastic::GetResidualStretch(MPMBase *mptr,double &dresStretch,ResidualStrains *res) const
{
	// total residual stretch exp(alpha dT(total) + beta csat dC(total))
	// incremental residual stretch (1 + alpha dT + beta csat dC)
	double dTemp=mptr->pPreviousTemperature-thermal.reference;
	double resStretch = CTE1*dTemp;
	dresStretch = CTE1*res->dT;
	if(DiffusionTask::active)
	{	double dConc=mptr->pPreviousConcentration-DiffusionTask::reference;
		resStretch += CME1*dConc;
		dresStretch += CTE1*res->dC;
	}
    dresStretch = exp(dresStretch);
	return exp(resStretch);
}

// Get current relative volume change = J = det F = lam1 lam2 lam3
// Need to have this call in material classes to allow small and large deformation material laws
//  to handle it differently. It is used on archiving to convert Kirchoff Stress/rho0 to Cauchy stress
double HyperElastic::GetCurrentRelativeVolume(MPMBase *mptr) const
{   return mptr->GetRelativeVolume();
}

// Return normal stress term (due to bulk modulus) and the pressure term (i.e. U(J)) for strain energy.
// Each block of lines is for a different U(J).
// Any change here must also be made in 2D MPMConstitutiveLaw for the numerical solution to find B.zz in plane stress
// Kse is strain energy, but no longer used
double HyperElastic::GetVolumetricTerms(double J,double Kred) const
{
    double Kterm;
    
    switch(UofJOption)
    {   case J_MINUS_1_SQUARED:
            // This is for *Kse = U(J) = (K/2)(J-1)^2
            Kterm = Kred*(J-1.);
            //*Kse = 0.5*Kterm*(J-1);
            break;
        
        case LN_J_SQUARED:
        {   // This is for for *Kse = U(J) = (K/2)(ln J)^2
            // Zienkiewicz & Taylor recommend not using this one
            double lj = log(J);
            Kterm = Kred*lj;
            //*Kse = 0.5*Kterm*lj;
            Kterm /= J;           // = Kred*(ln J)/J
            break;
        }
        
        case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
        default:
            // This is for *Kse = U(J) = (K/2)((1/2)(J^2-1) - ln J)
            // Zienkiewicz & Taylor note that stress is infinite as J->0 and J->infinity for this function, while others are not
            // Simo and Hughes also use this form (see Eq. 9.2.3)
            //*Kse = 0.5*Kred*(0.5*(J*J-1.)-log(J));
            Kterm = 0.5*Kred*(J - 1./J);      // = (Kred/2)*(J - 1/J)
            break;
    }
    
    // return result
    return Kterm;
}

// Return -J^2 P(J) and -d(J^2 P(J))/dJ
// These terms are used in Newton's method to zero out stress (such and in plane stress for Mooney
//    or in membranes for compressible materials)
void HyperElastic::GetNewtonPressureTerms(double J,double Kred,double &mJ2P,double &mdJ2PdJ) const
{
    switch(UofJOption)
    {   case J_MINUS_1_SQUARED:
            // This is for -J^2 P(J) = K(J^3-J^2)
            mJ2P = Kred*J*J*(J-1.);
            // K d(J^3-J^2) = K (3 J^2 - 2 J)
            mdJ2PdJ = Kred*J*(3.*J-2.);
            break;
            
        case LN_J_SQUARED:
            // This is for -J^2 P(J) = J K (ln J)
            mJ2P = Kred*J*log(J);
            // K d(J ln(J)) = K (ln(J) + 1)
            mdJ2PdJ = Kred*(log(J)+1.);
            break;
            
        case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
        default:
            // This is for -J^2 P(J) = (K/2) (J^3-J)
            mJ2P = 0.5*Kred*J*(J*J-1.);
            // (K/2) d(J^3 - J) = (K/2)*(3*J^2-1)
            mdJ2PdJ = 0.5*Kred*(3.*J*J-1.);
            break;
    }
}


// From thermodyanamics Cp-Cv = 9 K a^2 T/rho
// Ka2sp in nJ/(g-K^2) so output in nJ/(g-K)
// Here using K0 and rho0 - could modify if needed
double HyperElastic::GetCpMinusCv(MPMBase *mptr) const
{   return mptr!=NULL ? Ka2sp*mptr->pPreviousTemperature : Ka2sp*thermal.reference;
}

// store elastic B in alt strain
int HyperElastic::AltStrainContains(void) const { return LEFT_CAUCHY_TOTAL_B_STRAIN; }

