/********************************************************************************
    HEIsotropic.cpp
    NairnMPM
    
    Created by John Nairn, Sept 27, 2011.
    Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		HyperElastic.hpp (MaterialBase.hpp)
********************************************************************************/

#include "HEIsotropic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
// JAN: for hardenling law
#include "Materials/LinearHardening.hpp"

#pragma mark HEIsotropic::Constructors and Destructors

// Constructors
HEIsotropic::HEIsotropic() {}

// Constructors
HEIsotropic::HEIsotropic(char *matName) : HyperElastic(matName)
{
   	G1 = -1.;			// required
	G2 = 0.;			// zero is neo-Hookean
    
	// JAN: hard-code linear hardening - future can set to other laws
	plasticLaw = new LinearHardening(this);
}


#pragma mark HEIsotropic::Initialization

// Read material properties
char *HEIsotropic::InputMat(char *xName,int &input)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"G1")==0)
        return((char *)&G1);
    
    else if(strcmp(xName,"G2")==0)
        return((char *)&G2);

    // JAN: Move yielding properties to the hardening law (yield, Ep, Khard)
	// check plastic law
    char *ptr = plasticLaw->InputMat(xName,input);
    if(ptr != NULL) return ptr;
    
    return(HyperElastic::InputMat(xName,input));
}

// verify settings and maybe some initial calculations
const char *HEIsotropic::VerifyProperties(int np)
{
	// JAN: ask hardening law to verify properties too
	// call plastic law first
    const char *ptr = plasticLaw->VerifyProperties(np);
    if(ptr != NULL) return ptr;
    
    // user must enter G1 and Kbulk
    if(G1<0. || Kbulk < 0. )
		return "HEIsotropic Material needs non-negative G1 and K";
    
	// must call super class
	return HyperElastic::VerifyProperties(np);
}

// plane stress not allowed in viscoelasticity
void HEIsotropic::ValidateForUse(int np)
{	if(np==PLANE_STRESS_MPM || np==AXISYMMETRIC_MPM)
    {	throw CommonException("HEIsotropic materials require 2D plane strain or 3D MPM analysis",
                                "HEIsotropic::ValidateForUse");
    }
	
	//call super class (why can't call super class?)
	return HyperElastic::ValidateForUse(np);
}

// Constant properties used in constitutive law
void HEIsotropic::InitialLoadMechProps(int makeSpecific,int np)
{
	// JAN: hardening law gets reduced yielding properties
    // hardening law properties
    plasticLaw->InitialLoadMechProps(makeSpecific,np);
	
	// MU in specific units using initial rho
	// for MPM (units N/m^2 cm^3/g)
	G1sp=G1*1.0e+06/rho;
    G2sp=0.;
	
	// reduced properties
	Gred = G1*1.e6/rho;
    
    // call super class for the rest
    HyperElastic::InitialLoadMechProps(makeSpecific,np);
}

// print mechanical properties to the results
void HEIsotropic::PrintMechanicalProperties(void)
{	
	// call superclass here if it is not Material base
	HyperElastic::PrintMechanicalProperties();
	
    PrintProperty("G1",G1,"");
    PrintProperty("K",Kbulk,"");
    cout << endl;
        
    PrintProperty("a",aI,"");
    cout << endl;
    
	// JAN: print hardening law properties
	plasticLaw->PrintYieldProperties();
}

// JAN: put hardening law properties first and may need more than 1
// First ones for hardening law. Particle J appended at the end
char *HEIsotropic::MaterialData(void)
{	J_history = plasticLaw->HistoryDoublesNeeded();
	double *p = CreateAndZeroDoubles(J_history+1);
	p[J_history]=1.;                // J
	return (char *)p;
}

#pragma mark HEIsotropic:Methods

//===============================================================================================================
/*	Apply 2D constitutive law updating all needed terms for material type. Required updates are:
		stress, strain, plastic strain (all components) (stress should be a specific stress)
		rotation strain (single angle)
		strain energy, plastic energy, and dissipated energy (dissipated needed if want to couple to conduction)
	To support thermal and solvent expansion, include their effect on strains
    If there are material-related data on the particle, update them too
	dvij are (gradient rates X time increment) to give deformation gradient change
    For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
*/
//===============================================================================================================
void HEIsotropic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,double dvzz,
                                        double delTime,int np)
{
    // Compute Elastic Predictor
    // ============================================
    
	// Update total deformation gradient, and calculate trial B
    Tensor Btrial;
	double detdF = IncrementDeformation(mptr,dvxx,dvyy,dvxy,dvyx,dvzz,&Btrial);
    
	// Deformation gradients and Cauchy tensor differ in plane stress and plane strain
	double J2;
	if(np==PLANE_STRAIN_MPM)
	{
        //J2 = Btrial.xx*Btrial.yy - Btrial.xy*Btrial.xy;
      
        //Incremental calculation of J det(F)=det(dF)det(pF)
        J2 = detdF*mptr->GetHistoryDble(J_history);
	}
	else
	{	
		// first we don't consider a 2D plane stress
	}
    
    // save new J
    mptr->SetHistoryDble(J_history,J2);  // Stocking J
    
    // J as determinant of F (or sqrt root of determinant of B) normalized to residual stretch
	double resStretch = GetResidualStretch(mptr);
    double J = J2/(resStretch*resStretch*resStretch);
    //cout << "    #J  =   "<< J<<  endl;
	
	// JAN: Get hydrostatic stress component in subrountine
	double Kse;
	double sHydro = GetDilationalTerms(mptr,J,np,Kse);
    
    // Others constants
    double J23 = pow(J, 2./3.);
    
    // find Trial Specific Kirchoff stress Tensor (Trial_Tau/rho0)
    Tensor stk = GetTrialStressTensor2D(&Btrial,J23);
   
    // Checking for plastic loading
    // ============================================
    
    // Get magnitude of the deviatoric stress tensor
    // ||s|| = sqrt(2J2) where J2 = (1/6)((sx-sy)^2+(sy-sz)^2+(sx-sz)^2) + txy^2+txz^2+tyz^2
    // In 2D 2J2 = (2/3)(sx^2 + sy^2 - sx*sy + sz*(sz-sx-sy)) + 2*txy^2
    
	// JAN: set alpint for particle (not sure where done before)
	plasticLaw->UpdateTrialAlpha(mptr,np);			// initialize to last value and zero plastic strain rate
	// JAN: this not used because found below instead
    //magnitude_strial = sqrt(0.6666667*(stk.xx*stk.xx+stk.yy*stk.yy+stk.zz*stk.zz-stk.xx*stk.yy-stk.xx*stk.zz-stk.yy*stk.zz)+2*stk.xy*stk.xy);
    double magnitude_strial = GetMagnitudeS(&stk,np);
    double gyld=plasticLaw->GetYield(mptr,np,delTime);
    double ftrial = magnitude_strial-SQRT_TWOTHIRDS*gyld;
    //cout << "  #magnitude_strial =   "<< magnitude_strial<< "  GetYield =   "<< gyld<< "  ftrial =   "<< ftrial<< endl;
    //cout << "  #yldred =   "<< yldred << "  Epred =   "<< Epred << "  gyld =   "<< gyld <<"  alpint =   "<< alpint<< "  ftrial =   "<< ftrial<< endl;
        
    //============================
    //  TEST
    //============================
    //ftrial=-1.;
    if(ftrial<=0.)
        
        // if elastic
        //============================
    {   Tensor B = Btrial;
        
        // save on particle

        Tensor *pB = mptr->GetElasticLeftCauchyTensor();
        pB->xx = B.xx;
        pB->yy = B.yy;
        pB->zz = B.zz;
        pB->xy = B.xy;
        //cout << "# in elastic  B.yy  =    " << Btrial.yy <<"     B.zz  =    " << Btrial.zz << endl;
        
        // Get specifique stress i.e. (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
        Tensor *sp=mptr->GetStressTensor();
		
		// JAN: Use sHuydro calculated in subroutine
		// JAN: No G2sp?
        sp->xx = sHydro + stk.xx;		//1./3.*G1sp*(2.*B.xx-B.yy-B.zz)/J23;
        sp->yy = sHydro + stk.yy;		//1./3.*G1sp*(2.*B.yy-B.xx-B.zz)/J23;
        sp->zz = sHydro + stk.zz;		//1./3.*G1sp*(2.*B.zz-B.xx-B.yy)/J23;
        sp->xy = stk.xy;				//G1sp*B.xy/J23;
        
        //  J = 1
        //sp->xx = 1./3.*G1sp*(2.*B.xx-B.yy-B.zz);
        //sp->yy = 1./3.*G1sp*(2.*B.yy-B.xx-B.zz);
        //sp->zz = 1./3.*G1sp*(2.*B.zz-B.xx-B.yy);
        //sp->xy = G1sp*B.xy;
        
        // strain energy per unit mass (U/(rho0 V0)) and we are using
        // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
		// JAN: May need new energy methods
		// JAN: Note that you Ksp*(1/2*(J*J-1)-log(J)) was a bug to or operature precedence
		//       and gave negative energies
        double I1bar = (B.xx+B.yy+B.zz)/J23;
        mptr->SetStrainEnergy(0.5*(G1sp*(I1bar-3.) + Kse));
        
        return;
    }
    // Plastic behavior - Return Mapping algorithm 
    //=====================================================
    // if plastic
    
	// JAN: Use hardening law method (will need future expansion for other laws)
    Tensor B = Btrial;
    double Ie1bar = (1./3.)*(B.xx+B.yy+B.zz)/J23;
    double MUbar = Gred*Ie1bar;
    
    // Find  lambda for this plastic state
	double dlambda = plasticLaw->HESolveForLambdaBracketed(mptr,np,magnitude_strial,Gred,Ie1bar);

	Tensor nk = GetNormalTensor2D(&stk,magnitude_strial);
    //cout << "nk.xx  =    " << nk.xx << "nk.xy  =    " << nk.xy << endl;
    
    Tensor *sp=mptr->GetStressTensor();
    sp->xx = sHydro + (stk.xx-2.*MUbar*dlambda*nk.xx);
    sp->yy = sHydro + (stk.yy-2.*MUbar*dlambda*nk.yy);
    sp->zz = sHydro + (stk.zz-2.*MUbar*dlambda*nk.zz);
    sp->xy = (stk.xy-2.*MUbar*dlambda*nk.xy);
    
    
    //cout << "EXT B.xx  =    " << B.xx << endl;
    // save on particle
    
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
     pB->xx = ((stk.xx-2.*MUbar*dlambda*nk.xx)/Gred+Ie1bar)*J23;
     pB->yy = ((stk.yy-2.*MUbar*dlambda*nk.yy)/Gred+Ie1bar)*J23;
     pB->zz = ((stk.zz-2.*MUbar*dlambda*nk.zz)/Gred+Ie1bar)*J23;
     pB->xy = (stk.xy-2.*MUbar*dlambda*nk.xy)/Gred*J23;
     //cout << "# in pB plastic  nk.xx  =    " << nk->xx <<"     nk.yy  =    " << nk.yy<<"     nk.zz  =    " << nk.zzß << endl;
    
        
    // strain energy per unit mass (U/(rho0 V0)) and we are using
    // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
	// JAN: will need to get elastic and plastic energy
    double I1bar = (B.xx+B.yy+B.zz)/J23;
    mptr->SetStrainEnergy(0.5*(G1sp*(I1bar-3.) + Kse));

	// JAN: need call to update hardening law properties. Might revise to have in done in the solve method instead
	// update internal variables
	plasticLaw->UpdatePlasticInternal(mptr,np);
}

// 2D implementation of the strial stress tensor


Tensor HEIsotropic::GetTrialStressTensor2D(Tensor *B,double J23)
{
    // Trial specific Kirchhoff Stress Tensor 
    Tensor strial;
    ZeroTensor(&strial);
    strial.xx = 1./3.*Gred*(2.*B->xx-B->yy-B->zz)/J23;
    strial.yy = 1./3.*Gred*(2.*B->yy-B->xx-B->zz)/J23;
    strial.zz = 1./3.*Gred*(2.*B->zz-B->xx-B->yy)/J23;
    strial.xy = Gred*B->xy/J23;
    return strial;
}

double HEIsotropic::GetMagnitudeS(Tensor *st,int np)
{
	double s,t;
	
	switch(np)
	{	case PLANE_STRAIN_MPM:
			s = fmax(st->xx*st->xx + st->yy*st->yy - st->xx*st->yy + st->zz*(st->zz-st->xx-st->yy),0.);
			t = st->xy*st->xy;
			break;
		case PLANE_STRESS_MPM:
			s = fmax(st->xx*st->xx + st->yy*st->yy - st->xx*st->yy,0.);
			t = st->xy*st->xy;
			break;
		case THREED_MPM:
		default:
			s = fmax(st->xx*st->xx + st->yy*st->yy - st->xx*st->yy + st->zz*(st->zz-st->xx-st->yy),0.);
			t = st->xy*st->xy + st->xz*st->xz + st->yz*st->yz;
			break;
	}
	return sqrt(2.*(s/3. + t));
}
 

// 2D implementation of the normal to the yield surface

Tensor HEIsotropic::GetNormalTensor2D(Tensor *strial,double magnitude_strial)
{
    // Trial specific Kirchhoff Stress Tensor
    Tensor n;
    ZeroTensor(&n);
    n.xx = strial->xx/magnitude_strial;
    n.yy = strial->yy/magnitude_strial;
    n.xy = strial->xy/magnitude_strial;
    n.zz = strial->zz/magnitude_strial;
    //cout << "strial.yy  =    " << strial->yy << "      strial.zz  =    " << strial->zz << "magnitude_strial  =    " << magnitude_strial << endl;
    //cout << "n.xx  =    " << n.xx << "n.xy  =    " << n.xy << endl;
    return n;
}


//===============================================================================================================
/* Apply 3D constitutive law updating all needed terms for material type. Required updates are:
		stress, strain, plastic strain (all components) (stress should be a specific stress)
		rotation strain
		strain energy, plastic energy, and dissipated energy (dissipated needed if want to couple to conduction)
	To support thermal and solvent expansion, include their effect on strains
    If there are material-related data on the particle, update them too
    dvij are (gradient rates X time increment) to give deformation gradient change
*/
//===============================================================================================================


// 3D Stress Calculation New Law
void HEIsotropic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
                              double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)

{
	// Update strains and rotations and Left Cauchy strain
    Tensor Btrial;
	double detdF = IncrementDeformation(mptr,dvxx,dvyy,dvzz,dvxy,dvyx,dvxz,dvzx,dvyz,dvzy,&Btrial);
    
    // Incremental calculation of J det(F)=det(dF)det(pF)
    double J2 = detdF*mptr->GetHistoryDble(J_history);
    mptr->SetHistoryDble(J_history,J2);  // Stocking J
    
    // J as determinant of F (or sqrt root of determinant of B) normalized to residual stretch
	double resStretch = GetResidualStretch(mptr);
    double J = J2/(resStretch*resStretch*resStretch);
    //cout << " # J  =    " << J <<  endl;
    
	// JAN: Get hydrostatic stress component in subrountine
	double Kse;
	double sHydro = GetDilationalTerms(mptr,J,np,Kse);
    
    double J23=pow(J,2./3.);
    
    // find Trial Specific Kirchoff stress Tensor (Trial_Tau/rho0)
    
    Tensor stk = GetTrialStressTensor3D(&Btrial,J23);
    
    // Checking for plastic loading
    // ============================================
    
    // Get magnitude of the deviatoric stress tensor
  
    // magnitude_strial = sqrt(0.6666667*(stk.xx*stk.xx+stk.yy*stk.yy+stk.zz*stk.zz-stk.xx*stk.yy-stk.xx*stk.zz-stk.yy*stk.zz)+2*stk.xy*stk.xy+2*stk.xz*stk.xz+2*stk.yz*stk.yz);
    
	// JAN: set alpint for particle (not sure where done before)
	plasticLaw->UpdateTrialAlpha(mptr,np);			// initialize to last value and zero plastic strain rate
    double magnitude_strial = GetMagnitudeS(&stk,np);
    double gyld=plasticLaw->GetYield(mptr,np,delTime);
    double ftrial = magnitude_strial-SQRT_TWOTHIRDS*gyld;
    
    //cout << "  #yldred =   "<< yldred << "  Epred =   "<< Epred << "  gyld =   "<< gyld <<"  alpint =   "<< alpint<< "  ftrial =   "<< ftrial<< endl;
    
    //============================
    //  TEST
    //============================
    //ftrial=-1.;
    //cout << " # ftrial  =    " << ftrial <<  endl;
    if(ftrial<=0.)

        // if elastic
        //============================
    {   Tensor B = Btrial;
        
        // save on particle
        
        Tensor *pB = mptr->GetElasticLeftCauchyTensor();
        pB->xx = B.xx;
        pB->yy = B.yy;
        pB->zz = B.zz;
        pB->xy = B.xy;
        pB->xz = B.xz;
        pB->yz = B.yz;
        //cout << "# in elastic  B.yy  =    " << Btrial.yy <<"     B.zz  =    " << Btrial.zz << endl;
        
        // Get specifique stress i.e. (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
        Tensor *sp=mptr->GetStressTensor();
        
		// JAN: will want to move hydrostatic stress to subroutine
		// JAN: G2sp:
        sp->xx = sHydro + stk.xx;		//1./3.*G1sp*(2.*B.xx-B.yy-B.zz)/J23;
        sp->yy = sHydro + stk.yy;		//1./3.*G1sp*(2.*B.yy-B.xx-B.zz)/J23;
        sp->zz = sHydro + stk.zz;		//1./3.*G1sp*(2.*B.zz-B.xx-B.yy)/J23;
        sp->xy = stk.xy;				//G1sp*B.xy/J23;
        sp->xz = stk.xz;				//G1sp*B.xy/J23;
        sp->yz = stk.yz;				//G1sp*B.xy/J23;
        
        //  J = 1   for Simple Shear test
        
        //sp->xx = 1./3.*G1sp*(2.*B.xx-B.yy-B.zz);
        //sp->yy = 1./3.*G1sp*(2.*B.yy-B.xx-B.zz);
        //sp->zz = 1./3.*G1sp*(2.*B.zz-B.xx-B.yy);
        //sp->xy = G1sp*B.xy;
        //sp->xz = G1sp*B.xz;
        //sp->yz = G1sp*B.yz;
        
        // strain energy per unit mass (U/(rho0 V0)) and we are using
        // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
		// JAN: may need new energy methods
        double I1bar = (B.xx+B.yy+B.zz)/J23;
        mptr->SetStrainEnergy(0.5*(G1sp*(I1bar-3.) + Kse));
        
        return;
    }
    
    // Plastic behavior - Return Mapping algorithm
    //=====================================================
    // if plastic
    
	// JAN: Use hardening law method (will need future expansion for other laws)
    Tensor B = Btrial;
    double Ie1bar = (1./3.)*(B.xx+B.yy+B.zz)/J23;
    double MUbar = Gred*Ie1bar;
    
    // Find  lambda for this plastic state ueing hardening law
	double dlambda = plasticLaw->HESolveForLambdaBracketed(mptr,np,magnitude_strial,Gred,Ie1bar);
    
    Tensor nk = GetNormalTensor3D(&stk,magnitude_strial);
    //cout << "nk.xx  =    " << nk.xx << "nk.xy  =    " << nk.xy << endl;
    
    Tensor *sp=mptr->GetStressTensor();
    sp->xx = sHydro + (stk.xx-2.*MUbar*dlambda*nk.xx);
    sp->yy = sHydro + (stk.yy-2.*MUbar*dlambda*nk.yy);
    sp->zz = sHydro + (stk.zz-2.*MUbar*dlambda*nk.zz);
    sp->xy = (stk.xy-2.*MUbar*dlambda*nk.xy);
    sp->xz = (stk.xz-2.*MUbar*dlambda*nk.xz);
    sp->yz = (stk.yz-2.*MUbar*dlambda*nk.yz);
    
    
    //cout << "EXT B.xx  =    " << B.xx << endl;
    // save on particle
    
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    pB->xx = ((stk.xx-2.*MUbar*dlambda*nk.xx)/Gred+Ie1bar)*J23;
    pB->yy = ((stk.yy-2.*MUbar*dlambda*nk.yy)/Gred+Ie1bar)*J23;
    pB->zz = ((stk.zz-2.*MUbar*dlambda*nk.zz)/Gred+Ie1bar)*J23;
    pB->xy = (stk.xy-2.*MUbar*dlambda*nk.xy)/Gred*J23;
    pB->xz = (stk.xz-2.*MUbar*dlambda*nk.xz)/Gred*J23;
    pB->yz = (stk.yz-2.*MUbar*dlambda*nk.yz)/Gred*J23;
    
    //cout << "# in pB plastic  nk.xx  =    " << nk->xx <<"     nk.yy  =    " << nk.yy<<"     nk.zz  =    " << nk.zzß << endl;
    
    
    // strain energy per unit mass (U/(rho0 V0)) and we are using
    // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
    double I1bar = (B.xx+B.yy+B.zz)/J23;
    mptr->SetStrainEnergy(0.5*(G1sp*(I1bar-3.) + Kse));
    
	// JAN: need call to update hardening law properties. Might revise to have in done in the solve method instead
	// update internal variables
	plasticLaw->UpdatePlasticInternal(mptr,np);
}

// Return normal stress term (due to bulk modulus) and twice the pressure term for strain energy.
// Each block of lines is for a different U(J).
// Any change here must also be made in 2D MPMConstLaw for the numerical solution to find B.zz in plane stress



#pragma mark HEIsotropic::Custom Methods

// JAN: To allow better subclassing it is better to separate out calculations
//  of dilaiton energy. This version finds 1/3 new pressure (i.e. dilational
//  contribution to normal stress) and an energy term. This approach might
//  change with experience in using alternative laws
// Return normal Kirchoff stress/rho_0 and 2U(J)
// Note: This isolates Ksp to only be used here so subclass can replace wit
//  with something else
double HEIsotropic::GetDilationalTerms(MPMBase *mptr,double J,int np,double &Kse)
{
	Kse = Ksp*(0.5*(J*J-1)-log(J));
	return 0.5*Ksp*(J*J-1.);
}


// 2D implementation of the strial stress tensor


Tensor HEIsotropic::GetTrialStressTensor3D(Tensor *B,double J23)
{
    // Trial specific Kirchhoff Stress Tensor
    Tensor strial;
    ZeroTensor(&strial);
    strial.xx = 1./3.*Gred*(2.*B->xx-B->yy-B->zz)/J23;
    strial.yy = 1./3.*Gred*(2.*B->yy-B->xx-B->zz)/J23;
    strial.zz = 1./3.*Gred*(2.*B->zz-B->xx-B->yy)/J23;
    strial.xy = Gred*B->xy/J23;
    strial.xz = Gred*B->xz/J23;
    strial.yz = Gred*B->yz/J23;
   
    return strial;
}

// 2D implementation of the normal to the yield surface

Tensor HEIsotropic::GetNormalTensor3D(Tensor *strial,double magnitude_strial)
{
    // Trial specific Kirchhoff Stress Tensor
    Tensor n;
    ZeroTensor(&n);
    n.xx = strial->xx/magnitude_strial;
    n.yy = strial->yy/magnitude_strial;
    n.zz = strial->zz/magnitude_strial;
    n.xy = strial->xy/magnitude_strial;
    n.xz = strial->xz/magnitude_strial;
    n.yz = strial->yz/magnitude_strial;
    //cout << "strial.yy  =    " << strial->yy << "      strial.zz  =    " << strial->zz << "magnitude_strial  =    " << magnitude_strial << endl;
    //cout << "n.xx  =    " << n.xx << "n.xy  =    " << n.xy << endl;
    return n;
}


#pragma mark HEIsotropic::Accessors

// Return the material tag
int HEIsotropic::MaterialTag(void) { return HEISOTROPIC; }

// return unique, short name for this material
const char *HEIsotropic::MaterialType(void) { return "Hyperelastic Isotropic"; }

// calculate wave speed in mm/sec - needs to be implemented
double HEIsotropic::WaveSpeed(bool threeD,MPMBase *mptr)
{
    return sqrt(1.e9*(Kbulk+4.*G1/3.)/rho);
}

// this material has two history variables
double HEIsotropic::GetHistory(int num,char *historyPtr)
{
    double history=0.;
	if(num==1 || num==2)
	{	double *p=(double *)historyPtr;
		history=p[num-1];
	}
    return history;
}


