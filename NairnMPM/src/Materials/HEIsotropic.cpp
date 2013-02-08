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
#include "Materials/LinearHardening.hpp"
#include "Materials/NonlinearHardening.hpp"
#include "Materials/JohnsonCook.hpp"
#include "Materials/SCGLHardening.hpp"
#include "Materials/SLMaterial.hpp"
#include "Materials/Nonlinear2Hardening.hpp"

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

    // look for different plastic law
    if(strcmp(xName,"Hardening")==0)
    {	input = HARDENING_LAW_SELECTION;
        return (char *)this;
    }
    
    // JAN: Move yielding properties to the hardening law (yield, Ep, Khard)
	// check plastic law
    char *ptr = plasticLaw->InputMat(xName,input);
    if(ptr != NULL) return ptr;
    
    return(HyperElastic::InputMat(xName,input));
}

// change hardening law
void HEIsotropic::SetHardeningLaw(char *lawName)
{
    // delete old one
    delete plasticLaw;
    plasticLaw = NULL;
    
    // check options
    if(strcmp(lawName,"Linear")==0 || strcmp(lawName,"1")==0)
        plasticLaw = new LinearHardening(this);
    
    else if(strcmp(lawName,"Nonlinear")==0 || strcmp(lawName,"2")==0)
        plasticLaw = new NonlinearHardening(this);
    
    else if(strcmp(lawName,"Nonlinear2")==0 || strcmp(lawName,"6")==0)
        plasticLaw = new Nonlinear2Hardening(this);
    
    else if(strcmp(lawName,"JohnsonCook")==0 || strcmp(lawName,"3")==0)
        plasticLaw = new JohnsonCook(this);
    
    else if(strcmp(lawName,"SCGL")==0 || strcmp(lawName,"4")==0)
        plasticLaw = new SCGLHardening(this);
    
    else if(strcmp(lawName,"SL")==0 || strcmp(lawName,"5")==0)
        plasticLaw = new SLMaterial(this);
    
    // did it work
    if(plasticLaw == NULL)
    {   char errMsg[250];
        strcpy(errMsg,"The hardening law '");
        strcat(errMsg,lawName);
        strcat(errMsg,"' is not valid");
        throw CommonException(errMsg,"IsoPlasticity::SetHardeningLaw");
    }
    
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
{	if(np==PLANE_STRESS_MPM)
    {	throw CommonException("HEIsotropic materials cannot be used in plane stress MPM yet",
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
    switch(UofJOption)
    {   case J_MINUS_1_SQUARED:
            PrintProperty("U(J)",UofJOption,"[ = (K/2)(J-1)^2 ]");
            break;
            
        case LN_J_SQUARED:
            PrintProperty("U(J)",UofJOption,"[ = (K/2)(ln J)^2 ]");
            break;
            
        case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
        default:
            PrintProperty("U(J)",UofJOption,"[ = (K/2)((1/2)(J^2-1) - ln J) ]");
            break;
    }
    cout << endl;
    
	// JAN: print hardening law properties
	plasticLaw->PrintYieldProperties();
}

// JAN: put hardening law properties first and may need more than 1
// First ones for hardening law. Particle J appended at the end
char *HEIsotropic::InitHistoryData(void)
{	J_history = plasticLaw->HistoryDoublesNeeded();
	double *p = CreateAndZeroDoubles(J_history+1);
	p[J_history]=1.;					// J
	return (char *)p;
}

#pragma mark HEIsotropic:Methods

/* Take increments in strain and calculate new Particle: strains, rotation strain,
        stresses, strain energy,
    dvij are (gradient rates X time increment) to give deformation gradient change
    For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
    This material tracks pressure and stores deviatoric stress only in particle stress
        tensor
*/
void HEIsotropic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np)
{
    // store initial stress
    Tensor *sp = mptr->GetStressTensor();
    Tensor st0 = *sp;
    
    // Compute Elastic Predictor
    // ============================================
    
	// Update total deformation gradient, and calculate trial B
    Tensor B;
	double detdF = IncrementDeformation(mptr,du,&B,np);
    
	// Deformation gradients and Cauchy tensor differ in plane stress and plane strain
    // This code handles plane strain, axisymmetric, and 3D - Plane stress is blocked
	double J2 = detdF * mptr->GetHistoryDble(J_history);
    
    // save new J
    mptr->SetHistoryDble(J_history,J2);  // Stocking J
    
    // J as determinant of F (or sqrt root of determinant of B) normalized to residual stretch
	double resStretch = GetResidualStretch(mptr);
    double resStretch3 = resStretch*resStretch*resStretch;
    double J = J2/resStretch3;
    //cout << "    #J  =   "<< J<<  endl;
	
	// JAN: Get hydrostatic stress component in subroutine
    UpdatePressure(mptr,J,detdF,np);
    
    // Others constants
    double J23 = pow(J, 2./3.);
    
    // find Trial Specific Kirchoff stress Tensor (Trial_Tau/rho0)
    Tensor stk = GetTrialDevStressTensor(&B,J23,np);
    
    // Checking for plastic loading
    // ============================================
    
    // Get magnitude of the deviatoric stress tensor
    // ||s|| = sqrt(s.s)
    
	// JAN: set alpint for particle (not sure where done before)
	plasticLaw->UpdateTrialAlpha(mptr,np);                  // initialize to last value and zero plastic strain rate
	// JAN: this not used because found below instead
    //magnitude_strial = sqrt(0.6666667*(stk.xx*stk.xx+stk.yy*stk.yy+stk.zz*stk.zz-stk.xx*stk.yy-stk.xx*stk.zz-stk.yy*stk.zz)+2*stk.xy*stk.xy);
    double magnitude_strial = GetMagnitudeS(&stk,np);
    double gyld = plasticLaw->GetYield(mptr,np,delTime);
    double ftrial = magnitude_strial-SQRT_TWOTHIRDS*gyld;
    //cout << "  #magnitude_strial =   "<< magnitude_strial<< "  GetYield =   "<< gyld<< "  ftrial =   "<< ftrial<< endl;
    //cout << "  #yldred =   "<< yldred << "  Epred =   "<< Epred << "  gyld =   "<< gyld <<"  alpint =   "<< alpint<< "  ftrial =   "<< ftrial<< endl;
    
    // these will be needed for elastic or plasti
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    
    //============================
    //  TEST
    //============================
    if(ftrial<=0.)
	{	// if elastic
        //============================
		
		// save on particle
        *pB = B;
        //cout << "# in elastic  B.yy  =    " << Btrial.yy <<"     B.zz  =    " << Btrial.zz << endl;
        
        // Get specifique stress i.e. (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
		// JAN: Just store deviatoric stress
        *sp = stk;
        
        // strain energy per unit mass (U/(rho0 V0)) and we are using
        // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
		// JAN: May need new energy methods
        mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*du(0,0)
                                   + (st0.yy+sp->yy)*du(1,1)
                                   + (st0.zz+sp->zz)*du(2,2)
                                   + (st0.xy+sp->xy)*(du(1,0)+du(0,1)))/resStretch);
        if(np==THREED_MPM)
        {   mptr->AddStrainEnergy(0.5*((st0.yz+sp->yz)*(du(2,1)+du(1,2))
                                       + (st0.xz+sp->xz)*(du(2,0)+du(0,2)))/resStretch);
        }
        //double I1bar = (B.xx+B.yy+B.zz)/J23;
		//double shearEnergyFinal = 0.5*(Gred*(I1bar-3.));
        //mptr->AddStrainEnergy(shearEnergyFinal);
        
        return;
    }
    
    // Plastic behavior - Return Mapping algorithm 
    //=====================================================
    // if plastic
    
	// JAN: Use hardening law method (which can now ue other laws too)
    double Ie1bar = (1./3.)*(B.xx+B.yy+B.zz)/J23;
    double MUbar = Gred*Ie1bar;
    
    // Find  lambda for this plastic state
	double dlambda = plasticLaw->SolveForLambdaBracketed(mptr,np,magnitude_strial,&stk,MUbar,1.,1.,delTime);
    
	Tensor nk = GetNormalTensor(&stk,magnitude_strial,np);
    //cout << "nk.xx  =    " << nk.xx << "nk.xy  =    " << nk.xy << endl;
    
    // update deviatoric stress
    sp->xx = stk.xx - 2.*MUbar*dlambda*nk.xx;
    sp->yy = stk.yy - 2.*MUbar*dlambda*nk.yy;
    sp->zz = stk.zz - 2.*MUbar*dlambda*nk.zz;
    sp->xy = stk.xy - 2.*MUbar*dlambda*nk.xy;
    if(np==THREED_MPM)
    {   sp->xz = stk.xz - 2.*MUbar*dlambda*nk.xz;
        sp->xz = stk.yz - 2.*MUbar*dlambda*nk.yz;
    }
    
    // save on particle
    // JAN: reuse stress rather than calculate again
	pB->xx = (sp->xx/Gred+Ie1bar)*J23;
	pB->yy = (sp->yy/Gred+Ie1bar)*J23;
	pB->zz = (sp->zz/Gred+Ie1bar)*J23;
	pB->xy = sp->xy/Gred*J23;
    if(np==THREED_MPM)
    {   pB->xz = sp->xz/Gred*J23;
        pB->yz = sp->yz/Gred*J23;
    }
    
    // strain energy per unit mass (U/(rho0 V0)) and we are using
    // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
	// JAN: I1bar = 3 Ie1bar so no need to recalculate
    //double I1bar = (pB->xx+pB->yy+pB->zz)/J23;
    mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*du(0,0)
                               + (st0.yy+sp->yy)*du(1,1)
                               + (st0.zz+sp->zz)*du(2,2)
                               + (st0.xy+sp->xy)*(du(1,0)+du(0,1)))/resStretch);
    if(np==THREED_MPM)
    {   mptr->AddStrainEnergy(0.5*((st0.yz+sp->yz)*(du(2,1)+du(1,2))
                                   + (st0.xz+sp->xz)*(du(2,0)+du(0,2)))/resStretch);
    }
    //double shearEnergyFinal = 1.5*Gred*(Ie1bar-1.);
    //mptr->AddStrainEnergy(shearEnergyFinal);
    
    // Plastic or dissipated energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
    double dispEnergy = dlambda*(sp->xx*nk.xx + sp->yy*nk.yy + sp->zz*nk.zz + 2.*sp->xy*nk.xy);
    if(np==THREED_MPM)  dispEnergy += 2.*dlambda*(sp->xz*nk.xz + sp->yz*nk.yz);
    dispEnergy -= dlambda*SQRT_TWOTHIRDS*plasticLaw->GetYieldIncrement(mptr,np,delTime);;
    
    // increment energies
    mptr->AddStrainEnergy(-dispEnergy);
    mptr->AddPlastEnergy(dispEnergy);
	mptr->AddDispEnergy(dispEnergy);
	
	// JAN: need call to update hardening law properties. Might revise to have in done in the solve method instead
	// update internal variables
	plasticLaw->UpdatePlasticInternal(mptr,np);
}

#pragma mark HEIsotropic::Custom Methods

// JAN: To allow better subclassing it is better to separate out calculations
//  of dilaiton energy. This version finds 1/3 new pressure (i.e. dilational
//  contribution to normal stress) and an energy term. This approach might
//  change with experience in using alternative laws
// Return normal Kirchoff stress/rho_0 and 2U(J)
// Note: This isolates Ksp to only be used here so subclass can replace wit
//  with something else
void HEIsotropic::UpdatePressure(MPMBase *mptr,double J,double dJ,int np)
{
    double Kse;
	double Kterm = J*GetVolumetricTerms(J,&Kse);       // times J to get Kirchoff stress
    double P0 = mptr->GetPressure();
    mptr->SetPressure(-Kterm);
    
    // Elastic energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
    double avgP = 0.5*(P0-Kterm);
    double delV = 1. - 1./dJ;
    mptr->AddStrainEnergy(-avgP*delV);
    //mptr->SetStrainEnergy(Kse);
}

// get trial deviatoric stress tensor based on trial B
Tensor HEIsotropic::GetTrialDevStressTensor(Tensor *B,double J23,int np)
{
    // Trial specific Kirchhoff Stress Tensor
    Tensor strial;
    ZeroTensor(&strial);
    strial.xx = 1./3.*Gred*(2.*B->xx-B->yy-B->zz)/J23;
    strial.yy = 1./3.*Gred*(2.*B->yy-B->xx-B->zz)/J23;
    strial.zz = 1./3.*Gred*(2.*B->zz-B->xx-B->yy)/J23;
    strial.xy = Gred*B->xy/J23;
    if(np==THREED_MPM)
    {   strial.xz = Gred*B->xz/J23;
        strial.yz = Gred*B->yz/J23;
    }
    
    return strial;
}

// Get magnitude of s = sqrt(s.s) when s is a deviatoric stress
double HEIsotropic::GetMagnitudeS(Tensor *st,int np)
{
	double s,t;
	
	switch(np)
    {   case THREED_MPM:
            s = st->xx*st->xx + st->yy*st->yy + st->zz*st->zz;
            t = st->xy*st->xy + st->xz*st->xz + st->yz*st->yz;
            break;
            
		default:
			s = st->xx*st->xx + st->yy*st->yy + st->zz*st->zz;
			t = st->xy*st->xy;
			break;
	}
	return sqrt(s+t+t);
}

// Implementation of the normal to the yield surface
Tensor HEIsotropic::GetNormalTensor(Tensor *strial,double magnitude_strial,int np)
{
    // Trial specific Kirchhoff Stress Tensor
    Tensor n;
    ZeroTensor(&n);
    n.xx = strial->xx/magnitude_strial;
    n.yy = strial->yy/magnitude_strial;
    n.zz = strial->zz/magnitude_strial;
    n.xy = strial->xy/magnitude_strial;
    if(np==THREED_MPM)
    {   n.xz = strial->xz/magnitude_strial;
        n.yz = strial->yz/magnitude_strial;
    }
    //cout << "strial.yy  =    " << strial->yy << "      strial.zz  =    " << strial->zz << "magnitude_strial  =    " << magnitude_strial << endl;
    //cout << "n.xx  =    " << n.xx << "n.xy  =    " << n.xy << endl;
    return n;
}

#pragma mark HEIsotropic::Accessors

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor HEIsotropic::GetStress(Tensor *sp,double pressure)
{   Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}

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
	if(num>0 && num<=J_history+1)
	{	double *p=(double *)historyPtr;
		history=p[num-1];
	}
    return history;
}


