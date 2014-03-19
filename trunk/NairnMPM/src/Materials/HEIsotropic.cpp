/********************************************************************************
    HEIsotropic.cpp
    nairn-mpm-fea
    
    Created by John Nairn, Sept 27, 2011.
    Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		HyperElastic.hpp (MaterialBase.hpp)
********************************************************************************/

#include "HEIsotropic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Materials/HardeningLawBase.hpp"
#include "Materials/LinearHardening.hpp"
#include "Custom_Tasks/ConductionTask.hpp"

#pragma mark HEIsotropic::Constructors and Destructors

// Constructors
HEIsotropic::HEIsotropic() {}

// Constructors
HEIsotropic::HEIsotropic(char *matName) : HyperElastic(matName)
{
   	G1 = -1.;			// required
    
	// JAN: hard-code linear hardening - future can set to other laws
	plasticLaw = new LinearHardening(this);
}


#pragma mark HEIsotropic::Initialization

// Read material properties
char *HEIsotropic::InputMat(char *xName,int &input)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"G1")==0 || strcmp(xName,"G")==0)
        return((char *)&G1);
    
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

// Allows any hardening law
bool HEIsotropic::AcceptHardeningLaw(HardeningLawBase *pLaw,int lawID)
{   delete plasticLaw;
    plasticLaw = pLaw;
    return TRUE;
}

// verify settings and maybe some initial calculations
const char *HEIsotropic::VerifyAndLoadProperties(int np)
{
	// call plastic law first
    const char *ptr = plasticLaw->VerifyAndLoadProperties(np);
    if(ptr != NULL) return ptr;
    
    // user must enter G1 and Kbulk
    if(G1<0. || Kbulk < 0. )
		return "HEIsotropic Material needs non-negative G1 and K";
    
	// MU in specific units using initial rho
	// for MPM (units N/m^2 cm^3/g)
	G1sp = G1*1.0e+06/rho;
	
	// must call super class
	return HyperElastic::VerifyAndLoadProperties(np);
}

// plane stress not allowed in viscoelasticity
void HEIsotropic::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
    {	throw CommonException("HEIsotropic materials cannot be used in plane stress MPM yet",
                                "HEIsotropic::ValidateForUse");
    }
	
	//call super class (why can't call super class?)
	HyperElastic::ValidateForUse(np);
}

// print mechanical properties to the results
void HEIsotropic::PrintMechanicalProperties(void) const
{	
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
	
	// call superclass here if it is not Material base
	HyperElastic::PrintMechanicalProperties();
}

// First ones for hardening law. Particle J appended at the end
char *HEIsotropic::InitHistoryData(void)
{	J_history = plasticLaw->HistoryDoublesNeeded();
	double *p = CreateAndZeroDoubles(J_history+1);
    plasticLaw->InitPlasticHistoryData(p);
	p[J_history]=1.;					// J
	return (char *)p;
}

#pragma mark HEIsotropic:Methods

// buffer size for mechanical properties
int HEIsotropic::SizeOfMechanicalProperties(int &altBufferSize) const
{   altBufferSize = plasticLaw->SizeOfHardeningProps();
    return sizeof(HEPlasticProperties);
}

// Get elastic and plastic properties, return null on error
void *HEIsotropic::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
{
	HEPlasticProperties *p = (HEPlasticProperties *)matBuffer;
 	p->hardProps = plasticLaw->GetCopyOfHardeningProps(mptr,np,altBuffer);
	double Gratio = plasticLaw->GetShearRatio(mptr,mptr->GetPressure(),1.,p->hardProps);
	p->Gred = G1sp*Gratio;
	p->Kred = Ksp;
	return p;
}

/* Take increments in strain and calculate new Particle: strains, rotation strain,
        stresses, strain energy,
    dvij are (gradient rates X time increment) to give deformation gradient change
    For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
    This material tracks pressure and stores deviatoric stress only in particle stress
        tensor
*/
void HEIsotropic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	HEPlasticProperties *p = (HEPlasticProperties *)properties;

    // store initial stress
    Tensor *sp = mptr->GetStressTensor();
    Tensor st0 = *sp;
    
    // Compute Elastic Predictor
    // ============================================
    
	// Update total deformation gradient, and calculate trial B
    Tensor Btrial;
	double detdF = IncrementDeformation(mptr,du,&Btrial,np);
    
	// Deformation gradients and Cauchy tensor differ in plane stress and plane strain
    // This code handles plane strain, axisymmetric, and 3D - Plane stress is blocked
	double J = detdF * mptr->GetHistoryDble(J_history);
    
    // save new J
    mptr->SetHistoryDble(J_history,J);  // Stocking J
    
    // J is determinant of F (or sqrt root of determinant of B), Jeff is normalized to residual stretch
	double dresStretch,resStretch = GetResidualStretch(mptr,dresStretch,res);
    double resStretch2 = resStretch*resStretch;
    double Jres = resStretch2*resStretch;
    double Jeff = J/Jres;
	
	// Get hydrostatic stress component in subroutine
	double detdFres = dresStretch*dresStretch*dresStretch;		// increment volumetric stretch
    UpdatePressure(mptr,J,detdF,np,Jeff,delTime,p,res,detdFres);
    
    // Others constants
    double J23 = pow(J, 2./3.)/Jres;
    
    // find Trial Specific Kirchoff stress Tensor (Trial_Tau/rho0)
    Tensor stk = GetTrialDevStressTensor(&Btrial,J23,np,p->Gred);
    
    // Checking for plastic loading
    // ============================================
    
    // Get magnitude of the deviatoric stress tensor
    // ||s|| = sqrt(s.s)
    
	// Set alpint for particle
	HardeningAlpha alpha;
	plasticLaw->UpdateTrialAlpha(mptr,np,&alpha);
	
	// Trial stress state
    double magnitude_strial = GetMagnitudeS(&stk,np);
    double gyld = plasticLaw->GetYield(mptr,np,delTime,&alpha,p->hardProps);
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
        *pB = Btrial;
        
        // Get specifique stress i.e. (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
		// The deviatoric stress was calculated as (Kirchoff Stress)/rho0, so just save it now
        *sp = stk;
        
        // work energy per unit mass (U/(rho0 V0)) and we are using
        // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
        mptr->AddWorkEnergy(0.5*((st0.xx+sp->xx)*du(0,0)
                                   + (st0.yy+sp->yy)*du(1,1)
                                   + (st0.zz+sp->zz)*du(2,2)
                                   + (st0.xy+sp->xy)*(du(1,0)+du(0,1))));
        if(np==THREED_MPM)
        {   mptr->AddWorkEnergy(0.5*((st0.yz+sp->yz)*(du(2,1)+du(1,2))
                                       + (st0.xz+sp->xz)*(du(2,0)+du(0,2))));
        }
		
		// residual energy or sigma.deres - it is zero here for isotropic material
		// because deviatoric stress is traceless and deres has zero shear terms
		// residual energy due to pressure was added in the pressure update
		
		// heat energy is Cv(dT-dTq0) - dPhi, but dPhi is zero here
        // and Cv(dT-dTq0) was done in pressure update
        
        return;
    }
    
    // Plastic behavior - Return Mapping algorithm 
    //=====================================================
    // if plastic
    
	// JAN: Use hardening law method (which can now use other laws too)
    double Ie1bar = (Btrial.xx+Btrial.yy+Btrial.zz)/(3.*J23);
    double MUbar = Jres*p->Gred*Ie1bar;
    
    // Find  lambda for this plastic state
	double dlambda = plasticLaw->SolveForLambdaBracketed(mptr,np,magnitude_strial,&stk,MUbar,1.,1.,delTime,&alpha,p->hardProps);
    
    // update deviatoric stress
	Tensor nk = GetNormalTensor(&stk,magnitude_strial,np);
    //cout << "nk.xx  =    " << nk.xx << "nk.xy  =    " << nk.xy << endl;
	double twoMuLam = 2.*MUbar*dlambda;
    sp->xx = stk.xx - twoMuLam*nk.xx;
    sp->yy = stk.yy - twoMuLam*nk.yy;
    sp->zz = stk.zz - twoMuLam*nk.zz;
    sp->xy = stk.xy - twoMuLam*nk.xy;
    if(np == THREED_MPM)
    {   sp->xz = stk.xz - twoMuLam*nk.xz;
        sp->yz = stk.yz - twoMuLam*nk.yz;
    }
    
    // save on particle
	double twoThirdsLamI1bar = 2.*dlambda*Ie1bar;
    pB->xx = Btrial.xx - twoThirdsLamI1bar*nk.xx;
    pB->yy = Btrial.yy - twoThirdsLamI1bar*nk.yy;
    pB->zz = Btrial.zz - twoThirdsLamI1bar*nk.zz;
	pB->xy = Btrial.xy - twoThirdsLamI1bar*nk.xy;
    if(np == THREED_MPM)
    {   pB->xz = Btrial.xz - twoThirdsLamI1bar*nk.xz;
        pB->yz = Btrial.yz - twoThirdsLamI1bar*nk.yz;
    }
    /* Old method collecting B from stresses
     pB->xx = (sp->xx/p->Gred+Ie1bar)*J23;
     pB->yy = (sp->yy/p->Gred+Ie1bar)*J23;
     pB->zz = (sp->zz/p->Gred+Ie1bar)*J23;
     pB->xy = sp->xy*J23/p->Gred;
     if(np == THREED_MPM)
     {   pB->xz = sp->xz*J23/p->Gred;
     pB->yz = sp->yz*J23/p->Gred;
     }
     */
    
    // strain energy per unit mass (U/(rho0 V0)) and we are using
    double workEnergy = 0.5*((st0.xx+sp->xx)*du(0,0)
                               + (st0.yy+sp->yy)*du(1,1)
                               + (st0.zz+sp->zz)*du(2,2)
                               + (st0.xy+sp->xy)*(du(1,0)+du(0,1)));
    if(np==THREED_MPM)
    {   workEnergy += 0.5*((st0.yz+sp->yz)*(du(2,1)+du(1,2))
                                   + (st0.xz+sp->xz)*(du(2,0)+du(0,2)));
    }
    
    // Plastic work increment per unit mass (dw/(rho0 V0)) (uJ/g)
    double dispEnergy = dlambda*(sp->xx*nk.xx + sp->yy*nk.yy + sp->zz*nk.zz + 2.*sp->xy*nk.xy);
    if(np==THREED_MPM)  dispEnergy += 2.*dlambda*(sp->xz*nk.xz + sp->yz*nk.yz);
	
    // total work
    mptr->AddWorkEnergy(dispEnergy + workEnergy);
    
    // residual energy or sigma.deres - it is zero here for isotropic material
    // because deviatoric stress is traceless and deres has zero shear terms
    // residual energy due to pressure was added in the pressure update
    
    // disispated energy per unit mass (dPhi/(rho0 V0)) (uJ/g)
	// Subtract q.dalpha
    dispEnergy -= dlambda*SQRT_TWOTHIRDS*plasticLaw->GetYieldIncrement(mptr,np,delTime,&alpha,p->hardProps);
    
    // heat energy is Cv(dT-dTq0) - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adibatic heating occurs)
    // The Cv(dT-dTq0) was done in update pressure
    IncrementHeatEnergy(mptr,0.,0.,dispEnergy);
    
	// The cumulative dissipated energy is tracked in plastic energy
    mptr->AddPlastEnergy(dispEnergy);
    
	// update internal variables in the plastic law
	plasticLaw->UpdatePlasticInternal(mptr,np,&alpha);
}

#pragma mark HEIsotropic::Custom Methods

// To allow better subclassing it is better to separate out calculations
//  of dilation energy. This version updates pressure (i.e. dilational
//  contribution to normal stress) and adds incremental energy to strain energy
// J = V(T,c)/V0(T,c), Jeff = V(T,c)/V0(Tref,cref), Jres = V0(T,c)/V0(Tref,cref) = J/Jeff
// Jn+1 = (detdF/detdFres) Jn, Jresn+1 = detdFres Jresn, Jtot = detdF Jtotn
// detdFres = (1+dres)^3 (approximately)
// Here Tref and cref are starting conditions and T and c are current temperature and moisture
void HEIsotropic::UpdatePressure(MPMBase *mptr,double J,double detdF,int np,double Jeff,
								 double delTime,HEPlasticProperties *p,ResidualStrains *res,double detdFres) const
{
	double Kterm = J*GetVolumetricTerms(Jeff,p->Kred);       // times J to get Kirchoff stress
    double P0 = mptr->GetPressure();
    
    // artifical viscosity
	// delV is total incremental volumetric strain = total Delta(V)/V, here it is (Vn+1-Vn)/Vn+1
	double delV = 1. - 1./detdF;
    double QAVred = 0.,AVEnergy=0.;
    if(delV<0. && artificialViscosity)
	{	double c = sqrt(p->Kred/1000.);           // m/sec
        QAVred = GetArtificalViscosity(delV/delTime,c);
        if(ConductionTask::AVHeating) AVEnergy = fabs(QAVred*delV);
    }
    double Pfinal = -Kterm + QAVred;
    mptr->SetPressure(Pfinal);
    
    // work energy is dU = -P dV + s.de(total)
	// Here do hydrostatic term
    // Internal energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
	// Also get residual energy from -P (Vresn+1-Vresn)/Vresn+1 = -P delVres
    double avgP = 0.5*(P0+Pfinal);
	double delVres = 1. - 1./detdFres;
    mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-avgP*delVres);
	
    // heat energy is Cv dT  - dPhi
	// Here do Cv dT term and dPhi is done later
    IncrementHeatEnergy(mptr,res->dT,0.,AVEnergy);
}

// get trial deviatoric stress tensor based on trial B
// Note that input J23 is actually J^(2/3)/Jres
Tensor HEIsotropic::GetTrialDevStressTensor(Tensor *B,double J23,int np,double Gred) const
{
    // Trial Kirchhoff Stress Tensor / rho0 = J Cauchy/rho0 = Cauchy/rho
    Tensor strial;
	double J23normal = 3.*J23;
    strial.xx = Gred*(2.*B->xx-B->yy-B->zz)/J23normal;
    strial.yy = Gred*(2.*B->yy-B->xx-B->zz)/J23normal;
    //strial.zz = Gred*(2.*B->zz-B->xx-B->yy)/J23normal;
	strial.zz = -(strial.xx+strial.yy);			// slightly faster and making use of traceless
    strial.xy = Gred*B->xy/J23;
    if(np==THREED_MPM)
    {   strial.xz = Gred*B->xz/J23;
        strial.yz = Gred*B->yz/J23;
    }
	else
	{	strial.xz = 0.;
		strial.yz = 0.;
	}
    
    return strial;
}

// Get magnitude of s = sqrt(s.s) when s is a deviatoric stress
double HEIsotropic::GetMagnitudeS(Tensor *st,int np) const
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
Tensor HEIsotropic::GetNormalTensor(Tensor *strial,double magnitude_strial,int np) const
{
    // Trial specific Kirchhoff Stress Tensor
    Tensor n;
    n.xx = strial->xx/magnitude_strial;
    n.yy = strial->yy/magnitude_strial;
    n.zz = strial->zz/magnitude_strial;
    n.xy = strial->xy/magnitude_strial;
    if(np==THREED_MPM)
    {   n.xz = strial->xz/magnitude_strial;
        n.yz = strial->yz/magnitude_strial;
    }
	else
	{	n.xz = 0.;
		n.yz = 0.;
	}
    //cout << "strial.yy  =    " << strial->yy << "      strial.zz  =    " << strial->zz << "magnitude_strial  =    " << magnitude_strial << endl;
    //cout << "n.xx  =    " << n.xx << "n.xy  =    " << n.xy << endl;
    return n;
}

#pragma mark HEIsotropic::Accessors

// convert J to K using isotropic method
Vector HEIsotropic::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{	double nuLS = (3.*Kbulk-2.*G1)/(6.*Kbulk+2.*G1);
	return IsotropicJToK(d,C,J0,np,nuLS,G1);
}

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor HEIsotropic::GetStress(Tensor *sp,double pressure) const
{   Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}

// Return the material tag
int HEIsotropic::MaterialTag(void) const { return HEISOTROPIC; }

// return unique, short name for this material
const char *HEIsotropic::MaterialType(void) const { return "Hyperelastic Isotropic"; }

// calculate wave speed in mm/sec
double HEIsotropic::WaveSpeed(bool threeD,MPMBase *mptr) const
{
    return sqrt(1.e9*(Kbulk+4.*G1/3.)/rho);
}

// this material has two history variables
double HEIsotropic::GetHistory(int num,char *historyPtr) const
{
    double history=0.;
	if(num>0 && num<=J_history+1)
	{	double *p=(double *)historyPtr;
		history=p[num-1];
	}
    return history;
}

// if a subclass material supports artificial viscosity, override this and return TRUE
bool HEIsotropic::SupportsArtificialViscosity(void) const { return TRUE; }

