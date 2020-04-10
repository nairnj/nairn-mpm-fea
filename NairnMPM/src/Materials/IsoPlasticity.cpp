/********************************************************************************
    IsoPlasticity.cpp
    nairn-mpm-fea
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Orthotropic.hpp (TranIsotropic.hpp, MaterialBase.hpp)
********************************************************************************/

#include "stdafx.h"
#include "IsoPlasticity.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Materials/HardeningLawBase.hpp"
#include "Materials/LinearHardening.hpp"

#pragma mark IsoPlasticity::Constructors and Destructors

// Constructors
// throws std::bad_alloc
IsoPlasticity::IsoPlasticity(char *matName,int matID) : IsotropicMat(matName,matID)
{
    plasticLaw = new LinearHardening(this);
}

#pragma mark IsoPlasticity::Initialization

// Read material properties
char *IsoPlasticity::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    // look for different plastic law
    if(strcmp(xName,"Hardening")==0)
    {	input = HARDENING_LAW_SELECTION;
        return (char *)this;
    }
   
    // check plastic law
    char *ptr = plasticLaw->InputHardeningProperty(xName,input,gScaling);
    if(ptr != NULL) return ptr;
    
    // otherwise get material properties
    return(IsotropicMat::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *IsoPlasticity::VerifyAndLoadProperties(int np)
{
	// call plastic law that is used
    const char *ptr = plasticLaw->VerifyAndLoadProperties(np);
    if(ptr != NULL) return ptr;
	
	// check in superclass (along with its initialization)
	ptr = IsotropicMat::VerifyAndLoadProperties(np);
	
	// reduced prooperties
	G0red = C66/rho;
	pr.Gred = G0red;
	pr.Kred = C33/rho - 4.*G0red/3.;							// from C33 = lambda + 2G = K + 4G/3
	
	return ptr;
}

// Allows any hardening law
bool IsoPlasticity::AcceptHardeningLaw(HardeningLawBase *pLaw,int lawID)
{   delete plasticLaw;
    plasticLaw = pLaw;
    return true;
}

// return plastic law ID (or 0 if none)
HardeningLawBase *IsoPlasticity::GetPlasticLaw(void) const { return plasticLaw; }

// print mechanical properties to the results
void IsoPlasticity::PrintMechanicalProperties(void) const
{	
    IsotropicMat::PrintMechanicalProperties();
	plasticLaw->PrintYieldProperties();
}

#pragma mark IsoPlasticity::History Data Methods

// The IsoPlasticity has not history data, but its plasticity law might
char *IsoPlasticity::InitHistoryData(char *pchr,MPMBase *mptr)
{
	int num = plasticLaw->HistoryDoublesNeeded();
	if(num==0) return NULL;
	double *p = CreateAndZeroDoubles(pchr,num);
	plasticLaw->InitPlasticHistoryData(p);
    return (char *)p;
}

// Number of history variables - only the plastic law
int IsoPlasticity::NumberOfHistoryDoubles(void) const
{	return plasticLaw->HistoryDoublesNeeded();
}

#pragma mark IsoPlasticity::Methods

/* Take increments in strain and calculate new Particle: strains, rotation strain, plastic strain,
		stresses, strain energy, plastic energy, dissipated energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
	This is general analysis for isotropic plastic material. Subclass must define
		GetYield() and GetKPrime() and optionally can override more. Those methods
		require history dependent variables and rates (e.g. cum. plastic strain (alpint) and
		plastic strain rate (dalpha/delTime) to be set before they are called.
	This material tracks pressure and stores deviatoric stress only in particle stress
		tensor
*/
void IsoPlasticity::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{
	// Common calculations
    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
    //  (note: using unreduced terms in CTE3 and CME3)
    double delV,eres=CTE3*res->dT;
    if(DiffusionTask::HasFluidTransport())
        eres+=CME3*res->dC;
    
    // access properties
    PlasticProperties *p = (PlasticProperties *)properties;

    if(useLargeRotation)
	{	// get incremental strain and rotation
        Matrix3 dR;
        Matrix3 de = LRGetStrainIncrement(CURRENT_CONFIGURATION,mptr,du,&dR,NULL,NULL,NULL);
		
        // 3D
        if(np==THREED_MPM)
        {   delV = de.trace() - 3.*eres;
			PlasticityConstLaw(mptr,de,delTime,np,delV,eres,p,res,&dR,false);
            return;
        }
        
        // 2D
        if(np==PLANE_STRESS_MPM)
            delV = p->psRed*(de(0,0)+de(1,1)-2.*eres);
        else
            delV = de.trace() - 3.*eres;
        PlasticityConstLaw(mptr,de,delTime,np,delV,eres,p,res,&dR,true);
    }
	else
    {	// increment deformation gradient
        HypoIncrementDeformation(mptr,du);
    
        // 3D
        if(np==THREED_MPM)
        {   delV = du.trace() - 3.*eres;
            PlasticityConstLaw(mptr,du,delTime,np,delV,eres,p,res,NULL,false);
            return;
        }
        
        // 2D
        if(np==PLANE_STRESS_MPM)
            delV = p->psRed*(du(0,0)+du(1,1)-2.*eres);
        else
            delV = du.trace() - 3.*eres;
        PlasticityConstLaw(mptr,du,delTime,np,delV,eres,p,res,NULL,true);
    }
}

// To allow some subclasses to support large deformations, the initial calculation for incremental
// deformation gradient (the dvij), volume change (delV)
// handled first by the subclass. This method then finishes the constitutive law
void IsoPlasticity::PlasticityConstLaw(MPMBase *mptr,Matrix3 de,double delTime,int np,double delV,double eres,
									   PlasticProperties *p,ResidualStrains *res,Matrix3 *dR,bool is2D) const
{
#pragma mark ... Code to get trial stress and check for yielding
	// Task 1: Convert input form to effective strains
	//------------------------------------------------
	// engineering shear strains
	double dgxy = de(0,1)+de(1,0),dgxz,dgyz;
	if(!is2D)
	{	dgxz = de(0,2)+de(2,0);
		dgyz = de(1,2)+de(2,1);
	}
	
	// here deij is total strain increment, dexxr is relative strain by subtracting off eres
	double dexxr = de(0,0)-eres;
	double deyyr = de(1,1)-eres;
	double dezzr = de(2,2)-eres;			// In plane strain trial de(2,2)=0, but not in axisymmetric
	
	// Task 2: Update Pressure
	//--------------------------
	double P0 = mptr->GetPressure();
	double dTq0 = 0.,dispEnergy = 0.;
	UpdatePressure(mptr,delV,np,p,res,eres,dTq0,dispEnergy);
	double Pfinal = mptr->GetPressure();
	
	// Task 3: Rotate state n-1 to configuration n
	//----------------------------------------------
	Tensor *eplast=mptr->GetAltStrainTensor();
	Tensor *sp=mptr->GetStressTensor();
	Tensor st0;
	if(useLargeRotation)
	{	// plastic strain (stored on particle)
		Tensor etr = dR->RVoightRT(eplast,false,is2D);
		*eplast = etr;
		
		// prior stress (not stored on particle until later)
		st0 = dR->RVoightRT(sp,true,is2D);
	}
	else
	{	// plastic strain (stored on particle)
		double dwrotxy = de(1,0)-de(0,1),dwrotxz,dwrotyz;
		Tensor ep = *eplast;
		if(is2D)
		{	double dnorm = 0.5*dwrotxy*ep.xy;
			eplast->xx -= dnorm;
			eplast->yy += dnorm;
			eplast->xy += dwrotxy*(ep.xx-ep.yy);
		}
		else
		{	dwrotxz = de(2,0)-de(0,2);
			dwrotyz = de(2,1)-de(1,2);
			double dxy = 0.5*dwrotxy*ep.xy;
			double dxz = 0.5*dwrotxz*ep.xz;
			double dyz = 0.5*dwrotyz*ep.yz;
			eplast->xx += -dxy - dxz;
			eplast->yy += dxy - dyz;
			eplast->zz += dxz + dyz;
			eplast->yz += dwrotyz*(ep.yy-ep.zz) + 0.5*(dwrotxz*ep.xy+dwrotxy*ep.xz);
			eplast->xz += dwrotxz*(ep.xx-ep.zz) + 0.5*(dwrotyz*ep.xy-dwrotxy*ep.yz);
			eplast->xy += dwrotxy*(ep.xx-ep.yy) - 0.5*(dwrotyz*ep.xz+dwrotxz*ep.yz);
		}
		
		// prior stress (not stored on particle until later)
		st0=*sp;
		if(is2D)
		{	double dnorm = dwrotxy*sp->xy;
			st0.xx -= dnorm;
			st0.yy += dnorm;
			st0.xy += 0.5*dwrotxy*(ep.xx-ep.yy);
		}
		else
		{	double dxy = dwrotxy*sp->xy;
			double dxz = dwrotxz*sp->xz;
			double dyz = dwrotyz*sp->yz;
			st0.xx += -dxy - dxz;
			st0.yy += dxy - dyz;
			st0.zz += dxz + dyz;
			st0.yz += 0.5*(dwrotyz*(sp->yy-sp->zz) + dwrotxz*sp->xy + dwrotxy*sp->xz);
			st0.xz += 0.5*(dwrotxz*(sp->xx-sp->zz) + dwrotyz*sp->xy - dwrotxy*sp->yz);
			st0.xy += 0.5*(dwrotxy*(sp->xx-sp->yy) - dwrotyz*sp->xz - dwrotxz*sp->yz);
		}
	}
	
	// Task 3: Get deviatoric stress increment and trial value
	//--------------------------------------------------------
	Tensor dels;
	double thirdDelV = delV/3.;
	dels.xx = 2.*p->Gred*(dexxr-thirdDelV);
	dels.yy = 2.*p->Gred*(deyyr-thirdDelV);
	dels.xy = p->Gred*dgxy;
	if(is2D)
	{	if(np==PLANE_STRESS_MPM)
			dels.zz = Pfinal-P0;
		else
			dels.zz = 2.*p->Gred*(dezzr-thirdDelV);
	}
	else
	{	dels.zz = 2.*p->Gred*(dezzr-thirdDelV);
		dels.yz = p->Gred*dgyz;
		dels.xz = p->Gred*dgxz;
	}
	
	// trial deviatoric stress
	Tensor strial = st0;
	AddTensor(&strial,&dels);
	
	// Task 4: Get plastic potential with strial stress
	//-------------------------------------------------
	// initialize last alpla and zero plastic strain rate
	HardeningAlpha alpha;
	plasticLaw->UpdateTrialAlpha(mptr,np,&alpha,0);
	
	// check potential
	double ftrial = GetPlasticPotential(&strial,mptr,np,delTime,&alpha,p->hardProps);
	
#pragma mark ... Code for elastic update
	// Task 5: If ftrial<0 then elastic update and finish up
	//------------------------------------------------------
	if(ftrial<0.)
	{	// copy trial to particle
		*sp = strial;
		
		// work energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
		// energy units are also Pa mm^3/g, i.e., same as stress units
		if(!is2D)
		{	mptr->AddWorkEnergy(sp->xx*de(0,0) + sp->yy*de(1,1) + sp->zz*de(2,2)
								+ sp->yz*dgyz + sp->xz*dgxz + sp->xy*dgxy);
		}
		else if(np==AXISYMMETRIC_MPM)
		{	mptr->AddWorkEnergy(sp->xx*de(0,0) + sp->yy*de(1,1) + sp->xy*dgxy + sp->zz*de(2,2));
		}
		else if(np==PLANE_STRESS_MPM)
		{	// zz deformation
			de.set(2,2,-p->psLr2G*(dexxr+deyyr) + eres);
			mptr->IncrementDeformationGradientZZ(de(2,2));
			mptr->AddWorkEnergy(sp->xx*de(0,0) + sp->yy*de(1,1) + sp->xy*dgxy);
		}
		else
		{   mptr->AddWorkEnergy(sp->xx*de(0,0) + sp->yy*de(1,1) + sp->xy*dgxy);
		}
		
		// give subclass material chance to update history variables that change in elastic updates
		plasticLaw->ElasticUpdateFinished(mptr,np,delTime,0);
		
		// heat energy from pressure and any pressure method items
		IncrementHeatEnergy(mptr,dTq0,dispEnergy);
		
		return;
	}
	
#pragma mark ... Code to find plastic increment
	// Task 6: Use Radial Return to find plastic strain increment
	//-------------------------------------------------------------
	double lambdak = RRPlasticIncrement(&strial,mptr,np,Pfinal,delTime,&alpha,p);
	
	// Task 7: Find dep = lambdak df/dsigma
	//--------------------------------------
	Tensor dfds;
	if(np==PLANE_STRESS_MPM)
	{	// Note this plane stress codes assumes J2 plasticity, should generalize when write a subclass
		double d1 = (1 + p->psKred*lambdak);
		double d2 = (1.+2.*p->Gred*lambdak);
		double n1 = (strial.xx+strial.yy-2.*Pfinal)/d1;
		double n2 = (-strial.xx+strial.yy)/d2;
		double sxx = (n1-n2)/2.;
		double syy = (n1+n2)/2.;
		double txy = strial.xy/d2;
		
		// find increment in deviatoric stress
		dels.xx = sxx+Pfinal-st0.xx;
		dels.yy = syy+Pfinal-st0.yy;
		dels.xy = txy-st0.xy;
		
		// get final direction
		dfds.xx = (2.*sxx-syy)/3.;
		dfds.yy = (2.*syy-sxx)/3.;
		dfds.zz = -(dfds.xx+dfds.yy);
		dfds.xy = txy;				// tensorial shear strain
		
		// zz deformation
		de.set(2,2,-p->psLr2G*(dexxr+deyyr - lambdak*(dfds.xx+dfds.yy)) + eres + lambdak*dfds.zz);
		mptr->IncrementDeformationGradientZZ(de(2,2));
	}
	else
	{   // get final direction
		GetDfDsigma(alpha.strialmag,&strial,np,&dfds);
	}
	
	// Plastic strain increments on particle
	Tensor dep;
	dep.xx = lambdak*dfds.xx;
	dep.yy = lambdak*dfds.yy;
	dep.zz = lambdak*dfds.zz;
	dep.xy = 2.*lambdak*dfds.xy;     // 2 for engineering plastic shear strain
	if(!is2D)
	{	dep.xz = 2.*lambdak*dfds.xz;
		dep.yz = 2.*lambdak*dfds.yz;
	}
	AddTensor(eplast,&dep);
	
	// Task 8: Increment particle deviatoric stress
	//-----------------------------------------------
	if(!is2D)
	{	sp->xx = strial.xx - 2.*p->Gred*dep.xx;
		sp->yy = strial.yy - 2.*p->Gred*dep.yy;
		sp->zz = strial.zz - 2.*p->Gred*dep.zz;
		sp->yz = strial.yz - p->Gred*dep.yz;
		sp->xz = strial.xz - p->Gred*dep.xz;
		sp->xy = strial.xy - p->Gred*dep.xy;
	}
	else if(np!=PLANE_STRESS_MPM)
	{	sp->xx = strial.xx - 2.*p->Gred*dep.xx;
		sp->yy = strial.yy - 2.*p->Gred*dep.yy;
		sp->zz = strial.zz - 2.*p->Gred*dep.zz;
		sp->xy = strial.xy - p->Gred*dep.xy;
	}
	else
	{	// note that dels is change in code above specific to J2, plane stress plasticity
		sp->xx = st0.xx + dels.xx;
		sp->yy = st0.yy + dels.yy;
		sp->xy = st0.xy + dels.xy;
		sp->zz = Pfinal;          // now equal to Pfinal
	}
	
	// Task 9: Increment all energies
	//-------------------------------
	// Elastic work increment per unit mass
	double workEnergy = sp->xx*de(0,0) + sp->yy*de(1,1) + sp->xy*dgxy;
	if(np==AXISYMMETRIC_MPM)
		workEnergy += sp->zz*de(2,2);
	else if(!is2D)
		workEnergy += sp->zz*de(2,2) + sp->yz*dgyz + sp->xz*dgxz;
	mptr->AddWorkEnergy(workEnergy);
	
	// plastic strain work
	double plastEnergy = sp->xx*dep.xx + sp->yy*dep.yy + sp->zz*dep.zz + sp->xy*dep.xy;
	if(!is2D) plastEnergy += sp->xz*dep.xz + sp->yz*dep.yz;
	
	// and subtract q dalpha to get dissipated energy per unit mass
	double qdalphaTerm = lambdak*SQRT_TWOTHIRDS*plasticLaw->GetYieldIncrement(mptr,np,delTime,&alpha,p->hardProps);
	dispEnergy += plastEnergy - qdalphaTerm;
	
	// The cumulative dissipated energy is tracked in plastic energy
	mptr->AddPlastEnergy(dispEnergy);
	
	// heat energy
	IncrementHeatEnergy(mptr,dTq0,dispEnergy);
	
	// Task 9: Update Internal Variables
	//----------------------------------
	plasticLaw->UpdatePlasticInternal(mptr,np,&alpha,0);
}

#pragma mark IsoPlasticity::Custom Methods

// This method handles the pressure equation of state. Its tasks are
// 1. Calculate the new pressure
// 2. Update particle pressure
// 3. Increment the particle energy
// 4. Call plasticLaw to see if it wants to change the shear modulus
// 5. Optionally change delV (which is passed by reference)
// Notes:
//  delV is incremental volume change on this step.
void IsoPlasticity::UpdatePressure(MPMBase *mptr,double delV,int np,PlasticProperties *p,ResidualStrains *res,double eres,
								   double &dTq0,double &dispEnergy) const
{   // pressure change
    double dP = -p->Kred*delV;
    mptr->IncrementPressure(dP);
    
	// get total dV
	double dVoverV;
	if(np==PLANE_STRESS_MPM)
		dVoverV = delV + 2.*p->psRed*eres;
	else
		dVoverV = delV + 3.*eres;
	
    // work energy is dU = -P dVtot + s.de(total)
	// Here do hydrostatic term
    // Work energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
    double avgP = mptr->GetPressure()-0.5*dP;
    mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-3.*avgP*eres);
	
	// Isoentropic temperature rise = -(K 3 alpha T)/(rho Cv) (dV/V) = - gamma0 T (dV/V)
	dTq0 = -gamma0*mptr->pPreviousTemperature*dVoverV;
}

// return derivatives of the yield function wrt to components of deviatoric stress
// which for isotropic material with f = ||s|| - sqrt(2/3)*sy
// reduces to s/||s|| (written for tensorial plastic strains)
// But do not call for plane stress, which must be found in special case
void IsoPlasticity::GetDfDsigma(double smag,Tensor *st0,int np,Tensor *dfds) const
{
    // s/||s|| = n
    dfds->xx = st0->xx/smag;
    dfds->yy = st0->yy/smag;
    dfds->zz = st0->zz/smag;
    dfds->xy = st0->xy/smag;		// tensorial shear strain
    if(np==THREED_MPM)
    {	dfds->xz = st0->xz/smag;
        dfds->yz = st0->yz/smag;
    }
}

// Get plastic potential. This material is using J2 plasticity
// Subclasses can override for a different potential
double IsoPlasticity::GetPlasticPotential(Tensor *strial,MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	a->strialmag = GetMagnitudeSFromDev(strial,np);
	return a->strialmag - SQRT_TWOTHIRDS*plasticLaw->GetYield(mptr,np,delTime,a,properties);
}

// Radial return to find plastic strain increment
// This class uses J2 plasticity, subclass might override
double IsoPlasticity::RRPlasticIncrement(Tensor *strial,MPMBase *mptr,int np,double Pfinal,double delTime,
												   HardeningAlpha *a,PlasticProperties *p) const
{
	// Find  lambdak for this plastic state and current hardening law
	return plasticLaw->SolveForLambdaBracketed(mptr,np,a->strialmag,strial,p->Gred,p->psKred,Pfinal,delTime,a,p->hardProps,0);
	
}

#pragma mark IsoPlasticity::Accessors

// store plastic strain in alt strain
int IsoPlasticity::AltStrainContains(void) const
{	return ENG_BIOT_PLASTIC_STRAIN;
}

// buffer size for mechanical properties
int IsoPlasticity::SizeOfMechanicalProperties(int &altBufferSize) const
{   altBufferSize = plasticLaw->SizeOfHardeningProps();
    return sizeof(PlasticProperties);
}

// Isotropic material can use read-only initial properties
void *IsoPlasticity::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{
	PlasticProperties *p = (PlasticProperties *)matBuffer;
	*p = pr;
 	p->hardProps = plasticLaw->GetCopyOfHardeningProps(mptr,np,altBuffer,offset);
	double Gratio = plasticLaw->GetShearRatio(mptr,mptr->GetPressure(),1.,p->hardProps,offset);
	p->Gred = G0red*Gratio;
	
	if(np==PLANE_STRESS_MPM)
	{	// these are terms for plane stress calculations only
		p->psRed = 1./(p->Kred/(2.*p->Gred) + 2./3.);					// (1-2nu)/(1-nu) for plane stress
		p->psLr2G = (p->Kred/(2.*p->Gred) - 1./3.)*p->psRed;			// nu/(1-nu) to find ezz
		p->psKred = p->Kred*p->psRed;									// E/(3(1-v)) to find lambda
	}
	
	return p;
}

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor IsoPlasticity::GetStress(Tensor *sp,double pressure,MPMBase *mptr) const
{	return GetStressPandDev(sp,pressure,mptr);
}

// store a new total stress on a particle's stress and pressure variables
void IsoPlasticity::SetStress(Tensor *spnew,MPMBase *mptr) const
{	SetStressPandDev(spnew,mptr);
}

// Increment thickness (zz) stress through deviatoric stress and pressure
void IsoPlasticity::IncrementThicknessStress(double dszz,MPMBase *mptr) const
{	IncrementThicknessStressPandDev(dszz,mptr);
}

// return unique, short name for this material
const char *IsoPlasticity::MaterialType(void) const { return "Isotropic Elastic-Plastic"; }



