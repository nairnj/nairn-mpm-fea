/********************************************************************************
    IdealGas.cpp
    nairn-mpm-fea

    Created by Edward Le and Peter Mackenzie on May 26, 2012.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/IdealGas.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "System/UnitsController.hpp"

#pragma mark IdealGas::Constructors and Destructors

// Constructor 
IdealGas::IdealGas(char *matName,int matID) : HyperElastic(matName,matID)
{
	P0   = -1.;			// required initial pressure
	rho  = -1.;			// required density (override default of 1)
	T0   = -1.;			// required initial temperature in Kelvin
}

#pragma mark IdealGas::Initialization

// Read material properties
char *IdealGas::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"P0")==0)
        return UnitsController::ScaledPtr((char *)&P0,gScaling,1.e6);
    
    else if(strcmp(xName,"T0")==0)
        return((char *)&T0);
	
    return(HyperElastic::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *IdealGas::VerifyAndLoadProperties(int np)
{
	// make sure all were set
    if(P0 <= 0. || rho <= 0.0 || T0 <= 0.0 )
		return "Ideal gas material model needs positive parameters P0, rho, and T0";
	
	// Find ideal gas has Cv heat capacity in F-L/(mass-K)
	// For monotonic Ideal Gas, Cv = 1.5R for diatomic gas is 2.5R
	// If set to >1 is diatomic, otherwise monotonic (which is for not set too)
	if(heatCapacity>UnitsController::Scaling(1.e6))
	{	heatCapacity = 2.5*P0/(T0*rho);
		gammaAdiabatic = 7./5.;
	}
	else
	{	heatCapacity = 1.5*P0/(T0*rho);
		gammaAdiabatic = 5./3.;
	}
	CpMinusCv= P0/(T0*rho);
	
	// P0 in specific units (F/L^2 L^3/mass)
	P0sp=P0/rho;
	
    // call super class
    return HyperElastic::VerifyAndLoadProperties(np);
}

// if analysis not allowed, throw an exception
// throws CommonException()
void IdealGas::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
	{	throw CommonException("IdealGas material cannot do 2D plane stress analysis",
							  "IdealGas::ValidateForUse");
	}
	
	if(thermal.reference<=0)
	{	throw CommonException("IdealGas material requires the simulation to set the stress free temperature in degrees K",
							  "IdealGas::ValidateForUse");
	}
	
	// call super class (why can't call super class?)
	HyperElastic::ValidateForUse(np);
}

// print mechanical properties output window
void IdealGas::PrintMechanicalProperties(void) const
{
	PrintProperty("P0",P0*UnitsController::Scaling(1.e-6),"");
	PrintProperty("T0",T0,"");
	cout << endl;
}

// If needed, a material can initialize particle state
// For example, ideal gas initializes to base line pressure
void IdealGas::SetInitialParticleState(MPMBase *mptr,int np,int offset) const
{
    // The initial state has particle mass mp = Vref * rho0 at T = thermal.reference
    // Imagine heating from T0 to T holding volume constant and find Kirchoff stress / rho0
    double Psp = -P0sp * (thermal.reference/T0);
    
    // set the particle pressure
	Tensor *sp=mptr->GetStressTensor();
	sp->xx = Psp;
	sp->yy = Psp;
	sp->zz = Psp;
    
    // Initial particle strains are zero or define Jref=1 for undeformed particle
    
    // call super class for Cauchy Green strain
    HyperElastic::SetInitialParticleState(mptr,np,offset);
}

#pragma mark Mooney::History Data Methods

// return number of bytes needed for history data
int IdealGas::SizeOfHistoryData(void) const { return sizeof(double); }

// Store J, which is calculated incrementally, and available for archiving
// initialize to 1
char *IdealGas::InitHistoryData(char *pchr,MPMBase *mptr)
{	double *p = CreateAndZeroDoubles(pchr,2);
	p[0] = 1.;
	return (char *)p;
}

// Number of history variables
int IdealGas::NumberOfHistoryDoubles(void) const { return 1; }

#pragma mark IdealGas::Methods

// In unit nJ/(g-K)
double IdealGas::GetCpMinusCv(MPMBase *) const { return CpMinusCv; }

/* Take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    du are (gradient rates X time increment) to give deformation gradient change
*/
void IdealGas::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{
    // Update strains and rotations and Left Cauchy strain
    // get determinent of incremental deformation gradient
    double detf = IncrementDeformation(mptr,du,NULL,np);
	double delV = 1. - 1./detf;
	
	// Get new J and save result on the particle (not really needed, by tracked for phase transitions)
	double Jprev = mptr->GetHistoryDble(J_History,historyOffset);
	double J = detf * Jprev;
	mptr->SetHistoryDble(J_History,J,historyOffset);
	
	// store pressure strain as elastic B (only needed for phase transition materials)
	Tensor *pB = mptr->GetAltStrainTensor() ;
	if(np==THREED_MPM || np==AXISYMMETRIC_MPM)
	{	double J23 = pow(J,2./3.);
		pB->xx = J23;
		pB->yy = J23;
		pB->zz = J23;
	}
	else
	{	pB->xx = J;
		pB->yy = J;
	}
	
    // update stress (which is -P)
	Tensor *sp=mptr->GetStressTensor();
    double mPnsp = sp->xx;

//#define INCREMENTAL_PSP
	// compute specific pressure as p/rho = P0sp * (T/T0)
#ifdef INCREMENTAL_PSP
	// incrementally p(n+1)/rho = P0sp * (T(n+1)/T0)
	//							= P0sp * (Tn/T0) * (T(n+1)/Tn)
	//							= (pn/rho) * (T(n+1)/Tn) = pn * (1+dT/Tn)
	// (note: pPreviousTemperature is set in particle update, which will differ
	//			in update before and after time step)
	double mPsp = mPnsp*(1.+res->dT/mptr->pPreviousTemperature);
#else
	// Calculate based on current temperature
	//   (note: CpminusCv = P0/(T0 rho))
	double mPsp = -CpMinusCv*mptr->pPreviousTemperature;
#endif

	// artificial viscosity
	double QAVred = 0.;
	double AVEnergy = 0.;
	if(delV<0. && artificialViscosity)
	{	QAVred = GetArtificalViscosity(delV/delTime,sqrt(fabs(gammaAdiabatic*mPnsp)),mptr);
		AVEnergy += fabs(QAVred*delV);
		mPsp -= QAVred;
	}
	
	// store in stress (which is minus the pressure)
	sp->xx = mPsp;
	sp->yy = mPsp;
	sp->zz = mPsp;
	
	// find the -P dV energy per unit mass dW/(rho0 V0) (nJ/g) as follows
    // dW/(rho0 V0) = - 0.5 * (pn+p(n+1))/rho0 * (V(n+1)-Vn)/V0, which simplifies to
    double dW = 0.5*(mPnsp*detf + mPsp)*delV;
    
    // this energy is tracked in work energy and no residual energy is tracked
    mptr->AddWorkEnergy(dW);
    
    // heat energy (treated as reversible)
	double dTq0 = dW/GetHeatCapacity(mptr);
    IncrementHeatEnergy(mptr,dTq0-AVEnergy,AVEnergy);
}

#pragma mark IdealGas::Accessors

// return material type
const char *IdealGas::MaterialType(void) const { return "Ideal Gas (Hyperelastic)"; }

// Calculate current wave speed in L/sec.
double IdealGas::WaveSpeed(bool threeD,MPMBase *mptr) const
{	double Pspcurr;
	if(mptr!=NULL)
	{	Tensor *sp=mptr->GetStressTensor();
		Pspcurr = fmax(-sp->xx,0.);
	}
	else
		Pspcurr = P0sp;
	return sqrt(gammaAdiabatic*Pspcurr);
}

// not supported because moisture swell may change MGEOS law (doe not work in Diffusion either)
bool IdealGas::SupportsDiffusion(void) const { return false; }

// if a subclass material supports artificial viscosity, override this and return TRUE
bool IdealGas::SupportsArtificialViscosity(void) const { return true; }



