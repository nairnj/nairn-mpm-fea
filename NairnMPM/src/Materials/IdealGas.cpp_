/********************************************************************************
    IdealGas.cpp
    nairn-mpm-fea

    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Materials/IdealGas.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "System/UnitsController.hpp"

// does ideal gas (not van der Waal yets) on pressure scale where 0 is 1 atm
#define REL_PRESSURE

#pragma mark IdealGas::Constructors and Destructors

// Constructor 
IdealGas::IdealGas(char *matName,int matID) : HyperElastic(matName,matID)
{
	P0   = -1.;			// required initial pressure
	rho  = -1.;			// required density (override default of 1)
	T0   = -1.;			// required initial temperature in Kelvin
    vdwa = -1.;
    vdwb = -1.;
    vanderWaalsGas = false;
    mu0 = -1.;          // reference viscoty at T0
    cSuth = -1.;        // Sutherlands constant
    realPressure = false;
}

#pragma mark IdealGas::Initialization

// Read material properties
char *IdealGas::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"P0")==0)
    {   realPressure = true;
        return UnitsController::ScaledPtr((char *)&P0,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"T0")==0)
        return((char *)&T0);
	
    else if(strcmp(xName,"vdwa")==0)
    {   // units of L^6*pressure/mol^2
        // Legacy input as mm^6*MPa/mol^2, scale by 1e6 to get mm^6 Pa/mol^2
        vanderWaalsGas = true;
        return UnitsController::ScaledPtr((char *)&vdwa,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"vdwb")==0)
    {   // units L^3/mol -  no scaling in Legcay units (expecting mm^3/mol)
        vanderWaalsGas = true;
        return((char *)&vdwb);
    }
    
    else if(strcmp(xName,"mu0")==0)
    {   // viscosity
        input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&mu0,gScaling,1.e-3);
    }
    
    else if(strcmp(xName,"Smu")==0)
    {   // Sutherlands constant (effectively degrees)
        input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&cSuth,gScaling,1.e-3);
    }
    
    return(HyperElastic::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *IdealGas::VerifyAndLoadProperties(int np)
{
    // make sure all were set
    if(realPressure)
    {   if(P0 <= 0. || rho <= 0.0 || T0 <= 0.0 )
            return "Gas material model needs positive parameters P0, rho, and T0";
    }
    else
    {   if(rho <= 0.0 || T0 <= 0.0)
            return "Gas material model needs positive parameters rho and T0";
        P0 = UnitsController::GetOneAtmosphere();
    }
    
    // is visocity used?
    if(mu0>0.)
    {   if(cSuth<0.)
            return "Modeling gas viscosity needs Sutherland's constant (Smu) >= 0";
        twoMu0sp = 2.*mu0/rho;
    }
    
    // Find ideal gas (and non-ideal) Cv heat capacity in F-L/(mass-K)
    // For monotonic Ideal Gas, Cv = 1.5R for diatomic gas is 2.5R
    // If set to >1 is diatomic, otherwise monotonic (which is for when not set too)
    if(heatCapacity>UnitsController::Scaling(1.e6))
        gammaAdiabatic = 7./5.;
    else
        gammaAdiabatic = 5./3.;
    heatCapacity = P0/((gammaAdiabatic-1.)*T0*rho);
    
    // Molar volume from ideal law (V0/n) = R T0/P0 in reference conditions
    double R = UnitsController::GetGasConstant(1.);
    double V0 = R*T0/P0;

    // P0 in specific units (F/L^2 L^3/mass)
    P0sp = P0/rho;
    
    if(vanderWaalsGas)
    {   if(vdwa<0. || vdwb<0.)
            return "Van der Waals gas material model needs positive parameters vdwa and vdwb";
        
        // solve cubic equation by Newton's formula
        double abar = vdwa/(V0*V0),bbar=vdwb/V0;
        double xn = 1.;
        for(int i=0;i<10;i++)
        {   // get f(x0) and f'(x0)
            double axbar = abar/(xn*xn);
            double fxn = (P0+axbar)*(xn-bbar)-P0;
            double fpxn = P0 - axbar + 2.*axbar*bbar/xn;

            // increment
            double dx = fxn/fpxn;
            xn -= dx;
            //cout << "i=" << i << ", xn=" << xn << ", fx0=" << fxn << ", dx=" << fxn/fpxn << endl;
            if(fabs(dx)<1.e-12) break;
        }
        
        V0n = V0*xn;
        aprime = vdwa/(V0n*V0n);
        bprime = vdwb/V0n;
        
        // P0', and aprime in specific units
        P0prime = (P0+aprime)*(1-bprime)/rho;
        aprime /= rho;
    }
    
    else
    {   // use ideal molar volume
        V0n = V0;
        
        // Cp-Cv is constant and = (gamma-1)*Cv
        CpMinusCv = P0/(T0*rho);
    }
	
    // call super class, but skip Hyperelastic that calculates things not used by ideal gas
    return MaterialBase::VerifyAndLoadProperties(np);
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
    char units[100];
    
    // ideal gas (V0/n is correct in van dewr Waals gas
    if(realPressure)
        PrintProperty("P0",P0*UnitsController::Scaling(1.e-6),"");
	PrintProperty("T0",T0,"");
    strcpy(units,UnitsController::Label(CULENGTH_UNITS));
    strcat(units,"^3/mol");
    PrintProperty("V0/n",V0n,units);
    strcpy(units,UnitsController::Label(CUMASS_UNITS));
    strcat(units,"/mol");
    PrintProperty("MW",rho*V0n,units);
    cout << endl;

    if(vanderWaalsGas)
    {   // a
        strcpy(units,UnitsController::Label(CULENGTH_UNITS));
        strcat(units,"^6 ");
        strcat(units,UnitsController::Label(PRESSURE_UNITS));
        strcat(units,"/mol^2");
        PrintProperty("a",vdwa*UnitsController::Scaling(1.e-6),units);
        
        // b
        strcpy(units,UnitsController::Label(CULENGTH_UNITS));
        strcat(units,"^3/mol");
        PrintProperty("b",vdwb,units);
        cout << endl;
        
        // a'
        PrintProperty("a'",aprime*rho*UnitsController::Scaling(1.e-6),UnitsController::Label(PRESSURE_UNITS));
        
        // b'
        PrintProperty("b'",bprime,"");
        cout << endl;
        
        // Tc
        cout << "Valid for T>Tc = " << (8.*vdwa)/(27.*vdwb*UnitsController::GetGasConstant(1.)) << endl;
        
        // Pc
        PrintProperty("Pc",vdwa*UnitsController::Scaling(1.e-6)/(27.*vdwb*vdwb), UnitsController::Label(PRESSURE_UNITS));
        
        // Jc = (Vc/n)/(V0/n)
        PrintProperty("Jc",3.*vdwb/V0n,"");
        cout << endl;
    }
    
    if(mu0>0.)
    {   PrintProperty("mu0",mu0*UnitsController::Scaling(1.e3),UnitsController::Label(VISCOSITY_UNITS));
        PrintProperty("Smu",cSuth,"");
        cout << endl;
    }

}

// If needed, a material can initialize particle state
// For example, ideal gas initializes to base line pressure
void IdealGas::SetInitialParticleState(MPMBase *mptr,int np,int offset) const
{
    // The initial state has particle mass mp = Vref * rho0 at T = thermal.references
    // Imagine heating from T0 to T holding volume constant and find Kirchoff stress / rho0
    // Not that if initial particle temperature differs from reference, its pressure
    //      will jump on the first time step
    double Psp;
    if(realPressure)
    {   Psp = vanderWaalsGas ? (P0sp+aprime)*(thermal.reference/T0)-aprime :
                                        P0sp*(thermal.reference/T0);
    }
    else
    {   Psp = vanderWaalsGas ? (P0sp+aprime)*((thermal.reference/T0) - 1.) :
                                        P0sp*((thermal.reference/T0) - 1.);
    }
    
    // set the particle pressure
    mptr->SetPressure(Psp);
    
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
{	double *p = CreateAndZeroDoubles(pchr,NumberOfHistoryDoubles());
	p[0] = 1.;
	return (char *)p;
}

// reset history data
void IdealGas::ResetHistoryData(char *pchr,MPMBase *mptr)
{	double *p = (double *)pchr;
	p[0] = 1.;
}

// Number of history variables
int IdealGas::NumberOfHistoryDoubles(void) const { return 1; }

#pragma mark IdealGas::Methods

// Get different - currently never called with mptr!=NULL, if was, need offset to handles phase transitions
double IdealGas::GetCpMinusCv(MPMBase *mptr) const
{
    if(!vanderWaalsGas)
    {   // precalculated (gamma-1)*Cv
        return CpMinusCv;
    }
    
    // for van der Waals gas
    double J = mptr!=NULL ? mptr->GetHistoryDble(J_History,0) : 1. ;
    double T = mptr!=NULL ? mptr->pPreviousTemperature : thermal.reference ;
    double oneb = 1.-bprime/J;
    double arg = 2.*aprime*oneb*oneb*T0/(P0prime*T*J);
    return (gammaAdiabatic-1.)*heatCapacity/(1.-arg);
}

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
	
    // save previous pressure
    double Pnsp = mptr->GetPressure();

    // get new pressure (a Kirchoff pressure)
    double Psp;
    if(vanderWaalsGas)
    {   if(realPressure)
            Psp = P0prime*(mptr->pPreviousTemperature)/(T0*(1.-bprime/J)) - aprime/J;
        else
            Psp = P0prime*(mptr->pPreviousTemperature)/(T0*(1.-bprime/J)) - (P0sp*J+aprime/J);
    }
    else
    {   if(realPressure)
        {   // (note: CpminusCv = P0/(T0 rho))
            Psp = CpMinusCv*mptr->pPreviousTemperature;
        }
        else
            Psp = P0sp*(mptr->pPreviousTemperature/T0 - J);
    }

	// artificial viscosity
	double QAVred = 0.;
	double AVEnergy = 0.;
	if(delV<0. && artificialViscosity)
	{	QAVred = GetArtificialViscosity(delV/delTime,sqrt(fabs(gammaAdiabatic*Pnsp)),mptr);
		AVEnergy += fabs(QAVred*delV);
		Psp += QAVred;
	}
	
	// store in pressure
    mptr->SetPressure(Psp);
	
	// find the -P dV energy per unit mass dW/(rho0 V0) (nJ/g) as follows
    // dW/(rho0 V0) = - 0.5 * (pn+p(n+1))/rho0 * (V(n+1)-Vn)/V0, which simplifies to
    double dW;
    if(realPressure)
        dW = -0.5*(Pnsp*detf + Psp)*delV;
    else
        dW = -0.5*(Pnsp*detf + Psp + 2.*P0sp*J)*delV;
    
    // this energy is tracked in work energy and no residual energy is tracked
    mptr->AddWorkEnergy(dW);
    
    // heat energy (reversible and irreversible)
    double dV = 1.-1./detf;
	double dTq0 = vanderWaalsGas ? (dW-aprime*dV/J-AVEnergy)/GetHeatCapacity(mptr) :
                                    (dW-AVEnergy)/GetHeatCapacity(mptr);

    // gas viscosity
    if(mu0>0.)
    {   // viscosity term = 2 eta (0.5(grad v) + 0.5*(grad V)^T - (1/3) tr(grad v) I) = 2 eta dev(grad v)
        // (i.e., deviatoric part of the symmetric strain tensor, 2 is for conversion to engineering shear strain)
        // simple shear rate = |2 dev(grad v)|
        Matrix3 shear;
        double c[3][3];
        c[0][0] = (2.*du(0,0)-du(1,1)-du(2,2))/3.;
        c[1][1] = (2.*du(1,1)-du(0,0)-du(2,2))/3.;
        c[2][2] = (2.*du(2,2)-du(0,0)-du(1,1))/3.;
        c[0][1] = 0.5*(du(0,1)+du(1,0));
        c[1][0] = c[0][1];
        if(np==THREED_MPM)
        {   c[0][2] = 0.5*(du(0,2)+du(2,0));
            c[2][0] = c[0][2];
            c[1][2] = 0.5*(du(1,2)+du(2,1));
            c[2][1] = c[1][2];
            shear.set(c);
        }
        else
            shear.set(c[0][0],c[0][1],c[1][0],c[1][1],c[2][2]);
        
        // Get 2 * effective viscosity
        double twoetasp = twoMu0sp*((T0+cSuth)/(mptr->pPreviousTemperature+cSuth))
                    *pow(mptr->pPreviousTemperature/T0,1.5);
        
        // Get Kirchoff shear stress (over rho0)
        shear.Scale(J*twoetasp/delTime);
        
        // update deviatoric stress
        Tensor *sp=mptr->GetStressTensor();
        sp->xx = shear(0,0);
        sp->yy = shear(1,1);
        sp->zz = shear(2,2);
        sp->xy = shear(0,1);
        if(np==THREED_MPM)
        {   sp->xz = shear(0,2);
            sp->yz = shear(1,2);
        }
        
        // shear work per unit mass = tau.du = tau.tau*delTime/twoetasp
        double shearWork = sp->xx*sp->xx + sp->yy*sp->yy + sp->zz*sp->zz + 2.*sp->xy*sp->xy;
        if(np==THREED_MPM) shearWork += 2.*(sp->xz*sp->xz + sp->yz*sp->yz);
        shearWork *= delTime/twoetasp;
        mptr->AddWorkEnergy(shearWork);
        
        // add to AVEnergy
        AVEnergy += shearWork;
    }
    
    // final heat energy
    IncrementHeatEnergy(mptr,dTq0,AVEnergy);
}

#pragma mark IdealGas::Accessors

// Copy stress to a read-only tensor variable after combining deviatoric and pressure
Tensor IdealGas::GetStress(Tensor *sp,double pressure,MPMBase *mptr) const
{    return GetStressPandDev(sp,pressure,mptr);
}

// store a new total stress on a particle's stress and pressure variables
void IdealGas::SetStress(Tensor *spnew,MPMBase *mptr) const
{    SetStressPandDev(spnew,mptr);
}

// Increment thickness (zz) stress through deviatoric stress and pressure
void IdealGas::IncrementThicknessStress(double dszz,MPMBase *mptr) const
{    IncrementThicknessStressPandDev(dszz,mptr);
}

// return material type
const char *IdealGas::MaterialType(void) const
{   return vanderWaalsGas ? "Van der Waals Gas (Hyperelastic)" : "Ideal Gas (Hyperelastic)";
}

// Reference condition wave speed in L/sec.
double IdealGas::WaveSpeed(bool threeD,MPMBase *mptr) const
{   if(vanderWaalsGas)
    {   double Pspinit = (P0sp+aprime)*(thermal.reference/(T0*(1-bprime)));
        return sqrt(gammaAdiabatic*Pspinit - 2.*aprime);
    }
    else
    {   double Pspinit = P0sp * (thermal.reference/T0);
        return sqrt(gammaAdiabatic*Pspinit);
    }
}

// Calculate current wave speed in L/sec. (need offset if in phase transition)
double IdealGas::CurrentWaveSpeed(bool threeD,MPMBase *mptr,int offset) const
{   double Pspcurr = mptr->GetPressure();
    if(realPressure)
    {   if(vanderWaalsGas)
        {   double J = mptr->GetHistoryDble(J_History,offset);
            return sqrt(gammaAdiabatic*(Pspcurr+aprime/J)/(1.-bprime/J) - 2*aprime/J);
        }
        else
            return sqrt(gammaAdiabatic*Pspcurr);
    }
    
    // gauge pressure
    if(vanderWaalsGas)
    {   double J = mptr->GetHistoryDble(J_History,offset);
        return sqrt(gammaAdiabatic*(Pspcurr+J*P0sp+aprime/J)/(1.-bprime/J) - 2*aprime/J);
    }
    else
    {   double J = mptr->GetHistoryDble(J_History,offset);
        return sqrt(gammaAdiabatic*(Pspcurr+J*P0sp));
    }
}

// not supported because moisture swell may change MGEOS law (doe not work in Diffusion either)
bool IdealGas::SupportsDiffusion(void) const { return false; }

// if a subclass material supports artificial viscosity, override this and return TRUE
bool IdealGas::SupportsArtificialViscosity(void) const { return true; }



