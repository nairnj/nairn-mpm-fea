/********************************************************************************
    TaitLiquid.cpp
    nairn-mpm-fea

    Created by John Nairn, Dec 4, 2013.
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    See Create_MPM_Material for details in google code web site wiki
 ********************************************************************************/

#include "stdafx.h"
#include "Materials/TaitLiquid.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "System/UnitsController.hpp"
#include "Read_XML/Expression.hpp"
#include "Materials/HardeningLawBase.hpp"
#include "Materials/ClampedNeohookean.hpp"

// Find pressure by incremental equation in constitutive law
// Code was removed after revsion 2715
//#define INCREMENTAL_PRESSURE

#pragma mark TaitLiquid::Constructors and Destructors

// Constructor
TaitLiquid::TaitLiquid(char *matName,int matID) : HyperElastic(matName,matID)
{
	viscosity.clear();
	logShearRate.clear();
	numViscosity = 0;
	function = NULL;
}

#pragma mark TaitLiquid::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *TaitLiquid::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	// read properties for this material
    if(strcmp(xName,"viscosity")==0)
    {	viscosity.push_back(0.);
		numViscosity = (int)viscosity.size();
		input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&viscosity[numViscosity-1],gScaling,1.e-3);
    }
	
    else if(strcmp(xName,"logshearrate")==0)
    {	logShearRate.push_back(0.);
		int numRates = (int)logShearRate.size();
		input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&logShearRate[numRates-1],gScaling,1.);
    }

	else if(strcmp(xName,"InitialPressure")==0)
	{	input=PRESSURE_FUNCTION_BLOCK;
		return((char *)this);
	}
	
    return HyperElastic::InputMaterialProperty(xName,input,gScaling);
}

// Verify input properties do calculations
const char *TaitLiquid::VerifyAndLoadProperties(int np)
{
	// make sure all were set
    if(Kbulk <= 0. || numViscosity<1)
		return "TaitLiquid material model needs positive K and at least one viscosity";
	
	// check number
	if(numViscosity>1 && numViscosity!=logShearRate.size())
		return "TaitLiquid material model needs same number of viscosities and log shear rates";
	
	// verify sorted
	int i;
	for(i=1;i<numViscosity;i++)
	{	if(logShearRate[i]<logShearRate[i-1])
			return "TaitLiquid entered shear rates must monotonically increase";
	}
    
	// Viscosity in Specific units using initial rho (times 2)
	// Units mass/(L sec) L^3/mass = L^2/sec
	TwoEtasp = new (nothrow) double[numViscosity];
	if(TwoEtasp == NULL)
		return "TaitLiquid material out of memory";
	for(i=0;i<numViscosity;i++)
		TwoEtasp[i] = 2.*viscosity[i]/rho;
    
    // heating gamma0
    double alphaV = 3.e-6*aI;
    gamma0 = Kbulk*alphaV/(rho*heatCapacity);
	   
	// must call super class
	return HyperElastic::VerifyAndLoadProperties(np);
}

// print mechanical properties to the results
void TaitLiquid::PrintMechanicalProperties(void) const
{
	// print properties
	PrintProperty("K",Kbulk*UnitsController::Scaling(1.e-6),"");
    PrintProperty("gam0",gamma0,"");
	PrintProperty("a",aI,"");
    cout << endl;
	
	int i;
	for(i=0;i<numViscosity;i++)
	{	PrintProperty("eta",viscosity[i]*UnitsController::Scaling(1.e3),UnitsController::Label(VISCOSITY_UNITS));
		if(numViscosity>1)
		{	char hline[100];
            size_t hsize=100;
			snprintf(hline,hsize,"1/(%s)",UnitsController::Label(TIME_UNITS));
			PrintProperty("log(rate)",logShearRate[i],hline);
		}
		cout << endl;
	}
	
    if(function!=NULL)
    {	cout << "Initial Pressure  = " << function->GetString() << " " << UnitsController::Label(PRESSURE_UNITS) << " (";
		cout << "mass adjusted to match)" << endl;
    }
}

// check if analysis is allowed
// throws CommonException()
void TaitLiquid::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
    {	throw CommonException("TaitLiquid material cannot do 2D plane stress analysis",
                          "TaitLiquid::ValidateForUse");
    }

    return HyperElastic::ValidateForUse(np);
}

// If needed, a material can initialize particle state
void TaitLiquid::SetInitialParticleState(MPMBase *mptr,int np,int offset) const
{
	// is a pressure being set?
	if(function!=NULL)
	{	// function should be for pressure in MPa (e.g., rho*g*h)
		unordered_map<string,double> vars;
		vars["x"] = mptr->pos.x;
		vars["y"] = mptr->pos.y;
		vars["z"] = mptr->pos.z;
		
		// convert to internal specific pressure units of N/m^2 mm^3/g
		// Divide by rho0, which cancels with rho0 in Ksp when getting Jinit
		double Psp = UnitsController::Scaling(1.e6)*function->EvaluateFunction(vars)/rho;
		
		// Find initial Jinit
		// Note that an initial temperature will cause change in pressure
		double Jinit = 1. - TAIT_C*log(1+Psp/(TAIT_C*Ksp));
		mptr->SetHistoryDble(J_History,Jinit,offset);

		// set the particle pressure (which needs to be  Kirchoff pressure/rho0 = J P0/rho0)
		mptr->SetPressure(Jinit*Psp);
		
		// change mass to match new relative volume
		// Tracked particle deformation is relative to this initial pressurized state
		mptr->mp /= Jinit;
	}
    
    // call super class for Cauchy Green strain
    HyperElastic::SetInitialParticleState(mptr,np,offset);
}

#pragma mark TaitLiquid:History Data Methods

// return number of bytes needed for history data
int TaitLiquid::SizeOfHistoryData(void) const { return 3*sizeof(double); }

// Particle J
char *TaitLiquid::InitHistoryData(char *pchr,MPMBase *mptr)
{	double *p = CreateAndZeroDoubles(pchr,3);
	p[J_History]=1.;					// J
	p[J_History+1]=1.;					// Jres
	// J_History+2 is shear rate, initially zero
	return (char *)p;
}

// reset history data
void TaitLiquid::ResetHistoryData(char *pchr,MPMBase *mptr)
{	double *p = (double *)pchr;
	p[J_History]=1.;					// J
	p[J_History+1]=1.;					// Jres
	p[J_History+2]=0.;
}

// Number of history variables
int TaitLiquid::NumberOfHistoryDoubles(void) const { return 3; }

#pragma mark TaitLiquid:Step Methods

// Apply Constitutive law, check np to know what type
void TaitLiquid::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,
									ResidualStrains *res,int historyOffset) const
{
	// Update total deformation gradient, and calculate trial B
	double detdF = IncrementDeformation(mptr,du,NULL,np);
	
	// Get new J and save result on the particle
	double Jprev = mptr->GetHistoryDble(J_History,historyOffset);
	double J = detdF * Jprev;
	mptr->SetHistoryDble(J_History,J,historyOffset);
	
    // account for residual stresses
    double Jresprev = mptr->GetHistoryDble(J_History+1,historyOffset);
	double dJres = GetIncrementalResJ(mptr,res,Jresprev);
	double Jres = dJres * Jresprev;
	mptr->SetHistoryDble(J_History+1,Jres,historyOffset);
    double Jeff = J/Jres;

    // new Kirchhoff pressure (over rho0) from Tait equation
	double p0=mptr->GetPressure();
    double pressure = J*TAIT_C*Ksp*(exp((1.-Jeff)/TAIT_C)-1.);
	mptr->SetPressure(pressure);
	
	// volume change
	double delV = 1. - 1./detdF;
	
	// artificial viscosity
	double QAVred = 0.;
	double AVEnergy = 0.;
	if(delV<0. && artificialViscosity)
	{	double Kratio = Jeff*(1.+pressure/(TAIT_C*Ksp));
		QAVred = GetArtificialViscosity(delV/delTime,sqrt((Kbulk*Kratio)/rho),mptr);
		AVEnergy += fabs(QAVred*delV);
		pressure += QAVred;
		mptr->IncrementPressure(QAVred);
	}

	// incremental energy per unit mass - dilational part
    double avgP = 0.5*(p0+pressure);
    double workEnergy = -avgP*delV;
    
	// incremental residual energy per unit mass
	double delVres = 1. - 1./dJres;
	double resEnergy = -avgP*delVres;
	
    // viscosity term = 2 eta (0.5(grad v) + 0.5*(grad V)^T - (1/3) tr(grad v) I) = 2 eta dev(grad v)
    // (i.e., deviatoric part of the symmetric strain tensor, 2 is for conversion to engineering shear strain)
	// simple shear rate = |2 dev(grad v)|
    Matrix3 shear;
    double c[3][3];
	double shearRate;
    c[0][0] = (2.*du(0,0)-du(1,1)-du(2,2))/3.;
    c[1][1] = (2.*du(1,1)-du(0,0)-du(2,2))/3.;
    c[2][2] = (2.*du(2,2)-du(0,0)-du(1,1))/3.;
    c[0][1] = 0.5*(du(0,1)+du(1,0));
    c[1][0] = c[0][1];
	shearRate = c[0][0]*c[0][0] + c[1][1]*c[1][1] + c[2][2]*c[2][2]
	+ 2.*c[0][1]*c[0][1];
    if(np==THREED_MPM)
    {   c[0][2] = 0.5*(du(0,2)+du(2,0));
        c[2][0] = c[0][2];
        c[1][2] = 0.5*(du(1,2)+du(2,1));
        c[2][1] = c[1][2];
		shearRate += 2.*(c[0][2]*c[0][2] + c[1][2]*c[1][2]);
        shear.set(c);
    }
    else
        shear.set(c[0][0],c[0][1],c[1][0],c[1][1],c[2][2]);
	shearRate = 2.*sqrt(shearRate)/delTime;
	
	// Store shear rate
	mptr->SetHistoryDble(J_History+2,shearRate,historyOffset);
	
	// Get effective viscosity
	double twoetaspRate = GetTwoEtaOverRho(shearRate);
	
    // Get Kirchoff shear stress (over rho0)
    shear.Scale(J*twoetaspRate/delTime);
    
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
    
    // shear work per unit mass = tau.du = tau.tau*delTime/twoetaspRate
    double shearWork = sp->xx*sp->xx + sp->yy*sp->yy + sp->zz*sp->zz + 2.*sp->xy*sp->xy;
    if(np==THREED_MPM) shearWork += 2.*(sp->xz*sp->xz + sp->yz*sp->yz);
    shearWork *= delTime/twoetaspRate;
    mptr->AddWorkEnergyAndResidualEnergy(workEnergy+shearWork,resEnergy);
    
    // particle isentropic temperature increment dT/T = - J (K/K0) gamma0 Delta(V)/V
    // Delta(V)/V = 1. - 1/detdF (total volume)
	double Kratio = Jeff*(1.+pressure/(TAIT_C*Ksp));
	double dTq0 = -J*Kratio*gamma0*mptr->pPreviousTemperature*delV;
    
    // heat energy
    // dPhi = shearWork is lost due to shear term
    IncrementHeatEnergy(mptr,dTq0,shearWork+AVEnergy);
}

// get 2 eta/rho for use in constitutive law (<=0 to get zero shear rate viscosity)
double TaitLiquid::GetTwoEtaOverRho(double shearRate) const
{
	double twoetaspRate = 0.;
	if(numViscosity==1 || shearRate<=0.)
	{	twoetaspRate = TwoEtasp[0];
	}
	else
	{	shearRate = log10(shearRate);
		if(shearRate < logShearRate[0])
			twoetaspRate = TwoEtasp[0];
		else if(shearRate > logShearRate[numViscosity-1])
			twoetaspRate = TwoEtasp[numViscosity-1];
		else
		{	// interpolate
			for(int i=1;i<numViscosity;i++)
			{	if(shearRate <= logShearRate[i])
				{	// between i-1 and i
					double fract = (logShearRate[i]-shearRate)/(logShearRate[i]-logShearRate[i-1]);
					twoetaspRate = fract*TwoEtasp[i-1] + (1.-fract)*TwoEtasp[i];
					break;
				}
			}
		}
	}
	return twoetaspRate;
}

// This method has several options:
//	1. Solve x(1+k*eta(x*gmaxdot)) - 1 = 0 on interval 0 < x < 1 and return eta(g(dot))
//		Note: x = g(dot)/gmaxdot
//  2. If solution not possible, bracket the solution to y = x(1+k*eta(x*gmaxdot)) - 1 = 0
//		Such that y1(x1)<0 and y(x2)>0
//	3. If can't help, set brackets to 0 and 1
double TaitLiquid::BracketContactLawShearRate(double gmaxdot,double k,double &x1,double &y1,double &x2,double &y2) const
{
	// if no shear dependence, return constant viscosity
	if(numViscosity==1) return viscosity[0];
	
	// get log(gmax(dot)) and exit if below first break point
	double logsmax = log10(gmaxdot);
	if(logsmax <= logShearRate[0]) return viscosity[0];
	
	// start at first point (and done if already positive)
	x1 = pow(10,logShearRate[0]-logsmax);			// always < 1
	y1 = x1*(1.+k*viscosity[0]) - 1.;
	if(y1 >= 0) return viscosity[0];
	
	// find interval were y(x2)>0
	for(int i=1;i<numViscosity;i++)
	{	// done if past the max, return bracket x1 and x2=1
		if(logShearRate[i] > logsmax)
		{	x2 = 1.;
			y2 = k*GetViscosity(gmaxdot);
			return -1.;
		}
		
		// get new x2 (which must be < 1), y2, exit if y2>0
		x2 = pow(10,logShearRate[i]-logsmax);
		y2 = x2*(1.+k*viscosity[i]) - 1.;
		if(y2 >= 0.) return -1;
		
		// reset starting point
		x1 = x2;
		y1 = y2;
	}
	
	// answer is between last shear rate and max shear rate where this material uses a constant
	return viscosity[numViscosity-1];
}

#pragma mark TaitLiquid::Accessors

// return unique, short name for this material
const char *TaitLiquid::MaterialType(void) const { return "Tait Liquid"; }

// Calculate wave speed for material in L/sec.
double TaitLiquid::WaveSpeed(bool threeD,MPMBase *mptr) const
{	return sqrt(Kbulk/rho);
}

// Calculate current wave speed in L/sec for a deformed particle
double TaitLiquid::CurrentWaveSpeed(bool threeD,MPMBase *mptr,int offset) const
{
    double J = mptr->GetHistoryDble(J_History,offset);;
	double dTemp=mptr->pPreviousTemperature-thermal.reference;
	double resStretch = CTE1*dTemp;
	if(DiffusionTask::HasDiffusion())
	{	double dConc = diffusion->GetDeltaConcentration(mptr);
		resStretch += exp(CME1*dConc);
	}
	double Jres = exp(3.*resStretch);
    double Kratio = (J/Jres)*(1.+mptr->GetPressure()/(TAIT_C*Ksp));
    return sqrt((Kbulk*Kratio)/rho);
}

// Copy stress to a read-only tensor variable after combininng deviatoric and pressure
Tensor TaitLiquid::GetStress(Tensor *sp,double pressure,MPMBase *mptr) const
{	return GetStressPandDev(sp,pressure,mptr);
}

// store a new total stress on a particle's stress and pressure variables
void TaitLiquid::SetStress(Tensor *spnew,MPMBase *mptr) const
{	SetStressPandDev(spnew,mptr);
}

// Increment thickness (zz) stress through deviatoric stress and pressure
void TaitLiquid::IncrementThicknessStress(double dszz,MPMBase *mptr) const
{	IncrementThicknessStressPandDev(dszz,mptr);
}

// Get current relative volume change = J (which this material tracks)
double TaitLiquid::GetCurrentRelativeVolume(MPMBase *mptr,int offset) const
{   return mptr->GetHistoryDble(J_History,offset);
}

// setting initial pressure function if needed
// Fuunction should evaulate to pressure
// For gravity, P0 = rho*g*depth
// throws std::bad_alloc, SAXException()
void TaitLiquid::SetPressureFunction(char *pFunction)
{
	// NULL or empty is an error
	if(pFunction==NULL)
		ThrowSAXException("Initial pressure function of position is missing");
	
	// duplicate is error
	if(function!=NULL)
		ThrowSAXException("Duplicate initial pressure function");

	function = Expression::CreateExpression(pFunction,"Initial pressure function is not valid");
}

// get viscosity eta = 0.5*rho*(2 eta/rho) for use in constitutive law
double TaitLiquid::GetViscosity(double shearRate) const { return 0.5*rho*GetTwoEtaOverRho(shearRate); }

// if a subclass material supports artificial viscosity, override this and return TRUE
bool TaitLiquid::SupportsArtificialViscosity(void) const { return true; }
