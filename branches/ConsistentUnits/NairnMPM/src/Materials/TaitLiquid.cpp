/********************************************************************************
    TaitLiquid.cpp
    nairn-mpm-fea

    Created by John Nairn, Dec 4, 2013.
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    See Create_MPM_Material for details in google code web site wiki
 ********************************************************************************/

#include "Materials/TaitLiquid.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Read_XML/mathexpr.hpp"
#include "System/UnitsController.hpp"

// This model tracks on volume check. It does prevent particles degenerating into
// needles, but it is probably not correct thing to do for correct shape functions.
// The better approach might be just to not plot the transformed particles
// Note that CPDI will soon fail if particle become needles, but uGIMP can procees
//#define NO_SHEAR_MODEL

// global expression variables
double TaitLiquid::xPos=0.;
double TaitLiquid::yPos=0.;
double TaitLiquid::zPos=0.;
PRVar tlTimeArray[4] = { NULL,NULL,NULL };

#pragma mark TaitLiquid::Constructors and Destructors

// Constructors
TaitLiquid::TaitLiquid() {}

// The default contructor should call a parent class constructor and
// then fill in any new initialization.
TaitLiquid::TaitLiquid(char *matName) : HyperElastic(matName)
{
    viscosity = -1.;
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
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&viscosity,gScaling,1.e-3);
    }
	
	else if(strcmp(xName,"InitialPressure")==0)
	{	input=PRESSURE_FUNCTION_BLOCK;
		return((char *)this);
	}
	
    return HyperElastic::InputMaterialProperty(xName,input,gScaling);
}

// Verify input properties do calculations; if problem return string with an error message
// If OK, MUST pass on to super class. This is called just before PrintMaterial
// (see also ValidateForUse() for checks that depend on MPM calculation mode)
const char *TaitLiquid::VerifyAndLoadProperties(int np)
{
	// make sure all were set
    if(Kbulk <= 0. || viscosity<0.)
		return "TaitLiquid material model needs positive K and zero or positive viscosity";
    
	// Viscosity in Specific units using initial rho (times 2)
	// Units mass/(L sec) L^3/mass = L^2/sec
	TwoEtasp = 2.*viscosity/rho;
    
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
	PrintProperty("eta",viscosity*UnitsController::Scaling(1.e3),UnitsController::Label(VISCOSITY_UNITS));
	PrintProperty("a",aI,"");
    cout << endl;
    PrintProperty("gam0",gamma0,"");
    cout << endl;
    if(function!=NULL)
    {   char *expr=function->Expr('#');
        cout << "Initial Pressure  = " << expr << " " << UnitsController::Label(PRESSURE_UNITS) << " (";
        delete [] expr;
		cout << "mass adjusted to match)" << endl;
    }
}

// if analysis not allowed, throw a CommonException
void TaitLiquid::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
    {	throw CommonException("TaitLiquid material cannot do 2D plane stress analysis",
                          "TaitLiquid::ValidateForUse");
    }

    return HyperElastic::ValidateForUse(np);
}

// If needed, a material can initialize particle state
void TaitLiquid::SetInitialParticleState(MPMBase *mptr,int np) const
{
	// is a pressure being set?
	if(function!=NULL)
	{	// function should be for pressure in MPa (e.g., rho*g*h)
        xPos=mptr->pos.x;
        yPos=mptr->pos.y;
        zPos=mptr->pos.z;
		
		// convert to internal specific pressure units of N/m^2 mm^3/g
		// Divide by rho0, which cancels with rho0 in Ksp when getting Jinit
		double Psp = UnitsController::Scaling(1.e6)*function->Val()/rho;
		
		// Find initial Jinit
		// Note that an initial temperature will cause change in pressure
		double Jinit = 1. - TAIT_C*log(1+Psp/(TAIT_C*Ksp));
		mptr->SetHistoryDble(J_history,Jinit);

		// set the particle pressure (which needs to be  Kirchoff pressure/rho0 = J P0/rho0)
		mptr->SetPressure(Jinit*Psp);
		
		// change mass to match new relative volume
		// Tracked particle deformation is relative to this initial pressurized state
		mptr->mp /= Jinit;
	}
    
    // call super class for Cauchy Green strain
    HyperElastic::SetInitialParticleState(mptr,np);
}

#pragma mark TaitLiquid:HistoryVariables

// Particle J
char *TaitLiquid::InitHistoryData(void)
{	J_history = 0;
	double *p = CreateAndZeroDoubles(J_history+1);
	p[J_history]=1.;					// J
	return (char *)p;
}

// this material has one
double TaitLiquid::GetHistory(int num,char *historyPtr) const
{
    double history=0.;
	if(num>0 && num<=J_history+1)
	{	double *p=(double *)historyPtr;
		history=p[num-1];
	}
    return history;
}

#pragma mark TaitLiquid:Step Methods

// Apply Constitutive law, check np to know what type
void TaitLiquid::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
#ifdef NO_SHEAR_MODEL
    // get incremental deformation gradient
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
    
    // decompose dF into dR and dU
    Matrix3 dR;
    Matrix3 dU = dF.RightDecompose(&dR,NULL);
	
	// current deformation gradient
    double detdF = dF.determinant();
	Matrix3 pF = mptr->GetDeformationGradientMatrix();
    Matrix3 F = dR*pF;
    if(np==THREED_MPM)
        F.Scale(pow(detdF,1./3.));
    else
        F.Scale2D(sqrt(detdF));
	
	// new deformation matrix with volume change onle
    mptr->SetDeformationGradientMatrix(F);
#else

	// Update total deformation gradient, and calculate trial B
	double detdF = IncrementDeformation(mptr,du,NULL,np);
#endif
    
    // Get new J and save result on the particle
	double J = detdF * mptr->GetHistoryDble(J_history);
    mptr->SetHistoryDble(J_history,J);
    
    // account for residual stresses
	double dresStretch,resStretch = GetResidualStretch(mptr,dresStretch,res);
	double Jres = resStretch*resStretch*resStretch;
    double Jeff = J/Jres;

    // new Kirchhoff pressure (over rho0) from Tait equation
	double p0=mptr->GetPressure();
    double pressure = J*TAIT_C*Ksp*(exp((1.-Jeff)/TAIT_C)-1.);
    mptr->SetPressure(pressure);
    
	// incremental energy per unit mass - dilational part
    double avgP = 0.5*(p0+pressure);
    double delV = 1. - 1./detdF;
    double workEnergy = -avgP*delV;
    
	// incremental residual energy per unit mass
	double delVres = 1. - 1./(dresStretch*dresStretch*dresStretch);
	double resEnergy = -avgP*delVres;
	
    // viscosity term = 2 eta (0.5(grad v) + 0.5*(grad V)^T - (1/3) tr(grad v) I)
    // (i.e., divatoric part of the symmetric strain tensor, 2 is for conversion to engineering shear strain)
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
    
    // Get Kirchoff shear stress (over rho0)
    shear.Scale(J*TwoEtasp/delTime);
    
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
    
    // shear work per unit mass = tau.du = tau.tau*delTime/TwoEtasp
    double shearWork = sp->xx*sp->xx + sp->yy*sp->yy + sp->zz*sp->zz + 2.*sp->xy*sp->xy;
    if(np==THREED_MPM) shearWork += 2.*(sp->xz*sp->xz + sp->yz*sp->yz);
    shearWork *= delTime/TwoEtasp;
    mptr->AddWorkEnergyAndResidualEnergy(workEnergy+shearWork,resEnergy);
    
    // particle isentropic temperature increment dT/T = - J (K/K0) gamma0 Delta(V)/V
    // Delta(V)/V = 1. - 1/detdF (total volume)
	double Kratio = Jeff*(1.+pressure/(TAIT_C*Ksp));
	double dTq0 = -J*Kratio*gamma0*mptr->pPreviousTemperature*delV;
    
    // heat energy is Cv (dT - dTq0) -dPhi
	// Here do Cv (dT - dTq0)
    // dPhi = shearWork is lost due to shear term
    IncrementHeatEnergy(mptr,res->dT,dTq0,shearWork);
}

#pragma mark TaitLiquid::Accessors

// return unique, short name for this material
const char *TaitLiquid::MaterialType(void) const { return "Tait Liquid"; }

// Return the material tag
int TaitLiquid::MaterialTag(void) const { return TAITLIQUID; }

// Calculate wave speed for material in L/sec.
double TaitLiquid::WaveSpeed(bool threeD,MPMBase *mptr) const
{	return sqrt(Kbulk/rho);
}

// Calculate current wave speed in L/sec for a deformed particle
double TaitLiquid::CurrentWaveSpeed(bool threeD,MPMBase *mptr) const
{
    double J = mptr->GetRelativeVolume();
	double dTemp=mptr->pPreviousTemperature-thermal.reference;
	double resStretch = CTE1*dTemp;
	if(DiffusionTask::active)
	{	double dConc=mptr->pPreviousConcentration-DiffusionTask::reference;
		resStretch += exp(CME1*dConc);
	}
	double Jres = exp(3.*resStretch);
    double Kratio = (J/Jres)*(1.+mptr->GetPressure()/(TAIT_C*Ksp));
    return sqrt((Kbulk*Kratio)/rho);
}

// Copy stress to a read-only tensor variable after combininng deviatoric and pressure
Tensor TaitLiquid::GetStress(Tensor *sp,double pressure) const
{   Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}

// Get current relative volume change = J (which this material tracks)
double TaitLiquid::GetCurrentRelativeVolume(MPMBase *mptr) const
{   return mptr->GetHistoryDble(J_history);
}

// setting initial pressure function if needed
// Fuunction should evaulate to pressure
// For gravity, P0 = rho*g*depth
void TaitLiquid::SetPressureFunction(char *pFunction)
{
	// NULL or empty is an error
	if(pFunction==NULL)
		ThrowSAXException("Initial pressure function of position is missing");
	
	// create time variable if needed
	if(tlTimeArray[0]==NULL)
	{	tlTimeArray[0]=new RVar("x",&xPos);
		tlTimeArray[1]=new RVar("y",&yPos);
		tlTimeArray[2]=new RVar("z",&zPos);
	}
	
	// set the function
	if(function!=NULL)
		ThrowSAXException("Duplicate initial pressure function");
	function=new ROperation(pFunction,3,tlTimeArray);
	if(function->HasError())
		ThrowSAXException("Initial pressure function is not valid");
}


