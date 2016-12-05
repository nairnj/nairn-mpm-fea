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
#include "Global_Quantities/ThermalRamp.hpp"
#include "Read_XML/mathexpr.hpp"
#include "System/UnitsController.hpp"

// This model tracks only volume change. It does prevent particles degenerating into
// needles, but it is probably not correct thing to do for correct shape functions.
// The better approach might be just to not plot the transformed particles
// Note that CPDI will soon fail if particle become needles, but uGIMP can procees
//#define NO_SHEAR_MODEL

// This tracks full deformation, but B matrix is set to volume change only so B
// is Diagonal with Bii = J^(2/3) axisymmetric and 3D or Bxx=Byy=J for plane strain
#define ELASTIC_B_MODEL

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

// Verify input properties do calculations; if problem return string with an error message
// If OK, MUST pass on to super class. This is called just before PrintMaterial
// (see also ValidateForUse() for checks that depend on MPM calculation mode)
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
			sprintf(hline,"1/(%s)",UnitsController::Label(TIME_UNITS));
			PrintProperty("log(rate)",logShearRate[i],hline);
		}
		cout << endl;
	}
	
    if(function!=NULL)
    {   char *expr=function->Expr('#');
        cout << "Initial Pressure  = " << expr << " " << UnitsController::Label(PRESSURE_UNITS) << " (";
        delete [] expr;
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
// (Don't call from parallel code due to function)
void TaitLiquid::SetInitialParticleState(MPMBase *mptr,int np,int offset) const
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

// this material has two
double TaitLiquid::GetHistory(int num,char *historyPtr) const
{
    double history=0.;
	if(num>0 && num<=3)
	{	double *p=(double *)historyPtr;
		history=p[num-1];
	}
    return history;
}

#pragma mark TaitLiquid:Step Methods

// Apply Constitutive law, check np to know what type
void TaitLiquid::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,
									ResidualStrains *res,int historyOffset) const
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
#ifdef ELASTIC_B_MODEL
    // get incremental deformation gradient
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
	double detdF = dF.determinant();
	
	// current deformation gradient
	Matrix3 pF = mptr->GetDeformationGradientMatrix();
	
	// new deformation matrix
	const Matrix3 F = dF*pF;
    mptr->SetDeformationGradientMatrix(F);
#else

	// Update total deformation gradient, and calculate trial B
	double detdF = IncrementDeformation(mptr,du,NULL,np);
#endif
#endif
    
    // Get new J and save result on the particle
	double J = detdF * mptr->GetHistoryDble(J_History,historyOffset);
    mptr->SetHistoryDble(J_History,J,historyOffset);

#ifdef ELASTIC_B_MODEL
	// store pressure strain as elastic B
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
#endif
    
    // account for residual stresses
	double dJres = GetIncrementalResJ(mptr,res);
	double Jres = dJres * mptr->GetHistoryDble(J_History+1,historyOffset);
	mptr->SetHistoryDble(J_History+1,Jres,historyOffset);
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
	
	// Get effective visocisy
	double twoetaspRate = 0.;
	if(numViscosity==1)
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
    
    // heat energy is Cv (dT - dTq0) -dPhi
	// Here do Cv (dT - dTq0)
    // dPhi = shearWork is lost due to shear term
    IncrementHeatEnergy(mptr,res->dT,dTq0,shearWork);
}

// When becomes active update J, set B elastic, set pressure, and zero deviatorix stress
void TaitLiquid::BeginActivePhase(MPMBase *mptr,int np,int historyOffset) const
{	double J = mptr->GetRelativeVolume();
	mptr->SetHistoryDble(J_History,J,historyOffset);
#ifdef ELASTIC_B_MODEL
	double Jres = mptr->GetHistoryDble(J_History+1,historyOffset);
    double Jeff = J/Jres;
	
	// store pressure strain as elastic B
	// no need here, will happen in next time step
	
    // new Kirchhoff pressure (over rho0) from Tait equation
    double pressure = J*TAIT_C*Ksp*(exp((1.-Jeff)/TAIT_C)-1.);
    mptr->SetPressure(pressure);
	
	// set deviatoric stress to zero
	ZeroTensor(mptr->GetStressTensor());
#endif
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
	if(DiffusionTask::active)
	{	double dConc=mptr->pPreviousConcentration-DiffusionTask::reference;
		resStretch += exp(CME1*dConc);
	}
	double Jres = exp(3.*resStretch);
    double Kratio = (J/Jres)*(1.+mptr->GetPressure()/(TAIT_C*Ksp));
    return sqrt((Kbulk*Kratio)/rho);
}

// Copy stress to a read-only tensor variable after combininng deviatoric and pressure
Tensor TaitLiquid::GetStress(Tensor *sp,double pressure,MPMBase *mptr) const
{   Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
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

