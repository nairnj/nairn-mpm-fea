/********************************************************************************
    SCGLHardening.hpp
    nairn-mpm-fea

    Created by John Nairn, 1/18/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Hardening Law is
********************************************************************************/

#include "Materials/SCGLHardening.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "System/UnitsController.hpp"

#pragma mark NonlinearHardening::Constructors and Destructors

// Constructors
SCGLHardening::SCGLHardening() {}

// Constructors
SCGLHardening::SCGLHardening(MaterialBase *pair) : HardeningLawBase(pair)
{
	// defaults are some Tungsten properties
	yield=2200.*UnitsController::Scaling(1.e6);
	yieldMax=0.;
	beta=0.0;
	nhard=1.0;
	
	GPp=0.01e-3*UnitsController::Scaling(1.e-6);
	GTp=-2.2e-4;
}

#pragma mark LinearHardening::Initialize

// Read hardening law properties
char *SCGLHardening::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"GPpG0")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&GPp,gScaling,1.e-6);
    }
	
    else if(strcmp(xName,"GTpG0")==0)
    {	input=DOUBLE_NUM;
        return((char *)&GTp);
    }
	
    else if(strcmp(xName,"betahard")==0)
    {	input=DOUBLE_NUM;
        return((char *)&beta);
    }
	
    else if(strcmp(xName,"nhard")==0)
    {	input=DOUBLE_NUM;
        return((char *)&nhard);
    }
	
    else if(strcmp(xName,"yieldMax")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&yieldMax,gScaling,1.e6);
    }
	
    return HardeningLawBase::InputMaterialProperty(xName,input,gScaling);
}

// verify settings and some initial calculations
const char *SCGLHardening::VerifyAndLoadProperties(int np)
{
    // but fails is not order correctly
    if(yieldMax < yield)
    {   return "The maximum yield stress in SCGLHardening law is less than the yield stress";
    }
    
    // reduced yield stress in base class
	HardeningLawBase::VerifyAndLoadProperties(np);
    
	// reduced maximum yield stress
	yldMaxred = yieldMax/parent->rho;
    
    // reduced shear modulus pressure dependence
    GPpred = GPp*parent->rho;
	
	// base class never has an error
    return NULL;
}

// print just yield properties to output window
void SCGLHardening::PrintYieldProperties(void) const
{
    cout << GetHardeningLawName() << endl;
    
    // yield
    MaterialBase::PrintProperty("yld",yield*UnitsController::Scaling(1.e-6),"");
    MaterialBase::PrintProperty("beta",beta,"");
    MaterialBase::PrintProperty("nhard",nhard,"");
    MaterialBase::PrintProperty("yMax",yieldMax*UnitsController::Scaling(1.e-6),"");
    cout << endl;

	// shear temperature and pressure dependence
	char glabel[20];
	strcpy(glabel,UnitsController::Label(PRESSURE_UNITS));
	strcat(glabel,"^-1");
	MaterialBase::PrintProperty("Gp'/G0",GPp*UnitsController::Scaling(1.e6),glabel);
	MaterialBase::PrintProperty("GT'/G0",GTp,"K^-1");
	cout << endl;
}


#pragma mark NonlinearHardening:Methods

// size of hardening law properties needed in strain updates
int SCGLHardening::SizeOfHardeningProps(void) const { return sizeof(SCGLProperties); }

// Get particle-state dependent properties (filled by Get Shear Ratio)
void *SCGLHardening::GetCopyOfHardeningProps(MPMBase *mptr,int np,void *altBuffer)
{
	SCGLProperties *p = (SCGLProperties *)altBuffer;
	return p;
}
	
// Cast void * to correct pointer and delete it
void SCGLHardening::DeleteCopyOfHardeningProps(void *properties,int np) const
{
	SCGLProperties *p = (SCGLProperties *)properties;
	delete p;
}

// handle prressure and temperture depence of the shear modulus
// Find ratio of current shear modulus to initial shear modulus including factor of J
//  because dealing in Kirchoff stress
// Store results needed later in hardenling law properties
double SCGLHardening::GetShearRatio(MPMBase *mptr,double pressure,double J,void *properties) const
{
    double dTemp = mptr->pPreviousTemperature - thermal.reference;
    double neta = pow(1./J,ONETHIRD);
    double Gratio = J * (1. + GPpred*pressure/neta + GTp*dTemp);
	if(Gratio < 0.) Gratio = 0.;
	if(properties!=NULL)
	{	SCGLProperties *p = (SCGLProperties *)properties;
		p->Gratio = Gratio;
	}
    return Gratio;
}

#pragma mark NonlinearHardening::Law Methods

// Return yield stress for current conditions and it is specific Cauchy stress
//   (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
// yield = yldred*(1 + beta ep)^n * Gratio, where ep=alpint
// but leading term is limited to yldMaxred
double SCGLHardening::GetYield(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	SCGLProperties *p = (SCGLProperties *)properties;
    return fmin(yldred*pow(1.+beta*a->alpint,nhard),yldMaxred)*p->Gratio;
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
// ... and as specfic Cauchy stress
double SCGLHardening::GetKPrime(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{	
    // slope zero if in constant max yield condition
    if(yldred*pow(1.+beta*a->alpint,nhard)>=yldMaxred) return 0.;

    // return slope
	SCGLProperties *p = (SCGLProperties *)properties;
    double factor=yldred*p->Gratio;
    double bfactor = DbleEqual(nhard,1.) ? beta :
    beta*nhard*pow(1.+beta*a->alpint,nhard-1.) ;
    return TWOTHIRDS*factor*bfactor;
}

// this material does not support plane stress calculations
double SCGLHardening::GetK2Prime(MPMBase *mptr,double fnp1,double delTime,HardeningAlpha *a,void *properties) const
{
    // slope zero if in constant max yield condition
    if(yldred*pow(1.+beta*a->alpint,nhard)>=yldMaxred) return 0.;

	SCGLProperties *p = (SCGLProperties *)properties;
    double factor=yldred*p->Gratio;
    return SQRT_EIGHT27THS*factor*factor*beta*nhard*pow(1.+beta*a->alpint,2.*nhard-1)*fnp1;
}

// Return (K(alpha)-K(0)), which is used in dissipated energy calculation
double SCGLHardening::GetYieldIncrement(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{	SCGLProperties *p = (SCGLProperties *)properties;
	return (fmin(yldred*pow(1.+beta*a->alpint,nhard),yldMaxred)-yldred)*p->Gratio;
}

#pragma mark NonlinearHardening::Accessors

// hardening law name
const char *SCGLHardening::GetHardeningLawName(void) const { return "SCGL hardening"; }

