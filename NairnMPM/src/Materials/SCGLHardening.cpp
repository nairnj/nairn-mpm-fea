/********************************************************************************
    SCGLHardening.hpp
    NairnMPM

    Created by John Nairn, 1/18/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Hardening Law is
********************************************************************************/

#include "Materials/SCGLHardening.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "MPM_Classes/MPMBase.hpp"

#pragma mark NonlinearHardening::Constructors and Destructors

// Constructors
SCGLHardening::SCGLHardening() {}

// Constructors
SCGLHardening::SCGLHardening(MaterialBase *pair) : HardeningLawBase(pair)
{
	// defaults are some Tungsten properties
	yield=2200.;		// MPa
	yieldMax=0.;		// MPa
	beta=0.0;			// dimensionless
	nhard=1.0;			// dimensionless
	
	GPp=0.01e-3;		// MPa^-1
	GTp=-2.2e-4;		// K^-1
}

//if(yieldMax<yield) return "The maximum yield stress is less than the initial yield stress";

#pragma mark LinearHardening::Initialize

// Read hardening law properties
char *SCGLHardening::InputMat(char *xName,int &input)
{
    if(strcmp(xName,"GPpG0")==0)
    {	input=DOUBLE_NUM;
        return((char *)&GPp);
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
        return((char *)&yieldMax);
    }
	
    return HardeningLawBase::InputMat(xName,input);
}

// verify settings and some initial calculations
const char *SCGLHardening::VerifyProperties(int np)
{
    // but fails is not order correctly
    if(yieldMax < yield)
    {   return "The maximum yield stress in SCGLHardening law is less than the yield stress";
    }
    
    return NULL;
}

// get reduced stress than done
void SCGLHardening::InitialLoadMechProps(int makeSpecific,int np)
{   
    // reduced yield stress in base class
	HardeningLawBase::InitialLoadMechProps(makeSpecific,np);
    
	// reduced maximum yield stress
	yldMaxred = yieldMax*1.e6/parent->rho;
    
    // reduce shear modulus pressure dependence
    GPpred = GPp*parent->rho*1.e-6;
    Gratio = 1.;
    
}

// print just yield properties to output window
void SCGLHardening::PrintYieldProperties(void)
{
    cout << GetHardeningLawName() << endl;
    
    // yield
    MaterialBase::PrintProperty("yld",yield,"");
    MaterialBase::PrintProperty("beta",beta,"");
    MaterialBase::PrintProperty("nhard",nhard,"");
    MaterialBase::PrintProperty("yMax",yieldMax,"");
    cout << endl;

	// shear temperature and pressure dependence
	MaterialBase::PrintProperty("Gp'/G0",GPp,"MPa^-1");
	MaterialBase::PrintProperty("GT'/G0",GTp,"K^-1");
	cout << endl;
}



#pragma mark NonlinearHardening::Law Methods

// handle prressure and temperture depence of the shear modulus
// Find ratio of current shear modulus to G0red including factor of J
//  because dealing in Kirchoff stress
double SCGLHardening::GetShearRatio(MPMBase *mptr,double pressure,double J)
{
    double dTemp = mptr->pPreviousTemperature - thermal.reference;
    double neta = pow(1./J,ONETHIRD);
    Gratio = J * (1. + GPpred*pressure/neta + GTp*dTemp);
	if(Gratio < 0.) Gratio = 0.;
    return Gratio;
}

// Return yield stress for current conditions and it is specific Cauchy stress
//   (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
// yield = yldred*(1 + beta ep)^n * Gred/G0red, where ep=alpint
// but leading term is limited to yldMaxred
double SCGLHardening::GetYield(MPMBase *mptr,int np,double delTime)
{
    return fmin(yldred*pow(1.+beta*alpint,nhard),yldMaxred)*Gratio;
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
// ... and as specfic Cauchy stress
double SCGLHardening::GetKPrime(MPMBase *mptr,int np,double delTime)
{	
    // slope zero if in constant max yield condition
    if(yldred*pow(1.+beta*alpint,nhard)>=yldMaxred) return 0.;

    // return slope
    double factor=yldred*Gratio;
    double bfactor = DbleEqual(nhard,1.) ? beta :
    beta*nhard*pow(1.+beta*alpint,nhard-1.) ;
    return TWOTHIRDS*factor*bfactor;
}

// this material does not support plane stress calculations
double SCGLHardening::GetK2Prime(MPMBase *mptr,double fnp1,double delTime)
{
    // slope zero if in constant max yield condition
    if(yldred*pow(1.+beta*alpint,nhard)>=yldMaxred) return 0.;

    double factor=yldred*Gratio;
    return SQRT_EIGHT27THS*factor*factor*beta*nhard*pow(1.+beta*alpint,2.*nhard-1)*fnp1;
}

#pragma mark NonlinearHardening::Accessors

// hardening law name
const char *SCGLHardening::GetHardeningLawName(void) { return "SCGL hardening"; }

