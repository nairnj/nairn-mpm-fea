/********************************************************************************
    JohnsonCook.hpp
    nairn-mpm-fea
    
    Created by John Nairn, August 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
 
    The Johnson-Cook hardening law is
        (Ajc + Bjc εpnjc) [1 + Cjc ln(dεp/ep0jc) ] (1 - Trmjc)
********************************************************************************/

#include "Materials/JohnsonCook.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Exceptions/CommonException.hpp"
#include "System/UnitsController.hpp"

#pragma mark JohnsonCook::Constructors and Destructors

JohnsonCook::JohnsonCook() {}

JohnsonCook::JohnsonCook(MaterialBase *pair) : HardeningLawBase(pair)
{
    // Ajc is in the yield stress in HardenLawBase class
	Bjc = 0.;             // pressure
	njc = 1.;             // dimensionless
	Cjc = 0.;             // dimensionless
	ep0jc = 1.;           // time^-1
	Tmjc = 1000.;         // Melting point in K relative to thermal.reference
    mjc = 1.;             // dimensionless
}

#pragma mark JohnsonCook::Initialization

// Read material properties
char *JohnsonCook::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"Bjc")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&Bjc,gScaling,1.e6);
    }
    else if(strcmp(xName,"Ajc")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&yield,gScaling,1.e6);
	}
    else if(strcmp(xName,"njc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&njc);
    }
    else if(strcmp(xName,"Cjc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&Cjc);
    }
    else if(strcmp(xName,"ep0jc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&ep0jc);
    }
    else if(strcmp(xName,"Tmjc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&Tmjc);
    }
    else if(strcmp(xName,"mjc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&mjc);
    }
    
    return HardeningLawBase::InputMaterialProperty(xName,input,gScaling);
}

// print just yield properties to output window
void JohnsonCook::PrintYieldProperties(void) const
{
    cout << GetHardeningLawName() << endl;
    MaterialBase::PrintProperty("A",yield*UnitsController::Scaling(1.e-6),"");
	MaterialBase::PrintProperty("B",Bjc*UnitsController::Scaling(1.e-6),"");
	MaterialBase::PrintProperty("n",njc,"");
    cout << endl;
	MaterialBase::PrintProperty("C",Cjc,"");
	char glabel[20];
	strcpy(glabel,UnitsController::Label(TIME_UNITS));
	strcat(glabel,"^-1");
	MaterialBase::PrintProperty("ep0",ep0jc,glabel);
    cout << endl;
	MaterialBase::PrintProperty("Tm",Tmjc,"K");
    MaterialBase::PrintProperty("T0",thermal.reference,"K");
	MaterialBase::PrintProperty("m",mjc,"");
    cout << endl;
}

// Private properties used in hardening law
const char *JohnsonCook::VerifyAndLoadProperties(int np)
{	
    if(Tmjc <= thermal.reference)
        return "The melting temperature must be >= the reference temperature";
    
	// reduced prooperties
    Bred = Bjc/parent->rho;
	
    // reduced yield stress or Ajc
	HardeningLawBase::VerifyAndLoadProperties(np);
    
    // ignore strain rates below this
    edotMin = Cjc!=0. ? exp(-0.5/Cjc) : 1.e-20 ;
    eminTerm = 1. + Cjc*log(edotMin) ;
	
	// base class never has error
	return NULL;
}

#pragma mark JohnsonCook:Methods

// size of hardening law properties needed in strain updates
int JohnsonCook::SizeOfHardeningProps(void) const { return sizeof(JCProperties); }

// Get copy of particle-state dependent properties
void *JohnsonCook::GetCopyOfHardeningProps(MPMBase *mptr,int np,void *altBuffer)
{
	JCProperties *p = (JCProperties *)altBuffer;
	
	// homologous temperature (as needed by Johnson and Cook)
	p->hmlgTemp=(mptr->pPreviousTemperature - thermal.reference) /
                (Tmjc - thermal.reference);
    
    if(p->hmlgTemp>1.)
    {   // above the melting point
        p->TjcTerm = 0.;
    }
    else if(p->hmlgTemp>0.)
    {   // between T ref and melting and TjcTerm between 1 (at Tref) and 0 (at T melt)
        p->TjcTerm = 1. - pow(p->hmlgTemp,mjc);
    }
    else
    {   // below T ref or out of range. Pick some number >= 1
        p->TjcTerm = 1.;
    }
    
    // nothing needed from superclass (HardenLawBase)
	return p;
}

// Cast void * to correct pointer and delete it
void JohnsonCook::DeleteCopyOfHardeningProps(void *properties,int np) const
{
	JCProperties *p = (JCProperties *)properties;
	delete p;
}

#pragma mark JohnsonCook::Law Methods

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
// yield = (A + B ep^n + n epdot), where ep=alpint, epdot=dalpha/delTime
double JohnsonCook::GetYield(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	JCProperties *p = (JCProperties *)properties;
    if(p->hmlgTemp>=1.) return 0.;
    double term1 = yldred + Bred*pow(a->alpint,njc);
    double ep = a->dalpha/(delTime*ep0jc);
    double term2 = ep>edotMin ? 1. + Cjc*log(ep) : eminTerm ;
    return term1 * term2 * p->TjcTerm ;
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lambda or depdot/dlambda = sqrt(2./3.)/delTime
double JohnsonCook::GetKPrime(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	JCProperties *p = (JCProperties *)properties;
    if(p->hmlgTemp>=1.) return 0.;
    double ep = a->dalpha/(delTime*ep0jc);
    if(ep>edotMin)
    {   double term1 = yldred + Bred*pow(a->alpint,njc);
        double term2 = 1. + Cjc*log(ep) ;
        return TWOTHIRDS * p->TjcTerm * (Bred*njc*pow(a->alpint,njc-1.)*term2 + Cjc*term1/a->dalpha ) ;
    }
    else
        return TWOTHIRDS * p->TjcTerm * Bred*njc*pow(a->alpint,njc-1.) * eminTerm ;
}

// Get derivative of (1./3.)*yield^2 with respect to lambda for plane stress only
// ... and using dep/dlambda = sqrt(2./3.)*fnp1 where ep=alpint
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)*lambda*fnp1 or depdot/dlambda = sqrt(2./3.)*fnp1/delTime
// Also equal to sqrt(2./3.)*GetYield()*GetKPrime()*fnp1, but in separate call for efficiency
double JohnsonCook::GetK2Prime(MPMBase *mptr,double fnp1,double delTime,HardeningAlpha *a,void *properties) const
{   
    if(DbleEqual(a->alpint,0.)) return 0.;
	JCProperties *p = (JCProperties *)properties;
	if(p->hmlgTemp>=1.) return 0.;
    double term1 = yldred + Bred*pow(a->alpint,njc);
    double ep = a->dalpha/(delTime*ep0jc);
    if(ep>edotMin)
    {   double term2 = 1. + Cjc*log(ep) ;
        return SQRT_EIGHT27THS * term1 * term2 * fnp1 * p->TjcTerm * p->TjcTerm *
                        (Bred*njc*pow(a->alpint,njc-1.)*term2 + Cjc*term1/a->dalpha ) ;
    }
    else
    {   return SQRT_EIGHT27THS * term1 * fnp1 * p->TjcTerm * p->TjcTerm * eminTerm * eminTerm *
                    (Bred*njc*pow(a->alpint,njc-1.)) ;
    }
}

// Return (K(alpha)-K(0)), which is used in dissipated energy calculation
// If K(0) in current particle state differs from yldred, will need to override
double JohnsonCook::GetYieldIncrement(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	JCProperties *p = (JCProperties *)properties;
    if(p->hmlgTemp>=1.) return 0.;
    double ep = a->dalpha/(delTime*ep0jc);
    double term2 = ep>edotMin ? 1. + Cjc*log(ep) : eminTerm ;
	return Bred*pow(a->alpint,njc) * term2 * p->TjcTerm ;
}

// watch for temperature above the melting point and zero out the deviatoric stress
double JohnsonCook::SolveForLambdaBracketed(MPMBase *mptr,int np,double strial,Tensor *stk,
							double Gred,double psKred,double Pfinal,double delTime,
							HardeningAlpha *a,void *properties) const
{
    // if melted, return for zero deviatoric stress
	JCProperties *p = (JCProperties *)properties;
    if(p->hmlgTemp>=1.)
    {   return strial/(2.*Gred);
    }
    
    // assume error in bracking is because near melting, convert error to zero deviatoric stress
    return HardeningLawBase::SolveForLambdaBracketed(mptr,np,strial,stk,Gred,psKred,Pfinal,delTime,a,properties);
}

#pragma mark JohnsonCook::Accessors

// hardening law name
const char *JohnsonCook::GetHardeningLawName(void) const { return "Johnson-Cook hardening"; }

