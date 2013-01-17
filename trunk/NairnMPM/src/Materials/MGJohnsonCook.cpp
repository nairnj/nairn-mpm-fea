/********************************************************************************
    MGJohnsonCook.hpp
    NairnMPM

    Created by John Nairn, 1/17/2012.
    Copyright (c) 2013 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/MGJohnsonCook.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark MGJohnsonCook::Constructors and Destructors

// Constructors
MGJohnsonCook::MGJohnsonCook() {}
MGJohnsonCook::MGJohnsonCook(char *matName) : MGSCGLMaterial(matName)
{
    // Ajc is in the yield stress, which is required propertie
	Bjc=0.;             // in MPa
	njc=1.;             // dimensionless
	Cjc=0.;             // dimensionless
	ep0jc=1.;           // sec^-1
	Tmjc=1000.;         // Melting point in K relative to thermal.reference
    mjc=1.;             // dimensionless
    
    // to avvoid errors
    yieldMax = yield+1.;
    readYield = TRUE;
}

#pragma mark MGJohnsonCook::Initialization

// Read material properties
char *MGJohnsonCook::InputMat(char *xName,int &input)
{
    if(strcmp(xName,"Bjc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&Bjc);
    }
    if(strcmp(xName,"Ajc")==0)
    {	input=DOUBLE_NUM;
        readYield=TRUE;
        return((char *)&yield);
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
    
    return(MGSCGLMaterial::InputMat(xName,input));
}

// print just yield properties to output window
void MGJohnsonCook::PrintYieldProperties(void)
{
	PrintProperty("A",yield,"");
	PrintProperty("B",Bjc,"");
	PrintProperty("n",njc,"");
    cout << endl;
	PrintProperty("C",Cjc,"");
	PrintProperty("ep0",ep0jc,"s^-1");
    cout << endl;
	PrintProperty("Tm",Tmjc,"K");
    PrintProperty("T0",thermal.reference,"K");
	PrintProperty("m",mjc,"");
    cout << endl;
}

// Private properties used in constitutive law
// For variable shear and bulk moduli, subclass can overrive
//		LoadMechanicalProps(MPMBase *mptr,int np) and set new
//		Gred and Kred
void MGJohnsonCook::InitialLoadMechProps(int makeSpecific,int np)
{	
	// reduced prooperties (Units Pa - cm^3/g)
    Bred = Bjc*1.e6/rho;
	
	IsoPlasticity::InitialLoadMechProps(makeSpecific,np);
    
    // ignore strain rates below this
    edotMin = Cjc!=0. ? exp(-0.5/Cjc) : 1.e-20 ;
    eminTerm = 1. + Cjc*log(edotMin) ;
    
	// call superclass
	MGSCGLMaterial::InitialLoadMechProps(makeSpecific,np);
}

#pragma mark JohnsonCook:Methods

// State dependent material properties
void MGJohnsonCook::LoadMechanicalProps(MPMBase *mptr,int np)
{
	// homologous temperature (as named by Johnson and Cook)
	double hmlgTemp=(mptr->pPreviousTemperature - thermal.reference) / 
    (Tmjc - thermal.reference);
    
    TjcTerm = hmlgTemp < 0. ? 1. - hmlgTemp : 1. - pow(hmlgTemp,mjc) ;
}

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
// yield = (A + B ep^n + n epdot), where ep=alpint, epdot=dalpha/delTime
double MGJohnsonCook::GetYield(MPMBase *mptr,int np,double delTime)
{
    double term1 = yldred + Bred*pow(alpint,njc);
    double ep = dalpha/(delTime*ep0jc);
    double term2 = ep>edotMin ? 1. + Cjc*log(ep) : eminTerm ;
    return term1 * term2 * TjcTerm ;
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lambda or depdot/dlambda = sqrt(2./3.)/delTime
double MGJohnsonCook::GetKPrime(MPMBase *mptr,int np,double delTime)
{
    double ep = dalpha/(delTime*ep0jc);
    if(ep>edotMin)
    {   double term1 = yldred + Bred*pow(alpint,njc);
        double term2 = 1. + Cjc*log(ep) ;
        return TWOTHIRDS * TjcTerm * (Bred*njc*pow(alpint,njc-1.)*term2 + Cjc*term1/dalpha ) ;
    }
    else
        return TWOTHIRDS * TjcTerm * Bred*njc*pow(alpint,njc-1.) * eminTerm ;
}

// Get derivative of (1./3.)*yield^2 with respect to lambda for plane stress only
// ... and using dep/dlambda = sqrt(2./3.)*fnp1 where ep=alpint
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)*lambda*fnp1 or depdot/dlambda = sqrt(2./3.)*fnp1/delTime
// Also equal to sqrt(2./3.)*GetYield()*GetKPrime()*fnp1, but in separate call for efficiency
double MGJohnsonCook::GetK2Prime(MPMBase *mptr,double fnp1,double delTime)
{
    double term1 = yldred + Bred*pow(alpint,njc);
    double ep = dalpha/(delTime*ep0jc);
    if(ep>edotMin)
    {   double term2 = 1. + Cjc*log(ep) ;
        return SQRT_EIGHT27THS * term1 * term2 * fnp1 * TjcTerm * TjcTerm *
        (Bred*njc*pow(alpint,njc-1.)*term2 + Cjc*term1/dalpha ) ;
    }
    else
    {   return SQRT_EIGHT27THS * term1 * fnp1 * TjcTerm * TjcTerm * eminTerm * eminTerm *
        (Bred*njc*pow(alpint,njc-1.)) ;
    }
}


#pragma mark NewMaterial::Accessors

// Return the material tag
int MGJohnsonCook::MaterialTag(void) { return MGJOHNSONCOOK; }

// return unique, short name for this material
const char *MGJohnsonCook::MaterialType(void) { return "MG-Johnson-Cook Material"; }

