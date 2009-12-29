/********************************************************************************
    JohnsonCook.hpp
    NairnMPM
    
    Created by John Nairn, August 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/JohnsonCook.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark JohnsonCook::Constructors and Destructors

// Constructors
JohnsonCook::JohnsonCook() {}
JohnsonCook::JohnsonCook(char *matName) : IsoPlasticity(matName)
{
	Bjc=0.;
	njc=1.;
	Cjc=0.;
	ep0jc=1.;
}

#pragma mark JohnsonCook::Initialization

/* If material has new property types, it must override this method and
	1. Define XML tag in the DTD file
	2. If xName matches a new property tag, set input to the type
		of variable (DOUBLE_NUM or INT_NUM) and return pointer
		to the class variable to be set.
	c. If no match, call InputMat() of superclass
*/
// Read material properties
char *JohnsonCook::InputMat(char *xName,int &input)
{
    if(strcmp(xName,"Bjc")==0)
    {	input=DOUBLE_NUM;
        return((char *)&Bjc);
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
    
    return(IsoPlasticity::InputMat(xName,input));
}

// print mechanical properties to the results
void JohnsonCook::PrintMechanicalProperties(void)
{	
    IsotropicMat::PrintMechanicalProperties();
	PrintYieldProperties();
}

// print just yield properties to output window
void JohnsonCook::PrintYieldProperties(void)
{
	PrintProperty("A",yield,"");
	PrintProperty("B",Bjc,"");
	PrintProperty("n",njc,"");
    cout << endl;
	PrintProperty("C",Cjc,"");
	PrintProperty("ep0",ep0jc,"");
    cout << endl;
}

// Private properties used in constitutive law
// For variable shear and bulk moduli, subclass can overrive
//		LoadMechanicalProps(MPMBase *mptr,int np) and set new
//		Gred and Kred
void JohnsonCook::InitialLoadMechProps(int makeSpecific,int np)
{	
	// reduced prooperties
    Bred = Bjc*1.e6/rho;
	etared = Cjc*1.e6/rho;
	
	IsoPlasticity::InitialLoadMechProps(makeSpecific,np);
}

#pragma mark JohnsonCook:Methods

// State dependent material properties
void JohnsonCook::LoadMechanicalProps(MPMBase *mptr,int np)
{
	// calculate new bulk and shear moduli properties needed by IsoPlasticity
}

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
// yield = (A + B ep^n + n epdot), where ep=alpha, epdot=dalpha/dt
double JohnsonCook::GetYield(MPMBase *mptr,int np,double delTime)
{
	return yldred + Bred*pow(alpint,njc) + etared*dalpha/delTime;
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
double JohnsonCook::GetKPrime(MPMBase *mptr,int np,double delTime)
{
	double hard = DbleEqual(njc,1.) ? Bred : njc*Bred*pow(alpint,njc-1) ;
	return (2./3.)*(hard + etared/delTime);
}

// Get derivative of (1./3.)*yield^2 with respect to lambda for plane stress only
// ... and using dep/dlambda = sqrt(2./3.)*fnp1 where ep=alpint
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)*lambda*fnp1 or depdot/dlambda = sqrt(2./3.)*fnp1/delTime
double JohnsonCook::GetK2Prime(MPMBase *mptr,double fnp1,double delTime)
{
	double hard = DbleEqual(njc,1.) ? Bred : njc*Bred*pow(alpint,njc-1) ;
	return sqrt(8./27.)*(yldred + Bred*pow(alpint,njc) + etared*dalpha/delTime)*(hard+etared/delTime)*fnp1;
}

#pragma mark NewMaterial::Accessors

// Return the material tag
int JohnsonCook::MaterialTag(void) { return JOHNSONCOOK; }

// return unique, short name for this material
const char *JohnsonCook::MaterialType(void) { return "Johnson-Cook Material"; }

