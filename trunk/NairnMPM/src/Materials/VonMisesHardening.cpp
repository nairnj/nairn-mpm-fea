/********************************************************************************
    VonMisesHardening.cpp
    NairnMPM
    
    Created by Yajun Guo in Jan 2005.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/VonMisesHardening.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Read_XML/CommonReadHandler.hpp"

#pragma mark VonMisesHardening::Constructors and Destructors

// Constructors
VonMisesHardening::VonMisesHardening()
{
}

// Constructors
VonMisesHardening::VonMisesHardening(char *matName) : IsoPlasticity(matName)
{
    // default value of plastic modulus (ideal elastic-plasticity)
    Ep=Epred=0.0;  
    ET=-1.;
	linearHardening=TRUE;
	beta=0.;
	npow=1.;
}

#pragma mark VonMisesHardening::Initialization

// print to output window
void VonMisesHardening::PrintMechanicalProperties(void)
{	
    IsotropicMat::PrintMechanicalProperties();
	PrintYieldProperties();
}

// print just yield properties to output window
void VonMisesHardening::PrintYieldProperties(void)
{
	PrintProperty("yld",yield,"");
	if(linearHardening)
		PrintProperty("Ep",Ep,"");
	else
	{	PrintProperty("K",beta,"");
		PrintProperty("n",npow,"");
	}
    cout << endl;
}

// Read material properties
char *VonMisesHardening::InputMat(char *xName,int &input)
{
    // Ep: plastic modulus, i.e. slope of unidirectional stress-plastic strain curve
    if(strcmp(xName,"Ep")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Ep);
    }

    // ET: Tangential modulus of unidirectional stress-plastic strain curve
    else if(strcmp(xName,"ET")==0)
    {   input=DOUBLE_NUM;
        return((char *)&ET);
    }

    // Khard: coefficient of plastic strains for non-linear hardening (beta)
    else if(strcmp(xName,"Khard")==0)
    {   input=DOUBLE_NUM;
		linearHardening=FALSE;
        return((char *)&beta);
    }
	
    // mhard: power in non-linear hardening
    else if(strcmp(xName,"nhard")==0)
    {   input=DOUBLE_NUM;
		linearHardening=FALSE;
        return((char *)&npow);
    }
	
    return(IsoPlasticity::InputMat(xName,input));
}

// Constant reduced properties used in constitutive law
void VonMisesHardening::InitialLoadMechProps(int makeSpecific,int np)
{
    // Get plastic modulus if needed
	// if got ET, then override/calculate Ep
    if(ET>=0.) Ep=E*ET/(E-ET);
	
	// reduced properties
    Epred=Ep*1.e6/rho;
	
	// beta and npow (if used) are dimensionless
	// if either is entered, Ep and ET are ignored
	
	IsoPlasticity::InitialLoadMechProps(makeSpecific,np);
}

#pragma mark VonMisesHardening::Methods

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
double VonMisesHardening::GetYield(MPMBase *mptr,int np,double delTime)
{
	return linearHardening ? yldred + Epred*alpint : yldred*pow(1.+beta*alpint,npow) ;
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
double VonMisesHardening::GetKPrime(MPMBase *mptr,int np,double delTime)
{
	return linearHardening ? 2.*Epred/3. : 2.*yldred*beta*npow*pow(1.+beta*alpint,npow-1)/3. ;
}

// Get derivative of (1./3.)*yield^2 with respect to lambda for plane stress only
// ... and using dep/dlambda = sqrt(2./3.)*fnp1 where ep=alpint
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)*lambda*fnp1 or depdot/dlambda = sqrt(2./3.)*fnp1/delTime
double VonMisesHardening::GetK2Prime(MPMBase *mptr,double fnp1,double delTime)
{
	return linearHardening ? sqrt(8./27.)*(yldred + Epred*alpint)*Epred*fnp1 : 
								sqrt(8./27.)*yldred*yldred*beta*npow*pow(1.+beta*alpint,2.*npow-1)*fnp1;
}

// Solve this linear model in closed form for lambdak and update trial alpha
double VonMisesHardening::SolveForLambda(MPMBase *mptr,int np,double strial,Tensor *stk,double delTime)
{
	// plane strain is numerical
	if(np==PLANE_STRESS_MPM || !linearHardening)
		return IsoPlasticity::SolveForLambda(mptr,np,strial,stk,delTime);
		
	// closed form for plane strain and 3D and linear hardening
	double lambdak = (strial - sqrt(2./3.)*(yldred + Epred*alpint))/(2.*(Gred + Epred/3.));
	UpdateTrialAlpha(mptr,np,lambdak,(double)1.);
	return lambdak;
}

#pragma mark VonMises::Accessors

// Return the material tag
int VonMisesHardening::MaterialTag(void) { return VONMISESHARDENING; }

// return material type
const char *VonMisesHardening::MaterialType(void) { return "Von Mises Hardening Plastic"; }

