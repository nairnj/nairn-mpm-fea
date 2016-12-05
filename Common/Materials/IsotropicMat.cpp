/********************************************************************************
    IsotropicMat.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Jan 31 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/IsotropicMat.hpp"
#include "System/UnitsController.hpp"

#pragma mark IsotropicMat::Constructors and Destructors

// Constructors
IsotropicMat::IsotropicMat() {}

// Constructors
IsotropicMat::IsotropicMat(char *matName) : Elastic(matName)
{
    int i;
    
    aI = 40;
    E = -1;
    G = -1.;
    nu = -1.;
    
    for(i=0;i<ISO_PROPS;i++)
        read[i]=0;
}

#pragma mark IsotropicMat::Initialization

// Read material properties
char *IsotropicMat::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;

#ifdef MPM_CODE
    if(strcmp(xName,"E")==0)
    {	read[E_PROP]=1;
		return UnitsController::ScaledPtr((char *)&E,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"G")==0)
    {	read[G_PROP]=1;
		return UnitsController::ScaledPtr((char *)&G,gScaling,1.e6);
    }
#else
    if(strcmp(xName,"E")==0)
    {	read[E_PROP]=1;
		return (char *)&E;
    }
    
    else if(strcmp(xName,"G")==0)
    {	read[G_PROP]=1;
		return (char *)&G;
    }
#endif
    else if(strcmp(xName,"nu")==0)
    {	read[NU_PROP]=1;
        return((char *)&nu);
    }
    
    else if(strcmp(xName,"alpha")==0)
        return((char *)&aI);

    return Elastic::InputMaterialProperty(xName,input,gScaling);
}

// print to output window
void IsotropicMat::PrintMechanicalProperties(void) const
{
#ifdef MPM_CODE
	PrintProperty("E",E*UnitsController::Scaling(1.e-6),"");
	PrintProperty("v",nu,"");
	PrintProperty("G",G*UnitsController::Scaling(1.e-6),"");
#else
	PrintProperty("E",E,"");
	PrintProperty("v",nu,"");
	PrintProperty("G",G,"");
#endif
	cout << endl;
	
	PrintProperty("a",aI,"");
    cout << endl;
}

// calculate properties used in analyses
const char *IsotropicMat::VerifyAndLoadProperties(int np)
{
    // finish input and verify all there
    if(!read[G_PROP])
    {	G = E/(2.*(1.+nu));
        read[G_PROP]=1;
    }
    else if(!read[E_PROP])
    {	E = 2.*G*(1.+nu);
        read[E_PROP]=1;
    }
    else if(!read[NU_PROP])
    {	nu = E/(2.*G)-1.;
        read[NU_PROP]=1;
    }
    else
		return "E, nu, and G all specified. Only two allowed";
		
    for(int i=0;i<ISO_PROPS;i++)
    {	if(!read[i])
			return "A required material property is missing";
    }
    
    // analysis properties
    const char *err=SetAnalysisProps(np,E,E,E,nu,nu,nu,G,G,G,
							1.e-6*aI,1.e-6*aI,1.e-6*aI,
							betaI*concSaturation,betaI*concSaturation,betaI*concSaturation);
	if(err!=NULL) return err;
	
	// load elastic properties with constant values
	FillUnrotatedElasticProperties(&pr,np);
	
#ifdef MPM_CODE
	// heating gamma0 (dimensionless) (K 3 alpha)/(rho Cv)
    double alphaV = 3.e-6*aI;
	double Kbulk = E/(3.*(1-2*nu));
    gamma0 = Kbulk*alphaV/(rho*heatCapacity);
#endif
	
	// superclass call (skip Elastic, which has no VerifyAndLoadProperties()
	return MaterialBase::VerifyAndLoadProperties(np);
}

#pragma mark IsotropicMat::Accessors

// return material type
const char *IsotropicMat::MaterialType(void) const { return "Isotropic"; }

