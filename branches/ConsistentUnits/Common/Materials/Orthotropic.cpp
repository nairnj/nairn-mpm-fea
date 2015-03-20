/********************************************************************************
    Orthotropic.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Tue Jan 28 2003.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/Orthotropic.hpp"
#ifdef MPM_CODE
	#include "Custom_Tasks/ConductionTask.hpp"
	#include "Custom_Tasks/DiffusionTask.hpp"
#endif
#include "System/UnitsController.hpp"
#include "System/UnitsController.hpp"

#pragma mark Orthotropic::Constructors and Destructors

// Constructors
Orthotropic::Orthotropic() {}

// Constructors
Orthotropic::Orthotropic(char *matName) : TransIsotropic(matName,ORTHO)
{
#ifdef MPM_CODE
	Dz=0.;
	kCondz=0.;
#endif
	betax=0.;
	betay=0.;
	betaz=0.;
}

#pragma mark Orthotropic::Initialization

// print mechanical properties to output window
void Orthotropic::PrintMechanicalProperties(void) const
{	
    PrintProperty("E1",Ex*UnitsController::Scaling(1.e-6),"");
	PrintProperty("E2",Ey*UnitsController::Scaling(1.e-6),"");
	PrintProperty("E3",Ez*UnitsController::Scaling(1.e-6),"");
	PrintProperty("v12",nuxy,"");
    cout << endl;
    
	PrintProperty("v13",nuxz,"");
	PrintProperty("v23",nuyz,"");
    PrintProperty("G12",Gxy*UnitsController::Scaling(1.e-6),"");
	PrintProperty("G13",Gxz*UnitsController::Scaling(1.e-6),"");
    cout << endl;
    
	PrintProperty("G23",Gyz*UnitsController::Scaling(1.e-6),"");
    cout << endl;
	
    PrintProperty("a1",ax,"");
	PrintProperty("a2",ay,"");
	PrintProperty("a3",az,"");
    cout << endl;
}
    
#ifdef MPM_CODE
// print transport properties to output window
void Orthotropic::PrintTransportProperties(void) const
{
    char mline[200];
	
	// Diffusion constants
	if(DiffusionTask::active)
	{   sprintf(mline,"D1 =%12.3g   D2 =%12.3g   D3 =%12.3g mm^2/sec  csat = %9.5lf",diffT,diffA,Dz,concSaturation);
		cout << mline << endl;
	    sprintf(mline,"b1 =%12.6g   b2 =%12.6g   b3 =%12.6g 1/wt fr",betax,betay,betaz);
		cout << mline << endl;
	}
	// Conductivity constants
	if(ConductionTask::active)
	{   sprintf(mline,"k1 =%12.3g   k2 =%12.3g   k3 =%12.3g %s\nC   =%12.3g %s",
                    rho*kCondT*UnitsController::Scaling(1.e-6),
					rho*kCondA*UnitsController::Scaling(1.e-6),
					rho*kCondz*UnitsController::Scaling(1.e-6),UnitsController::Label(CONDUCTIVITY_UNITS),
					heatCapacity*UnitsController::Scaling(1.e-6),UnitsController::Label(HEATCAPACITY_UNITS));
		cout << mline << endl;
	}
}
#endif

// Read material properties
char *Orthotropic::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"Ex")==0 || strcmp(xName,"ER")==0)
    {	read[EX_PROP]=1;
        return UnitsController::ScaledPtr((char *)&Ex,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"Ey")==0 || strcmp(xName,"EZ")==0)
    {	read[EY_PROP]=1;
        return UnitsController::ScaledPtr((char *)&Ey,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"Ez")==0 || strcmp(xName,"ET")==0)
    {	read[EZ_PROP]=1;
        return UnitsController::ScaledPtr((char *)&Ez,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"Gxy")==0 || strcmp(xName,"Gyx")==0 || strcmp(xName,"GRZ")==0 || strcmp(xName,"GZR")==0)
    {	read[GXY_PROP]=1;
        return UnitsController::ScaledPtr((char *)&Gxy,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"Gxz")==0 || strcmp(xName,"Gzx")==0 || strcmp(xName,"GRT")==0 || strcmp(xName,"GTR")==0)
    {	read[GXZ_PROP]=1;
        return UnitsController::ScaledPtr((char *)&Gxz,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"Gyz")==0 || strcmp(xName,"Gzy")==0 || strcmp(xName,"GZT")==0 || strcmp(xName,"GTZ")==0)
    {	read[GYZ_PROP]=1;
        return UnitsController::ScaledPtr((char *)&Gyz,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"nuyx")==0 || strcmp(xName,"nuZR")==0)
    {	read[NUYX_PROP]=1;
        return((char *)&nuyx);
    }
    
    else if(strcmp(xName,"nuxy")==0 || strcmp(xName,"nuRZ")==0)
    {	read[NUXY_PROP]=1;
        return((char *)&nuxy);
    }
    
    else if(strcmp(xName,"nuyz")==0 || strcmp(xName,"nuZT")==0)
    {	read[NUYZ_PROP]=1;
        return((char *)&nuyz);
    }
    
    else if(strcmp(xName,"nuzy")==0 || strcmp(xName,"nuTZ")==0)
    {	read[NUZY_PROP]=1;
        return((char *)&nuzy);
    }
    
    else if(strcmp(xName,"nuxz")==0 || strcmp(xName,"nuRT")==0)
    {	read[NUXZ_PROP]=1;
        return((char *)&nuxz);
    }
    
    else if(strcmp(xName,"nuzx")==0 || strcmp(xName,"nuTR")==0)
    {	read[NUZX_PROP]=1;
        return((char *)&nuzx);
    }
    
    else if(strcmp(xName,"alphax")==0 || strcmp(xName,"alphaR")==0)
    {	read[AX_PROP]=1;
        return((char *)&ax);
    }
    
    else if(strcmp(xName,"alphay")==0 || strcmp(xName,"alphaZ")==0)
    {	read[AY_PROP]=1;
        return((char *)&ay);
    }
    
    else if(strcmp(xName,"alphaz")==0 || strcmp(xName,"alphaT")==0)
    {	read[AZ_PROP]=1;
        return((char *)&az);
    }

#ifdef MPM_CODE
    else if(strcmp(xName,"betax")==0 || strcmp(xName,"betaR")==0)
        return((char *)&betax);
    
    else if(strcmp(xName,"betay")==0 || strcmp(xName,"betaZ")==0)
        return((char *)&betay);
    
    else if(strcmp(xName,"betaz")==0 || strcmp(xName,"betaT")==0)
        return((char *)&betaz);

    else if(strcmp(xName,"Dx")==0 || strcmp(xName,"DR")==0)
        return((char *)&diffT);
		
    else if(strcmp(xName,"Dy")==0 || strcmp(xName,"DZ")==0)
        return((char *)&diffA);
		
    else if(strcmp(xName,"Dz")==0 || strcmp(xName,"DT")==0)
        return((char *)&Dz);
		
    else if(strcmp(xName,"kCondx")==0 || strcmp(xName,"kCondR")==0)
		return UnitsController::ScaledPtr((char *)&kCondT,gScaling,1.e6);
		
    else if(strcmp(xName,"kCondy")==0 || strcmp(xName,"kCondZ")==0)
		return UnitsController::ScaledPtr((char *)&kCondA,gScaling,1.e6);
		
    else if(strcmp(xName,"kCondz")==0 || strcmp(xName,"kCondT")==0)
		return UnitsController::ScaledPtr((char *)&kCondz,gScaling,1.e6);
#endif
		
	return(MaterialBase::InputMaterialProperty(xName,input,gScaling));
}

// calculate properties used in analyses
const char *Orthotropic::VerifyAndLoadProperties(int np)
{
    // finish input
    if(!read[NUXY_PROP])
    {	nuxy=nuyx*Ex/Ey;
        read[NUXY_PROP]=1;
    }
    else if(!read[NUYX_PROP])
    {	nuyx=nuxy*Ey/Ex;
        read[NUYX_PROP]=1;
    }
    else
		return "nuxy and nuyx both given. Only one is allowed";
    
    if(!read[NUXZ_PROP])
    {	nuxz=nuzx*Ex/Ez;
        read[NUXZ_PROP]=1;
    }
    else if(!read[NUZX_PROP])
    {	nuzx=nuxz*Ez/Ex;
        read[NUZX_PROP]=1;
    }
    else
		return "nuxz and nuzx both given. Only one is allowed";
    
    if(!read[NUYZ_PROP])
    {	nuyz=nuzy*Ey/Ez;
        read[NUYZ_PROP]=1;
    }
    else if(!read[NUZY_PROP])
    {	nuzy=nuyz*Ez/Ey;
        read[NUZY_PROP]=1;
    }
    else
		return "nuyz and nuzy both given. Only one is allowed";
    
    int i;
    for(i=0;i<ORTHO_PROPS;i++)
    {	if(!read[i])
			return "A required material property is missing";
    }

#ifdef MPM_CODE
    // make conductivty specific (N mm^3/(sec-K-g))
    kCondz /= rho;

    // set properties
    const char *err=SetAnalysisProps(np,Ex,Ey,Ez,nuxy,nuxz,nuyz,
						Gxy,Gxz,Gyz,1.e-6*ax,1.e-6*ay,1.e-6*az,
						betax*concSaturation,betay*concSaturation,betaz*concSaturation);
#else
	double rescale = UnitsController::Scaling(1.e-6);
    const char *err=SetAnalysisProps(np,rescale*Ex,rescale*Ey,rescale*Ez,nuxy,nuxz,nuyz,
						rescale*Gxy,rescale*Gxz,rescale*Gyz,1.e-6*ax,1.e-6*ay,1.e-6*az,
						betax*concSaturation,betay*concSaturation,betaz*concSaturation);
#endif
	if(err!=NULL) return err;
	
	// superclass call
	return MaterialBase::VerifyAndLoadProperties(np);
}

#pragma mark Orthotropic::Accessors

// Return the material tag
int Orthotropic::MaterialTag(void) const { return ORTHO; }

// return material type
const char *Orthotropic::MaterialType(void) const { return "Orthotropic (3 axis normal to x-y plane)"; }

#ifdef MPM_CODE

// calculate maximum wave speed in mm/sec (moduli in MPa, rho in g/mm^3)
double Orthotropic::WaveSpeed(bool threeD,MPMBase *mptr) const
{
    double xx,cnorm,cshear;
    
    xx=1.-nuxy*nuyx-nuxz*nuzx-nuyz*nuzy-2.*nuxy*nuyz*nuzx;
    cnorm=fmax(Ex*(1.-nuyz*nuzy),Ey*(1.-nuxz*nuzx))/xx;
	if(threeD) cnorm=fmax(cnorm,Ez*(1-nuxy*nuyx)/xx);
	cshear = threeD ? fmax(Gxy,fmax(Gxz,Gyz)) : Gxy;
    return sqrt(fmax(cnorm,cshear)/rho);
}

// diffusion and conductivity in the z direction
double Orthotropic::GetDiffZ(void) const { return Dz; }
double Orthotropic::GetKcondZ(void) const { return kCondz; }

#endif





