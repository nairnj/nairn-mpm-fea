/********************************************************************************
	OrthoSoftening.cpp
	nairn-mpm-fea
 
	Created by John Nairn on 30 JAN 2018
	Copyright (c) 2018 John A. Nairn, All rights reserved.
 
 	Next Tasks
  		Java NFMViz - Support orthosoftening, softeningXX, etc
********************************************************************************/

#include "stdafx.h"
#include "Materials/OrthoSoftening.hpp"
#include "Materials/OrthoFailureSurface.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "System/UnitsController.hpp"
#include "Materials/LinearSoftening.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"

#pragma mark OrthoSoftening::Constructors and Destructors

// Constructor
// throws std::bad_alloc
OrthoSoftening::OrthoSoftening(char *matName,int matID) : TransIsoSoftening(matName,matID)
{
	Dz=0.;
	kCondz=0.;
	betax=0.;
	betay=0.;
	betaz=0.;
	
	// We need to change initiation law
	delete initiationLaw;
	initiationLaw = (FailureSurface *)new OrthoFailureSurface(this);

	// We define nine softening laws
	//   XX, YY, ZZ, XYX, XYY, XZX, XZZ, YZY, and YZZ
	// We store XX and YY and AI and I of parent class
	//	 and store XYX and XYY in AII and TII of parent class
	//   and store YZY in II of parent class
	// ZZ, YZZ, XZX, and XZZ are four new law
	softeningZZ = new LinearSoftening();
	softeningYZZ = new LinearSoftening();
	softeningXZX = new LinearSoftening();
	softeningXZZ = new LinearSoftening();
}

// Read material properties
// Because this material inherits from TransIsoSoftening and does not have Orthotropic in its
// 	super classes, it needs to repeat some parts of the Orthotropic material code
char *OrthoSoftening::InputMaterialProperty(char *xName,int &input,double &gScaling)
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
	
	// four new softening laws
	else if(strcmp(xName,"SofteningZZ")==0)
	{	input = SOFTZZ_LAW_SELECTION;
		return (char *)this;
	}
	else if(strcmp(xName,"SofteningYZZ")==0)
	{	input = SOFTYZZ_LAW_SELECTION;
		return (char *)this;
	}
	else if(strcmp(xName,"SofteningXZX")==0)
	{	input = SOFTXZX_LAW_SELECTION;
		return (char *)this;
	}
	else if(strcmp(xName,"SofteningXZZ")==0)
	{	input = SOFTXZZ_LAW_SELECTION;
		return (char *)this;
	}
	
	// softening law properties
	if(strlen(xName)>3)
	{	char *ptr;
		if(xName[0]=='Z' && xName[1]=='Z' && xName[2]=='-')
		{	ptr = softeningZZ->InputSofteningProperty(&xName[3],input,gScaling);
			if(ptr != NULL) return ptr;
		}
		else if(xName[0]=='X' && xName[1]=='X' && xName[2]=='-')
		{	ptr = softeningXX()->InputSofteningProperty(&xName[3],input,gScaling);
			if(ptr != NULL) return ptr;
		}
		else if(xName[0]=='Y' && xName[1]=='Y' && xName[2]=='-')
		{	ptr = softeningYY()->InputSofteningProperty(&xName[3],input,gScaling);
			if(ptr != NULL) return ptr;
		}
		if(strlen(xName)>4)
		{	if(xName[0]=='X' && xName[1]=='Z' && xName[2]=='X' && xName[3]=='-')
			{	ptr = softeningXZX->InputSofteningProperty(&xName[4],input,gScaling);
				if(ptr != NULL) return ptr;
			}
			else if(xName[0]=='X' && xName[1]=='Z' && xName[2]=='Z' && xName[3]=='-')
			{	ptr = softeningXZZ->InputSofteningProperty(&xName[4],input,gScaling);
				if(ptr != NULL) return ptr;
			}
			else if(xName[0]=='X' && xName[1]=='Y' && xName[2]=='X' && xName[3]=='-')
			{	ptr = softeningXYX()->InputSofteningProperty(&xName[4],input,gScaling);
				if(ptr != NULL) return ptr;
			}
			else if(xName[0]=='X' && xName[1]=='Y' && xName[2]=='Y' && xName[3]=='-')
			{	ptr = softeningXYY()->InputSofteningProperty(&xName[4],input,gScaling);
				if(ptr != NULL) return ptr;
			}
			else if(xName[0]=='Y' && xName[1]=='Z' && xName[2]=='Y' && xName[3]=='-')
			{	ptr = softeningYZY()->InputSofteningProperty(&xName[4],input,gScaling);
				if(ptr != NULL) return ptr;
			}
			else if(xName[0]=='Y' && xName[1]=='Z' && xName[2]=='Z' && xName[3]=='-')
			{	ptr = softeningYZZ->InputSofteningProperty(&xName[4],input,gScaling);
				if(ptr != NULL) return ptr;
			}
		}
	}
    
    // block unused transversely isotropic properties
    // some "T" properties used above to mean tangential in cylindrical orthotropy
    if(strcmp(xName,"EA")==0  || strcmp(xName,"GT")==0 || strcmp(xName,"GA")==0 || strcmp(xName,"nuT")==0 ||
       strcmp(xName,"nuA")==0 || strcmp(xName,"alphaA")==0 || strcmp(xName,"betaA")==0 || strcmp(xName,"DA")==0 ||
       strcmp(xName,"kCondA")==0)
    {   return (char *)NULL;
    }
	
	return TransIsoSoftening::InputMaterialProperty(xName,input,gScaling);
}

// Allows selected initiation laws
bool OrthoSoftening::AcceptInitiationLaw(FailureSurface *iLaw,int lawID)
{
	// only allows max principle stress
	if(lawID != ORTHOFAILURESURFACE) return false;
	
	delete initiationLaw;
	initiationLaw = iLaw;
	return true;
}

// Softening laws in parent call
//   Fxx for XX (stored in softeningAI) with strength in sigmaXXc (in sigmaAc)
//   Fyy for YY (stored in softeningI) with strength in sigmaYYc (in sigmac)
//   Fxyx for XYX (stored in softeningAII) with strength in tauXY-X (in tauAc) (normal in y direction)
//   Fxyy for XYY (stored in softeningTII) with strength in tauXY-Y (in tauTc) (normal in x direction)
//   Fyzy for YZY (stored in softeningII) with strength in tauYZY (in tauc) (normal in z direction)
// Additional laws for 3D
//   Fzz for ZZ with strength in sigmaZZc
//   Fyzz for YZZ with strength in tauYZ-Z (normal in y direction)
//   Fxzx for XZX with strength in tauXZ-X (normal in z direction)
//   Fxzz for XZZ with strength in tauXZ-Z (normal in 1 direction)
bool OrthoSoftening::AcceptSofteningLaw(SofteningLaw *sLaw,int mode)
{
	if(mode==SOFTZZ_LAW_SELECTION)
	{	delete softeningZZ;
		softeningZZ = sLaw;
	}
	else if(mode==SOFTYZZ_LAW_SELECTION)
	{	delete softeningYZZ;
		softeningYZZ = sLaw;
	}
	else if(mode==SOFTXZX_LAW_SELECTION)
	{	delete softeningXZX;
		softeningXZX = sLaw;
	}
	else if(mode==SOFTXZZ_LAW_SELECTION)
	{	delete softeningXZZ;
		softeningXZZ = sLaw;
	}
	else
	{	return TransIsoSoftening::AcceptSofteningLaw(sLaw,mode);
	}
	return true;
}

// verify settings and some initial calculations
const char *OrthoSoftening::VerifyAndLoadProperties(int np)
{
	if(!useLargeRotation)
		return "OrthoSoftening materials require activation of large rotation option";
	
	if(np==PLANE_STRESS_MPM)
		return "OrthoSoftening materials cannot yet be used in plane stress calculations";
    
    // Remap softening laws
    if(swapz==1)
    {   SwapLaws(&softeningAI,&softeningZZ);        // XX with ZZ
        SwapLaws(&softeningXZX,&softeningXZZ);      // XZX with ZXZ/XZZ
        SwapLaws(&softeningTII,&softeningII);       // XYY with ZYY/YZY
        SwapLaws(&softeningAII,&softeningYZZ);      // XYX with ZYZ/YZZ
    }
    else if(swapz>1)
    {   SwapLaws(&softeningI,&softeningZZ);         // YY with ZZ
        SwapLaws(&softeningII,&softeningYZZ);       // YZY with ZYZ/YZZ
        SwapLaws(&softeningTII,&softeningXZZ);      // XYY with XZZ
        SwapLaws(&softeningAII,&softeningXZX);      // XYX with XZX
    }
    
    // sway y and z of initiation law
    initiationLaw->RemapProperties(swapz);
	
	// call initiation law that is used
	const char *ptr = initiationLaw->VerifyAndLoadProperties(np);
	if(ptr != NULL) return ptr;
	
	// check in superclass (along with its initialization)
	const char *err = VerifyAndLoadOrthoProperties(np);
	if(err!=NULL) return err;
	
	if(softenStatsMode<VARY_STRENGTH || softenStatsMode>VARY_STRENGTH_AND_TOUGHNESS)
		return "Invalid option selected for which softening properties to vary";
	
	// verify distibution parameters
	if(distributionMode==SOFTDIST_NORMAL)
	{	if(softenCV<0.)
			return "Coefficient of variation cannot be negative";
		else if(softenCV==0.)
			distributionMode=SOFTDIST_NONE;
	}
	else if(distributionMode==SOFTDIST_WEIBULL)
	{	if(wAlpha<=0. || wV0<=0.)
			return "Weibull scale and reference volume must be greater than 0";
		wGam1A = gamma_fxn(1.+1./wAlpha);
		softenCV = sqrt(gamma_fxn(1.+2./wAlpha)/(wGam1A*wGam1A) - 1.);
	}

    // make sure have a valid setting
    if(tractionFailureSurface==CYLINDER_SURFACE)
    {    return "Cylindrical failure surface not allowed for OrthoSoftening material";
    }
	else if(tractionFailureSurface==COUPLED_CUBOID_SURFACE)
	{	// treat same as OVOID
		tractionFailureSurface = OVOID_SURFACE;
	}

#ifdef POROELASTICITY
    // poroelasticity conversions
    if(DiffusionTask::HasPoroelasticity())
    {    return "OrthoSoftening does not implement poroelasticity";
    }
    else
    {    // diffusion CT is 1
        diffusionCT = 1.;
    }
#else
    // diffusion CT is 1
    //diffusionCT = 1.;
#endif
    
	// no parent class actions used
	return NULL;
}

// calculate properties used in analyses
const char *OrthoSoftening::VerifyAndLoadOrthoProperties(int np)
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
	
    // reorder all terms if needed
    if(swapz==1)
    {   // swap x amd z
        SwapProperties(Ex,Ez,Gyz,Gxy);
        SwapProperties(nuxz,nuzx,nuxy,nuzy);
        SwapProperties(nuyz,nuyx,ax,az);
        SwapProperties(betax,betaz,diffT,Dz);
        SwapProperties(kCondT,kCondz);
#ifdef POROELASTICITY
        // not supported in this material
        //SwapProperties(alphaTPE,alphazPE,DarcyT,Darcyz);
#endif
    }
    else if(swapz>1)
    {   SwapProperties(Ey,Ez,Gxz,Gxy);
        SwapProperties(nuxz,nuxy,nuzx,nuyx);
        SwapProperties(nuyz,nuzy,ay,az);
        SwapProperties(betay,betaz,diffA,Dz);
        SwapProperties(kCondA,kCondz);
#ifdef POROELASTICITY
        // not supported in this material
        //SwapProperties(alphaAPE,alphazPE,DarcyA,Darcyz);
#endif
    }

    // make conductivty specific (N mm^3/(sec-K-g))
	kCondA /= rho;
	kCondT /= rho;
	kCondz /= rho;

#ifdef POROELASTICITY
	// poroelasticity conversions
	if(DiffusionTask::HasPoroelasticity())
	{	return "OrthoSoftening does not implement poroelasticity";
	}
	else
	{	// diffusion CT is 1
		diffusionCT = 1.;
	}
#else
	// diffusion CT is 1
	//diffusionCT = 1.;
#endif
	
	// set properties
	const char *err=SetAnalysisProps(np,Ex,Ey,Ez,nuxy,nuxz,nuyz,
									 Gxy,Gxz,Gyz,1.e-6*ax,1.e-6*ay,1.e-6*az,
									 betax*concSaturation,betay*concSaturation,betaz*concSaturation);
	if(err!=NULL) return err;
	
	// load elastic properties with constant values
	FillUnrotatedElasticProperties(&pr,np);
	
	// superclass call (but skip over TransIsotropic)
	return MaterialBase::VerifyAndLoadProperties(np);
}

// print mechanical properties to the results
void OrthoSoftening::PrintMechanicalProperties(void) const
{
	// Orthotropic material properties
	PrintProperty("Exx",Ex*UnitsController::Scaling(1.e-6),"");
	PrintProperty("Eyy",Ey*UnitsController::Scaling(1.e-6),"");
	PrintProperty("Ezz",Ez*UnitsController::Scaling(1.e-6),"");
	PrintProperty("vxy",nuxy,"");
	cout << endl;
	
	PrintProperty("vxz",nuxz,"");
	PrintProperty("vyz",nuyz,"");
	PrintProperty("Gxy",Gxy*UnitsController::Scaling(1.e-6),"");
	PrintProperty("Gxz",Gxz*UnitsController::Scaling(1.e-6),"");
	cout << endl;
	
	PrintProperty("Gyz",Gyz*UnitsController::Scaling(1.e-6),"");
    PrintProperty("vyx",nuyx,"");
    PrintProperty("vzx",nuzx,"");
    PrintProperty("vzy",nuzy,"");
	cout << endl;
	
	PrintProperty("ax",ax,"");
	PrintProperty("ay",ay,"");
	PrintProperty("az",az,"");
	cout << endl;

	cout << "Failure Surface: ";
	initiationLaw->PrintInitiationProperties();
	cout << "XX Mode I Softening: ";
	softeningAI->PrintSofteningProperties(rho*((OrthoFailureSurface *)initiationLaw)->sigmaXX());
	cout << "YY Mode I Softening: ";
	softeningI->PrintSofteningProperties(rho*((OrthoFailureSurface *)initiationLaw)->sigmaYY());
	cout << "XY-X Mode II Softening: ";
	softeningAII->PrintSofteningProperties(rho*((OrthoFailureSurface *)initiationLaw)->tauXYX());
	cout << "XY-Y Mode II Softening: ";
	softeningTII->PrintSofteningProperties(rho*((OrthoFailureSurface *)initiationLaw)->tauXYY());
	if(fmobj->IsThreeD())
	{	cout << "ZZ Mode I Softening: ";
		softeningZZ->PrintSofteningProperties(rho*((OrthoFailureSurface *)initiationLaw)->sigmaZZ());
		cout << "XZ-X Mode II Softening: ";
		softeningXZX->PrintSofteningProperties(rho*((OrthoFailureSurface *)initiationLaw)->tauXZX());
		cout << "XZ-Z Mode II Softening: ";
		softeningXZZ->PrintSofteningProperties(rho*((OrthoFailureSurface *)initiationLaw)->tauXZZ());
		cout << "YZ-Y Mode II Softening: ";
		softeningII->PrintSofteningProperties(rho*((OrthoFailureSurface *)initiationLaw)->tauYZY());
		cout << "YZ-Z Mode II Softening: ";
		softeningII->PrintSofteningProperties(rho*((OrthoFailureSurface *)initiationLaw)->tauYZZ());
	}
	
	if(distributionMode==SOFTDIST_WEIBULL)
	{	cout << "Weibull scale = " << wAlpha << " and V0 = " << wV0 <<
                " " << UnitsController::Label(CULENGTH_UNITS) << "^3" << endl;
	}
	cout << "Coefficient of Variation = " << softenCV << endl;
	if(distributionMode!=SOFTDIST_NONE)
	{   if(softenStatsMode==VARY_STRENGTH_AND_TOUGHNESS)
			cout << "    Vary strength and toughness" << endl;
		else if(softenStatsMode==VARY_TOUGHNESS)
			cout << "    Vary toughness" << endl;
		else
			cout << "    Vary strength" << endl;
	}
	
	// currently only cuboid, but save for future
	initiationLaw->SetFailureSurface(tractionFailureSurface);
	cout << "Traction Failure Surface: ";
	if(tractionFailureSurface == CYLINDER_SURFACE)
		cout << "cylindrical" << endl;
	else if(tractionFailureSurface == CUBOID_SURFACE)
		cout << "cuboid" << endl;
	else if(tractionFailureSurface == OVOID_SURFACE)
		cout << "coupled cuboid" << endl;
	else
		cout << "unsupported failure surface" << endl;
	
	cout << "Post-failure contact: ";
	if(frictionCoeff>0.)
		cout << "Coulomb coefficient of friction = " << frictionCoeff << endl;
	else
		cout << "frictionless" << endl;
}

// print transport properties to output window
void OrthoSoftening::PrintTransportProperties(void) const
{
	char mline[200];
	size_t msize=200;
	
	// Diffusion constants
	if(DiffusionTask::HasFluidTransport())
	{   snprintf(mline,msize,"Dx =%12.3g   Dy =%12.3g   Dz =%12.3g mm^2/sec  csat = %9.5lf",diffT,diffA,Dz,concSaturation);
		cout << mline << endl;
		snprintf(mline,msize,"bx =%12.6g   by =%12.6g   bz =%12.6g 1/wt fr",betax,betay,betaz);
		cout << mline << endl;
	}
	// Conductivity constants
	if(ConductionTask::active)
	{   snprintf(mline,msize,"kx =%12.3g   ky =%12.3g   kz =%12.3g %s\nC   =%12.3g %s",
				rho*kCondT*UnitsController::Scaling(1.e-6),
				rho*kCondA*UnitsController::Scaling(1.e-6),
				rho*kCondz*UnitsController::Scaling(1.e-6),UnitsController::Label(CONDUCTIVITY_UNITS),
				heatCapacity*UnitsController::Scaling(1.e-6),UnitsController::Label(HEATCAPACITY_UNITS));
		cout << mline << endl;
	}
}

#pragma mark TransIsoSoftening::History Data Methods

// Create history variables needed for softening behavior
char *OrthoSoftening::InitHistoryData(char *pchr,MPMBase *mptr)
{
	// Validate this law (in use) is stable on current grid
	//	This validation uses average values. Large stochastic variations might
	//	cause some particles to be unstable. Enhancement could be to validate
	//	using minimum Gc and maximum strength instead
	// Does for all particles, but this makes sure softening material
	//   is being used and works with variable size particles
	
	// particle size
	double delx = mptr->GetMinParticleLength();
	
	// XX mode I failure (sigmaXX/C11)
	double sigMax = ((OrthoFailureSurface *)initiationLaw)->sigmaXX();
	double estr0 = fmobj->IsThreeD() ? sigMax/pr.C[0][0] : sigMax/pr.C[1][1] ;
	VerifyStability(sigMax,estr0,delx,softeningXX(),"XX","sigmaXXc");
	
	// YY mode I failure (sigmaYY/C22)
	sigMax = ((OrthoFailureSurface *)initiationLaw)->sigmaYY();
	estr0 = fmobj->IsThreeD() ? sigMax/pr.C[1][1] : sigMax/pr.C[2][2] ;
	VerifyStability(sigMax,estr0,delx,softeningYY(),"YY","sigmaYYc");
	
	// XY-X mode II failure (tauXYX/C66)
	sigMax = ((OrthoFailureSurface *)initiationLaw)->tauXYX();
	estr0 = fmobj->IsThreeD() ? sigMax/pr.C[5][5] : sigMax/pr.C[3][3] ;
	VerifyStability(sigMax,estr0,delx,softeningXYX(),"XY-X","tauXY-Xc");
	
	// XY-Y mode II failure (tauXYY/C66)
	sigMax = ((OrthoFailureSurface *)initiationLaw)->tauXYY();
	estr0 = fmobj->IsThreeD() ? sigMax/pr.C[5][5] : sigMax/pr.C[3][3] ;
	VerifyStability(sigMax,estr0,delx,softeningXYY(),"XY-Y","tauXY-Yc");
	
	if(fmobj->IsThreeD())
	{	// ZZ mode I failure (sigmaZZ/C33)
		sigMax = ((OrthoFailureSurface *)initiationLaw)->sigmaZZ();
		estr0 = sigMax/pr.C[2][2];
		VerifyStability(sigMax,estr0,delx,softeningZZ,"ZZ","sigmaZZc");
		
		// XZX mode II failure (tauXZX/C55)
		sigMax = ((OrthoFailureSurface *)initiationLaw)->tauXZX();
		estr0 = sigMax/pr.C[4][4];
		VerifyStability(sigMax,estr0,delx,softeningXZX,"XZ-X","tauXZ-Xc");
		
		// XZZ mode II failure (tauXZZ/C55)
		sigMax = ((OrthoFailureSurface *)initiationLaw)->tauXZZ();
		estr0 = sigMax/pr.C[4][4];
		VerifyStability(sigMax,estr0,delx,softeningXZZ,"XZ-Z","tauXZ-Zc");
		
		// YZY mode II failure (tauYZY/C44)
		sigMax = ((OrthoFailureSurface *)initiationLaw)->tauYZY();
		estr0 = sigMax/pr.C[3][3];
		VerifyStability(sigMax,estr0,delx,softeningYZY(),"YZ-Y","tauYZ-Yc");
		
		// YZZ mode II failure (tauYZZ/C44)
		sigMax = ((OrthoFailureSurface *)initiationLaw)->tauYZZ();
		estr0 = sigMax/pr.C[3][3];
		VerifyStability(sigMax,estr0,delx,softeningYZZ,"YZ-Z","tauYZ-Zc");
	}

	// historu variables
	double *p = CreateAndZeroDoubles(pchr,NumberOfHistoryDoubles());
	
	// set damage state to 0.1 instead of zero for help in plotting on a grid
	p[SOFT_DAMAGE_STATE] = 0.1;
	
	// set relative strength and toughness
	// Note that both set to see if used. Future idea would be to use different
	//		stats for strength and toughness
	p[RELATIVE_STRENGTH] = 1.;
	p[RELATIVE_TOUGHNESS] = 1.;
	if(distributionMode==SOFTDIST_NORMAL)
	{	double fract = RandomRange(0.001,0.999);		// .001 to .999, max +/- 3.09 std dev
		double relValue = fmax(1. + softenCV*NormalCDFInverse(fract),0.1);
		if(softenStatsMode&VARY_STRENGTH)
			p[RELATIVE_STRENGTH] = relValue;
		if(softenStatsMode&VARY_TOUGHNESS)
			p[RELATIVE_TOUGHNESS] = relValue;
	}
	else if(distributionMode==SOFTDIST_WEIBULL)
	{	double fract = RandomRange(0.001,0.999);		// .001 to .999 limited range
		double vp = mptr->GetUnscaledVolume();
		double relValue = fmax(pow(-wV0*log(1.-fract)/vp,1./wAlpha)/wGam1A,0.05);
		if(softenStatsMode&VARY_STRENGTH)
			p[RELATIVE_STRENGTH] = relValue;
		if(softenStatsMode&VARY_TOUGHNESS)
			p[RELATIVE_TOUGHNESS] = relValue;
	}
	
	return (char *)p;
}

// reset history data
void OrthoSoftening::ResetHistoryData(char *pchr,MPMBase *mptr)
{	double *p = (double *)pchr;
	double relStrength = p[RELATIVE_STRENGTH];
	double relToughness = p[RELATIVE_TOUGHNESS];
	ZeroDoubles(pchr,NumberOfHistoryDoubles());
	p[SOFT_DAMAGE_STATE] = 0.1;
	p[RELATIVE_STRENGTH] = relStrength;
	p[RELATIVE_TOUGHNESS] = relToughness;
}

#pragma mark OrthoSoftening::Methods

// Convert damage initiation mode and normal to stored setting with information
//      to later get the correct damage tensor
// 3D, norm = ZYZ rotation from material axes to crack normal. For TI softening either (0,pi/2,0)
//		to indicate crack axis normal in axial direction (DX_DAMAGE) (switches z and x) or (theta,0,0) to keep axial direction
//		in z direction (DZ_DAMAGE), but rotate around that axis
// 2D, norm = (cos(theta),sin(theta))
int OrthoSoftening::DecodeDamageInitiation(int np,Vector *norm,int failureMode,double *soft) const
{
	// initiate failure in material axis system with predamageState=0.5
	int initForm = DX_DAMAGE;
	
	if(failureMode == EA_FAILURE)
		soft[SOFT_DAMAGE_STATE] = 0.95;						// xx tension failure
	else if(failureMode == TENSILE_FAILURE)
	{	soft[SOFT_DAMAGE_STATE] = 1.15;						// yy tension failure
		initForm = DY_DAMAGE;
	}
	else if(failureMode == EZZ_FAILURE)
	{	soft[SOFT_DAMAGE_STATE] = 0.75;						// zz tension failure
		initForm = DZ_DAMAGE;
	}
	else if(failureMode == GA_FAILURE)
	{	if((np==THREED_MPM && norm->x>1.) || (np!=THREED_MPM && norm->y>0.99))
		{	soft[SOFT_DAMAGE_STATE] = 1.20;					// xy shear failure, y normal
			initForm = DY_DAMAGE;
		}
		else
			soft[SOFT_DAMAGE_STATE] = 1.0;					// xy shear failure, x normal
	}
	else if(failureMode == GXZ_FAILURE)
	{	if(norm->y>0.99)
		{	soft[SOFT_DAMAGE_STATE] = 0.8;					// xz shear failure, z normal
			initForm = DZ_DAMAGE;
		}
		else
			soft[SOFT_DAMAGE_STATE] = 1.05;					// xz shear failure, x normal
	}
	else
	{	if(norm->x>0.99)
		{	soft[SOFT_DAMAGE_STATE] = 1.25;					// yz shear failure, y normal
			initForm = DY_DAMAGE;
		}
		else
		{	soft[SOFT_DAMAGE_STATE] = 0.85;					// yz shear failure, z normal
			initForm = DZ_DAMAGE;
		}
	}

	return initForm;
}

// Fill d with properties in the crack axis system depending on the type of damage
// tensor as specified by DForm
void OrthoSoftening::LoadCrackAxisProperties(int np,CrackAxisProperties *d,int DForm,ElasticProperties *p) const
{
	// some for 3D only are cleared
	d->C44=0.;
	d->C55=0.;
	d->fXZLaw=NULL;
	d->vzxc=0.;
	d->vzyc=0.;
	
	if(DForm==DZ_DAMAGE)
	{	// crack axis normal in material z direction, swap x and z
		// 1 = Z, 2 = Y, 3 = X (only occurs in 3D)
		d->C11 = p->C[2][2];			// Czz
		d->C22 = p->C[1][1];			// Cyy
		d->C33 = p->C[0][0];			// Czz
		d->C12 = p->C[1][2];			// Cyz
		d->C13 = p->C[0][2];			// Cxz
		d->C23 = p->C[0][1];			// Cxy
		d->C44 = p->C[5][5];			// Gxy=C66
		
		// laws
		d->fnLaw = softeningZZ;
		d->sigNc = ((OrthoFailureSurface *)initiationLaw)->sigmaZZ();
		
		d->fXYLaw = softeningYZY();
		d->tauXYc = ((OrthoFailureSurface *)initiationLaw)->tauYZY();
		d->C66 = p->C[3][3];            // Gyz=C44
		
		if(np==THREED_MPM)
		{	d->fXZLaw = softeningXZX;
			d->tauXZc = ((OrthoFailureSurface *)initiationLaw)->tauXZX();
			d->C55 = p->C[4][4];        // Gxz=C55
		}
	}
	else if(DForm==DX_DAMAGE)
	{	// crack axis normal in material x direction, keep original form
		// 1 = X, 2 = Y, 3 = Z
		if(np==THREED_MPM)
		{	d->C11 = p->C[0][0];			// Cxx
			d->C22 = p->C[1][1];			// Cyy
			d->C33 = p->C[2][2];			// Czz
			d->C12 = p->C[0][1];			// Cxy
			d->C13 = p->C[0][2];			// Cxz
			d->C23 = p->C[1][2];			// Cyz
			d->C44 = p->C[3][3];			// Gyz=C44
		}
		else
		{	d->C11 = p->C[1][1];			// Cxx
			d->C22 = p->C[2][2];			// Cyy
			d->C12 = p->C[1][2];			// Cxy
			d->C13 = p->C[4][1];			// Cxz for sigmaxx in dsig.zz
			d->C23 = p->C[4][2];			// Cyz for sigmayy in dsig.zz
			d->C33 = p->C[4][4];			// Czz for ezzr in dsig.zz
			d->vzxc = p->alpha[5];			// for dsig.zz
			d->vzyc = p->alpha[6];			// for dsig.zz
		}
		
		// laws
		d->fnLaw = softeningXX();
		d->sigNc = ((OrthoFailureSurface *)initiationLaw)->sigmaXX();
		
		d->fXYLaw = softeningXYY();
		d->tauXYc = ((OrthoFailureSurface *)initiationLaw)->tauXYY();
		d->C66 = (np==THREED_MPM) ? p->C[5][5] : p->C[3][3];        // Gxy==C66

		if(np == THREED_MPM)
		{	d->fXZLaw = softeningXZZ;
			d->tauXZc = ((OrthoFailureSurface *)initiationLaw)->tauXZZ();
			d->C55 = p->C[4][4];			// Gxz=C55
		}
	}
	else
	{	// crack axis normal in material y direction, swap x and y
		// 1 = Y, 2 = X, 3 = Z
		if(np==THREED_MPM)
		{	d->C11 = p->C[1][1];			// Cyy
			d->C22 = p->C[0][0];			// Cxx
			d->C33 = p->C[2][2];			// Czz
			d->C12 = p->C[0][1];			// Cxy
			d->C13 = p->C[1][2];			// Cyz
			d->C23 = p->C[0][2];			// Cxz
			d->C44 = p->C[4][4];			// Gxz=C55
		}
		else
		{	d->C11 = p->C[2][2];			// Cyy
			d->C22 = p->C[1][1];			// Cxx
			d->C12 = p->C[1][2];			// Cxy
			d->C13 = p->C[4][2];			// Cxz for sigmaxx in dsig.zz
			d->C23 = p->C[4][1];			// Cyz for sigmayy in dsig.zz
			d->C33 = p->C[4][4];			// Czz for ezzr in dsig.zz
			d->vzxc = p->alpha[6];			// for dsig.zz
			d->vzyc = p->alpha[5];			// for dsig.zz
		}
		
		// laws
		d->fnLaw = softeningYY();
		d->sigNc = ((OrthoFailureSurface *)initiationLaw)->sigmaYY();
		
		d->fXYLaw = softeningXYX();
		d->tauXYc = ((OrthoFailureSurface *)initiationLaw)->tauXYX();
		d->C66 = (np==THREED_MPM) ? p->C[5][5] : p->C[3][3];    // Gxy=C66
		
		if(np == THREED_MPM)
		{	d->fXZLaw = softeningYZZ;
			d->tauXZc = ((OrthoFailureSurface *)initiationLaw)->tauYZZ();
			d->C55 = p->C[3][3];			// Gyz=C44
		}
	}
	
	// ratios
	d->C12C11 = d->C12/d->C11;
	d->C13C11 = d->C13/d->C11;
}

// Calculate rotation matrix from crack to material axes
bool OrthoSoftening::GetRToCrack(Matrix3 *R,double *soft, bool is2D, int Dstyle) const
{	// none if undamaged
	if(soft[SOFT_DAMAGE_STATE] < predamageState) return false;
	
	// 3D or 2D
	if(!is2D)
	{	if(Dstyle==DY_DAMAGE)
		{	// (pi/2,0,0) swap y and x
			R->set(0.,-1.,0.,  1.,0., 0.,  0.,0.,1.);
		}
		else if(Dstyle==DZ_DAMAGE)
		{	// (0,pi/2,0) swap x and z
			R->set(0.,0.,1.,  0.,1., 0.,  -1.,0.,0.);
		}
		else
		{	// no rotation
			R->set(1.,0.,0.,  0.,1., 0.,  0.,0.,1.);
		}
	}
	else
	{	// already cos and sin
		R->set(soft[NORMALDIR1], -soft[NORMALDIR2], soft[NORMALDIR2], soft[NORMALDIR1], 1);
	}
	return true;
}
#pragma mark OrthoSoftening::Accessors

SofteningLaw *OrthoSoftening::softeningXX(void) const { return softeningAI; }
SofteningLaw *OrthoSoftening::softeningYY(void) const { return softeningI; }
SofteningLaw *OrthoSoftening::softeningXYX(void) const { return softeningAII; }
SofteningLaw *OrthoSoftening::softeningXYY(void) const { return softeningTII; }
SofteningLaw *OrthoSoftening::softeningYZY(void) const { return softeningII; }

// return material type
const char *OrthoSoftening::MaterialType(void) const
{	return "Orthotropic softening";
}

// calculate maximum wave speed in mm/sec (moduli in MPa, rho in g/mm^3)
double OrthoSoftening::WaveSpeed(bool threeD,MPMBase *mptr) const
{
	double xx,cnorm,cshear;
	
	xx=1.-nuxy*nuyx-nuxz*nuzx-nuyz*nuzy-2.*nuxy*nuyz*nuzx;
	cnorm=fmax(Ex*(1.-nuyz*nuzy),Ey*(1.-nuxz*nuzx))/xx;
	if(threeD) cnorm=fmax(cnorm,Ez*(1-nuxy*nuyx)/xx);
	cshear = threeD ? fmax(Gxy,fmax(Gxz,Gyz)) : Gxy;
	return sqrt(fmax(cnorm,cshear)/rho);
}

// diffusion and conductivity in the z direction
double OrthoSoftening::GetDiffZ(void) const { return Dz; }
double OrthoSoftening::GetKcondZ(void) const { return kCondz; }

