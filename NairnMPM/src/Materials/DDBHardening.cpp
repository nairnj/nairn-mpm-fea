/********************************************************************************
    DDBHardening.cpp
	Written by Tim Fagan, June, 2013.
	
    nairn-mpm-fea
    
    Created by John Nairn, 11/12/2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
 
	A dislocation density based polycrystal plasticity model
	See Toth et. al., J. Eng. Mat. and Techn., v124, 71-77 (2002)
	"Strain hardening at large strains as predicted by 
	dislocation based polycrystal plasticity model"

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#include "Materials/DDBHardening.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "System/UnitsController.hpp"

#pragma mark DDBHardening::Constructors and Destructors

// Constructors
DDBHardening::DDBHardening() {}

// The default contructor should call a parent class constructor and
// then fill in any new initialization.
DDBHardening::DDBHardening(MaterialBase *pair) : HardeningLawBase(pair)
{
	// defaults are for copper under high strain rate. see "Grain refinement under 
	// high strain rate impact: A numerical approach" - V. Lemiale (2010)

	rhoW0=1e13;			// 1/L^2
	rhoC0=1e14;			// 1/L^2
	fLim=0.06;			// fraction (unitless)
	fo=0.25;			// fraction (unitless)
	fsto=3.2;			// fraction (unitless)
	sto=1e6;			// 1/T
	alp=0.25;			// (unitless)
	burg=2.56e-10;		// L
	K1=10;				// (unitless)
	esal=0.120;			// (unitless)
	esbe=0.006;			// (unitless)
	disk1=5.8;			// (unitless)
	tayM=3.06;			// (unitless)
	MMG=47.4e3;			// F/L^2
	Atd=30000;			// 1/K
	Btd=14900;			// 1/K
	SHM0=50;			// (unitless)
	N0=10;				// (unitless)
	tempDepend = 0;		// option for fixed m and n or temperature dependence

}

#pragma mark DDBHardening::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *DDBHardening::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	// tayM: Taylor Factor
    if(strcmp(xName,"tayM")==0)
    {   input=DOUBLE_NUM;
        return((char *)&tayM);
    }

	// rhoW: Initial dislocation density in cell wall (Legacy 1/m^2 to 1/mm^2)
    else if(strcmp(xName,"rhoW")==0)
    {   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&rhoW0,gScaling,1.e-6);
    }

	// rhoC: Initial dislocation density in cell (Legacy 1/m^2 to 1/mm^2)
    else if(strcmp(xName,"rhoC")==0)
    {   input=DOUBLE_NUM;
 		return UnitsController::ScaledPtr((char *)&rhoC0,gScaling,1.e-6);
    }

	// fo: Initial volume fraction
    else if(strcmp(xName,"fo")==0)
    {   input=DOUBLE_NUM;
        return((char *)&fo);
    }

	// fLim: limit of volume fraction
    else if(strcmp(xName,"flim")==0)
    {   input=DOUBLE_NUM;
        return((char *)&fLim);
    }

	// fsto: reference volume fraction
    else if(strcmp(xName,"fsto")==0)
    {   input=DOUBLE_NUM;
        return((char *)&fsto);
    }

	// sto: reference slip rate
    else if(strcmp(xName,"sto")==0)
    {   input=DOUBLE_NUM;
        return((char *)&sto);
    }

	// m: initial m
    else if(strcmp(xName,"m")==0)
    {   input=DOUBLE_NUM;
        return((char *)&SHM0);
    }

	// n: initial n
    else if(strcmp(xName,"n")==0)
    {   input=DOUBLE_NUM;
        return((char *)&N0);
    }

	// alp: alpha
    else if(strcmp(xName,"alp")==0)
    {   input=DOUBLE_NUM;
        return((char *)&alp);
    }

	// burg: burgers vector (Legacy m to mm)
    else if(strcmp(xName,"burg")==0)
    {   input=DOUBLE_NUM;
 		return UnitsController::ScaledPtr((char *)&burg,gScaling,1.e3);
    }

	// K1: K in average cell size
    else if(strcmp(xName,"K1")==0)
    {   input=DOUBLE_NUM;
        return((char *)&K1);
    }

	// alpstar: alpha*
    else if(strcmp(xName,"alpstar")==0)
    {   input=DOUBLE_NUM;
        return((char *)&esal);
    }

	// betastar: beta*
    else if(strcmp(xName,"betastar")==0)
    {   input=DOUBLE_NUM;
        return((char *)&esbe);
    }

	// ko: ko
    else if(strcmp(xName,"k0")==0)
    {   input=DOUBLE_NUM;
        return((char *)&disk1);
    }
	
	// Atd: A
    else if(strcmp(xName,"Atd")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Atd);
    }
	
	// Btd: B
    else if(strcmp(xName,"Btd")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Btd);
    }

    // true or false for temperature dependence
    else if(strcmp(xName,"tempDepend")==0)
    {   input=DOUBLE_NUM;
        return((char *)&tempDepend);
    }
	
	// MMG: Shear Modulus (Legacy MPa to Pa)
    else if(strcmp(xName,"MMG")==0)
    {   input=DOUBLE_NUM;
 		return UnitsController::ScaledPtr((char *)&MMG,gScaling,1.e6);
    }
	
    return HardeningLawBase::InputMaterialProperty(xName,input,gScaling);
}

// Verify input properties do calculations; if problem return string with an error message
// If OK, MUST pass on to super class. This is called just before PrintMaterial
// (see also ValidateForUse() for checks that depend on MPM calculation mode)
const char *DDBHardening::VerifyAndLoadProperties(int np)
{
	// check only isotropic material

	// check not plane stress
	if(np==PLANE_STRESS_MPM) return "The dislocation-density based hardening does not support plane stress calculations yet.";
    
	// check properties
	
	// reduced yield stress in base class
	HardeningLawBase::VerifyAndLoadProperties(np);
	
	// initial grain size for material
	rhoT0 = fo *rhoW0 +(1.-fo)*rhoC0;			// 1/L^2
	dSize0 = K1/sqrt(rhoT0);					// L
	
	// base class never has an error
    return NULL;
}

// print just yield properties to output window
void DDBHardening::PrintYieldProperties(void) const
{
    cout << GetHardeningLawName() << endl;
    
    // yield
	char ulab[50];
	sprintf(ulab,"1/%s^2",UnitsController::Label(CULENGTH_UNITS));
    MaterialBase::PrintProperty("Initial cell dislocation density, pc",rhoC0,ulab);
	MaterialBase::PrintProperty("Initial cell wall dislocation density, pw",rhoW0,ulab);
	MaterialBase::PrintProperty("Taylor Factor",tayM,"");
    MaterialBase::PrintProperty("fo",fo,"");
    MaterialBase::PrintProperty("flim",fLim,"");
	MaterialBase::PrintProperty("fsto",fsto,"");
	MaterialBase::PrintProperty("sto",sto,"");
	MaterialBase::PrintProperty("M",SHM0,"");
	MaterialBase::PrintProperty("N",N0,"");
	MaterialBase::PrintProperty("Alpha",alp,"");
	MaterialBase::PrintProperty("Burgers Vector",burg,UnitsController::Label(CULENGTH_UNITS));
	MaterialBase::PrintProperty("K",K1,"");
	MaterialBase::PrintProperty("Alpha*",esal,"");
	MaterialBase::PrintProperty("Beta*",esbe,"");
	MaterialBase::PrintProperty("ko",disk1,"");
	MaterialBase::PrintProperty("A",Atd,"K^-1");
	MaterialBase::PrintProperty("B",Btd,"K^-1");
	MaterialBase::PrintProperty("G",MMG*UnitsController::Scaling(1.e-6),UnitsController::Label(PRESSURE_UNITS));
	
    cout << endl;
}

#pragma mark DDBHardening:Methods


// size of hardening law properties needed in strain updates
int DDBHardening::SizeOfHardeningProps(void) const { return sizeof(DDBHProperties); }


// Get particle-state dependent properties (can also filled by Get Shear Ratio)
// Called in GetCopyofMechanicalProps()...
// UpdateStrainsFirstTask calls MechanicalProps, which is then passed to UpdateStrain,
// which belongs to the MatPoint3D class. This function calls MPMConstitutiveLaw().
// So this function below is called every timestep prior to constitutive calcs..
void *DDBHardening::GetCopyOfHardeningProps(MPMBase *mptr,int np,void *altBuffer)
{
	DDBHProperties *p = (DDBHProperties *)altBuffer;

	if(altBuffer==NULL) return p;
	
	// Calculation of m and n for temperature dependent strain-rate sensitivity
	if( tempDepend>0. && ConductionTask::active )
	{	p->SHM=Atd/mptr->pPreviousTemperature;
		p->N=Btd/mptr->pPreviousTemperature;
	}

	else
	{	p->SHM=SHM0;
		p->N=N0;
	}

	// update variables from stored history values
	p->rhoC = mptr->GetHistoryDble(RHOC_HISTORY);
	p->rhoW = mptr->GetHistoryDble(RHOW_HISTORY);
	p->rhoCDot = mptr->GetHistoryDble(RHOCDOT_HISTORY);
	p->rhoWDot = mptr->GetHistoryDble(RHOWDOT_HISTORY); 
	p->yieldC = mptr->GetHistoryDble(YLDC_HISTORY)*1.e6/parent->rho;

	// nothing needed from superclass (HardenLawBase)
	return p;
}

// Cast void * to correct pointer and delete it
void DDBHardening::DeleteCopyOfHardeningProps(void *properties,int np) const
{
	DDBHProperties *p = (DDBHProperties *)properties;
	delete p;
}


#pragma mark DDBHardening::Law Methods

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
double DDBHardening::GetYield(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	DDBHProperties *p = (DDBHProperties *)properties;

	p->rhoCtemp = p->rhoC;				// 1/L^2
	p->rhoWtemp = p->rhoW;				// 1/L^2

	// update rhoc/w from previous increment
	p->rhoWtemp += p->rhoWDot*delTime;			// 1/(L^2-T)
	p->rhoCtemp += p->rhoCDot*delTime;			// 1/(L^2-T)
		
	// update values for strain and strain rate
	double eqss = SQRT_THREE*a->alpint;			// unitless
	double rss = tayM*eqss;						// unitless
	double eqssra = SQRT_THREE*(a->dalpha/delTime);		// 1/T
	double rssra = tayM*eqssra;							// 1/T

	// update value of volume fraction
	p->fr = fLim + (fo-fLim)*exp(-1.*rss/fsto);			// unitless
		
	p->rhoT = p->fr*p->rhoWtemp+(1.-p->fr)*p->rhoCtemp;		// 1/L^2
	p->dSize = K1/sqrt(p->rhoT);							// L
		
	// update dislocation density evolution rate in cell wall and cell interior
	if(!DbleEqual(rssra,0.))
	{	double wAdd = (6.*esbe*rssra*pow(1-p->fr,TWOTHIRDS))/(burg*p->dSize*p->fr);				// 1/(L^2-T)
		double wRem = (SQRT_THREE*esbe*rssra*(1.-p->fr)*sqrt(p->rhoWtemp))/(p->fr*burg);		// 1/(L^2-T)
		double wDis = -disk1*pow(rssra/sto,-1./p->N)*rssra*p->rhoWtemp;							// 1/(L^2-T)
		p->rhoWDotTemp = wAdd+wRem+wDis;														// 1/(L^2-T)
		
		double cAdd = esal*SQRT_ONETHIRD*(sqrt(p->rhoWtemp)/burg)*rssra;						// 1/(L^2-T)
		double cRem = -esbe*((6.*rssra)/(burg*p->dSize*pow(1.-p->fr,ONETHIRD)));				// 1/(L^2-T)
		double cDis = -disk1*pow(rssra/sto,-1./p->N)*rssra*p->rhoCtemp;							// 1/(L^2-T)
		p->rhoCDotTemp = cAdd+cRem+cDis;														// 1/(L^2-T)
	}
		
	// update 
	double sigoc=alp*MMG*burg*sqrt(p->rhoCtemp);			// F/L^2
	double rstc=sigoc*pow(rssra/sto,1./p->SHM);				// F/L^2
	
	// update 
	double sigow=alp*MMG*burg*sqrt(p->rhoWtemp);			// F/L^2
	double rstw=sigow*pow(rssra/sto,1./p->SHM);				// F/L^2
	
	// update stress
	double rst=p->fr*rstw+(1.-p->fr)*rstc;					// F/L^2
	p->yieldP = p->yieldC;									// F/L^2
	p->yieldInc = tayM*rst*(SQRT_THREE/parent->rho);		// F/L^2
	p->yieldC = p->yieldInc + yldred;						// F/L^2
	
	// return new stress value
	return p->yieldC;										// F/L^2
		
}


// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lambda or depdot/dlambda = sqrt(2./3.)/delTime
double DDBHardening::GetKPrime(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	DDBHProperties *p = (DDBHProperties *)properties;
	if(!DbleEqual(a->dalpha,0.))
		//return (p->yieldC-p->yieldP)/(a->dalpha*parent->rho); 
		return (p->yieldC-p->yieldP)/a->dalpha; // Removed the rho as already reduced in GetYield.
	else 
		return 0;
		
}

double DDBHardening::GetK2Prime(MPMBase *mptr,double fnp1,double delTime,HardeningAlpha *a,void *properties) const
{ return 0.;
}

// Return (K(alpha)-K(0)), which is used in dissipated energy calculation
// If K(0) in current particle state differs from yldred, will need to override
/*double DDBHardening::GetYieldIncrement(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	DDBHProperties *p = (DDBHProperties *)properties;
	// check that this is correct value. Can yieldInc be called before being initiated?
    return p->yieldInc;
}*/

#pragma mark DDBHardening::Return Mapping

double DDBHardening::SolveForLambdaBracketed(MPMBase *mptr,int np,double strial,Tensor *stk,double Gred,double psKred,double Pfinal,double delTime, HardeningAlpha *a,void *p) const
{
	// overwrite for saving of internal variables after...

	double lambdak = HardeningLawBase::SolveForLambdaBracketed(mptr, np, strial, stk, Gred, psKred, Pfinal, delTime, a, p);

	// set all history variables when done:
	// update internal variables from temporary values
	DDBHProperties *prop = (DDBHProperties *)p;
		
	// update history variables
	mptr->SetHistoryDble(EP_HISTORY,a->alpint);
	mptr->SetHistoryDble(EPDOT_HISTORY,a->dalpha/delTime);
	mptr->SetHistoryDble(YLDC_HISTORY,prop->yieldC*parent->rho/1.e6);
	mptr->SetHistoryDble(GRAIN_HISTORY,prop->dSize);
	mptr->SetHistoryDble(RHOC_HISTORY,prop->rhoCtemp);
	mptr->SetHistoryDble(RHOW_HISTORY,prop->rhoWtemp);
	mptr->SetHistoryDble(RHOCDOT_HISTORY,prop->rhoCDotTemp);
	mptr->SetHistoryDble(RHOWDOT_HISTORY,prop->rhoWDotTemp);

    // return final answer
	return lambdak;
}



#pragma mark DDBHardening:HistoryVariables

// The base class hardening law has cumulative equivalent plastic strain
//		(defined as dalpha = sqrt((2/3)||dep||))
// This class has a further 5 history variables
int DDBHardening::HistoryDoublesNeeded(void) const 
{ 
	/* See header file for current numbered list, but generally: 
	1. cumulative plastic strain,  
	2. equivalent von mises yield stress (MPa)
	3. strain rate (s^-1),
	4. grain size (m)
	x. volume fraction of dislocation density in cell wall
	x. total dislocation density (m^-2)  */

	return NUMBER_HISTORY;
}

void DDBHardening::InitPlasticHistoryData(double *p) const
{
	// initial values of dislocation densities
	p[GRAIN_HISTORY] = dSize0; // calculated in VerifyProps above
	p[RHOC_HISTORY] = rhoC0;
	p[RHOW_HISTORY] = rhoW0;

}

void DDBHardening::ElasticUpdateFinished(MPMBase *mptr,int np,double delTime) const
{	mptr->SetHistoryDble(EP_HISTORY,0.);
	mptr->SetHistoryDble(EPDOT_HISTORY,0.);
	mptr->SetHistoryDble(YLDC_HISTORY,0.);
	mptr->SetHistoryDble(GRAIN_HISTORY,dSize0);
	mptr->SetHistoryDble(RHOC_HISTORY,rhoC0);
	mptr->SetHistoryDble(RHOW_HISTORY,rhoW0);
	mptr->SetHistoryDble(RHOCDOT_HISTORY,0.);
	mptr->SetHistoryDble(RHOWDOT_HISTORY,0.);
}


// over-riding base class
void DDBHardening::UpdatePlasticInternal(MPMBase *mptr,int np,HardeningAlpha *a) const
{		// all variables updated in DDBHardening::SolveForLambdaBracketed
}

// Return history data for this material type when requested (this material has 8)
double DDBHardening::GetHistory(int num,char *historyPtr) const
{
    double history=0.;
	if(num>0 && num<9)
	{	double *hist=(double *)historyPtr;
		history=hist[num-1];
	}
	return history;
}


#pragma mark DDBHardening:Step Methods

#pragma mark DDBHardening::Custom Methods

#pragma mark DDBHardening::Accessors

// hardening law name
const char *DDBHardening::GetHardeningLawName(void) const { return "Dislocation-Density Based Hardening"; }



