/********************************************************************************
    SLMaterial.cpp
    NairnMPM
    
    Created by John Nairn, 11/12/2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
	
	Add rate dependence to yield stress of the MGSCGLMaterial
	See D. J. Steinberg and C. M. Lund, J. Appl. Phys., v64, 1528-1533 (1989)
	"A constitutive model for strain rates from 10^-4 to 10^6 s^-1
********************************************************************************/

#include "Materials/SLMaterial.hpp"
#include "MPM_Classes/MPMBase.hpp"

#pragma mark SLMaterial::Constructors and Destructors

// Constructors
SLMaterial::SLMaterial() {}

/* The default contructor should call a parent class constructor and
	then fill in any new initialization.
	*/
// Constructors
SLMaterial::SLMaterial(char *matName) : MGSCGLMaterial(matName)
{
	// defaults are from Steinberg and Lund
	Uk=0.31;			// eV
	YP=1000.;			// MPa
	C1=0.71e6;			// 1/s
	C2=0.012;			// MPa s
}

#pragma mark SLMaterial::Initialization

// Read material properties
char *SLMaterial::InputMat(char *xName,int &input)
{
	// unique properties here
	
	// here are the rest
    if(strcmp(xName,"Uk")==0)
    {	input=DOUBLE_NUM;
        return((char *)&Uk);
    }
    
    else if(strcmp(xName,"YP")==0)
    {	input=DOUBLE_NUM;
        return((char *)&YP);
    }
    
    else if(strcmp(xName,"C1")==0)
    {	input=DOUBLE_NUM;
        return((char *)&C1);
    }
	
    else if(strcmp(xName,"C2")==0)
    {	input=DOUBLE_NUM;
        return((char *)&C2);
    }
	
	return(MGSCGLMaterial::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *SLMaterial::VerifyProperties(int np)
{	
	// check properties
	if(np==PLANE_STRESS_MPM) return "The Steinberg-Lund material does not support plane stress calculations yet.";
    
	// call super class
	return MGSCGLMaterial::VerifyProperties(np);
}

// Constant properties used in constitutive law
void SLMaterial::InitialLoadMechProps(int makeSpecific,int np)
{	
	YPred = YP*1.e6/rho;					// reduced stress units
	C2red = C2*1.e6/rho;					// reduced stress units
	YTmin = YPred/PRECISION_FACTOR;			// below this, strain dependent yield stress is zero
	YTprecision = YPred/PRECISION_FACTOR;	// precision as ratio to YPred
	
	// call superclass
	MGSCGLMaterial::InitialLoadMechProps(makeSpecific,np);
}

// print just yield properties to output window
void SLMaterial::PrintYieldProperties(void)
{
	MGSCGLMaterial::PrintYieldProperties();
	PrintProperty("Uk",Uk,"eV");
	PrintProperty("YP",YP,"");
	PrintProperty("C1",C1,"s^-1");
	PrintProperty("C2",C2,"MPa-s");
	cout << endl;
}

// The base class history variable is cummulative equivalent plastic strain
//		(defined as dalpha = sqrt((2/3)||dep||))
// 1: YT (unreduced in MPA), 2: plastic strain rate in sec^-1
char *SLMaterial::MaterialData(void)
{
	double *p=new double[3];
	p[0]=0.;
	p[1]=0.;
	p[2]=0.;
	return (char *)p;
}

#pragma mark SLMaterial::Custom Methods

// Get temperature changed needed in T-dependent yield stress latter
// Then pass on to super class for pressure calculation
// Also get shear modulus here since it is not needed until after this is called
double SLMaterial::GetPressureChange(MPMBase *mptr,double &delV,double J,int np)
{
	// thermal term in yield stress
	TwoUkkT = 2.*Uk/(8.61734e-5*mptr->pPreviousTemperature);			// 2Uk/kT (add 3 digit at end)
	epdotmin=GetEpdot(YTmin);											// rate at small fraction of YP
	epdotmax=GetEpdot(YPred);											// rate to get YP
	double YTlast=mptr->GetHistoryDble(YT_HISTORY)*1.e6/rho;
	currentYTred=fmax(YTmin,YTlast);

	// super class handles pressure calculation
	return MGSCGLMaterial::GetPressureChange(mptr,delV,J,np);
}

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
// yield = (yldred*(1 + beta ep)^n + YT(epdot)) * Gred/G0red, where ep=alpint, epdot=dalpha/delTime
// but lyldred*(1 + beta ep)^n is limited to yldMaxred and YT(epdot) is limited to YPred
double SLMaterial::GetYield(MPMBase *mptr,int np,double delTime)
{
	// aThermal term
	double aThermal=fmin(yldred*pow(1.+beta*alpint,nhard),yldMaxred);

	// find rate-dependent term by Newton's method
	double YTred;
	if(isConstantYT)
	{	YTred=constantYT;
		currentYTred=YTred;
	}
	else
	{	double epdot=dalpha/delTime;
		YTred=YTmin;
		if(epdot>epdotmax)
			YTred=YPred;
		else if(epdot>epdotmin)
		{	// use Newton's method to find YTred(epdot)
			double YTi=currentYTred;				// start at previous value
			double lnepgoal=log(epdot);
			double epdoti,arg,slope;
			int iter=1;
			while(iter<MAX_ITERATIONS)
			{	epdoti = GetEpdot(YTi);
				arg = C2red/YTi;
				slope = epdoti*arg/YTi + (1. - epdoti*arg)*2.*TwoUkkT*(1.-YTi/YPred)/YPred;
				YTred = YTi + (lnepgoal-log(epdoti))/slope;
				if(YTred<YTmin)
				{	YTred=YTmin;
					break;
				}
				else if(YTred>YPred)
				{	YTred=YPred;
					break;
				}
				else if(fabs(YTred-YTi)<=YTprecision)
					break;
				iter++;
				YTi=YTred;
			}
			if(iter>=MAX_ITERATIONS)
				cout << "#max iterations exceeded with " << (YTi*rho/1.e6) << " and " << (YTred*rho/1.e6) << " at epdot= " << epdot << endl;
		}
		
		// save the result
		currentYTred=YTred;
	}
	
	// combine and return
	return (aThermal + YTred)*Gred/G0red;
}

// Get derivative of sqrt(2./3)*yield wrt ln(epdot), but if constantYT,
//  get derivative of sqrt(2./3.)*yield wrt lambda or (2./3)*yield wrt alpha for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
double SLMaterial::GetKPrime(MPMBase *mptr,int np,double delTime)
{	
	// aThermal terms
	
	// slope zero if in constant max yield condition
	double aThermal=0.;
	if(yldred*pow(1.+beta*alpint,nhard)<yldMaxred)
	{	aThermal = DbleEqual(nhard,1.) ? yldred*beta :
				yldred*beta*nhard*pow(1.+beta*alpint,nhard-1.) ;
	}
	
	// if constantYT, then not rate dependent terms and want derivative wrt lambda
	if(isConstantYT)
	{	return TWOTHIRDS*aThermal*Gred/G0red;
	}
	
	// Get derivative of sqrt(2./3.)*yield with respect to epdot times epdot (i.e., wrt ln(epdot)) for plane strain and 3D
	double YTslope=0.;
	double epdot=dalpha/delTime;
	if(epdot>0. && epdot<=epdotmax)
	{	// Find slope from current values
		double YTred=currentYTred;				// recall current
		double arg = C2red/YTred;
		double slope = epdot*arg/YTred + (1. - epdot*arg)*2.*TwoUkkT*(1.-YTred/YPred)/YPred;
		YTslope = 1./slope;
	}
	
	// combine and return derivative wrt ln(epdot)
	return SQRT_TWOTHIRDS*(aThermal*dalpha + YTslope)*Gred/G0red;
}

// find epdot from Steinberg-Lund expression for give reduced YT
double SLMaterial::GetEpdot(double YT)
{	double arg=1.-YT/YPred;
	return 1./(exp(TwoUkkT*arg*arg)/C1 + C2red/YT);
}

// Solve numerically for lambda
// Subclasses can override for analytical solution if possible
double SLMaterial::SolveForLambdaBracketed(MPMBase *mptr,int np,double strial,Tensor *stk,double delTime)
{
	if(np==PLANE_STRESS_MPM)
	{	// not allowed
		return 0.;
	}
	else
	{	// solve - sqrt(2/3)GetYield(alpha+dalpha,dalpha) + strial - 2 GRed sqrt(3/2)dalpha  = 0
	
		// test lower limit
		dalpha = delTime*epdotmin;
		alpint = mptr->GetHistoryDble() + dalpha;
		isConstantYT=true;
		constantYT=YTmin;
		double gmin = strial - 2.*Gred*dalpha/SQRT_TWOTHIRDS - SQRT_TWOTHIRDS*GetYield(mptr,np,delTime);
		if(gmin<0.)
		{	// low strain rate answer between 0 and epdotmin
			double lambdak=IsoPlasticity::SolveForLambda(mptr,np,strial,stk,delTime);
			//cout << "# low strain rate condition " << lambdak << " should be below " << delTime*epdotmin/SQRT_TWOTHIRDS << endl;
			isConstantYT=false;
			mptr->SetHistoryDble(YT_HISTORY,currentYTred*rho*1.e-6);
			mptr->SetHistoryDble(EPDOT_HISTORY,SQRT_TWOTHIRDS*lambdak/delTime);
			return lambdak;
		}
		
		// test upper limit
		dalpha = delTime*epdotmax;
		alpint = mptr->GetHistoryDble() + dalpha;
		constantYT=YPred;
		double gmax =  strial - 2.*Gred*dalpha/SQRT_TWOTHIRDS - SQRT_TWOTHIRDS*GetYield(mptr,np,delTime);
		if(gmax>0.)
		{	// high string rate answer for rate higher than epmax
			double lambdak=IsoPlasticity::SolveForLambda(mptr,np,strial,stk,delTime);
			//cout << "# high strain rate condition " << lambdak << " should be above " << delTime*epdotmax/SQRT_TWOTHIRDS << endl;
			isConstantYT=false;
			mptr->SetHistoryDble(YT_HISTORY,currentYTred*rho*1.e-6);
			mptr->SetHistoryDble(EPDOT_HISTORY,SQRT_TWOTHIRDS*lambdak/delTime);
			return lambdak;
		}
		isConstantYT=false;
		
		// Newton method in ln epdot space
		currentYTred=fmax(YTmin,mptr->GetHistoryDble(YT_HISTORY)*1.e6/rho);
		double epdot=GetEpdot(currentYTred);
		double logepdot = log(epdot);
		dalpha = epdot*delTime;
		alpint = mptr->GetHistoryDble() + dalpha;
		int step=1;
		while(true)
		{	// update iterative variables (alpha, dalpha)
			double glam = -SQRT_TWOTHIRDS*GetYield(mptr,np,delTime) + strial - 2.*Gred*dalpha/SQRT_TWOTHIRDS;
			double slope = -2.*Gred*dalpha/SQRT_TWOTHIRDS - GetKPrime(mptr,np,delTime);
			double delLogepdot = -glam/slope;
			logepdot += delLogepdot;
			
			// check for convergence
			dalpha = exp(logepdot)*delTime;
			alpint = mptr->GetHistoryDble() + dalpha;
			if(step>20 || fabs(delLogepdot)<0.0001) break;
			step++;
		}
	
		// set history when done
		mptr->SetHistoryDble(YT_HISTORY,currentYTred*rho*1.e-6);
		mptr->SetHistoryDble(EPDOT_HISTORY,dalpha/delTime);
		return dalpha/SQRT_TWOTHIRDS;
		
	}
	
}

// Elastic means zero strain rate to YT is zero (or YTmin)
void SLMaterial::ElasticUpdateFinished(MPMBase *mptr,int np,double delTime)
{	mptr->SetHistoryDble(YT_HISTORY,(double)0.0);
	mptr->SetHistoryDble(EPDOT_HISTORY,(double)0.0);
}

#pragma mark SLMaterial::Accessors

// Return the material tag
int SLMaterial::MaterialTag(void) { return SLMATERIAL; }

// return unique, short name for this material
const char *SLMaterial::MaterialType(void) { return "Steinberg-Lund Material"; }

#pragma mark IsoPlaticity::Accessors

// this material has two history variables
double SLMaterial::GetHistory(int num,char *historyPtr)
{
    double history=0.;
	if(num==1 || num==2 || num==3)
	{	double *cumStrain=(double *)historyPtr;
		history=cumStrain[num-1];
	}
    return history;
}

