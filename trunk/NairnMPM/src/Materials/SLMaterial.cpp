/********************************************************************************
    SLMaterial.cpp
    nairn-mpm-fea
    
    Created by John Nairn, 11/12/2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
	
	Add rate dependence to yield stress of the MGEOS + SCGL Material
	See D. J. Steinberg and C. M. Lund, J. Appl. Phys., v64, 1528-1533 (1989)
	"A constitutive model for strain rates from 10^-4 to 10^6 s^-1
********************************************************************************/

#include "Materials/SLMaterial.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "System/UnitsController.hpp"

#pragma mark SLMaterial::Constructors and Destructors

// Constructors
SLMaterial::SLMaterial() {}

// Constructors
SLMaterial::SLMaterial(MaterialBase *pair) : SCGLHardening(pair)
{
	// defaults are from Steinberg and Lund
	UkOverk = 0.31/8.617332478e-5;						// eV/k in K
	YP = 1000.*UnitsController::Scaling(1.e6);
	C1 = 0.71e6;
	C2 = 0.012*UnitsController::Scaling(1.e6);
}

#pragma mark SLMaterial::Initialization

// Read material properties
char *SLMaterial::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	// unique properties here
	
	// here are the rest
    if(strcmp(xName,"Uk")==0)
    {	// Legacy units in eV and scaling converts to K
		// but need UkOverk when using consisten units
		input=DOUBLE_NUM;
		gScaling = UnitsController::Scaling(1./8.617332478e-5);
        return	gScaling>2. ? (char *)&UkOverk : NULL ;
    }
    
    else if(strcmp(xName,"UkOverk")==0)
    {	input=DOUBLE_NUM;
        return (char *)&UkOverk;
    }
    
    else if(strcmp(xName,"YP")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&YP,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"C1SL")==0)
    {	input=DOUBLE_NUM;
        return (char *)&C1;
    }
	
    else if(strcmp(xName,"C2SL")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&YP,gScaling,1.e6);
    }
	
	return(SCGLHardening::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *SLMaterial::VerifyAndLoadProperties(int np)
{	
	// check properties
	if(np==PLANE_STRESS_MPM) return "The Steinberg-Lund hardening does not support plane stress calculations yet.";
    
	YPred = YP/parent->rho;							// reduced stress units
	C2red = C2/parent->rho;					// reduced stress units
	YTmin = YPred/PRECISION_FACTOR;					// below this, strain dependent yield stress is zero
	YTprecision = YPred/PRECISION_FACTOR;			// precision as ratio to YPred
	
	// call super class
	return SCGLHardening::VerifyAndLoadProperties(np);
}

// print just yield properties to output window
void SLMaterial::PrintYieldProperties(void) const
{
    cout << GetHardeningLawName() << endl;
    
    // yield
    MaterialBase::PrintProperty("yld",yield*UnitsController::Scaling(1.e-6),"");
    MaterialBase::PrintProperty("beta",beta,"");
    MaterialBase::PrintProperty("nhard",nhard,"");
    MaterialBase::PrintProperty("yMax",yieldMax*UnitsController::Scaling(1.e-6),"");
    cout << endl;
    
	// shear temperature and pressure dependence
	char glabel[20];
	strcpy(glabel,UnitsController::Label(PRESSURE_UNITS));
	strcat(glabel,"^-1");
	MaterialBase::PrintProperty("Gp'/G0",GPp*UnitsController::Scaling(1.e6),glabel);
	MaterialBase::PrintProperty("GT'/G0",GTp,"K^-1");
	cout << endl;
    
    // Steinberg-Lund additinos
	MaterialBase::PrintProperty("Uk/k",UkOverk,"K");
	MaterialBase::PrintProperty("YP",YP*UnitsController::Scaling(1.e-6),"");
	strcpy(glabel,UnitsController::Label(TIME_UNITS));
	strcat(glabel,"^-1");
	MaterialBase::PrintProperty("C1",C1,glabel);
	strcpy(glabel,UnitsController::Label(PRESSURE_UNITS));
	strcat(glabel,"-");
	strcat(glabel,UnitsController::Label(TIME_UNITS));
	MaterialBase::PrintProperty("C2",C2*UnitsController::Scaling(1.e-6),glabel);
	cout << endl;
}

// 0: The base class history variable is cumulative equivalent plastic strain
//		(defined as dalpha = sqrt((2/3)||dep||))
// 1: YT (unreduced in MPA), 2: plastic strain rate in sec^-1
int SLMaterial::HistoryDoublesNeeded(void) const { return 3; }

#pragma mark SLMaterial:Methods

// size of hardening law properties needed in strain updates
int SLMaterial::SizeOfHardeningProps(void) const { return sizeof(SLProperties); }

// Get particle-state dependent properties (filled by Get Shear Ratio)
void *SLMaterial::GetCopyOfHardeningProps(MPMBase *mptr,int np,void *altBuffer)
{
	SLProperties *p = (SLProperties *)altBuffer;
	return p;
}

// Cast void * to correct pointer and delete it
void SLMaterial::DeleteCopyOfHardeningProps(void *properties,int np) const
{
	SLProperties *p = (SLProperties *)properties;
	delete p;
}

// Get temperature changed needed in T-dependent yield stress latter
// Then pass on to super class for pressure calculation
// Also get shear modulus here since it is not needed until after this is called
// Store results needed later in hardening law properties available
double SLMaterial::GetShearRatio(MPMBase *mptr,double pressure,double J,void *properties) const
{
	// fetch just Gratio from parent
	double Gratio = SCGLHardening::GetShearRatio(mptr,pressure,J,NULL);
	if(properties==NULL) return Gratio;
	
	// thermal term in yield stress
	SLProperties *p = (SLProperties *)properties;
	p->Gratio = Gratio;
	p->TwoUkkT = 2.*UkOverk/mptr->pPreviousTemperature;						// 2Uk/kT
	double YTlast=mptr->GetHistoryDble(YT_HISTORY)*1.e6/parent->rho;
	p->currentYTred=fmax(YTmin,YTlast);
	p->epdotmin=GetEpdot(YTmin,p->TwoUkkT);								// rate at small fraction of YP
	p->epdotmax=GetEpdot(YPred,p->TwoUkkT);								// rate to get YP
	
	return Gratio;
}

#pragma mark SLMaterial::Law Methods

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
// yield = (yldred*(1 + beta ep)^n + YT(epdot)) * Gratio, where ep=alpint, epdot=dalpha/delTime
// but yldred*(1 + beta ep)^n is limited to yldMaxred and YT(epdot) is limited to YPred
double SLMaterial::GetYield(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	SLProperties *p = (SLProperties *)properties;
	
	// aThermal term
	double aThermal=fmin(yldred*pow(1.+beta*a->alpint,nhard),yldMaxred);

	// find rate-dependent term by Newton's method
	double YTred;
	if(p->isConstantYT)
	{	YTred = p->constantYT;
		p->currentYTred = YTred;
	}
	else
	{	double epdot=a->dalpha/delTime;
		YTred=YTmin;
		if(epdot>p->epdotmax)
			YTred=YPred;
		else if(epdot>p->epdotmin)
		{	// use Newton's method to find YTred(epdot)
			double YTi = p->currentYTred;				// start at previous value
			double lnepgoal=log(epdot);
			double epdoti,arg,slope;
			int iter=1;
			while(iter<MAX_ITERATIONS)
			{	epdoti = GetEpdot(YTi,p->TwoUkkT);
				arg = C2red/YTi;
				slope = epdoti*arg/YTi + (1. - epdoti*arg)*2.*p->TwoUkkT*(1.-YTi/YPred)/YPred;
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
            {   cout << "#max iterations exceeded with " << (YTi*parent->rho/1.e6) << " and " << (YTred*parent->rho/1.e6)
                        << " at epdot= " << epdot << endl;
            }
		}
		
		// save the result
		p->currentYTred = YTred;
	}
	
	// combine and return
	return (aThermal + YTred)*p->Gratio;
}

// Get derivative of sqrt(2./3)*yield wrt ln(epdot), but if constantYT,
//  get derivative of sqrt(2./3.)*yield wrt lambda or (2./3)*yield wrt alpha for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
double SLMaterial::GetKPrime(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	SLProperties *p = (SLProperties *)properties;
	
	// aThermal terms
	
	// slope zero if in constant max yield condition
	double aThermal=0.;
	if(yldred*pow(1.+beta*a->alpint,nhard)<yldMaxred)
	{	aThermal = DbleEqual(nhard,1.) ? yldred*beta :
				yldred*beta*nhard*pow(1.+beta*a->alpint,nhard-1.) ;
	}
	
	// if constantYT, then not rate dependent terms and want derivative wrt lambda
	if(p->isConstantYT)
	{	return TWOTHIRDS*aThermal*p->Gratio;
	}
	
	// Get derivative of sqrt(2./3.)*yield with respect to epdot times epdot (i.e., wrt ln(epdot)) for plane strain and 3D
	double YTslope = 0.;
	double epdot = a->dalpha/delTime;
	if(epdot>0. && epdot<=p->epdotmax)
	{	// Find slope from current values
		double YTred = p->currentYTred;				// recall current
		double arg = C2red/YTred;
		double slope = epdot*arg/YTred + (1. - epdot*arg)*2.*p->TwoUkkT*(1.-YTred/YPred)/YPred;
		YTslope = 1./slope;
	}
	
	// combine and return derivative wrt ln(epdot)
	return SQRT_TWOTHIRDS*(aThermal*a->dalpha + YTslope)*p->Gratio;
}

// place holder until plane stress is allowed
double SLMaterial::GetK2Prime(MPMBase *mptr,double fnp1,double delTime,HardeningAlpha *a,void *properties) const
{
	return 0.0;
}


// find epdot from Steinberg-Lund expression for give reduced YT
double SLMaterial::GetEpdot(double YT,double TwoUkkT) const
{	double arg = 1.-YT/YPred;
	return 1./(exp(TwoUkkT*arg*arg)/C1 + C2red/YT);
}

// Solve numerically for lambda
// Subclasses can override for analytical solution if possible
double SLMaterial::SolveForLambdaBracketed(MPMBase *mptr,int np,double strial,Tensor *stk,double Gred,
								double psKred,double Ptrial,double delTime,HardeningAlpha *a,void *properties) const
{
	if(np==PLANE_STRESS_MPM)
	{	// not allowed
		return 0.;
	}
	else
	{	// solve - sqrt(2/3)GetYield(alpha+dalpha,dalpha) + strial - 2 GRed sqrt(3/2)dalpha  = 0
		SLProperties *p = (SLProperties *)properties;
		
		// test lower limit
		a->dalpha = delTime*p->epdotmin;
		a->alpint = mptr->GetHistoryDble() + a->dalpha;
		p->isConstantYT = true;
		p->constantYT = YTmin;
		double gmin = strial - 2.*Gred*a->dalpha/SQRT_TWOTHIRDS - SQRT_TWOTHIRDS*GetYield(mptr,np,delTime,a,p);
		if(gmin<0.)
		{	// low strain rate answer between 0 and epdotmin
			double lambdak = HardeningLawBase::SolveForLambda(mptr,np,strial,stk,Gred,psKred,1.,delTime,a,p);
			//cout << "# low strain rate condition " << lambdak << " should be below " << delTime*epdotmin/SQRT_TWOTHIRDS << endl;
			p->isConstantYT = false;
			mptr->SetHistoryDble(YT_HISTORY,p->currentYTred*parent->rho*1.e-6);
			mptr->SetHistoryDble(EPDOT_HISTORY,SQRT_TWOTHIRDS*lambdak/delTime);
			return lambdak;
		}
		
		// test upper limit
		a->dalpha = delTime*p->epdotmax;
		a->alpint = mptr->GetHistoryDble() + a->dalpha;
		p->constantYT=YPred;
		double gmax =  strial - 2.*Gred*a->dalpha/SQRT_TWOTHIRDS - SQRT_TWOTHIRDS*GetYield(mptr,np,delTime,a,p);
		if(gmax>0.)
		{	// high string rate answer for rate higher than epmax
			double lambdak=HardeningLawBase::SolveForLambda(mptr,np,strial,stk,Gred,psKred,1.,delTime,a,p);
			//cout << "# high strain rate condition " << lambdak << " should be above " << delTime*epdotmax/SQRT_TWOTHIRDS << endl;
			p->isConstantYT = false;
			mptr->SetHistoryDble(YT_HISTORY,p->currentYTred*parent->rho*1.e-6);
			mptr->SetHistoryDble(EPDOT_HISTORY,SQRT_TWOTHIRDS*lambdak/delTime);
			return lambdak;
		}
		p->isConstantYT=false;
		
		// Newton method in ln epdot space
		p->currentYTred=fmax(YTmin,mptr->GetHistoryDble(YT_HISTORY)*1.e6/parent->rho);
		double epdot=GetEpdot(p->currentYTred,p->TwoUkkT);
		double logepdot = log(epdot);
		a->dalpha = epdot*delTime;
		a->alpint = mptr->GetHistoryDble() + a->dalpha;
		int step=1;
		while(true)
		{	// update iterative variables (alpha, dalpha)
			double glam = -SQRT_TWOTHIRDS*GetYield(mptr,np,delTime,a,p) + strial - 2.*Gred*a->dalpha/SQRT_TWOTHIRDS;
			double slope = -2.*Gred*a->dalpha/SQRT_TWOTHIRDS - GetKPrime(mptr,np,delTime,a,p);
			double delLogepdot = -glam/slope;
			logepdot += delLogepdot;
			
			// check for convergence
			a->dalpha = exp(logepdot)*delTime;
			a->alpint = mptr->GetHistoryDble() + a->dalpha;
			if(step>20 || fabs(delLogepdot)<0.0001) break;
			step++;
		}
	
		// set history when done
		mptr->SetHistoryDble(YT_HISTORY,p->currentYTred*parent->rho*1.e-6);
		mptr->SetHistoryDble(EPDOT_HISTORY,a->dalpha/delTime);
		return a->dalpha/SQRT_TWOTHIRDS;
		
	}
	
}

// Elastic means zero strain rate o YT is zero (or YTmin)
void SLMaterial::ElasticUpdateFinished(MPMBase *mptr,int np,double delTime) const
{	mptr->SetHistoryDble(YT_HISTORY,(double)0.0);
	mptr->SetHistoryDble(EPDOT_HISTORY,(double)0.0);
}

#pragma mark SLMaterial::Accessors

// hardening law name
const char *SLMaterial::GetHardeningLawName(void) const { return "SL hardening"; }

// this hardening law has three history variables
double SLMaterial::GetHistory(int num,char *historyPtr) const
{
    double history=0.;
	if(num==1 || num==2 || num==3)
	{	double *cumStrain=(double *)historyPtr;
		history=cumStrain[num-1];
	}
    return history;
}

