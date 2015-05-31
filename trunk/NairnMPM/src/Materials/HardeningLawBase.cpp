/********************************************************************************
    HardeningLawBase.cpp
    nairn-mpm-fea

    Created by John Nairn, 1/17/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    This is elastic plastic, and never assigned to a material
********************************************************************************/

#include "Materials/HardeningLawBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "System/UnitsController.hpp"

#pragma mark HardeningLawBase::Constructors and Destructors

// Constructors
HardeningLawBase::HardeningLawBase() {}

// Constructors
HardeningLawBase::HardeningLawBase(MaterialBase *pair)
{
    yield = 1.e50;
    parent = pair;
}

HardeningLawBase::~HardeningLawBase() {}

#pragma mark LinearHardening::Initialize

// Read hardening law properties
char *HardeningLawBase::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    // base yield stress
    if(strcmp(xName,"yield")==0)
    {   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&yield,gScaling,1.e6);
    }
    
    // is not a hardening law property
    return NULL;
}

// verify settings and some initial calculations
const char *HardeningLawBase::VerifyAndLoadProperties(int np)
{
	// reduced yield stress
    yldred = yield/parent->rho;
	
	return NULL;
}
// The base class hardening law has cumulative equivalent plastic strain
//		(defined as dalpha = sqrt((2/3)||dep||))
// This should be put in history variable 0
int HardeningLawBase::HistoryDoublesNeeded(void) const { return 1; }

// initialize any non-zero history data needed in this plastic law
void HardeningLawBase::InitPlasticHistoryData(double *p) const {}

#pragma mark HardeningLawBase::Methods

// size of hardening law properties needed in strain updates
int HardeningLawBase::SizeOfHardeningProps(void) const { return 0; }

// If needed, hardenling law can create properties that depend on particle state
void *HardeningLawBase::GetCopyOfHardeningProps(MPMBase *mptr,int np,void *altBuffer)
{	return NULL;
}

// If need, cast void * to correct pointer and delete it
void HardeningLawBase::DeleteCopyOfHardeningProps(void *properties,int np) const {}

// handle prressure and temperature depence of the shear modulus
double HardeningLawBase::GetShearRatio(MPMBase *mptr,double pressure,double J,void *p) const { return 1.; }

#pragma mark HardeningLawBase::Law Methods

// Return (K(alpha)-K(0)), which is used in dissipated energy calculation
// If K(0) in current particle state differs from yldred, will need to override
double HardeningLawBase::GetYieldIncrement(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	return GetYield(mptr,np,delTime,a,properties) - yldred;
}

#pragma mark HardeningLawBase::Plastic Strain Updates

/* Get internal variables while iterating to find lambda
    Here alpha (in alpint) is cumulative equivalent plastic strain and
        The increment in alpha is dalpha. dalpha/delTime is plastic strain rate
    For plane strain and 3D
        dalpha = lambda sqrt(2/3) df/dsigma::df/dsigma = lambda sqrt(2/3) since df/dsigma is unit vector
    For plane stress
        dalpha = lambda sqrt(2/3) fnp1
    Version with two arguments called at start to initialize alpint and dalpha
*/
void HardeningLawBase::UpdateTrialAlpha(MPMBase *mptr,int np,HardeningAlpha *a) const
{	a->alpint = mptr->GetHistoryDble();
	a->dalpha = 0.;
}
void HardeningLawBase::UpdateTrialAlpha(MPMBase *mptr,int np,double lambdak,double fnp1,HardeningAlpha *a) const
{	a->dalpha = (np==PLANE_STRESS_MPM) ? SQRT_TWOTHIRDS*lambdak*fnp1 : SQRT_TWOTHIRDS*lambdak ;
	a->alpint = mptr->GetHistoryDble() + a->dalpha;
}

// The prior soltution for lambda tracked the internal variable. Just move to the particle now
void HardeningLawBase::UpdatePlasticInternal(MPMBase *mptr,int np,HardeningAlpha *a) const
{	mptr->SetHistoryDble(a->alpint);
}

// Chance to update when in elastic regime
void HardeningLawBase::ElasticUpdateFinished(MPMBase *mptr,int np,double delTime) const
{
}

#pragma mark HardeningLawBase::Return Mapping

/* Solve numerically for lambda by Newton's method
    This method is not used by the by base hardening laws, but it may be used by subclasses
        by overriding SolveForLambdaBracketed() and calling this instead
    Uses Newton's law with initial guess being lambda = dalpha/sqrt(2/3). The solution is not
        bracketed, which means it may not be safe. Its faster than bracketing when it is safe, but'
        otherwise should not be used (currently used by LinearHardening and SLHardening)
    Set alpint and dalpha before calling
    Note: psKred and Pfinal only needed for plane stress, otherwise ignored
    Hyperelastic can adjust Gred to be mubar = Gred*Ielbar
*/
double HardeningLawBase::SolveForLambda(MPMBase *mptr,int np,double strial,Tensor *stk,
							double Gred,double psKred,double Pfinal,double delTime,
							HardeningAlpha *a,void *p) const
{
	// initial lambdk from dalpha set before call, often 0, but might be otherwise
	double lambdak=a->dalpha/SQRT_TWOTHIRDS;
	int step=1;
	
	if(np==PLANE_STRESS_MPM)
	{	double n2trial = -stk->xx+stk->yy;
		n2trial *= n2trial/2;
		n2trial += 2.*stk->xy*stk->xy;
		double n1trial = stk->xx+stk->yy-2.*Pfinal;
		n1trial *= n1trial/6.;
		while(true)
		{	// update iterative variables (lambda, alpha, stress)
			double d1 = (1 + psKred*lambdak);
			double d2 = (1.+2.*Gred*lambdak);
			double fnp12 = n1trial/(d1*d1) + n2trial/(d2*d2);
			double kyld = GetYield(mptr,np,delTime,a,p);
			double glam = 0.5*fnp12 - kyld*kyld/3.;
			double fnp1 = sqrt(fnp12);
			double slope = -(psKred*n1trial/(d1*d1*d1) + 2*Gred*n2trial/(d2*d2*d2)) - GetK2Prime(mptr,fnp1,delTime,a,p);
			double delLam = -glam/slope;
			lambdak += delLam;
			UpdateTrialAlpha(mptr,np,lambdak,fnp1,a);
			
			// check for convergence
			if(LambdaConverged(step++,lambdak,delLam)) break;
		}
	}
	else
	{	while(true)
        {	// update iterative variables (lambda, alpha)
            double glam = -SQRT_TWOTHIRDS*GetYield(mptr,np,delTime,a,p) + strial - 2*Gred*lambdak;
            double slope = -2.*Gred - GetKPrime(mptr,np,delTime,a,p);
            double delLam = -glam/slope;
            lambdak += delLam;
            UpdateTrialAlpha(mptr,np,lambdak,(double)0.,a);
            
            // check for convergence
            if(LambdaConverged(step++,lambdak,delLam)) break;
        }
	}
	return lambdak;
}

/* Solve numerically for lambda by safe Newton's method (i.e., with bracketing)
    Subclasses can override for analytical solution is possible or if more efficient method
        is available (e.g., non-bracketed method in SolveForLambda())
    The input ftrial is f function when lambda=0 (but not useful in in plane stress)
    Note: psKred and Pfinal only needed for plane stress, otherwise ignored
    Hyperelastic can adjust Gred to be mubar = Gred*Ielbar
*/
double HardeningLawBase::SolveForLambdaBracketed(MPMBase *mptr,int np,double strial,Tensor *stk,
								double Gred,double psKred,double Pfinal,double delTime,
								HardeningAlpha *a,void *p) const
{
    double xl,xh;
    BracketSolution(mptr,np,strial,stk,Gred,psKred,Pfinal,delTime,&xl,&xh,a,p);
    
    // if fails to bracket, convert to zero deviatoric stress and continue
    // This option does not happen in plane stress calculations
    if(xh>xl) return xh;
    
	// initial lambdk midpoint of the brackets
	double lambdak=0.5*(xl+xh);
    UpdateTrialAlpha(mptr,np,lambdak,(double)0.,a);
    double dxold=fabs(xh-xl);
    double dx=dxold;
	int step=1;
	
	if(np==PLANE_STRESS_MPM)
	{	double n2trial = -stk->xx+stk->yy;
		n2trial *= n2trial/2;
		n2trial += 2.*stk->xy*stk->xy;
		double n1trial = stk->xx+stk->yy-2.*Pfinal;
		n1trial *= n1trial/6.;
        while(true)
        {	// update iterative variables (lambda, alpha)
			double d1 = (1 + psKred*lambdak);
			double d2 = (1.+2.*Gred*lambdak);
			double fnp12 = n1trial/(d1*d1) + n2trial/(d2*d2);
			double kyld = GetYield(mptr,np,delTime,a,p);
			double glam = 0.5*fnp12 - kyld*kyld/3.;
			double fnp1 = sqrt(fnp12);
			double slope = -(psKred*n1trial/(d1*d1*d1) + 2*Gred*n2trial/(d2*d2*d2)) - GetK2Prime(mptr,fnp1,delTime,a,p);
            
            // bisect if Newton out of range
            if( ((lambdak-xh)*slope-glam) * ((lambdak-xl)*slope-glam) >= 0. ||
               fabs(2.*glam) > fabs(dxold*slope) )
            {   dxold = dx;
                dx = 0.5*(xh-xl);
                lambdak = xl+dx;
                if(xl == lambdak) break;    // change in root is negligible
            }
            else
            {   dxold = dx;
                dx = glam/slope;
                double temp = lambdak;
                lambdak -= dx;
                if(temp == lambdak) break;  // change in root is negligible
            }
            
            // update and check convergence
            UpdateTrialAlpha(mptr,np,lambdak,fnp1,a);
            if(LambdaConverged(step++,lambdak,dx)) break;
            
            // reset limits
            if(glam < 0.)
                xl = lambdak;
            else
                xh = lambdak;
        }
	}
	else
	{	while(true)
        {	// update iterative variables (lambda, alpha)
            double glam = strial - 2*Gred*lambdak - SQRT_TWOTHIRDS*GetYield(mptr,np,delTime,a,p);
            double slope = -2.*Gred - GetKPrime(mptr,np,delTime,a,p);
            
            // bisect if Newton out of range
            if( ((lambdak-xh)*slope-glam) * ((lambdak-xl)*slope-glam) >= 0. ||
               fabs(2.*glam) > fabs(dxold*slope) )
            {   dxold = dx;
                dx = 0.5*(xh-xl);
                lambdak = xl+dx;
 				if(xl == lambdak) break;    // change in root is negligible
            }
            else
            {   dxold = dx; 
                dx = glam/slope;
                double temp = lambdak;
                lambdak -= dx;
                if(temp == lambdak) break;  // change in root is negligible
            }
            
            // update and check convergence
            UpdateTrialAlpha(mptr,np,lambdak,(double)0.,a);
            if(LambdaConverged(step++,lambdak,dx)) break;
            
            // reset limits
            if(glam < 0.)
                xl = lambdak;
            else
                xh = lambdak;
        }
	}
    
    // return final answer
	return lambdak;
}

/* Bracket the solution for lambda for safe Newton's method
    Subclass can override if have faster way to bracket
    ftrial is 3D or plane strain result for lambda=0 and it is positive
    Return lamNeg for f<0 (higher lambda) and lamPos where f>0 (lower lambda)
    Note: psKred and Pfinal only needed for plane stress, otherwise ignored
    Hyperelastic can adjust Gred to be mubar = Gred*Ielbar
*/
void HardeningLawBase::BracketSolution(MPMBase *mptr,int np,double strial,Tensor *stk,double Gred,double psKred,double Pfinal,
                                       double delTime,double *lamNeg,double *lamPos,
									   HardeningAlpha *a,void *p) const
{
    double epdot = 1.,gmax;
    int step=0;
    
    // take lambda = 0 as positive limit (to start)
    *lamPos = 0.;
    
    if(np==PLANE_STRESS_MPM)
	{	double n2trial = -stk->xx+stk->yy;
		n2trial *= 0.5*n2trial;
		n2trial += 2.*stk->xy*stk->xy;
		double n1trial = stk->xx+stk->yy-2.*Pfinal;
		n1trial *= n1trial/6.;
		
        // find when plane stress term become negative
        while(step<20)
        {   // try above
            a->dalpha = epdot*delTime;
            // note that dalpha here is actually dalpha*fnp1
            a->alpint = mptr->GetHistoryDble() + a->dalpha;
            double lambdak = a->dalpha/SQRT_TWOTHIRDS;
			double d1 = (1 + psKred*lambdak);
			double d2 = (1.+2.*Gred*lambdak);
			double fnp12 = n1trial/(d1*d1) + n2trial/(d2*d2);
 			double kyld = GetYield(mptr,np,delTime,a,p);
			gmax = 0.5*fnp12 - kyld*kyld/3.;
            if(gmax<0.)
			{	// set upper limits
				*lamNeg = a->dalpha/SQRT_TWOTHIRDS;
				return;
			}
            
            // update positive limit and go to next order of magnitude
            *lamPos = lambdak;
            epdot *= 10.;
            step++;
        }
        
        // exception if did not find answer in 20 orders of magnitude in strain rate
		cout << "# Material point information that caused the exception:" << endl;
		mptr->Describe();
		char errMsg[250];
		strcpy(errMsg,"Plasticity plane stress solution for material type '");
		strcat(errMsg,GetHardeningLawName());
		strcat(errMsg,"' could not be bracketed in 20 steps");
		throw CommonException(errMsg,"IsoPlasticity::BracketSolution");
    }
    else
	{	// It must be negative when all plastic or when alphamax = strial/(2.*Gred);
		// (assuming yield stress  monotonically increases with alpha and epdot and is not zero)
		// It will usually be closer to zero
		// Test intervals dx, dx*r, ... dx*r^n where r (>1) is ratio of sizes and dx*r^n = alphamax
		double dalpha = strial/(2.*Gred);
		a->alpint = mptr->GetHistoryDble() + dalpha;
		if(GetYield(mptr,np,delTime,a,p) <= 0.)
		{	*lamNeg = 0.;
			*lamPos = dalpha/SQRT_TWOTHIRDS;
		}
		
		// take full range and new Newton's method with bracketing do the binary bisection
		*lamNeg = dalpha/SQRT_TWOTHIRDS;
		return;
	}
}

// decide if the numerical solution for lambda has converged
// subclass can override to change convergence rules
bool HardeningLawBase::LambdaConverged(int step,double lambda,double delLam) const
{
	if(step>20 || fabs(delLam/lambda)<0.0001) return true;
	return false;
}

#pragma mark HardeningLawBase::Accessors

// IsoPlasticity has no history data, but the hardening law might (and are all assumed doubles here)
double HardeningLawBase::GetHistory(int num,char *historyPtr) const
{	
    double history=0.;
    if(num>0 && num<=HistoryDoublesNeeded())
    {	double *p=(double *)historyPtr;
		history = p[num-1];
    }
    return history;
}


