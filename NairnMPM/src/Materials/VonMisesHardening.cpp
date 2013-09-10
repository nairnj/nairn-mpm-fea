/********************************************************************************
    VonMisesHardening.cpp
    NairnMPM
    
    Created by Yajun Guo in Jan 2005.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/VonMisesHardening.hpp"
#include "MPM_Classes/MPMBase.hpp"

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

// The base class history variable is cummulative equivalent plastic strain
//		(defined as dalpha = sqrt((2/3)||dep||))
// 1: Yield (in PA), 2: RhoC 3: RhoW 4: Cell Size 5: Total Dislocation Density 6: plastic strain rate in sec^-1 (empty) 7: volume fraction (empty)
char *VonMisesHardening::MaterialData(void)
{
	double *p=new double[2];
	p[0]=0.;
	p[1]=0.;
	p[2]=0.;
	//p[3]=0.;
	//p[4]=0.;
	//p[5]=0.;
	//p[6]=0.;
	//p[7]=0.;
	return (char *)p;
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
	
	// modiftf #temperature-dependent-yield-stress
	//Temperature Dependent Shear Stress:
		// Copper (C10200 Oxygen-free Copper) attempt 1	
		if(alpint>0)
		{
		double Tore = 0;
		double pTemp = mptr->pTemperature; // get the current temperature of the particle
		
		// create linear interpolations to determine the shear yield stress
		if (pTemp<=24)
			Tore = 242e6;	
		else if(pTemp<=100)
			Tore = -0.30387e6*(pTemp-24)+242e6;
		else if(pTemp<=150)
			Tore = -0.34641e6*(pTemp-100)+219.39e6;
		else if(pTemp<=200)
			Tore = -0.34641e6*(pTemp-150)+202.07e6;
		else if(pTemp<=250)
			Tore = -0.34641e6*(pTemp-200)+184.75e6;
		else if(pTemp<=300)
			Tore = -0.46188e6*(pTemp-250)+167.43e6;
		else if(pTemp<=350)
			Tore = -0.63509e6*(pTemp-300)+144.34e6;
		else if(pTemp<=400)
			Tore = -0.80829e6*(pTemp-350)+112.58e6;
		else if(pTemp<=450)
			Tore = -0.86603e6*(pTemp-400)+72.17e6;
		else if(pTemp<=500)
			Tore = -0.26558e6*(pTemp-450)+28.87e6;
		else if(pTemp<=550)
			Tore = -0.12702e6*(pTemp-500)+15.59e6;
		else if(pTemp<=600)
			Tore = -0.09238e6*(pTemp-550)+9.24e6;
		else if(pTemp<=625)
			Tore = -0.03464e6*(pTemp-600)+4.62e6;
		else if(pTemp<=650)
			Tore = -0.08083e6*(pTemp-625)+3.75e6;
		else if(pTemp<=675)
			Tore = -0.06928e6*(pTemp-650)+1.73e6;
		else 
			Tore=0;
			
			mptr->yieldC = sqrt(3)*Tore/rho; 
			return mptr->yieldC;  					//return the temperature dependent yield stress
		}	
		// modiftf #temperature-dependent-yield-stress	
			
	else 
		return linearHardening ? yldred + Epred*alpint : yldred*pow(1.+beta*alpint,npow) ;
	
	
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
double VonMisesHardening::GetKPrime(MPMBase *mptr,int np,double delTime)
{
	return linearHardening ? TWOTHIRDS*Epred : TWOTHIRDS*yldred*beta*npow*pow(1.+beta*alpint,npow-1) ;
}

// Get derivative of (1./3.)*yield^2 with respect to lambda for plane stress only
// ... and using dep/dlambda = sqrt(2./3.)*fnp1 where ep=alpint
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)*lambda*fnp1 or depdot/dlambda = sqrt(2./3.)*fnp1/delTime
// Also equal to sqrt(2./3.)*GetYield()*GetKPrime()*fnp1, but in separate call for efficiency
double VonMisesHardening::GetK2Prime(MPMBase *mptr,double fnp1,double delTime)
{
	return linearHardening ? SQRT_EIGHT27THS*(yldred + Epred*alpint)*Epred*fnp1 : 
								SQRT_EIGHT27THS*yldred*yldred*beta*npow*pow(1.+beta*alpint,2.*npow-1)*fnp1;
}

// Solve this linear model in closed form for lambdak and update trial alpha
double VonMisesHardening::SolveForLambdaBracketed(MPMBase *mptr,int np,double strial,Tensor *stk,double delTime)
{
	// plane strain is numerical
	if(np==PLANE_STRESS_MPM || !linearHardening)
    {   // The unbracketed one is faster and seems stable for this hardening law
        return IsoPlasticity::SolveForLambda(mptr,np,strial,stk,delTime);
        //return IsoPlasticity::SolveForLambdaBracketed(mptr,np,strial,stk,delTime);
    }
		
	// closed form for plane strain and 3D and linear hardening
	double lambdak = (strial - sqrt(2./3.)*(yldred + Epred*alpint))/(2.*(Gred + Epred/3.));
	UpdateTrialAlpha(mptr,np,lambdak,(double)1.);
	return lambdak;
}

void VonMisesHardening::UpdatePlasticInternal(MPMBase *mptr,int np)
{
		mptr->SetHistoryDble(0,alpint);
		mptr->SetHistoryDble(YT_HISTORY,mptr->yieldC*rho/1.e6);
		mptr->SetHistoryDble(EPDOT_HISTORY,1.); // no delTime variable to apply...
		//mptr->SetHistoryDble(RHOC,mptr->archiverhoC);
		//mptr->SetHistoryDble(RHOW,mptr->archiverhoW);
		//mptr->SetHistoryDble(DSIZE,mptr->archiveDSize);
		//mptr->SetHistoryDble(TDL,mptr->archiveTDL);
		//mptr->SetHistoryDble(FR,mptr->fr);
}

#pragma mark VonMises::Accessors

// Return the material tag
int VonMisesHardening::MaterialTag(void) { return VONMISESHARDENING; }

// return material type
const char *VonMisesHardening::MaterialType(void) { return "Von Mises Hardening Plastic"; }


// over-riding base class IsPlasticity
// this material has three additional history variables
double VonMisesHardening::GetHistory(int num,char *historyPtr)
{
    double history=0.;
	if(num==1 || num==2 || num==3 ) //|| num==4 || num==5 || num==6)
	{	double *cumStrain=(double *)historyPtr;
		history=cumStrain[num-1];
	}
	return history;
}
