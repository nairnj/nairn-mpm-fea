/********************************************************************************
    MicrostructureModel.cpp
    NairnMPM
    
    Created by Tim Fagan in May 2012.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/MicrostructureModel.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"

#pragma mark MicrostructureModel::Constructors and Destructors

// Constructors
MicrostructureModel::MicrostructureModel()
{
}

// Constructors
MicrostructureModel::MicrostructureModel(char *matName) : IsoPlasticity(matName)
{
    // default values for constants
    rhoW=1.;
	rhoC=1.;
	fLim=1.;
	fo=1.;
	fsto=1.;
	sto=1.;
	SHM=1.;
	N=1.;
	alp=1.;
	burg=1.;
	K1=1.;
	esal=1.;
	esbe=1.;
	disk1=1.;
	tayM=3.06;
	MMG=1.;
	// for first step calcs:
	//fr=-3;
	//tdl = fr*rhoW+(1-fr)*rhoC;
	//dSize = K1/sqrt(tdl);
	//yieldPrevious=1.0;

}

#pragma mark MicrostructureModel::Initialization

// print to output window
void MicrostructureModel::PrintMechanicalProperties(void)
{	
    IsotropicMat::PrintMechanicalProperties();
	PrintYieldProperties();
}

// print just yield properties to output window
void MicrostructureModel::PrintYieldProperties(void)
{
	PrintProperty("p_total_0",tdl0,"\n");
	PrintProperty("p_c_0",rhoC0,"\n");
	PrintProperty("p_w_0",rhoW0,"\n");
	PrintProperty("Initial Average Cell Size",dSize0,"");
    cout << endl;
}

// The base class history variable is cummulative equivalent plastic strain
//		(defined as dalpha = sqrt((2/3)||dep||))
// 1: Yield (in PA), 2: RhoC 3: RhoW 4: Cell Size 5: Total Dislocation Density 6: plastic strain rate in sec^-1 (empty) 7: volume fraction (empty)
char *MicrostructureModel::MaterialData(void)
{
	double *p=new double[6];
	p[0]=0.;
	p[1]=0.;
	p[2]=0.;
	p[3]=0.;
	p[4]=0.;
	p[5]=0.;
	//p[6]=0.;
	//p[7]=0.;
	return (char *)p;
}

// Read material properties
char *MicrostructureModel::InputMat(char *xName,int &input)
{
    // tayM: Taylor Factor
    if(strcmp(xName,"tayM")==0)
    {   input=DOUBLE_NUM;
        return((char *)&tayM);
    }

	// rhoW: Initial dislocation density in cell wall
    else if(strcmp(xName,"rhoW")==0)
    {   input=DOUBLE_NUM;
        return((char *)&rhoW0);
    }

	// rhoC: Initial dislocation density in cell
    else if(strcmp(xName,"rhoC")==0)
    {   input=DOUBLE_NUM;
        return((char *)&rhoC0);
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
        return((char *)&SHM);
    }

	// n: initial n
    else if(strcmp(xName,"n")==0)
    {   input=DOUBLE_NUM;
        return((char *)&N);
    }

	// alp: alpha
    else if(strcmp(xName,"alp")==0)
    {   input=DOUBLE_NUM;
        return((char *)&alp);
    }

	// burg: burgers vector
    else if(strcmp(xName,"burg")==0)
    {   input=DOUBLE_NUM;
        return((char *)&burg);
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

    return(IsoPlasticity::InputMat(xName,input));
}

// Constant reduced properties used in constitutive law
void MicrostructureModel::InitialLoadMechProps(int makeSpecific,int np)
{
	//fr0=fo;
	//tdl0 = fr0*rhoW0+(1-fr0)*rhoC0;
	//dSize0 = K1/sqrt(tdl0);
	MMG=G*1.e6;				// input to IsotropicMat. IsoPlasticity uses Ge6/rho ?
	

    IsoPlasticity::InitialLoadMechProps(makeSpecific,np);
	
	/* Called once at beginning (by VerifyProperties() in MaterialBase). For efficiency,
	use this method to calculate new terms that are independent of the particle
	state and thus will remain constant throughout the calculation. When done
	(or before), pass on to super class (but MaterialBase and Elastic do not need it)
*/
}

#pragma mark MicrostructureModel::Methods

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
double MicrostructureModel::GetYield(MPMBase *mptr,int np,double delTime)
{	
		// create initial values for variables at initial timestep
		if(mptr->initfr==-3.)
		{	mptr->initfr=0.;
		//	mptr->fr=fr0;
			mptr->archiverhoC=rhoC0;
			mptr->archiverhoW=rhoW0;
		//	mptr->archiveTDL=tdl0;
		//	mptr->archiveDSize=dSize0;
		}
		
		
		// Calculation of m for temperature dependent strain-rate sensitivity
		if(ConductionTask::active)
		{	SHM=30000/mptr->pTemperature;
			N=14900/mptr->pTemperature;
		}
		
		double rhoc = mptr->archiverhoC;
		double rhow = mptr->archiverhoW;
		
		// update rhoc/w from previous increment
			rhow += mptr->rhoWDot*delTime;
			rhoc += mptr->rhoCDot*delTime;
		
		// update values for strain and strain rate
		eqss = SQRT_THREE*alpint;
		rss = tayM*eqss;
		eqssra = SQRT_THREE*(dalpha/delTime);
		rssra = tayM*eqssra;
		
		// update value of volume fraction
		mptr->fr = fLim + (fo-fLim)*exp(-1.*rss/fsto);
		
		double tdl = mptr->fr*rhow+(1.-mptr->fr)*rhoc;
		double dSize = K1/sqrt(tdl);
		
		
		double rhowDot,rhocDot = 0;
		
		// update dislocation density evolution rate in cell wall and cell interior
		if(rssra!=0)
		{	wAdd = (6.*esbe*rssra*pow(1-mptr->fr,TWOTHIRDS))/(burg*dSize*mptr->fr);
			wRem = (SQRT_THREE*esbe*rssra*(1.-mptr->fr)*sqrt(rhow))/(mptr->fr*burg);
			wDis = -disk1*pow(rssra/sto,-1./N)*rssra*rhow;
			rhowDot = wAdd+wRem+wDis;
			
			cAdd = esal*SQRT_ONETHIRD*(sqrt(rhow)/burg)*rssra;
			cRem = -esbe*((6.*rssra)/(burg*dSize*pow(1.-mptr->fr,ONETHIRD)));
			cDis = -disk1*pow(rssra/sto,-1./N)*rssra*rhoc;
			rhocDot = cAdd+cRem+cDis;
		}
			

		// update 
		sigoc=alp*MMG*burg*sqrt(rhoc);
		rstc=sigoc*pow(rssra/sto,1./SHM);
		
		// update 
		sigow=alp*MMG*burg*sqrt(rhow);
		rstw=sigow*pow(rssra/sto,1./SHM);
		
		// update stress
		rst=mptr->fr*rstw+(1.-mptr->fr)*rstc;
		mptr->yieldP = mptr->yieldC;
		mptr->yieldC = tayM*rst*(SQRT_THREE/rho);
		//mptr->yieldC *= SQRT_THREE/rho;
		mptr->yieldC += yldred;
		
		// save data if saving:
		if(saving)
		{	mptr->archiverhoC=rhoc;
			mptr->archiverhoW=rhow;
			mptr->archiveTDL = tdl;
			mptr->archiveDSize = dSize;
			mptr->rhoWDot = rhowDot;
			mptr->rhoCDot = rhocDot;
		}
		
		
		// return new stress value
		return mptr->yieldC;
		
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
double MicrostructureModel::GetKPrime(MPMBase *mptr,int np,double delTime)
{ 
	//if(dalpha>0.00000001||dalpha<-0.00000001)
	if(dalpha!=0)
		return (mptr->yieldC-mptr->yieldP)/(dalpha*rho);
	else 
		return 0;
		
}

double MicrostructureModel::GetK2Prime(MPMBase *mptr,double fnp1,double delTime)
{
	return 0;
}

// Solve this linear model in closed form for lambdak and update trial alpha
/*double MicrostructureModel::SolveForLambdaBracketed(MPMBase *mptr,int np,double strial,Tensor *stk,double delTime)
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
}*/

// over-riding base class IsoPlasticity
void MicrostructureModel::UpdatePlasticInternal(MPMBase *mptr,int np)
{
		mptr->SetHistoryDble(0,alpint);
		mptr->SetHistoryDble(YT_HISTORY,mptr->yieldC);
		//mptr->SetHistoryDble(EPDOT_HISTORY,(double)0.0);
		mptr->SetHistoryDble(RHOC,mptr->archiverhoC);
		mptr->SetHistoryDble(RHOW,mptr->archiverhoW);
		mptr->SetHistoryDble(DSIZE,mptr->archiveDSize);
		mptr->SetHistoryDble(TDL,mptr->archiveTDL);
		//mptr->SetHistoryDble(FR,mptr->fr);
}

#pragma mark VonMises::Accessors

// Return the material tag
int MicrostructureModel::MaterialTag(void) { return MICROSTRUCTUREMODEL; }

// return material type
const char *MicrostructureModel::MaterialType(void) { return "Dislocation Density Based Microstructure Model"; }

// over-riding base class IsPlasticity
// this material has three additional history variables
double MicrostructureModel::GetHistory(int num,char *historyPtr)
{
    double history=0.;
	if(num==1 || num==2 || num==3 || num==4 || num==5 || num==6)
	{	double *cumStrain=(double *)historyPtr;
		history=cumStrain[num-1];
	}
	return history;
}



