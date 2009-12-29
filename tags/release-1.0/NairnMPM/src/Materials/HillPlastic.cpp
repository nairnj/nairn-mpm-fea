/********************************************************************************
    HillPlastic.cpp
    NairnMPM
    
    Created by John Nairn, June 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		AnisoPlasticity.hpp (Orthotropic.hpp, TranIsotropic.hpp, MaterialBase.hpp)
********************************************************************************/

#include "HillPlastic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Read_XML/CommonReadHandler.hpp"

#pragma mark HillPlastic::Constructors and Destructors

// Constructors
HillPlastic::HillPlastic()
{
}

// Constructors
HillPlastic::HillPlastic(char *matName) : AnisoPlasticity(matName)
{
    // default value of hardening (elastic plastic)
	Khard=0.;
	nhard=1.;
}

#pragma mark HillPlastic::Initialize

// Read material properties
char *HillPlastic::InputMat(char *xName,int &input)
{
    // Khard: hardening parameter
    if(strcmp(xName,"Khard")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Khard);
    }

    // nhard - hardening exponent
    else if(strcmp(xName,"nhard")==0)
    {   input=DOUBLE_NUM;
        return((char *)&nhard);
    }
	
    return(AnisoPlasticity::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *HillPlastic::VerifyProperties(int np)
{
	double rsxx=0.,rsyy=0.,rszz=0.;
	
	// check potential terms
	if(syxx>=0.)
		rsxx=1./(syxx*syxx);
	if(syyy>=0.)
		rsyy=1./(syyy*syyy);
	if(syzz>=0.)
		rszz=1./(syzz*syzz);
	double arg = rsxx*rsxx + rsyy*rsyy + rszz*rszz - rsyy*rsxx - rszz*rsxx - rsyy*rszz ;
	double fgh = 0.5*(rsxx+rsyy+rszz);
	if(arg<0.) return "Hill plastic potential is not postive semidefinite (1)";
	if(fgh-sqrt(arg)<0.) return "Hill plastic potential is not postive semidefinite (2)";
	
	// call super class
	return AnisoPlasticity::VerifyProperties(np);
}

// Private properties used in constitutive law
void HillPlastic::InitialLoadMechProps(int makeSpecific,int np)
{
	AnisoPlasticity::InitialLoadMechProps(makeSpecific,np);
	
	// combination terms
	fTerm=(syyyred2 + syzzred2 - syxxred2)/2.;
	gTerm=(syzzred2 + syxxred2 - syyyred2)/2.;
	hTerm=(syxxred2 + syyyred2 - syzzred2)/2.;
}

// print to output window
void HillPlastic::PrintYieldProperties(void)
{	
    AnisoPlasticity::PrintYieldProperties();
	PrintProperty("K",Khard,"");
	PrintProperty("n",nhard,"");
    cout << endl;
}

// history is cumulative strain
char *HillPlastic::MaterialData(void)
{
	double *p=new double;
	*p=0.;
	return (char *)p;
}

#pragma mark HillPlastic::Methods

// Load current internal variables into local alpha variables
void HillPlastic::UpdateTrialAlpha(MPMBase *mptr,int np)
{	aint = mptr->GetHistoryDble();
}

// Update alpha: Here dalpha = lamda R df = lambda dfrot
void HillPlastic::UpdateTrialAlpha(MPMBase *mptr,int np,double lambdak)
{	aint = mptr->GetHistoryDble() + lambdak*minush;
}

// return the yield function in a subroutine so sub class can define new yield functions
// For Von-Mises isotropic linear hardening material, plastic potential 
// f=s_bar^2-yield^2-2 Ep Wp, where s_bar - equivalent Von-Mises stress,
// Ep - plastic modulus; Wp - plastic work.
double HillPlastic::GetF(MPMBase *mptr,double sxx,double syy,double txy,double szz,int np)
{
	// clockwise rotation from analysis to material axes
	sxrot=sxx*cos2t + syy*sin2t - 2.*txy*costsint;
	syrot=sxx*sin2t + syy*cos2t + 2.*txy*costsint;
	txyrot=(sxx-syy)*costsint + txy*(cos2t-sin2t);
	
	// the potential
	//double sAs = sxrot*sxrot*syxxred2 + syrot*syrot*syyyred2 + txyrot*txyrot*tyxyred2 - 2.*hTerm*sxrot*syrot;
	//if(np==PLANE_STRAIN_MPM)
	//	sAs += szz*szz*syzzred2 - 2.*fTerm*syrot*szz - 2.*gTerm*sxrot*szz;
	double dxy=sxrot-syrot,sAs;
	if(np==PLANE_STRAIN_MPM)
	{	double dyz=syrot-szz;
		double dxz=sxrot-szz;
		sAs=fTerm*dyz*dyz + gTerm*dxz*dxz + hTerm*dxy*dxy + txyrot*txyrot*tyxyred2;
	}
	else
		sAs=fTerm*syrot*syrot + gTerm*sxrot*sxrot + hTerm*dxy*dxy + txyrot*txyrot*tyxyred2;
	
	// return sqrt(sAs) - pow(1. + Khard*aint,nhard);
	// check on negative sAs can happen due to round-off error when stresses near zero
	return sAs>0. ? sqrt(sAs) - (1. + Khard*pow(aint,nhard)) : -1.; 
}

// return derivatives of the yield function in a subroutine so sub class can define new yield functions
void HillPlastic::GetDfDsigma(MPMBase *mptr,Tensor *st0,int np)
{
	// clockwise rotation from analysis to material axes
	sxrot = st0->xx*cos2t + st0->yy*sin2t - 2.*st0->xy*costsint;
	syrot = st0->xx*sin2t + st0->yy*cos2t + 2.*st0->xy*costsint;
	txyrot = (st0->xx-st0->yy)*costsint + st0->xy*(cos2t-sin2t);
	
	// sqrt(s As)
	double sAs = sxrot*sxrot*syxxred2 + syrot*syrot*syyyred2 + txyrot*txyrot*tyxyred2 - 2.*hTerm*sxrot*syrot;
	if(np==PLANE_STRAIN_MPM)
		sAs += st0->zz*st0->zz*syzzred2 - 2.*fTerm*syrot*st0->zz - 2.*gTerm*sxrot*st0->zz;
	double rootSAS = sqrt(sAs);
	
	// the derivatives = dfrot = R df
	dfdsxxrot = syxxred2*sxrot - hTerm*syrot ;
	dfdsyyrot = -hTerm*sxrot + syyyred2*syrot;
	dfdtxyrot = tyxyred2*txyrot;
    if(np==PLANE_STRAIN_MPM)
	{	dfdsxxrot -= gTerm*st0->zz;
		dfdsyyrot -= fTerm*st0->zz;
		dfdszzrot = -gTerm*sxrot - fTerm*syrot + syzzred2*st0->zz;
		dfdszzrot/=rootSAS;
   }
	else
		dfdszzrot=0.;
	dfdsxxrot/=rootSAS;
	dfdsyyrot/=rootSAS;
	dfdtxyrot/=rootSAS;
		
	// rotate to analysis coordinates df = R^(-1) dfrot
	dfdsxx = dfdsxxrot*cos2t + dfdsyyrot*sin2t + dfdtxyrot*costsint;
	dfdsyy = dfdsxxrot*sin2t + dfdsyyrot*cos2t - dfdtxyrot*costsint;
	dfdtxy = -2*costsint*(dfdsxxrot - dfdsyyrot) + dfdtxyrot*(cos2t-sin2t);
	dfdszz = dfdszzrot;
	
	// for use in alpha upate
	minush = dfdsxxrot*dfdsxxrot + dfdsyyrot*dfdsyyrot + dfdszzrot*dfdszzrot + 0.5*dfdtxyrot*dfdtxyrot;
	minush = sqrt(minush/1.5);
}

// Non-linear hardening - df . h
double HillPlastic::GetDfAlphaDotH(MPMBase *mptr,int np,Tensor *st0)
{	//return DbleEqual(nhard,1.) ? Khard*minush :
	//			Khard*nhard*pow(1+Khard*aint,nhard-1)*minush ;
	return DbleEqual(nhard,1.) ? Khard*minush :
				Khard*nhard*pow(aint,nhard-1)*minush ;
}

// transfer final alpha variables to the material point
void HillPlastic::UpdatePlasticInternal(MPMBase *mptr,int np)
{	mptr->SetHistoryDble(aint);
}

#pragma mark HillPlastic:Accessors

// Return the material tag
int HillPlastic::MaterialTag(void) { return HILLPLASTIC; }

// return material type
const char *HillPlastic::MaterialType(void) { return "Elastic-Plastic Hill Material"; }

// hardening history - equivalent plastic strain (absolute strain)
double HillPlastic::GetHistory(int num,char *historyPtr)
{
    double history=0.;
	if(num==1)
	{	double *h=(double *)historyPtr;
		history=*h;
	}
    return history;
}

