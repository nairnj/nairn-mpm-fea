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
// 2D only. 3D in separate method for speed
double HillPlastic::GetF(MPMBase *mptr,Tensor *st0,int np)
{	// special case for 3D
	if(np==THREED_MPM) return GetF3D(mptr,st0);
	
	// clockwise rotation from analysis to material axes
	sxrot=st0->xx*cos2t + st0->yy*sin2t - 2.*st0->xy*costsint;
	syrot=st0->xx*sin2t + st0->yy*cos2t + 2.*st0->xy*costsint;
	txyrot=(st0->xx-st0->yy)*costsint + st0->xy*(cos2t-sin2t);
	
	// the potential
	//double sAs = sxrot*sxrot*syxxred2 + syrot*syrot*syyyred2 + txyrot*txyrot*tyxyred2 - 2.*hTerm*sxrot*syrot;
	//if(np==PLANE_STRAIN_MPM)
	//	sAs += st0->zz*st0->zz*syzzred2 - 2.*fTerm*syrot*st0->zz - 2.*gTerm*sxrot*st0->zz;
	double dxy=sxrot-syrot,sAs;
	if(np==PLANE_STRAIN_MPM)
	{	double dyz=syrot-st0->zz;
		double dxz=sxrot-st0->zz;
		sAs=fTerm*dyz*dyz + gTerm*dxz*dxz + hTerm*dxy*dxy + txyrot*txyrot*tyxyred2;
	}
	else
		sAs=fTerm*syrot*syrot + gTerm*sxrot*sxrot + hTerm*dxy*dxy + txyrot*txyrot*tyxyred2;
	
	// return sqrt(sAs) - pow(1. + Khard*aint,nhard);
	// check on negative sAs can happen due to round-off error when stresses near zero
	return sAs>0. ? sqrt(sAs) - (1. + Khard*pow(aint,nhard)) : -1.; 
}

// return the yield function in a subroutine so sub class can define new yield functions
double HillPlastic::GetF3D(MPMBase *mptr,Tensor *st0)
{
	// rotation from analysis to material axes
	Tensor srot;
	srot.xx = rzyx[0][0]*st0->xx+rzyx[0][1]*st0->yy+rzyx[0][2]*st0->zz+rzyx[0][3]*st0->yz+rzyx[0][4]*st0->xz+rzyx[0][5]*st0->xy;
	srot.yy = rzyx[1][0]*st0->xx+rzyx[1][1]*st0->yy+rzyx[1][2]*st0->zz+rzyx[1][3]*st0->yz+rzyx[1][4]*st0->xz+rzyx[1][5]*st0->xy;
	srot.zz = rzyx[2][0]*st0->xx+rzyx[2][1]*st0->yy+rzyx[2][2]*st0->zz+rzyx[2][3]*st0->yz+rzyx[2][4]*st0->xz+rzyx[2][5]*st0->xy;
	srot.yz = rzyx[3][0]*st0->xx+rzyx[3][1]*st0->yy+rzyx[3][2]*st0->zz+rzyx[3][3]*st0->yz+rzyx[3][4]*st0->xz+rzyx[3][5]*st0->xy;
	srot.xz = rzyx[4][0]*st0->xx+rzyx[4][1]*st0->yy+rzyx[4][2]*st0->zz+rzyx[4][3]*st0->yz+rzyx[4][4]*st0->xz+rzyx[4][5]*st0->xy;
	srot.xy = rzyx[5][0]*st0->xx+rzyx[5][1]*st0->yy+rzyx[5][2]*st0->zz+rzyx[5][3]*st0->yz+rzyx[5][4]*st0->xz+rzyx[5][5]*st0->xy;

	// the potential
	double dyz=srot.yy-srot.zz;
	double dxz=srot.xx-srot.zz;
	double dxy=srot.xx-srot.yy;
	double sAs=fTerm*dyz*dyz + gTerm*dxz*dxz + hTerm*dxy*dxy + srot.xy*srot.xy*tyxyred2
				 + srot.xz*srot.xz*tyxzred2 + srot.yz*srot.yz*tyyzred2;
	
	// return sqrt(sAs) - pow(1. + Khard*aint,nhard);
	// check on negative sAs can happen due to round-off error when stresses near zero
	return sAs>0. ? sqrt(sAs) - (1. + Khard*pow(aint,nhard)) : -1.; 
}

// return derivatives of the yield function in a subroutine so sub class can define new yield functions
void HillPlastic::GetDfDsigma(MPMBase *mptr,Tensor *st0,int np)
{	// special case for 3D
	if(np==THREED_MPM)
	{	GetDfDsigma3D(mptr,st0);
		return;
	}
	
	// clockwise rotation from analysis to material axes
	sxrot = st0->xx*cos2t + st0->yy*sin2t - 2.*st0->xy*costsint;
	syrot = st0->xx*sin2t + st0->yy*cos2t + 2.*st0->xy*costsint;
	txyrot = (st0->xx-st0->yy)*costsint + st0->xy*(cos2t-sin2t);
	
	// sqrt(s As)
	double sAs = sxrot*sxrot*syxxred2 + syrot*syrot*syyyred2 + txyrot*txyrot*tyxyred2 - 2.*hTerm*sxrot*syrot;
	if(np==PLANE_STRAIN_MPM)
		sAs += st0->zz*st0->zz*syzzred2 - 2.*fTerm*syrot*st0->zz - 2.*gTerm*sxrot*st0->zz;
	
	if(sAs>0.)
	{	double rootSAS=sqrt(sAs);
	
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
	
	else
	{	dfdsxxrot = dfdsyyrot = dfdszzrot = dfdtxyrot = 0.;
		dfdsxx = dfdsyy = dfdszz = dfdtxy = 0.;
		minush=0.;
	}
}

// return derivatives of the yield function in a subroutine so sub class can define new yield functions
void HillPlastic::GetDfDsigma3D(MPMBase *mptr,Tensor *st0)
{
	// rotation from analysis to material axes
	Tensor srot;
	srot.xx = rzyx[0][0]*st0->xx+rzyx[0][1]*st0->yy+rzyx[0][2]*st0->zz+rzyx[0][3]*st0->yz+rzyx[0][4]*st0->xz+rzyx[0][5]*st0->xy;
	srot.yy = rzyx[1][0]*st0->xx+rzyx[1][1]*st0->yy+rzyx[1][2]*st0->zz+rzyx[1][3]*st0->yz+rzyx[1][4]*st0->xz+rzyx[1][5]*st0->xy;
	srot.zz = rzyx[2][0]*st0->xx+rzyx[2][1]*st0->yy+rzyx[2][2]*st0->zz+rzyx[2][3]*st0->yz+rzyx[2][4]*st0->xz+rzyx[2][5]*st0->xy;
	srot.yz = rzyx[3][0]*st0->xx+rzyx[3][1]*st0->yy+rzyx[3][2]*st0->zz+rzyx[3][3]*st0->yz+rzyx[3][4]*st0->xz+rzyx[3][5]*st0->xy;
	srot.xz = rzyx[4][0]*st0->xx+rzyx[4][1]*st0->yy+rzyx[4][2]*st0->zz+rzyx[4][3]*st0->yz+rzyx[4][4]*st0->xz+rzyx[4][5]*st0->xy;
	srot.xy = rzyx[5][0]*st0->xx+rzyx[5][1]*st0->yy+rzyx[5][2]*st0->zz+rzyx[5][3]*st0->yz+rzyx[5][4]*st0->xz+rzyx[5][5]*st0->xy;
	
	// sqrt(s.As)
	double dyz=srot.yy-srot.zz;
	double dxz=srot.xx-srot.zz;
	double dxy=srot.xx-srot.yy;
	double sAs=fTerm*dyz*dyz + gTerm*dxz*dxz + hTerm*dxy*dxy + srot.xy*srot.xy*tyxyred2
						+ srot.xz*srot.xz*tyxzred2 + srot.yz*srot.yz*tyyzred2;
	if(sAs>0.)
	{	double rootSAS=sqrt(sAs);
	
		// the derivatives = dfrot = A srot/sqrt(sAS)
		dfdsxxrot = (syxxred2*srot.xx - hTerm*srot.yy - gTerm*srot.zz) / rootSAS;
		dfdsyyrot = (-hTerm*srot.xx + syyyred2*srot.yy - fTerm*srot.zz) / rootSAS;
		dfdszzrot = (-gTerm*srot.xx - fTerm*srot.yy + syzzred2*srot.zz) / rootSAS;
		dfdtyzrot = tyyzred2*srot.yz / rootSAS;
		dfdtxzrot = tyxzred2*srot.xz / rootSAS;
		dfdtxyrot = tyxyred2*srot.xy / rootSAS;
		
		// rotate to analysis coordinates df = R^T dfrot
		dfdsxx = rzyx[0][0]*dfdsxxrot+rzyx[1][0]*dfdsyyrot+rzyx[2][0]*dfdszzrot+rzyx[3][0]*dfdtyzrot+rzyx[4][0]*dfdtxzrot+rzyx[5][0]*dfdtxyrot;
		dfdsyy = rzyx[0][1]*dfdsxxrot+rzyx[1][1]*dfdsyyrot+rzyx[2][1]*dfdszzrot+rzyx[3][1]*dfdtyzrot+rzyx[4][1]*dfdtxzrot+rzyx[5][1]*dfdtxyrot;
		dfdszz = rzyx[0][2]*dfdsxxrot+rzyx[1][2]*dfdsyyrot+rzyx[2][2]*dfdszzrot+rzyx[3][2]*dfdtyzrot+rzyx[4][2]*dfdtxzrot+rzyx[5][2]*dfdtxyrot;
		dfdtyz = rzyx[0][3]*dfdsxxrot+rzyx[1][3]*dfdsyyrot+rzyx[2][3]*dfdszzrot+rzyx[3][3]*dfdtyzrot+rzyx[4][3]*dfdtxzrot+rzyx[5][3]*dfdtxyrot;
		dfdtxz = rzyx[0][4]*dfdsxxrot+rzyx[1][4]*dfdsyyrot+rzyx[2][4]*dfdszzrot+rzyx[3][4]*dfdtyzrot+rzyx[4][4]*dfdtxzrot+rzyx[5][4]*dfdtxyrot;
		dfdtxy = rzyx[0][5]*dfdsxxrot+rzyx[1][5]*dfdsyyrot+rzyx[2][5]*dfdszzrot+rzyx[3][5]*dfdtyzrot+rzyx[4][5]*dfdtxzrot+rzyx[5][5]*dfdtxyrot;
		
		// for use in alpha upate
		minush = dfdsxxrot*dfdsxxrot + dfdsyyrot*dfdsyyrot + dfdszzrot*dfdszzrot 
						+ 0.5*(dfdtyzrot*dfdtyzrot + dfdtxzrot*dfdtxzrot + dfdtxyrot*dfdtxyrot);
		minush = sqrt(minush/1.5);
	}
	
	else
	{	dfdsxxrot = dfdsyyrot = dfdszzrot = dfdtyzrot = dfdtxzrot = dfdtxyrot = 0.;
		dfdsxx = dfdsyy = dfdszz = dfdtyz = dfdtxz = dfdtxy = 0.;
		minush=0.;
	}
	
}

// Non-linear hardening : df^(alpha) . h and it assumes g(alpha) = 1 + Khard alpha^nhard
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

