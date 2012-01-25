/********************************************************************************
    AnisoPlasticity.cpp
    NairnMPM
    
    Created by John Nairn, June 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Orthotropic.hpp (TranIsotropic.hpp, MaterialBase.hpp)
********************************************************************************/

#include "AnisoPlasticity.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"

#include "NairnMPM_Class/NairnMPM.hpp"

#pragma mark AnisoPlasticity::Constructors and Destructors

// Constructors
AnisoPlasticity::AnisoPlasticity() {}

// Constructors
AnisoPlasticity::AnisoPlasticity(char *matName) : Orthotropic(matName)
{
	// negative yield stress implies not yielding in that direction
	syxx=-1.;
	syyy=-1.;
	syzz=-1.;
	tyxy=-1.;
	tyxz=-1.;
	tyyz=-1.;
}

#pragma mark AnisoPlasticity::Initialization

// Read material properties
char *AnisoPlasticity::InputMat(char *xName,int &input)
{
    if(strcmp(xName,"yldxx")==0)
    {	input=DOUBLE_NUM;
        return((char *)&syxx);
    }
    
    else if(strcmp(xName,"yldyy")==0)
    {	input=DOUBLE_NUM;
        return((char *)&syyy);
    }
	
    else if(strcmp(xName,"yldzz")==0)
    {	input=DOUBLE_NUM;
        return((char *)&syzz);
    }

    else if(strcmp(xName,"yldxy")==0)
    {	input=DOUBLE_NUM;
        return((char *)&tyxy);
    }
	
    else if(strcmp(xName,"yldxz")==0)
    {	input=DOUBLE_NUM;
        return((char *)&tyxz);
    }
	
    else if(strcmp(xName,"yldyz")==0)
    {	input=DOUBLE_NUM;
        return((char *)&tyyz);
    }
	
	return(Orthotropic::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *AnisoPlasticity::VerifyProperties(int np)
{
	// check at least some yielding
	if(syzz<0. && syxx<0. && syyy<0. && tyxy<0. && tyxz<0. && tyyz<0.)
		return "No yield stresses were defined";

	// call super class
	return Orthotropic::VerifyProperties(np);
}

// print to output window
void AnisoPlasticity::PrintMechanicalProperties(void)
{	
    Orthotropic::PrintMechanicalProperties();
	PrintYieldProperties();
}

// print just yield properties to output window
void AnisoPlasticity::PrintYieldProperties(void)
{
	if(syxx>=0.)
		PrintProperty("yld1",syxx,"");
	else
		PrintProperty("yld1= inf",false);
		
	if(syyy>=0.)
		PrintProperty("yld2",syyy,"");
	else
		PrintProperty("yld2= inf",false);
		
	if(syzz>=0.)
		PrintProperty("yld3",syzz,"");
	else
		PrintProperty("yld3= inf",false);
		
    cout << endl;

	if(tyyz>=0.)
		PrintProperty("yld23",tyyz,"");
	else
		PrintProperty("yld23= inf",false);

	if(tyxz>=0.)
		PrintProperty("yld13",tyxz,"");
	else
		PrintProperty("yld13= inf",false);
	
	if(tyxy>=0.)
		PrintProperty("yld12",tyxy,"");
	else
		PrintProperty("yld12= inf",false);
	
    cout << endl;
}

// Private properties used in constitutive law
void AnisoPlasticity::InitialLoadMechProps(int makeSpecific,int np)
{
	// reciprocals of reduced yield stresses
	if(syxx>=0.)
    {	syxxred2=rho/(syxx*1.e6);
		syxxred2*=syxxred2;
	}
	else
		syxxred2=0.;		// 1/inf^2
	if(syyy>=0.)
    {	syyyred2=rho/(syyy*1.e6); 
		syyyred2*=syyyred2;
	}
	else
		syyyred2=0.;		// 1/inf^2
	if(syzz>=0.)
	{	syzzred2=rho/(syzz*1.e6); 
		syzzred2*=syzzred2;
	}
	else
		syzzred2=0.;		// 1/inf^2
	if(tyxy)
	{	tyxyred2=rho/(tyxy*1.e6); 
		tyxyred2*=tyxyred2;
	}
	else
		tyxyred2=0.;		// 1/inf^2
	if(tyxz)
	{	tyxzred2=rho/(tyxz*1.e6); 
		tyxzred2*=tyxzred2;
	}
	else
		tyxzred2=0.;		// 1/inf^2
	if(tyyz)
	{	tyyzred2=rho/(tyyz*1.e6); 
		tyyzred2*=tyyzred2;
	}
	else
		tyyzred2=0.;		// 1/inf^2
	
	Orthotropic::InitialLoadMechProps(makeSpecific,np);
}

// if cannot be used in current analysis type throw MPMTermination()
void AnisoPlasticity::MPMConstLaw(int np)
{	if(np!=PLANE_STRAIN_MPM && np!=THREED_MPM)
		throw CommonException("Anisotropic plasticity materials require 3D or 2D plane strain MPM analysis","NairnMPM::ValidateOptions");
}

#pragma mark VonMisesHardening::Methods

/* For 2D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, plastic strain, stresses, strain energy, 
		plastic energy, dissipated energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
	This is general analysis for isotropic material. Subclass can define any
		desired plastic potential using methods GetF(), GetDfDsigma(), and GetDfDWp()
*/
void AnisoPlasticity::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
        double delTime,int np)
{
    // Effective strain by deducting thermal strain
	// (note me0[1] and mc0[1] are reduced in plane strain, but CTE3 and CME3 are not)
	double erxx=me0[1]*ConductionTask::dTemperature;
	double eryy=me0[2]*ConductionTask::dTemperature;
	double erxy=me0[3]*ConductionTask::dTemperature;
	double erzz=CTE3*ConductionTask::dTemperature;
	if(DiffusionTask::active)
	{	erxx+=mc0[1]*DiffusionTask::dConcentration;
		eryy+=mc0[2]*DiffusionTask::dConcentration;
		erxy+=mc0[3]*DiffusionTask::dConcentration;
		erzz+=CME3*DiffusionTask::dConcentration;
	}
    double dexx=dvxx-erxx;  
    double deyy=dvyy-eryy; 
    double dgxy=dvxy+dvyx-erxy;

    // Elastic stress increment
	Tensor *sp=mptr->GetStressTensor();
    Tensor stk=*sp;
    stk.xx += mdm[1][1]*dexx+mdm[1][2]*deyy+mdm[1][3]*dgxy;
    stk.yy += mdm[1][2]*dexx+mdm[2][2]*deyy+mdm[2][3]*dgxy;
    stk.xy += mdm[1][3]*dexx+mdm[2][3]*deyy+mdm[3][3]*dgxy;
    if(np==PLANE_STRAIN_MPM)
	{	stk.zz += mdm[4][1]*(dexx+me0[5]*erzz) + mdm[4][2]*(deyy+me0[6]*erzz) + mdm[4][3]*(dgxy+me0[7]*erzz) - mdm[4][4]*erzz;
	}
  
	// get rotation matrix elements
	double angle=mptr->GetRotationZ();
	double c=cos(angle);
	cos2t=c*c;
	double s=sin(angle);
	sin2t=s*s;
	costsint=c*s;

    // Calculate plastic potential f
	UpdateTrialAlpha(mptr,np);
	double ftrial = GetF(mptr,&stk,np);
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
		Elastic::MPMConstLaw(mptr,dvxx,dvyy,dvxy,dvyx,delTime,np); 
		return; 
    }
    
	// Find  lambda for this plastic state
	// This material overrides to have custom solution
	double lambda = SolveForLambdaAP(mptr,np,ftrial,&stk);
    
    // Plastic strain increments on particle
    double dexxp=lambda*dfdsxx;
    double deyyp=lambda*dfdsyy;
    double dgxyp=lambda*dfdtxy;
	Tensor *eplast=mptr->GetPlasticStrainTensor();
    eplast->xx+=dexxp;
    eplast->yy+=deyyp;
    eplast->xy+=dgxyp;
	double dezzp,dezz;
    if(np==PLANE_STRAIN_MPM)
	{	// elastic strain in plane strain is -dezzp
		dezzp=lambda*dfdszz;
		eplast->zz+=dezzp;
		dezz=-dezzp;
    }
	else
	{	// not correct, or used, now
		dezzp=0.;
		dezz=mdm[4][1]*(dexx-dexxp) + mdm[4][2]*(deyy-deyyp) + mdm[4][3]*(dgxy-dgxyp) + erzz;
	}
    
    // Elastic strain increments on particle
	Tensor *ep=mptr->GetStrainTensor();
    ep->xx+=(dvxx-dexxp);
    ep->yy+=(dvyy-deyyp);
    ep->xy+=(dvxy+dvyx-dgxyp);
    ep->zz+=dezz;
	
	// rotational strain
	double dwrotxy=dvyx-dvxy;
	
	// Effective elastic strain increment = de - deres
	dexx -= dexxp;
	deyy -= deyyp;
	dgxy -= dgxyp;
	dezz -= erzz;
	
	// increment particle stresses
    // Hypoelastic - Eliminate effect of rotation 
    Tensor st0=*sp;			// save previous stress for energy updates below
    double c1 = mdm[1][1]*dexx+mdm[1][2]*deyy+mdm[1][3]*dgxy;
    double c2 = mdm[1][2]*dexx+mdm[2][2]*deyy+mdm[2][3]*dgxy;
    double c3 = mdm[1][3]*dexx+mdm[2][3]*deyy+mdm[3][3]*dgxy;
	Hypo2DCalculations(mptr,-dwrotxy,c1,c2,c3);
	
	// out of plane stress
    if(np==PLANE_STRAIN_MPM)
	{	sp->zz += mdm[4][1]*(dexx+me0[5]*erzz) + mdm[4][2]*(deyy+me0[6]*erzz) + mdm[4][3]*(dgxy+me0[7]*erzz) + mdm[4][4]*dezz;
	}

	// Elastic energy density
    mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexx
                        + (st0.yy+sp->yy)*deyy
                        + (st0.xy+sp->xy)*dgxy));

    // Plastic energy increment
	double dispEnergy=0.5*((st0.xx+sp->xx)*dexxp
                        + (st0.yy+sp->yy)*deyyp
                        + (st0.xy+sp->xy)*dgxyp);

    // Extra term for plain strain when there are thermal stresses
    if(np==PLANE_STRAIN_MPM)
    {   dispEnergy+=0.5*(st0.zz+sp->zz)*dezzp;
        mptr->AddStrainEnergy(0.5*(st0.zz+sp->zz)*dezz);
    }
	
	// add dissipated and plastic energy to the particle
	mptr->AddDispEnergy(dispEnergy);
    mptr->AddPlastEnergy(dispEnergy);
	
	// update internal variables
	UpdatePlasticInternal(mptr,np);
}

/* For 3D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
	Assumes linear elastic, uses hypoelastic correction
*/
void AnisoPlasticity::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
								double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
    // Effective strain by deducting thermal strain
	double erxx=me0[0]*ConductionTask::dTemperature;
	double eryy=me0[1]*ConductionTask::dTemperature;
	double erzz=me0[2]*ConductionTask::dTemperature;
	double eryz=me0[3]*ConductionTask::dTemperature;
	double erxz=me0[4]*ConductionTask::dTemperature;
	double erxy=me0[5]*ConductionTask::dTemperature;
	if(DiffusionTask::active)
	{	erxx+=mc0[0]*DiffusionTask::dConcentration;
		eryy+=mc0[1]*DiffusionTask::dConcentration;
		erzz+=mc0[2]*DiffusionTask::dConcentration;
		eryz+=mc0[3]*DiffusionTask::dConcentration;
		erxz+=mc0[4]*DiffusionTask::dConcentration;
		erxy+=mc0[5]*DiffusionTask::dConcentration;
	}
    double dexx=dvxx-erxx;  
    double deyy=dvyy-eryy; 
	double dezz=dvzz-erzz;
    double dgyz=dvyz+dvzy-eryz;
    double dgxz=dvxz+dvzx-erxz;
    double dgxy=dvxy+dvyx-erxy;
	
    // Elastic stress increment
	Tensor *sp=mptr->GetStressTensor();
    Tensor stk=*sp;
    stk.xx += mdm[0][0]*dexx+mdm[0][1]*deyy+mdm[0][2]*dezz+mdm[0][3]*dgyz+mdm[0][4]*dgxz+mdm[0][5]*dgxy;
    stk.yy += mdm[1][0]*dexx+mdm[1][1]*deyy+mdm[1][2]*dezz+mdm[1][3]*dgyz+mdm[1][4]*dgxz+mdm[1][5]*dgxy;
    stk.zz += mdm[2][0]*dexx+mdm[2][1]*deyy+mdm[2][2]*dezz+mdm[2][3]*dgyz+mdm[2][4]*dgxz+mdm[2][5]*dgxy;
    stk.yz += mdm[3][0]*dexx+mdm[3][1]*deyy+mdm[3][2]*dezz+mdm[3][3]*dgyz+mdm[3][4]*dgxz+mdm[3][5]*dgxy;
    stk.xz += mdm[4][0]*dexx+mdm[4][1]*deyy+mdm[4][2]*dezz+mdm[4][3]*dgyz+mdm[4][4]*dgxz+mdm[4][5]*dgxy;
    stk.xy += mdm[5][0]*dexx+mdm[5][1]*deyy+mdm[5][2]*dezz+mdm[5][3]*dgyz+mdm[5][4]*dgxz+mdm[5][5]*dgxy;
	
	// get rotation matrix elements rzyx[i][j]
	double z=mptr->GetRotationZ();
	double cz=cos(z);
	double sz=sin(z);
	double cz2=cz*cz;
	double sz2=sz*sz;
	double c2z=cos(2.*z);
	double s2z=sin(2.*z);
	
	double y=mptr->GetRotationY();
	double cy=cos(y);
	double sy=sin(y);
	double cy2=cy*cy;
	double sy2=sy*sy;
	double s2y=sin(2.*y);
	double c2y=cos(2.*y);
	
	double x=mptr->GetRotationX();
	double cx=cos(x);
	double sx=sin(x);
	double cx2=cx*cx;
	double sx2=sx*sx;
	double c2x=cos(2.*x);
	double s2x=sin(2.*x);
	
	// Explicit equations form Mathematica for speed in evaluation
	// Trig terms done once above
	rzyx[0][0]=cy2*cz2;
	rzyx[0][1]=cy2*sz2;
	rzyx[0][2]=sy2;
	rzyx[0][3]=-2*cy*sy*sz;
	rzyx[0][4]=cz*s2y;
	rzyx[0][5]=-2*cy2*cz*sz;
	
	double arg1=cz*sx*sy + cx*sz;
	rzyx[1][0]=arg1*arg1;
	double arg2=cx*cz - sx*sy*sz;
	rzyx[1][1]=arg2*arg2;
	rzyx[1][2]=cy2*sx2;
	rzyx[1][3]=-2*cy*sx*arg2;
	rzyx[1][4]=-2*cy*sx*arg1;
	rzyx[1][5]=cx2*s2z + c2z*s2x*sy - 2*cz*sx2*sy2*sz;

	double arg3=cx*cz*sy - sx*sz;
	rzyx[2][0]=arg3*arg3;
	double arg4=cz*sx + cx*sy*sz;
	rzyx[2][1]=arg4*arg4;
	rzyx[2][2]=cx2*cy2;
	rzyx[2][3]=2*cx*cy*arg4;
	rzyx[2][4]=-2*cx*cy*arg3;
	rzyx[2][5]=s2z*sx2 - c2z*s2x*sy - 2*cx2*cz*sy2*sz;
	
	double arg5=c2x*s2z*sy/2.;
	rzyx[3][0]=-arg5 + cx*sx*(sz2 - cz2*sy2);
	rzyx[3][1]=arg5 + cx*sx*(cz2 - sy2*sz2);
	rzyx[3][2]=-cx*cy2*sx;
	rzyx[3][3]=cy*(c2x*cz - 2*cx*sx*sy*sz);
	rzyx[3][4]=cy*(cz*s2x*sy + c2x*sz);
	rzyx[3][5]=-c2x*c2z*sy + 2*cx*cz*sx*(1 + sy2)*sz;
	
	rzyx[4][0]=-cy*cz*arg3;
	rzyx[4][1]=-cy*sz*arg4;
	rzyx[4][2]=cx*cy*sy;
	rzyx[4][3]=cz*sx*sy - c2y*cx*sz;
	rzyx[4][4]=c2y*cx*cz + sx*sy*sz;
	rzyx[4][5]=cy*(c2z*sx + cx*s2z*sy);
	
	rzyx[5][0]=cy*cz*arg1;
	rzyx[5][1]=-cy*sz*arg2;
	rzyx[5][2]=-cy*sx*sy;
	rzyx[5][3]=cx*cz*sy + c2y*sx*sz;
	rzyx[5][4]=cx*sy*sz - c2y*cz*sx;
	rzyx[5][5]=cy*(c2z*cx - 2*cz*sx*sy*sz);
	
    // Calculate plastic potential f
	UpdateTrialAlpha(mptr,np);
	double ftrial = GetF(mptr,&stk,np);
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
		Elastic::MPMConstLaw(mptr,dvxx,dvyy,dvzz,dvxy,dvyx,dvxz,dvzx,dvyz,dvzy,delTime,np);
		return; 
    }
    
	//int keyParticle=-1;
	//if(MPMBase::currentParticleNum==keyParticle)
	//{	cout << "# STEP " << fmobj->mstep << ", f = " << ftrial << endl;
	//}
	
	// Find  lambda for this plastic state
	// This material overrides to custom solution
	double lambda = SolveForLambdaAP(mptr,np,ftrial,&stk);
    
    // Plastic strain increments on particle
    double dexxp=lambda*dfdsxx;
    double deyyp=lambda*dfdsyy;
	double dezzp=lambda*dfdszz;
    double dgyzp=lambda*dfdtyz;
    double dgxzp=lambda*dfdtxz;
    double dgxyp=lambda*dfdtxy;
	Tensor *eplast=mptr->GetPlasticStrainTensor();
    eplast->xx+=dexxp;
    eplast->yy+=deyyp;
	eplast->zz+=dezzp;
    eplast->yz+=dgyzp;
    eplast->xz+=dgxzp;
    eplast->xy+=dgxyp;
    
    // Elastic strain increments on particle
	Tensor *ep=mptr->GetStrainTensor();
    ep->xx+=(dvxx-dexxp);
    ep->yy+=(dvyy-deyyp);
    ep->zz+=(dvzz-dezzp);
    ep->xy+=(dvyz+dvzy-dgyzp);
    ep->xy+=(dvxz+dvzx-dgxzp);
    ep->xy+=(dvxy+dvyx-dgxyp);
	
	// rotational strain increments (particle updated by Hypo3D)
	double dwrotyz=dvzy-dvyz;
	double dwrotxz=dvzx-dvxz;
	double dwrotxy=dvyx-dvxy;
	
	// Effective elastic strain increment = de - deres
	dexx -= dexxp;
	deyy -= deyyp;
	dezz -= dezzp;
	dgyz -= dgyzp;
	dgxz -= dgxzp;
	dgxy -= dgxyp;
	
	// increment particle stresses
	Tensor st0=*sp;			// save previous stress for energy updates below
	double dsig[6];
    dsig[XX] = mdm[0][0]*dexx+mdm[0][1]*deyy+mdm[0][2]*dezz+mdm[0][3]*dgyz+mdm[0][4]*dgxz+mdm[0][5]*dgxy;
    dsig[YY] = mdm[1][0]*dexx+mdm[1][1]*deyy+mdm[1][2]*dezz+mdm[1][3]*dgyz+mdm[1][4]*dgxz+mdm[1][5]*dgxy;
    dsig[ZZ] = mdm[2][0]*dexx+mdm[2][1]*deyy+mdm[2][2]*dezz+mdm[2][3]*dgyz+mdm[2][4]*dgxz+mdm[2][5]*dgxy;
    dsig[YZ] = mdm[3][0]*dexx+mdm[3][1]*deyy+mdm[3][2]*dezz+mdm[3][3]*dgyz+mdm[3][4]*dgxz+mdm[3][5]*dgxy;
    dsig[XZ] = mdm[4][0]*dexx+mdm[4][1]*deyy+mdm[4][2]*dezz+mdm[4][3]*dgyz+mdm[4][4]*dgxz+mdm[4][5]*dgxy;
    dsig[XY] = mdm[5][0]*dexx+mdm[5][1]*deyy+mdm[5][2]*dezz+mdm[5][3]*dgyz+mdm[5][4]*dgxz+mdm[5][5]*dgxy;
	Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,dsig);

    // Elastic energy density
	mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexx
							   + (st0.yy+sp->yy)*deyy
							   + (st0.zz+sp->zz)*dezz
							   + (st0.yz+sp->yz)*dgyz
							   + (st0.xz+sp->xz)*dgxz
							   + (st0.xy+sp->xy)*dgxy));
	
    // Plastic energy increment
	double dispEnergy=0.5*(0.5*((st0.xx+sp->xx)*dexxp
								+ (st0.yy+sp->yy)*deyyp
								+ (st0.zz+sp->zz)*dezzp
								+ (st0.yz+sp->yz)*dgyzp
								+ (st0.xz+sp->xz)*dgxzp
								+ (st0.xy+sp->xy)*dgxyp));
	
	// add plastic energy to the particle
	mptr->AddDispEnergy(dispEnergy);
    mptr->AddPlastEnergy(dispEnergy);
	
	// update internal variables
	UpdatePlasticInternal(mptr,np);
	
	//if(MPMBase::currentParticleNum==keyParticle)
	//{	cout << "#   lam = " << lambda << ", cum ep = " << mptr->GetHistoryDble() << endl;
	//}
}

#pragma mark AnisoPlasticity::Custom Methods

// Solve numerically for lambda
// Ouptut is lambdak, final df, final alpha
double AnisoPlasticity::SolveForLambdaAP(MPMBase *mptr,int np,double ftrial,Tensor *strial)
{
	// step 0 = stress is strial, alpha is previous alpha
	GetDfCdf(mptr,strial,np);								// find df/dsigma
	double lambda2=ftrial/(dfCdf + GetDfAlphaDotH(mptr,np,strial));
	double lambdaInitial=lambda2;
	
	// find df/dsigma and -h at this initial guessed stress rather than at strial as above
	Tensor stk;
	UpdateStress(strial,&stk,lambda2,np);
	GetDfCdf(mptr,&stk,np);
	UpdateStress(strial,&stk,lambda2,np);			// second pass might help
	GetDfCdf(mptr,&stk,np);
	
	// Find f using new slopes at lambda1
	UpdateTrialAlpha(mptr,np,lambda2);
	double f2=GetF(mptr,&stk,np);
	
	// pick second lambda that is lower
	double lambda1=0.5*lambda2;
	double f1=GetFkFromLambdak(mptr,strial,&stk,lambda1,np);
	
	//if(MPMBase::currentParticleNum==5018) cout << "# initial (lambda1,f1) = (" << lambda1 << "," << f1 << ") and (lambda2,f2) = (" << lambda2 << "," << f2 << ")" << endl;
	
	// bracket the solution
	int step=1;
	while(true)
	{	if(f1*f2<0.) break;
		if(fabs(f1)<fabs(f2))
		{	lambda1+=1.6*(lambda1-lambda2);
			if(lambda1>0.)
				f1=GetFkFromLambdak(mptr,strial,&stk,lambda1,np);
			else
			{	f1=ftrial;
				lambda1=0.;
			}
		}
		else
		{	lambda2+=1.6*(lambda2-lambda1);
			f2=GetFkFromLambdak(mptr,strial,&stk,lambda2,np);
		}
		
		// if fails to bracket in 50 tries, return with single step solution
		step++;
		if(step>50)
		{	//if(MPMBase::currentParticleNum==5018) cout << "# Failed to bracket solution" << endl;
			GetDfCdf(mptr,strial,np);
			UpdateStress(strial,&stk,lambdaInitial,np);
			GetDfCdf(mptr,&stk,np);
			UpdateTrialAlpha(mptr,np,lambdaInitial);
			return lambdaInitial;
		}
	}
	
	//if(MPMBase::currentParticleNum==5018) cout << "# bracketed (lambda1,f1) = (" << lambda1 << "," << f1 << ") and (lambda2,f2) = (" << lambda2 << "," << f2 << "), steps = " << step << endl;
	
	// order solutions so f1<0
	if(f1>0.)
	{	double ftemp=f1;
		double lambdatemp=lambda1;
		f1=f2;
		lambda1=lambda2;
		f2=ftemp;
		lambda2=lambdatemp;
	}
	
	// set up
	double lambdak=0.5*(lambda1+lambda2);				// initial guess
	double dxold=fabs(lambda2-lambda1);					// the step size before last
	double dx=dxold;									// the last step sie
	
	// initial guess
	double fk=GetFkFromLambdak(mptr,strial,&stk,lambdak,np);
	double dfkdlam=-(dfCdf + GetDfAlphaDotH(mptr,np,&stk));
	
	// Loop
	step=1;
	while(true)
	{	// if out of range, or not decreasing fast enough, use bisection method
		if((((lambdak-lambda2)*dfkdlam-fk)*((lambdak-lambda1)*dfkdlam-fk) >= 0.0)
		   || (fabs(2.0*fk) > fabs(dxold*dfkdlam)))
		{	dxold=dx;
			dx=0.5*(lambda2-lambda1);
			lambdak=lambda1+dx;
		}
		else
		{	dxold=dx;
			dx=fk/dfkdlam;
			lambdak-=dx;
		}
		
		// convergence check
		if(step++>20 || fabs(dx/lambdak)<0.001) break;
		
		// next value
		fk=GetFkFromLambdak(mptr,strial,&stk,lambdak,np);
		dfkdlam=-(dfCdf + GetDfAlphaDotH(mptr,np,&stk));
		
		// maintain the bracket
		if(fk<0.0)
		{	lambda1=lambdak;
			f1=fk;
		}
		else
		{	lambda2=lambdak;
			f2=fk;
		}
		
		//if(MPMBase::currentParticleNum==5018) cout << "# ..." << step << "," << lambdak << "," << fk << endl;
	}
	
	// output df (from initial setting of GetDfCdf()), alpha (here with latest lambda), and lamda (the return vale)
	UpdateTrialAlpha(mptr,np,lambdak);
	
	//UpdateStress(strial,&stk,lambdak,np);
	//if(MPMBase::currentParticleNum==5018) cout << "# ..." << step << "," << lambdak << "," << GetF(mptr,&stk,np) << endl;
	
	return lambdak;
}

// Find C.df and df.C.df at given stress - stored in material variables active only during the loop
void AnisoPlasticity::GetDfCdf(MPMBase *mptr,Tensor *stk,int np)
{
	// get C df and df C df, which need df/dsig (Function of yield criteria, and normally only current stress in stk)
	GetDfDsigma(mptr,stk,np);
	if(np==THREED_MPM)
	{	Cdfxx = mdm[0][0]*dfdsxx+mdm[0][1]*dfdsyy+mdm[0][2]*dfdszz+mdm[0][3]*dfdtyz+mdm[0][4]*dfdtxz+mdm[0][5]*dfdtxy;
		Cdfyy = mdm[1][0]*dfdsxx+mdm[1][1]*dfdsyy+mdm[1][2]*dfdszz+mdm[1][3]*dfdtyz+mdm[1][4]*dfdtxz+mdm[1][5]*dfdtxy;
		Cdfzz = mdm[2][0]*dfdsxx+mdm[2][1]*dfdsyy+mdm[2][2]*dfdszz+mdm[2][3]*dfdtyz+mdm[2][4]*dfdtxz+mdm[2][5]*dfdtxy;
		Cdfyz = mdm[3][0]*dfdsxx+mdm[3][1]*dfdsyy+mdm[3][2]*dfdszz+mdm[3][3]*dfdtyz+mdm[3][4]*dfdtxz+mdm[3][5]*dfdtxy;
		Cdfxz = mdm[4][0]*dfdsxx+mdm[4][1]*dfdsyy+mdm[4][2]*dfdszz+mdm[4][3]*dfdtyz+mdm[4][4]*dfdtxz+mdm[4][5]*dfdtxy;
		Cdfxy = mdm[5][0]*dfdsxx+mdm[5][1]*dfdsyy+mdm[5][2]*dfdszz+mdm[5][3]*dfdtyz+mdm[5][4]*dfdtxz+mdm[5][5]*dfdtxy;
		dfCdf = dfdsxx*Cdfxx + dfdsyy*Cdfyy + dfdszz*Cdfzz + dfdtyz*Cdfyz + dfdtxz*Cdfxz + dfdtxy*Cdfxy;
	}
	else
	{	Cdfxx = mdm[1][1]*dfdsxx + mdm[1][2]*dfdsyy + mdm[1][3]*dfdtxy;
		Cdfyy = mdm[1][2]*dfdsxx + mdm[2][2]*dfdsyy + mdm[2][3]*dfdtxy;
		Cdfxy = mdm[1][3]*dfdsxx + mdm[2][3]*dfdsyy + mdm[3][3]*dfdtxy;
		dfCdf = dfdsxx*Cdfxx + dfdsyy*Cdfyy + dfdtxy*Cdfxy;
		if(np==PLANE_STRAIN_MPM)
		{	Cdfzz = mdm[4][1]*dfdsxx + mdm[4][2]*dfdsyy + mdm[4][3]*dfdtxy + mdm[4][4]*dfdszz;
			dfCdf += dfdszz*Cdfzz;
		}
		else
			Cdfzz=0.;
	}
}

// Update stress (assumes Cdfij calculated before)
void AnisoPlasticity::UpdateStress(Tensor *strial,Tensor *stk,double lambda,int np)
{
	stk->xx = strial->xx - lambda*Cdfxx;
	stk->yy = strial->yy - lambda*Cdfyy;
	stk->xy = strial->xy - lambda*Cdfxy;
	if(np==PLANE_STRAIN_MPM)
		stk->zz = strial->zz - lambda*Cdfzz;
	else if(np==THREED_MPM)
	{	stk->zz = strial->zz - lambda*Cdfzz;
		stk->yz = strial->yz - lambda*Cdfyz;
		stk->xz = strial->xz - lambda*Cdfxz;
	}
}

// Given stress and lambda, find f (also find the stress and return dfCdf)
double AnisoPlasticity::GetFkFromLambdak(MPMBase *mptr,Tensor *strial,Tensor *stk,double lambda,int np)
{
	// Change stress to new value based on strial, lambda, and previous slope
	UpdateStress(strial,stk,lambda,np);
	
	// update alpha using new lambda and -h from most recent GetDfDsigma()
	UpdateTrialAlpha(mptr,np,lambda);
	
	// update fk using new stress and alpha
	return GetF(mptr,stk,np);
}

