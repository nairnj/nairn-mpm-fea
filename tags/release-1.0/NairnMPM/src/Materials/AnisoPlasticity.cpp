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
#include "Read_XML/CommonReadHandler.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"

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
	
	return(Orthotropic::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *AnisoPlasticity::VerifyProperties(int np)
{
	// check at least some yielding
	if(syzz<0. && syxx<0. && syyy<0. && tyxy<0.)
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
	
	Orthotropic::InitialLoadMechProps(makeSpecific,np);
}

// if cannot be used in current analysis type throw MPMTermination()
void AnisoPlasticity::ValidateUse(int np)
{	if(np!=PLANE_STRAIN_MPM)
		throw CommonException("Anisotropic plasticity materials require plane strain MPM analysis","NairnMPM::ValidateOptions");
	Orthotropic::ValidateUse(np);
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
	double angle=mptr->GetRotation();
	double c=cos(angle);
	cos2t=c*c;
	double s=sin(angle);
	sin2t=s*s;
	costsint=c*s;

    // Calculate plastic potential f
	UpdateTrialAlpha(mptr,np);
	double ftrial = GetF(mptr,stk.xx,stk.yy,stk.xy,stk.zz,np);
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
		Elastic::MPMConstLaw(mptr,dvxx,dvyy,dvxy,dvyx,delTime,np); 
		return; 
    }
    
	// Find  lambda for this plastic state
	// Base class finds it numerically, subclass can override if solvable
	double lambda = SolveForLambda(mptr,np,ftrial,&stk);
    
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

#pragma mark AnisoPlasticity::Custom Methods

// Solve numerically for lambda
// Subclasses can override for alternative solution if possible
// Ouptut is lambdak, final stress, final df, final alpha
double AnisoPlasticity::SolveForLambda(MPMBase *mptr,int np,double fk,Tensor *stk)
{
	double lambdak=0.;
	int step=1;
	Tensor strial=*stk;
	
	// cout << "# initial f " << fk << endl;
	while(true)
	{	// For current values
		GetDfDsigma(mptr,stk,np);
		double dfdAlpha=GetDfAlphaDotH(mptr,np,stk);

		// get df/dlam = - ( df C df + dfa h )
		double Cdfxx = mdm[1][1]*dfdsxx + mdm[1][2]*dfdsyy + mdm[1][3]*dfdtxy;
		double Cdfyy = mdm[1][2]*dfdsxx + mdm[2][2]*dfdsyy + mdm[2][3]*dfdtxy;
		double Cdfxy = mdm[1][3]*dfdsxx + mdm[2][3]*dfdsyy + mdm[3][3]*dfdtxy;
		double Cdfzz = 0.;
		double dfCdf = dfdsxx*Cdfxx + dfdsyy*Cdfyy + dfdtxy*Cdfxy;
		if(np==PLANE_STRAIN_MPM)
		{	Cdfzz = mdm[4][1]*dfdsxx + mdm[4][2]*dfdsyy + mdm[4][3]*dfdtxy + mdm[4][4]*dfdszz;
			dfCdf += dfdszz*Cdfzz;
		}
		double dfdlam = - (dfCdf + dfdAlpha);
		
		// increment lambda and alpha
		double delLam = -fk/dfdlam;
		lambdak += delLam;
		UpdateTrialAlpha(mptr,np,lambdak);
		
		// update stk and get new fk
		stk->xx = strial.xx - lambdak*Cdfxx;
		stk->yy = strial.yy - lambdak*Cdfyy;
		stk->xy = strial.xy - lambdak*Cdfxy;
		if(np==PLANE_STRAIN_MPM)
			stk->zz = strial.zz - lambdak*Cdfzz;
		
		//cout << "# ......" << step << "," << lambdak << "," << GetF(mptr,stk->xx,stk->yy,stk->xy,stk->zz,np) << endl;
		
		// check for convergence
		if(LambdaConverged(step++,lambdak,delLam)) break;
		
		// new f
		fk=GetF(mptr,stk->xx,stk->yy,stk->xy,stk->zz,np);
	}
	
	// get final derivatives
	GetDfDsigma(mptr,stk,np);
	
	return lambdak;
}

// decide if the numerical solution for lambda had converged
// subclass can override to change convergence rules
bool AnisoPlasticity::LambdaConverged(int step,double lambda,double delLam)
{
	if(step>20 || fabs(delLam/lambda)<0.001) return true;
	return false;
}


