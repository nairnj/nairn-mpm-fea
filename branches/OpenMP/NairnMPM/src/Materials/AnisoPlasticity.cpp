/********************************************************************************
    AnisoPlasticity.cpp
    NairnMPM
    
    Created by John Nairn, June 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
 
	Yield criterion is
 
		f = sqrt(sigma.A sigma) - g(alpha)
 
	where sigma is stress after rotation into the material axes, A is
	matrix of yield properties by Hill method, and g(alpha) is a
	hardening law.

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
	
	return Orthotropic::InputMat(xName,input);
}

// verify settings and some initial calculations
const char *AnisoPlasticity::VerifyAndLoadProperties(int np)
{
	// check at least some yielding
	if(syzz<0. && syxx<0. && syyy<0. && tyxy<0. && tyxz<0. && tyyz<0.)
		return "No yield stresses were defined";

	// check A is positive semi definite
	double rsxx=0.,rsyy=0.,rszz=0.;
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
	
	// combination terms
	fTerm=(syyyred2 + syzzred2 - syxxred2)/2.;
	gTerm=(syzzred2 + syxxred2 - syyyred2)/2.;
	hTerm=(syxxred2 + syyyred2 - syzzred2)/2.;
	
	// call super class
	return Orthotropic::VerifyAndLoadProperties(np);
}

// if cannot be used in current analysis type throw CommonException()
void AnisoPlasticity::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
	{	throw CommonException("Anisotropic plasticity materials cannot use 2D plane stress MPM analysis",
							  "AnisoPlasticity::ValidateForUse");
	}
	
	// call super class (why can't call super class?)
	return Orthotropic::ValidateForUse(np);
}

// print to output window
void AnisoPlasticity::PrintMechanicalProperties(void) const
{	
    Orthotropic::PrintMechanicalProperties();
	PrintYieldProperties();
}

// print just yield properties to output window
void AnisoPlasticity::PrintYieldProperties(void) const
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

#pragma mark AnisoPlasticity::Methods

// Isotropic material can use read-only initial properties
void *AnisoPlasticity::GetCopyOfMechanicalProps(MPMBase *mptr,int np) const
{
	// full plastic properties
	AnisoPlasticProperties *p = (AnisoPlasticProperties *)malloc(sizeof(AnisoPlasticProperties));
	
	// create new elastic properties
	p->ep = (ElasticProperties *)malloc(sizeof(ElasticProperties));
	if(np!=THREED_MPM)
		FillElasticProperties2D(p->ep,TRUE,mptr->GetRotationZ(),np);
	else
		FillElasticProperties3D(mptr,p->ep,np);
	
	return p;
}

// If need, cast void * to correct pointer and delete it
void AnisoPlasticity::DeleteCopyOfMechanicalProps(void *properties,int np) const
{
	AnisoPlasticProperties *p = (AnisoPlasticProperties *)properties;
	delete p->ep;
	delete p;
}

/* For 2D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, plastic strain, stresses, strain energy, 
		plastic energy, dissipated energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
   For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
   This is general analysis for anisotropic material based on hill criterion. Various
		hardening laws can implement the hardening term.
*/
void AnisoPlasticity::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
        double dvzz,double delTime,int np,void *properties,ResidualStrains *res) const
{
	AnisoPlasticProperties *p = (AnisoPlasticProperties *)properties;
	
    // Effective strain by deducting thermal strain
	// (note ep->alpha[1] and ep->beta[1] are reduced in plane strain, but CTE3 and CME3 are not)
	ElasticProperties *r = p->ep;
	double erxx = r->alpha[1]*res->dT;
	double eryy = r->alpha[2]*res->dT;
	double erxy = r->alpha[3]*res->dT;
	double erzz = CTE3*res->dT;
	if(DiffusionTask::active)
	{	erxx += r->beta[1]*res->dC;
		eryy += r->beta[2]*res->dC;
		erxy += r->beta[3]*res->dC;
		erzz += CME3*res->dC;
	}
    double dexx = dvxx-erxx;  
    double deyy = dvyy-eryy;
    double dgxy = dvxy+dvyx-erxy;
	double dezz = dvzz-erzz;

    // Elastic stress increment
	Tensor *sp = mptr->GetStressTensor();
    Tensor stk = *sp;
    stk.xx += r->C[1][1]*dexx+r->C[1][2]*deyy+r->C[1][3]*dgxy;
    stk.yy += r->C[1][2]*dexx+r->C[2][2]*deyy+r->C[2][3]*dgxy;
    stk.xy += r->C[1][3]*dexx+r->C[2][3]*deyy+r->C[3][3]*dgxy;
	stk.zz += r->C[4][1]*dexx + r->C[4][2]*deyy + r->C[4][3]*dgxy + r->C[4][4]*dezz;
    if(np==PLANE_STRAIN_MPM)
	{	stk.zz += r->C[4][1]*r->alpha[5]*erzz + r->C[4][2]*r->alpha[6]*erzz
					+ r->C[4][3]*r->alpha[7]*erzz;
	}
  
	// get rotation matrix elements and evaluate sin() and cos() once per update
	double angle=mptr->GetRotationZ();
	double c=cos(angle);
	double s=sin(angle);
	p->rzyx[0][0] = p->rzyx[1][1] = c*c;
	p->rzyx[0][1] = p->rzyx[1][0] = s*s;
	p->rzyx[0][5] = -2.*c*s;
	p->rzyx[1][5] = -p->rzyx[0][5];
	p->rzyx[5][0] = c*s;
	p->rzyx[5][1] = -p->rzyx[5][0];
	p->rzyx[5][5] = p->rzyx[0][0] - p->rzyx[0][1];

    // Calculate plastic potential f
	UpdateTrialAlpha(mptr,np,p);
	Tensor srot;
	double sAstrial = GetMagnitudeRotatedHill(&stk,&srot,np,p);
	double ftrial = sAstrial - GetYield(p);
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
		Elastic::MPMConstLaw(mptr,dvxx,dvyy,dvxy,dvyx,dvzz,delTime,np,r,res);
		return; 
    }
    
	// Find  lambda for this plastic state
	// This material overrides to have custom solution
	double lambda = SolveForLambdaAP(mptr,np,ftrial,&stk,p);
    
    // Plastic strain increments on particle
    double dexxp = lambda*p->dfds.xx;
    double deyyp = lambda*p->dfds.yy;
    double dgxyp = lambda*p->dfds.xy;
	double dezzp = lambda*p->dfds.zz;
	Tensor *eplast = mptr->GetPlasticStrainTensor();
    eplast->xx += dexxp;
    eplast->yy += deyyp;
    eplast->xy += dgxyp;
	eplast->zz += dezzp;
    
    // Elastic strain increments on particle
	Tensor *ep=mptr->GetStrainTensor();
    ep->xx += (dvxx-dexxp);
    ep->yy += (dvyy-deyyp);
    ep->xy += (dvxy+dvyx-dgxyp);
    ep->zz += (dvzz-dezzp);
	
	// rotational strain
	double dwrotxy=dvyx-dvxy;
	
	// Effective elastic strain increment = de - deres
	dexx -= dexxp;
	deyy -= deyyp;
	dgxy -= dgxyp;
	dezz -= dezzp;
	
	// increment particle stresses
    // Hypoelastic - Eliminate effect of rotation 
    Tensor st0=*sp;			// save previous stress for energy updates below
    double c1 = r->C[1][1]*dexx+r->C[1][2]*deyy+r->C[1][3]*dgxy;
    double c2 = r->C[1][2]*dexx+r->C[2][2]*deyy+r->C[2][3]*dgxy;
    double c3 = r->C[1][3]*dexx+r->C[2][3]*deyy+r->C[3][3]*dgxy;
	Hypo2DCalculations(mptr,-dwrotxy,c1,c2,c3);
	
	// out of plane stress
    if(np==PLANE_STRAIN_MPM)
	{	sp->zz += r->C[4][1]*(dexx+r->alpha[5]*erzz) + r->C[4][2]*(deyy+r->alpha[6]*erzz)
					+ r->C[4][3]*(dgxy+r->alpha[7]*erzz) + r->C[4][4]*dezz;
	}

	// Elastic energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
	double strainEnergy = 0.5*((st0.xx+sp->xx)*dexx
							   + (st0.yy+sp->yy)*deyy
							   + (st0.xy+sp->xy)*dgxy
							   + (st0.zz+sp->zz)*dezz);

    // Plastic energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
	double dispEnergy = 0.5*((st0.xx+sp->xx)*dexxp
                        + (st0.yy+sp->yy)*deyyp
                        + (st0.xy+sp->xy)*dgxyp
						+ (st0.zz+sp->zz)*dezzp);

	// add now
	mptr->AddStrainEnergy(strainEnergy + dispEnergy);
	
	// add dissipated energy to plastic energy to the particle
    mptr->AddPlastEnergy(dispEnergy);
	
    // heat energy is Cv(dT-dTq0) - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adiabatic heating occurs)
    IncrementHeatEnergy(mptr,res->dT,0.,dispEnergy);

	// update internal variables
	UpdatePlasticInternal(mptr,np,p);
}

/* For 3D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
	Assumes linear elastic, uses hypoelastic correction
*/
void AnisoPlasticity::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
								double dvxz,double dvzx,double dvyz,double dvzy,
								  double delTime,int np,void *properties,ResidualStrains *res) const
{
	AnisoPlasticProperties *p = (AnisoPlasticProperties *)properties;
	
    // Effective strain by deducting thermal strain
	ElasticProperties *r = p->ep;
	double erxx=r->alpha[0]*res->dT;
	double eryy=r->alpha[1]*res->dT;
	double erzz=r->alpha[2]*res->dT;
	double eryz=r->alpha[3]*res->dT;
	double erxz=r->alpha[4]*res->dT;
	double erxy=r->alpha[5]*res->dT;
	if(DiffusionTask::active)
	{	erxx+=r->beta[0]*res->dC;
		eryy+=r->beta[1]*res->dC;
		erzz+=r->beta[2]*res->dC;
		eryz+=r->beta[3]*res->dC;
		erxz+=r->beta[4]*res->dC;
		erxy+=r->beta[5]*res->dC;
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
    stk.xx += r->C[0][0]*dexx+r->C[0][1]*deyy+r->C[0][2]*dezz+r->C[0][3]*dgyz+r->C[0][4]*dgxz+r->C[0][5]*dgxy;
    stk.yy += r->C[1][0]*dexx+r->C[1][1]*deyy+r->C[1][2]*dezz+r->C[1][3]*dgyz+r->C[1][4]*dgxz+r->C[1][5]*dgxy;
    stk.zz += r->C[2][0]*dexx+r->C[2][1]*deyy+r->C[2][2]*dezz+r->C[2][3]*dgyz+r->C[2][4]*dgxz+r->C[2][5]*dgxy;
    stk.yz += r->C[3][0]*dexx+r->C[3][1]*deyy+r->C[3][2]*dezz+r->C[3][3]*dgyz+r->C[3][4]*dgxz+r->C[3][5]*dgxy;
    stk.xz += r->C[4][0]*dexx+r->C[4][1]*deyy+r->C[4][2]*dezz+r->C[4][3]*dgyz+r->C[4][4]*dgxz+r->C[4][5]*dgxy;
    stk.xy += r->C[5][0]*dexx+r->C[5][1]*deyy+r->C[5][2]*dezz+r->C[5][3]*dgyz+r->C[5][4]*dgxz+r->C[5][5]*dgxy;
	
	// get rotation matrix elements p->rzyx[i][j]
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
	p->rzyx[0][0]=cy2*cz2;
	p->rzyx[0][1]=cy2*sz2;
	p->rzyx[0][2]=sy2;
	p->rzyx[0][3]=-2*cy*sy*sz;
	p->rzyx[0][4]=cz*s2y;
	p->rzyx[0][5]=-2*cy2*cz*sz;
	
	double arg1=cz*sx*sy + cx*sz;
	p->rzyx[1][0]=arg1*arg1;
	double arg2=cx*cz - sx*sy*sz;
	p->rzyx[1][1]=arg2*arg2;
	p->rzyx[1][2]=cy2*sx2;
	p->rzyx[1][3]=-2*cy*sx*arg2;
	p->rzyx[1][4]=-2*cy*sx*arg1;
	p->rzyx[1][5]=cx2*s2z + c2z*s2x*sy - 2*cz*sx2*sy2*sz;

	double arg3=cx*cz*sy - sx*sz;
	p->rzyx[2][0]=arg3*arg3;
	double arg4=cz*sx + cx*sy*sz;
	p->rzyx[2][1]=arg4*arg4;
	p->rzyx[2][2]=cx2*cy2;
	p->rzyx[2][3]=2*cx*cy*arg4;
	p->rzyx[2][4]=-2*cx*cy*arg3;
	p->rzyx[2][5]=s2z*sx2 - c2z*s2x*sy - 2*cx2*cz*sy2*sz;
	
	double arg5=c2x*s2z*sy/2.;
	p->rzyx[3][0]=-arg5 + cx*sx*(sz2 - cz2*sy2);
	p->rzyx[3][1]=arg5 + cx*sx*(cz2 - sy2*sz2);
	p->rzyx[3][2]=-cx*cy2*sx;
	p->rzyx[3][3]=cy*(c2x*cz - 2*cx*sx*sy*sz);
	p->rzyx[3][4]=cy*(cz*s2x*sy + c2x*sz);
	p->rzyx[3][5]=-c2x*c2z*sy + 2*cx*cz*sx*(1 + sy2)*sz;
	
	p->rzyx[4][0]=-cy*cz*arg3;
	p->rzyx[4][1]=-cy*sz*arg4;
	p->rzyx[4][2]=cx*cy*sy;
	p->rzyx[4][3]=cz*sx*sy - c2y*cx*sz;
	p->rzyx[4][4]=c2y*cx*cz + sx*sy*sz;
	p->rzyx[4][5]=cy*(c2z*sx + cx*s2z*sy);
	
	p->rzyx[5][0]=cy*cz*arg1;
	p->rzyx[5][1]=-cy*sz*arg2;
	p->rzyx[5][2]=-cy*sx*sy;
	p->rzyx[5][3]=cx*cz*sy + c2y*sx*sz;
	p->rzyx[5][4]=cx*sy*sz - c2y*cz*sx;
	p->rzyx[5][5]=cy*(c2z*cx - 2*cz*sx*sy*sz);
	
    // Calculate plastic potential f
	UpdateTrialAlpha(mptr,np,p);
	Tensor srot;
	double sAstrial = GetMagnitudeRotatedHill(&stk,&srot,np,p);
	double ftrial = sAstrial - GetYield(p);
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
		Elastic::MPMConstLaw(mptr,dvxx,dvyy,dvzz,dvxy,dvyx,dvxz,dvzx,dvyz,dvzy,delTime,np,r,res);
		return; 
    }
    
	// Find  lambda for this plastic state
	// This material overrides to custom solution
	double lambda = SolveForLambdaAP(mptr,np,ftrial,&stk,p);
    
    // Plastic strain increments on particle
    double dexxp=lambda*p->dfds.xx;
    double deyyp=lambda*p->dfds.yy;
	double dezzp=lambda*p->dfds.zz;
    double dgyzp=lambda*p->dfds.yz;
    double dgxzp=lambda*p->dfds.xz;
    double dgxyp=lambda*p->dfds.xy;
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
    dsig[XX] = r->C[0][0]*dexx+r->C[0][1]*deyy+r->C[0][2]*dezz+r->C[0][3]*dgyz+r->C[0][4]*dgxz+r->C[0][5]*dgxy;
    dsig[YY] = r->C[1][0]*dexx+r->C[1][1]*deyy+r->C[1][2]*dezz+r->C[1][3]*dgyz+r->C[1][4]*dgxz+r->C[1][5]*dgxy;
    dsig[ZZ] = r->C[2][0]*dexx+r->C[2][1]*deyy+r->C[2][2]*dezz+r->C[2][3]*dgyz+r->C[2][4]*dgxz+r->C[2][5]*dgxy;
    dsig[YZ] = r->C[3][0]*dexx+r->C[3][1]*deyy+r->C[3][2]*dezz+r->C[3][3]*dgyz+r->C[3][4]*dgxz+r->C[3][5]*dgxy;
    dsig[XZ] = r->C[4][0]*dexx+r->C[4][1]*deyy+r->C[4][2]*dezz+r->C[4][3]*dgyz+r->C[4][4]*dgxz+r->C[4][5]*dgxy;
    dsig[XY] = r->C[5][0]*dexx+r->C[5][1]*deyy+r->C[5][2]*dezz+r->C[5][3]*dgyz+r->C[5][4]*dgxz+r->C[5][5]*dgxy;
	Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,dsig);

    // Elastic energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
	double strainEnergy = 0.5*((st0.xx+sp->xx)*dexx
							   + (st0.yy+sp->yy)*deyy
							   + (st0.zz+sp->zz)*dezz
							   + (st0.yz+sp->yz)*dgyz
							   + (st0.xz+sp->xz)*dgxz
							   + (st0.xy+sp->xy)*dgxy);
	
    // Plastic energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
	double dispEnergy=0.5*(0.5*((st0.xx+sp->xx)*dexxp
								+ (st0.yy+sp->yy)*deyyp
								+ (st0.zz+sp->zz)*dezzp
								+ (st0.yz+sp->yz)*dgyzp
								+ (st0.xz+sp->xz)*dgxzp
								+ (st0.xy+sp->xy)*dgxyp));
	
	// add now
	mptr->AddStrainEnergy(strainEnergy + dispEnergy);
	
	// add dissipated energy to plastic energy to the particle
    mptr->AddPlastEnergy(dispEnergy);
	
    // heat energy is Cv(dT-dTq0) - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adiabatic heating occurs)
    IncrementHeatEnergy(mptr,res->dT,0.,dispEnergy);
	
	// update internal variables
	UpdatePlasticInternal(mptr,np,p);
}

#pragma mark AnisoPlasticity::Hill Terms

// Get sqrt(s As) and also return rotated stresses in case caller needs them
double AnisoPlasticity::GetMagnitudeRotatedHill(Tensor *st0,Tensor *srot,int np,AnisoPlasticProperties *p) const
{
	double sAs=0.;;
	
	if(np==THREED_MPM)
	{	// rotation from analysis to material axes
		srot->xx = p->rzyx[0][0]*st0->xx+p->rzyx[0][1]*st0->yy+p->rzyx[0][2]*st0->zz
					+p->rzyx[0][3]*st0->yz+p->rzyx[0][4]*st0->xz+p->rzyx[0][5]*st0->xy;
		srot->yy = p->rzyx[1][0]*st0->xx+p->rzyx[1][1]*st0->yy+p->rzyx[1][2]*st0->zz
					+p->rzyx[1][3]*st0->yz+p->rzyx[1][4]*st0->xz+p->rzyx[1][5]*st0->xy;
		srot->zz = p->rzyx[2][0]*st0->xx+p->rzyx[2][1]*st0->yy+p->rzyx[2][2]*st0->zz
					+p->rzyx[2][3]*st0->yz+p->rzyx[2][4]*st0->xz+p->rzyx[2][5]*st0->xy;
		srot->yz = p->rzyx[3][0]*st0->xx+p->rzyx[3][1]*st0->yy+p->rzyx[3][2]*st0->zz
					+p->rzyx[3][3]*st0->yz+p->rzyx[3][4]*st0->xz+p->rzyx[3][5]*st0->xy;
		srot->xz = p->rzyx[4][0]*st0->xx+p->rzyx[4][1]*st0->yy+p->rzyx[4][2]*st0->zz
					+p->rzyx[4][3]*st0->yz+p->rzyx[4][4]*st0->xz+p->rzyx[4][5]*st0->xy;
		srot->xy = p->rzyx[5][0]*st0->xx+p->rzyx[5][1]*st0->yy+p->rzyx[5][2]*st0->zz
					+p->rzyx[5][3]*st0->yz+p->rzyx[5][4]*st0->xz+p->rzyx[5][5]*st0->xy;
		
		sAs = srot->xz*srot->xz*tyxzred2 + srot->yz*srot->yz*tyyzred2;
	}
	else
	{	// clockwise rotation from analysis to material axes
		srot->xx = st0->xx*p->rzyx[0][0] + st0->yy*p->rzyx[0][1] + st0->xy*p->rzyx[0][5];
		srot->yy = st0->xx*p->rzyx[1][0] + st0->yy*p->rzyx[1][1] + st0->xy*p->rzyx[1][5];
		srot->xy = (st0->xx-st0->yy)*p->rzyx[5][0] + st0->xy*p->rzyx[5][5];
		srot->zz = st0->zz;
		srot->xz = 0.;
		srot->yz = 0.;
	}
	
	// return sqrt(sAs);
	// check on negative sAs can happen due to round-off error when stresses near zero
	double dyz=srot->yy-srot->zz;
	double dxz=srot->xx-srot->zz;
	double dxy=srot->xx-srot->yy;
	sAs += fTerm*dyz*dyz + gTerm*dxz*dxz + hTerm*dxy*dxy + srot->xy*srot->xy*tyxyred2;
	return sAs>0. ? sqrt(sAs) : 0 ;
}

// Find C.df and df.C.df at given stress - store results in plastic property variables active only during the loop
// and only for current material point
void AnisoPlasticity::GetDfCdf(Tensor *stk,int np,AnisoPlasticProperties *p) const
{
	// get C df and df C df, which need df/dsig (Function of yield criteria, and normally only current stress in stk)
	GetDfDsigma(stk,np,p);
	ElasticProperties *r = p->ep;
	if(np==THREED_MPM)
	{	p->Cdf.xx = r->C[0][0]*p->dfds.xx+r->C[0][1]*p->dfds.yy+r->C[0][2]*p->dfds.zz
						+r->C[0][3]*p->dfds.yz+r->C[0][4]*p->dfds.xz+r->C[0][5]*p->dfds.xy;
		p->Cdf.yy = r->C[1][0]*p->dfds.xx+r->C[1][1]*p->dfds.yy+r->C[1][2]*p->dfds.zz
						+r->C[1][3]*p->dfds.yz+r->C[1][4]*p->dfds.xz+r->C[1][5]*p->dfds.xy;
		p->Cdf.zz = r->C[2][0]*p->dfds.xx+r->C[2][1]*p->dfds.yy+r->C[2][2]*p->dfds.zz
						+r->C[2][3]*p->dfds.yz+r->C[2][4]*p->dfds.xz+r->C[2][5]*p->dfds.xy;
		p->Cdf.yz = r->C[3][0]*p->dfds.xx+r->C[3][1]*p->dfds.yy+r->C[3][2]*p->dfds.zz
						+r->C[3][3]*p->dfds.yz+r->C[3][4]*p->dfds.xz+r->C[3][5]*p->dfds.xy;
		p->Cdf.xz = r->C[4][0]*p->dfds.xx+r->C[4][1]*p->dfds.yy+r->C[4][2]*p->dfds.zz
						+r->C[4][3]*p->dfds.yz+r->C[4][4]*p->dfds.xz+r->C[4][5]*p->dfds.xy;
		p->Cdf.xy = r->C[5][0]*p->dfds.xx+r->C[5][1]*p->dfds.yy+r->C[5][2]*p->dfds.zz
						+r->C[5][3]*p->dfds.yz+r->C[5][4]*p->dfds.xz+r->C[5][5]*p->dfds.xy;
		p->dfCdf = p->dfds.xx*p->Cdf.xx + p->dfds.yy*p->Cdf.yy + p->dfds.zz*p->Cdf.zz
						+ p->dfds.yz*p->Cdf.yz + p->dfds.xz*p->Cdf.xz + p->dfds.xy*p->Cdf.xy;
	}
	else
	{	p->Cdf.xx = r->C[1][1]*p->dfds.xx + r->C[1][2]*p->dfds.yy + r->C[1][3]*p->dfds.xy;
		p->Cdf.yy = r->C[1][2]*p->dfds.xx + r->C[2][2]*p->dfds.yy + r->C[2][3]*p->dfds.xy;
		p->Cdf.xy = r->C[1][3]*p->dfds.xx + r->C[2][3]*p->dfds.yy + r->C[3][3]*p->dfds.xy;
		p->Cdf.zz = r->C[4][1]*p->dfds.xx + r->C[4][2]*p->dfds.yy + r->C[4][3]*p->dfds.xy + r->C[4][4]*p->dfds.zz;
		p->dfCdf = p->dfds.xx*p->Cdf.xx + p->dfds.yy*p->Cdf.yy + p->dfds.xy*p->Cdf.xy + p->dfds.zz*p->Cdf.zz;;
	}
}

// Find A srot/sqrt(s As) in rotated coordinates (in dfdsijrot) and then rotate the result to
//   get df term such that dep = lamda df in the analysis coordinates. Rotated term
//	 is stored in dfds tensor in plastic properties
// As find R df = -h and store in minush hardning property
void AnisoPlasticity::GetDfDsigma(Tensor *st0,int np,AnisoPlasticProperties *p) const
{
	// clockwise rotation from analysis to material axes
	Tensor srot;
	double rootSAS = GetMagnitudeRotatedHill(st0,&srot,np,p);
	
	if(rootSAS>0.)
	{	double dfdsxxrot,dfdsyyrot,dfdszzrot,dfdtxyrot;
		dfdsxxrot = (syxxred2*srot.xx - hTerm*srot.yy - gTerm*srot.zz) / rootSAS;
		dfdsyyrot = (-hTerm*srot.xx + syyyred2*srot.yy - fTerm*srot.zz) / rootSAS;
		dfdszzrot = (-gTerm*srot.xx - fTerm*srot.yy + syzzred2*srot.zz) / rootSAS;
		dfdtxyrot = tyxyred2*srot.xy / rootSAS;
		
		if(np==THREED_MPM)
		{	double dfdtyzrot,dfdtxzrot;
			dfdtyzrot = tyyzred2*srot.yz / rootSAS;
			dfdtxzrot = tyxzred2*srot.xz / rootSAS;
			
			// rotate to analysis coordinates df = R^T dfrot
			p->dfds.xx = p->rzyx[0][0]*dfdsxxrot+p->rzyx[1][0]*dfdsyyrot+p->rzyx[2][0]*dfdszzrot
							+p->rzyx[3][0]*dfdtyzrot+p->rzyx[4][0]*dfdtxzrot+p->rzyx[5][0]*dfdtxyrot;
			p->dfds.yy = p->rzyx[0][1]*dfdsxxrot+p->rzyx[1][1]*dfdsyyrot+p->rzyx[2][1]*dfdszzrot
							+p->rzyx[3][1]*dfdtyzrot+p->rzyx[4][1]*dfdtxzrot+p->rzyx[5][1]*dfdtxyrot;
			p->dfds.zz = p->rzyx[0][2]*dfdsxxrot+p->rzyx[1][2]*dfdsyyrot+p->rzyx[2][2]*dfdszzrot
							+p->rzyx[3][2]*dfdtyzrot+p->rzyx[4][2]*dfdtxzrot+p->rzyx[5][2]*dfdtxyrot;
			p->dfds.yz = p->rzyx[0][3]*dfdsxxrot+p->rzyx[1][3]*dfdsyyrot+p->rzyx[2][3]*dfdszzrot
							+p->rzyx[3][3]*dfdtyzrot+p->rzyx[4][3]*dfdtxzrot+p->rzyx[5][3]*dfdtxyrot;
			p->dfds.xz = p->rzyx[0][4]*dfdsxxrot+p->rzyx[1][4]*dfdsyyrot+p->rzyx[2][4]*dfdszzrot
							+p->rzyx[3][4]*dfdtyzrot+p->rzyx[4][4]*dfdtxzrot+p->rzyx[5][4]*dfdtxyrot;
			p->dfds.xy = p->rzyx[0][5]*dfdsxxrot+p->rzyx[1][5]*dfdsyyrot+p->rzyx[2][5]*dfdszzrot
							+p->rzyx[3][5]*dfdtyzrot+p->rzyx[4][5]*dfdtxzrot+p->rzyx[5][5]*dfdtxyrot;
			
			// for use in alpha upate
			p->minush = dfdsxxrot*dfdsxxrot + dfdsyyrot*dfdsyyrot + dfdszzrot*dfdszzrot
							+ 0.5*(dfdtyzrot*dfdtyzrot + dfdtxzrot*dfdtxzrot + dfdtxyrot*dfdtxyrot);
			p->minush = sqrt(p->minush/1.5);
		}
		else
		{	// rotate to analysis coordinates df = R^(-1) dfrot = R^T dfrot
			p->dfds.xx = dfdsxxrot*p->rzyx[0][0] + dfdsyyrot*p->rzyx[1][0] + dfdtxyrot*p->rzyx[5][0];
			p->dfds.yy = dfdsxxrot*p->rzyx[0][1] + dfdsyyrot*p->rzyx[1][1] + dfdtxyrot*p->rzyx[5][1];
			p->dfds.xy = (dfdsxxrot - dfdsyyrot)*p->rzyx[0][5] + dfdtxyrot*p->rzyx[5][5];
			p->dfds.zz = dfdszzrot;
			
			// for use in alpha upate
			p->minush = dfdsxxrot*dfdsxxrot + dfdsyyrot*dfdsyyrot + dfdszzrot*dfdszzrot + 0.5*dfdtxyrot*dfdtxyrot;
			p->minush = sqrt(p->minush/1.5);
		}
	}
	
	else
	{	// negative root implies zero stress with roundoff error
		p->dfds.xx = p->dfds.yy = p->dfds.zz = p->dfds.yz = p->dfds.xz = p->dfds.xy = 0.;
		p->minush=0.;
	}
}

#pragma mark AnisoPlasticity::Custom Methods

// Solve numerically for lambda
// Ouptut is lambdak, final df, final alpha
double AnisoPlasticity::SolveForLambdaAP(MPMBase *mptr,int np,double ftrial,Tensor *strial,AnisoPlasticProperties *p) const
{
	// step 0 = stress is strial, alpha is previous alpha
	GetDfCdf(strial,np,p);								// find df/dsigma
	double lambda2 = ftrial/(p->dfCdf + GetDfAlphaDotH(mptr,np,strial,p));
	double lambdaInitial=lambda2;
	
	// find df/dsigma and -h at this initial guessed stress rather than at strial as above
	Tensor stk;
	UpdateStress(strial,&stk,lambda2,np,p);
	GetDfCdf(&stk,np,p);
	UpdateStress(strial,&stk,lambda2,np,p);			// second pass might help
	GetDfCdf(&stk,np,p);

	// Find f using new slopes at lambda1
	UpdateTrialAlpha(mptr,np,lambda2,p);
	Tensor srot;
	double f2 = GetMagnitudeRotatedHill(&stk,&srot,np,p) - GetYield(p);
	
	// pick second lambda that is lower
	double lambda1 = 0.5*lambda2;
	double f1 = GetFkFromLambdak(mptr,strial,&stk,lambda1,np,p);
	
	// bracket the solution
	int step=1;
	while(true)
	{	if(f1*f2<0.) break;
		if(fabs(f1)<fabs(f2))
		{	lambda1+=1.6*(lambda1-lambda2);
			if(lambda1>0.)
				f1=GetFkFromLambdak(mptr,strial,&stk,lambda1,np,p);
			else
			{	f1=ftrial;
				lambda1=0.;
			}
		}
		else
		{	lambda2+=1.6*(lambda2-lambda1);
			f2=GetFkFromLambdak(mptr,strial,&stk,lambda2,np,p);
		}
		
		// if fails to bracket in 50 tries, return with single step solution
		step++;
		if(step>50)
		{
			GetDfCdf(strial,np,p);
			UpdateStress(strial,&stk,lambdaInitial,np,p);
			GetDfCdf(&stk,np,p);
			UpdateTrialAlpha(mptr,np,lambdaInitial,p);
			return lambdaInitial;
		}
	}
	
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
	double fk=GetFkFromLambdak(mptr,strial,&stk,lambdak,np,p);
	double dfkdlam=-(p->dfCdf + GetDfAlphaDotH(mptr,np,&stk,p));
	
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
		fk=GetFkFromLambdak(mptr,strial,&stk,lambdak,np,p);
		dfkdlam=-(p->dfCdf + GetDfAlphaDotH(mptr,np,&stk,p));
		
		// maintain the bracket
		if(fk<0.0)
		{	lambda1=lambdak;
			f1=fk;
		}
		else
		{	lambda2=lambdak;
			f2=fk;
		}
		
	}
	
	// output df (from initial setting of GetDfCdf()), alpha (here with latest lambda), and lamda (the return value)
	UpdateTrialAlpha(mptr,np,lambdak,p);
	
	return lambdak;
}

// Update stress (assumes Cdfij calculated before)
void AnisoPlasticity::UpdateStress(Tensor *strial,Tensor *stk,double lambda,int np,AnisoPlasticProperties *p) const
{
	stk->xx = strial->xx - lambda*p->Cdf.xx;
	stk->yy = strial->yy - lambda*p->Cdf.yy;
	stk->xy = strial->xy - lambda*p->Cdf.xy;
	stk->zz = strial->zz - lambda*p->Cdf.zz;
	if(np==THREED_MPM)
	{	stk->yz = strial->yz - lambda*p->Cdf.yz;
		stk->xz = strial->xz - lambda*p->Cdf.xz;
	}
}

// Given stress and lambda, find f (also find the stress and return dfCdf)
double AnisoPlasticity::GetFkFromLambdak(MPMBase *mptr,Tensor *strial,Tensor *stk,
										 double lambda,int np,AnisoPlasticProperties *p) const
{
	// Change stress to new value based on strial, lambda, and previous slope
	UpdateStress(strial,stk,lambda,np,p);
	
	// update alpha using new lambda and -h from most recent GetDfDsigma()
	UpdateTrialAlpha(mptr,np,lambda,p);
	
	// update fk using new stress and alpha
	Tensor srot;
	return GetMagnitudeRotatedHill(stk,&srot,np,p) - GetYield(p);
}

// plastic strain needed to get deformation gradient for this material class
bool AnisoPlasticity::PartitionsElasticAndPlasticStrain(void) { return TRUE; }

