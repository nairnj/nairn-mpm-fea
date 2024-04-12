/********************************************************************************
    AnisoPlasticity.cpp
    nairn-mpm-fea
    
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

#include "stdafx.h"
#include "Materials/AnisoPlasticity.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"		// only need to trial solution method

int hstepsMax=0;

// class statics
int AnisoPlasticity::warnNonconvergence;

double maxLambda = 0.;

#pragma mark AnisoPlasticity::Constructors and Destructors

// Constructors
AnisoPlasticity::AnisoPlasticity(char *matName,int matID) : Orthotropic(matName,matID)
{
	// negative yield stress implies no yielding in that direction
	syxx=-1.;
	syyy=-1.;
	syzz=-1.;
	tyxy=-1.;
	tyxz=-1.;
	tyyz=-1.;
	
	h.style = SQUARED_TERMS;
}

#pragma mark AnisoPlasticity::Initialization

// Read material properties
char *AnisoPlasticity::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"yldxx")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&syxx,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"yldyy")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&syyy,gScaling,1.e6);
    }
	
    else if(strcmp(xName,"yldzz")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&syzz,gScaling,1.e6);
    }

    else if(strcmp(xName,"yldxy")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&tyxy,gScaling,1.e6);
    }
	
    else if(strcmp(xName,"yldxz")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&tyxz,gScaling,1.e6);
    }
	
    else if(strcmp(xName,"yldyz")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&tyyz,gScaling,1.e6);
    }
	
	else if(strcmp(xName,"HillStyle")==0)
	{	input=INT_NUM;
		return (char *)(&h.style);
	}
	
	return Orthotropic::InputMaterialProperty(xName,input,gScaling);
}

// verify settings and some initial calculations
const char *AnisoPlasticity::VerifyAndLoadProperties(int np)
{
#ifdef POROELASTICITY
	if(DiffusionTask::HasPoroelasticity())
		return "Anisotropic plasticity materials cannot yet be used when poroelasticity is activated";
#endif
    
    if(swapz==1)
    {   // swap x and z, xx->zz, xy->yz, xz same, yz->xy
        SwapProperties(syxx,syzz,tyyz,tyxy);
    }
    else if(swapz>1)
    {   // swap y and z, yy->zz, xy->xz, xz->xy, yz same
        SwapProperties(syyy,syzz,tyxz,tyxy);
    }

	// requires large rotation mode
	if(useLargeRotation==0) useLargeRotation = 1;
	
	// check at least some yielding
	if(syzz<0. && syxx<0. && syyy<0. && tyxy<0. && tyxz<0. && tyyz<0.)
		return "No yield stresses were defined";
	
	// check non zeros
	if(syzz==0. || syxx==0. || syyy==0. || tyxy==0. || tyxz==0. || tyyz==0.)
		return "No yield stresses can be zero";

	// check A is positive semi definite (watching for round off error is all the same)
	if(syzz>0. && syxx>0. && syyy>0. && syxx==syyy && syxx==syzz)
	{	// OK if all the same
	}
	else
	{	// <0 means infinite so reciprocal is zero
		double rsxx=0.,rsyy=0.,rszz=0.;
		if(syxx>0.)
			rsxx=1./(syxx*syxx);
		if(syyy>0.)
			rsyy=1./(syyy*syyy);
		if(syzz>0.)
			rszz=1./(syzz*syzz);
		// check all the same
		double arg = rsxx*rsxx + rsyy*rsyy + rszz*rszz - rsyy*rsxx - rszz*rsxx - rsyy*rszz ;
		double fgh = 0.5*(rsxx+rsyy+rszz);
		if(arg<0.) return "Hill plastic potential is not postive semidefinite (1)";
		if(fgh-sqrt(arg)<0.) return "Hill plastic potential is not postive semidefinite (2)";
	}
	
	// reciprocals of reduced normal yield stresses
	if(syxx>0.)
    {	h.syxx2=rho/syxx;
		h.syxx2*=h.syxx2;
	}
	else
		h.syxx2=0.;		// 1/inf^2
	if(syyy>0.)
    {	h.syyy2=rho/syyy;
		h.syyy2*=h.syyy2;
	}
	else
		h.syyy2=0.;		// 1/inf^2
	if(syzz>0.)
    {	h.syzz2=rho/syzz;
		h.syzz2*=h.syzz2;
	}
	else
		h.syzz2=0.;		// 1/inf^2
	
	// reciprocals of reduced shear yield stresses (actually 2*standard Hill L, N, and M)
	if(tyxy>0.)
    {	h.N=rho/tyxy;
		h.N*=h.N;
	}
	else
		h.N=0.;		// 1/inf^2
	if(tyxz>0.)
    {	h.M=rho/tyxz;
		h.M*=h.M;
	}
	else
		h.M=0.;		// 1/inf^2
	if(tyyz>0.)
    {	h.L=rho/tyyz;
        h.L*=h.L;
	}
	else
		h.L=0.;		// 1/inf^2
	
	// combination terms
    h.F = 0.5*(h.syyy2 + h.syzz2 - h.syxx2);
    h.G = 0.5*(h.syzz2 + h.syxx2 - h.syyy2);
    h.H = 0.5*(h.syxx2 + h.syyy2 - h.syzz2);
    //cout << "F=" << h.F*1.e12/(rho*rho) << ", G=" << h.G*1.e12/(rho*rho) << ", H=" << h.H*1.e12/(rho*rho) << endl;
    //cout << "L=" << h.L*1.e12/(2.*rho*rho) << ", M=" << h.M*1.e12/(2.*rho*rho) << ", N=" << h.N*1.e12/(2.*rho*rho) << endl;
	
	// reference yield stress (it is actually sqrt(2/3)/sigma(Y,ref))
    double sumNormal = h.F+h.G+h.H;
	if(sumNormal<=0.) return "(F+G+H) for Hill plastic potential must be positive";
	
	// The materials.tex notes define this term as 1/sigmaYref = sqrt23OversigmaYref.
	sqrt23OversigmaYref = sqrt(2.*sumNormal/3.);
	
	// for convergence problems
	warnNonconvergence=warnings.CreateWarning("Anisotropic plastic algorithm failed to converge",-1,3);
    
	// call super class
	const char *errMsg = Orthotropic::VerifyAndLoadProperties(np);
	if(errMsg!=NULL) return errMsg;
	
	// set more hill properties
	FillHillStyleProperties(np,h,pr);

	return NULL;
}

// if cannot be used in current analysis type
// throws std::bad_alloc
void AnisoPlasticity::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
	{	throw CommonException("Anisotropic plasticity materials cannot use 2D plane stress MPM analysis",
							  "AnisoPlasticity::ValidateForUse");
	}
	
	// call super class (why can't call super class?)
	Orthotropic::ValidateForUse(np);
}

// print to output window
void AnisoPlasticity::PrintMechanicalProperties(void) const
{	
    Orthotropic::PrintMechanicalProperties();
	PrintYieldProperties();
}

#pragma mark AnisoPlasticity::Methods

// buffer size for mechanical properties
int AnisoPlasticity::SizeOfMechanicalProperties(int &altBufferSize) const
{   altBufferSize = 0;
    return sizeof(AnisoPlasticProperties) ;
}

// Get current anisotropic properties (NULL on memory error)
void *AnisoPlasticity::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{
	AnisoPlasticProperties *p = (AnisoPlasticProperties *)matBuffer;
	p->ep = (ElasticProperties *)&pr;
	return p;
}

// return pointer to elastic properties
ElasticProperties *AnisoPlasticity::GetElasticPropertiesPointer(void *properties) const
{
	AnisoPlasticProperties *p = (AnisoPlasticProperties *)properties;
	return p->ep;
}

#pragma mark AnisoPlasticity::Methods (Large Rotation)

// Once stress deformation has been decomposed, finish calculations in the material axis system
// When stress are found, they are rotated back to the global axes (using Rtot and dR)
// Similarly, strain increments are rotated back to find work energy (done in global system)
void AnisoPlasticity::LRElasticConstitutiveLaw(MPMBase *mptr,Matrix3 &de,Matrix3 &er,Matrix3 &Rtot,Matrix3 &dR,
											   Matrix3 &Rnm1tot,int np,void *properties,ResidualStrains *res) const
{
	// Properties
	AnisoPlasticProperties *p = (AnisoPlasticProperties *)properties;
	ElasticProperties *r = p->ep;
	
	// get stress increment in material axes by first rotating current stress
	// to material and then adding elastic increment
	Tensor *sp = mptr->GetStressTensor();
	
	// effective strains
	double dvxxeff = de(0,0)-er(0,0);
	double dvyyeff = de(1,1)-er(1,1);
	double dvzzeff = de(2,2)-er(2,2);
	double dgamxy = de(0,1)+de(1,0);

	// Step 1: Get trial stress in material axes assuming increment is elastic
	Tensor strial;
	Tensor stnm1,stnm1Mat;
	bool is2D;
	if(np==THREED_MPM)
	{	double dgamyz = de(1,2)+de(2,1);
		double dgamxz = de(0,2)+de(2,0);
		
		// rotate current stress (stnm1) to material axes (stnm1Mat) and increment with effective strain increment
		stnm1 = MakeTensor(sp->xx,sp->yy,sp->zz,sp->yz,sp->xz,sp->xy);
		stnm1Mat = Rnm1tot.RTVoightR(&stnm1,true,false);
		strial.xx = stnm1Mat.xx + r->C[0][0]*dvxxeff + r->C[0][1]*dvyyeff + r->C[0][2]*dvzzeff;
		strial.yy = stnm1Mat.yy + r->C[1][0]*dvxxeff + r->C[1][1]*dvyyeff + r->C[1][2]*dvzzeff;
		strial.zz = stnm1Mat.zz + r->C[2][0]*dvxxeff + r->C[2][1]*dvyyeff + r->C[2][2]*dvzzeff;
		strial.yz = stnm1Mat.yz + r->C[3][3]*dgamyz;
		strial.xz = stnm1Mat.xz + r->C[4][4]*dgamxz;
		strial.xy = stnm1Mat.xy + r->C[5][5]*dgamxy;
		
		is2D = false;
 	}
	else
	{	// rotate current stress (stnm1) to material axes (stnm1Mat) and increment with effective strain increment
		stnm1 = MakeTensor2D(sp->xx,sp->yy,sp->zz,sp->xy);
		stnm1Mat = Rnm1tot.RTVoightR(&stnm1,true,true);
		strial.xx = stnm1Mat.xx + r->C[1][1]*dvxxeff + r->C[1][2]*dvyyeff;
		strial.yy = stnm1Mat.yy + r->C[1][2]*dvxxeff + r->C[2][2]*dvyyeff;
		strial.xy = stnm1Mat.xy + r->C[3][3]*dgamxy;
		
		// sigma(zz)
		strial.zz = stnm1Mat.zz + r->C[4][1]*dvxxeff + r->C[4][2]*dvyyeff + r->C[4][4]*dvzzeff;
		if(np==PLANE_STRAIN_MPM)
		{	strial.zz += r->C[4][1]*r->alpha[5]*er(2,2) + r->C[4][2]*r->alpha[6]*er(2,2);
		}
		
		is2D = true;
	}
	
	// Step 2: rotate plastic strain incrementally
	// done if elastic, or add new plastic strain later if plastic
	Tensor *eplast=mptr->GetAltStrainTensor();		// in global axes
	*eplast = dR.RVoightRT(eplast,false,is2D);
	
    // Step 3: Calculate plastic potential f
	p->aint = mptr->GetHistoryDble(0,0);
	// NEW_HILL two surface options
	double halfsPs = GetHillMagnitude(strial,&h,np);
	double yield = GetYield(p);
	double ftrial;
	if(h.style==SQRT_TERMS)
	{	p->sAQsmag = halfsPs;
		ftrial = p->sAQsmag - yield;
	}
	else
		ftrial = halfsPs - yield*yield;
	
	// Step 4: Done if elastic
	if(ftrial<=0.)
	{	// elastic, update stress and strain energy as usual (lazy method is recalculate elastic update)
		Elastic::LRElasticConstitutiveLaw(mptr,de,er,Rtot,dR,Rnm1tot,np,properties,res);
		return;
    }
    
	// Step 5: Solve for plastic strain increment and new stress and new alpha
    double alpha0 = p->aint;
	Tensor dep = SolveForPlasticIncrement(mptr,np,ftrial,strial,p);
	
	// Step 8: update the particle
	
	// get energies in material axes
	double workEnergy,dispEnergy,resEnergy;
	if(np==THREED_MPM)
	{	// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		workEnergy = strial.xx*de(0,0) + strial.yy*de(1,1) + strial.zz*de(2,2)
						+ 2.*(strial.xy*de(0,1) + strial.xz*de(0,2) + strial.yz*de(1,2));
		
		// Plastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		dispEnergy = strial.xx*dep.xx + strial.yy*dep.yy + strial.zz*dep.zz
						+ strial.xy*dep.xy + strial.xz*dep.xz + strial.yz*dep.yz;
        
		// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		resEnergy = strial.xx*er(0,0) + strial.yy*er(1,1) + strial.zz*er(2,2)
						+ 2.*(strial.xy*er(0,1) + strial.xz*er(0,2) + strial.yz*er(1,2));
	}
	else
	{	// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		workEnergy = strial.xx*de(0,0) + strial.yy*de(1,1) + strial.zz*de(2,2) + 2.*strial.xy*de(0,1);
		
		// Plastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		dispEnergy = strial.xx*dep.xx + strial.yy*dep.yy + strial.zz*dep.zz + strial.xy*dep.xy;
		
		// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		resEnergy = strial.xx*er(0,0) + strial.yy*er(1,1) + strial.zz*er(2,2) + 2.*strial.xy*er(0,1);
	}
	
	// add now
	mptr->AddWorkEnergyAndResidualEnergy(workEnergy,resEnergy);
	
	// add dissipated energy to plastic energy to the particle after subtracting q*dalpha
    dispEnergy -= (GetYield(p)-1.)*(p->aint-alpha0)/sqrt23OversigmaYref;
	mptr->AddPlastEnergy(dispEnergy);
	
	// rotate strain increment to current config and add to eplast in current config
	dep = Rtot.RVoightRT(&dep,false,is2D);
	AddTensor(eplast,&dep);
	
	//  set stress in current configuration on the particle
	Tensor stk = Rtot.RVoightRT(&strial,true,is2D);
	sp->xx = stk.xx;
	sp->yy = stk.yy;
	sp->zz = stk.zz;
	sp->xy = stk.xy;
	if(np==THREED_MPM)
	{	sp->xz = stk.xz;
		sp->yz = stk.yz;
	}
	
    // heat energy is Cv(dT-dTq0) - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adiabatic heating occurs)
    IncrementHeatEnergy(mptr,0.,dispEnergy);
	
	// update hardening variable
	mptr->SetHistoryDble(0,p->aint,0);
	
}

// Solve numerically plastic strain increment by implicit radial return
// stk is trial stress state assuming elastic increment (lambda=0)
// ftrial>0 is in Phi(lambda=0) for elastic increment
// alphak in p->aint in alpha(lambda=0)
Tensor AnisoPlasticity::SolveForPlasticIncrement(MPMBase *mptr,int np,double fk,Tensor &stk,AnisoPlasticProperties *p) const
{
	Tensor dep;
	ZeroTensor(&dep);
	double cutoff = 1.e-10*GetYield(p);
	int nsteps=0;
	
	//int dstep=1000;
	//if(fmobj->mstep>dstep)
	//{	cout << "Start with f = " << fk << ", a = " << p->aint << ", step = "
	//		<< fmobj->mstep << ", cutoff = " << cutoff << endl;
	//}
	
	// NEW_HILL plastic increment
	if(h.style==SQUARED_TERMS)
	{	double lambdak = 0;
		Tensor strial = stk;
		double atrial = p->aint;
		
		// P strial = DPhi/dSigma, (2/3)Q strial, and sqrt((2/3)strial Q strial) stored in p->
		GetHillDfDsigmaPQsigma(strial,np,p,&h);
		
		// get (I + lambdak CP)^{-1} for first step when lambdak=0
		Matrix3 IlamCPInv = Matrix3::Identity();
		double lamInv[3][3];
		IlamCPInv.get(lamInv);
		
		while(nsteps<10)
		{	// stk = sigma(lambdak) = (I + lambda CP)^{-1}strial (from above=strial or end of previous step)
			
			// get P stk, (2/3)Q stk, and sqrt((2/3)stk Q stk) (from above using strial or end of previous step)

			// dsigma/dlambda = -(I + lamdak CP)^{-1) CP.sigma(lambdak)
			Tensor CPSigma,dSigdLam;
			CPSigma.xx = h.CP11*stk.xx + h.CP12*stk.yy + h.CP13*stk.zz;
			CPSigma.yy = h.CP21*stk.xx + h.CP22*stk.yy + h.CP23*stk.zz;
			CPSigma.zz = h.CP31*stk.xx + h.CP32*stk.yy + h.CP33*stk.zz;
			CPSigma.xy = h.CP66*stk.xy;
			if(np==THREED_MPM)
			{   CPSigma.yz = h.CP44*stk.yz;
				CPSigma.xz = h.CP55*stk.xz;
			}
			dSigdLam.xx = -lamInv[0][0]*CPSigma.xx - lamInv[0][1]*CPSigma.yy - lamInv[0][2]*CPSigma.zz;
			dSigdLam.yy = -lamInv[1][0]*CPSigma.xx - lamInv[1][1]*CPSigma.yy - lamInv[1][2]*CPSigma.zz;
			dSigdLam.zz = -lamInv[2][0]*CPSigma.xx - lamInv[2][1]*CPSigma.yy - lamInv[2][2]*CPSigma.zz;
			dSigdLam.xy = -CPSigma.xy/(1.+lambdak*h.CP66);
			if(np==THREED_MPM)
			{   dSigdLam.yz = -CPSigma.yz/(1.+lambdak*h.CP44);
                dSigdLam.xz = -CPSigma.xz/(1.+lambdak*h.CP55);
			}
			
			// hardining terms
			p->aint = atrial + lambdak*p->sAQsmag;
			double hardTerm = 2.*GetYield(p)*GetGPrime(p);
			
			// remaining terms
			double ppterm,a2term;
			if(np==THREED_MPM)
			{	ppterm = DotTensors(&(p->dfdsPsigma),&dSigdLam);
				a2term = p->sAQsmag + lambdak*DotTensors(&(p->CdfQsigma),&dSigdLam)/p->sAQsmag;
			}
			else
			{	ppterm = DotTensors2D(&(p->dfdsPsigma),&dSigdLam);
				a2term = p->sAQsmag + lambdak*DotTensors2D(&(p->CdfQsigma),&dSigdLam)/p->sAQsmag;
			}
			
			// the final scalar derivative
			double dPhidLambda = ppterm - hardTerm*a2term;
			
			// the update
			double dlambda = -fk/dPhidLambda;
			lambdak += dlambda;
			
			// get new stress (I + lambdak CP)^{-1}strial
			// get (I+lambdak CP) - 3X3 upper left, diagonal bottom right
			Matrix3 IlamCP = Matrix3(1.+lambdak*h.CP11,lambdak*h.CP12,lambdak*h.CP13,
									lambdak*h.CP21,1.+lambdak*h.CP22,lambdak*h.CP23,
									lambdak*h.CP31,lambdak*h.CP32,1.+lambdak*h.CP33);
			IlamCPInv = IlamCP.Inverse();
			IlamCPInv.get(lamInv);
			stk.xx = lamInv[0][0]*strial.xx + lamInv[0][1]*strial.yy + lamInv[0][2]*strial.zz;
			stk.yy = lamInv[1][0]*strial.xx + lamInv[1][1]*strial.yy + lamInv[1][2]*strial.zz;
			stk.zz = lamInv[2][0]*strial.xx + lamInv[2][1]*strial.yy + lamInv[2][2]*strial.zz;
			stk.xy = strial.xy/(1.+lambdak*h.CP66);
			if(np==THREED_MPM)
			{   stk.yz = strial.yz/(1.+lambdak*h.CP44);
				stk.xz = strial.xz/(1.+lambdak*h.CP55);
			}
			
			// P stk = DPhi/dSigma, (2/3)Q stk, and sqrt((2/3)stk Q stk)
			GetHillDfDsigmaPQsigma(stk,np,p,&h);

			// get new alpha
			p->aint = atrial + lambdak*p->sAQsmag;

			// check convergenve with new values
			double halfsPs = GetHillMagnitude(stk,&h,np);
			double yield = GetYield(p);
			fk = halfsPs - yield*yield;
			
			// check for convergence
			nsteps++;
			//if(fmobj->mstep>dstep)
			//{	cout << "   " << nsteps << ": f = " << fk << ", a = " << p->aint << ", dlambda = "
			//            << dlambda << ", lambda = " << lambdak << ", dfk/dlam = " << dPhidLambda << endl;
			//}
			if(fabs(fk)<cutoff)
			{	// get final plastic strain
				dep = p->dfdsPsigma;
				ScaleTensor(&dep,lambdak);
				break;
			}
		}
	}
	else
	{	// Newton's method to solve for plastic strain
		while(nsteps<10)
		{	// find dlambda = f/((df/dsig)C(df/dsig) + g'(alphak))
			GetDfCdf(stk,np,p);
			double dlambda = fk/(p->dfCdf + sqrt23OversigmaYref*GetGPrime(p));
			//cout << "     dfCdf=" << p->dfCdf << ", sqrt(2/3)g'(a)/sYRef="
			//            << sqrt23OversigmaYref*GetGPrime(p) << ", dlambda=" << dlambda << endl;

			// update variables
			p->aint += sqrt23OversigmaYref*dlambda;
			
			// next subincrement in plastic strain
			Tensor ddep = p->dfdsPsigma;
			ScaleTensor(&ddep,dlambda);
			
			// update stress
			stk.xx -= dlambda*p->CdfQsigma.xx;
			stk.yy -= dlambda*p->CdfQsigma.yy;
			stk.zz -= dlambda*p->CdfQsigma.zz;
			stk.xy -= dlambda*p->CdfQsigma.xy;
			if(np==THREED_MPM)
			{	stk.xz -= dlambda*p->CdfQsigma.xz;
				stk.yz -= dlambda*p->CdfQsigma.yz;
			}
			
			// total incremental plastic strain accumulated
			AddTensor(&dep,&ddep);
			
			// get new magniture and potential
			p->sAQsmag = GetHillMagnitude(stk,&h,np);
			fk = p->sAQsmag - GetYield(p);
			
			// check for convergence
			nsteps++;
			//cout << "   " << nsteps << ": f = " << fk << ", a = " << p->aint << ", dlambda = " << dlambda << endl;
			if(fabs(fk)<cutoff) break;
			
		}
	}

	// check number of steps
	if(nsteps>hstepsMax)
	{	hstepsMax = nsteps;
#pragma omp critical (output)
        {
            cout << "# AnisoPlasticity took " << nsteps << " steps to converge" << endl;
        }
	}
	
	// return full plastic strain increment
	return dep;
}

// Only called when h.style == SQRT_TERMS
// Find C.df and df.C.df at given stress - store results in plastic property variables active only during the loop
// and only for current material point
void AnisoPlasticity::GetDfCdf(Tensor &stk,int np,AnisoPlasticProperties *p) const
{
	// Get dfds (in dfdsPsigma)
	GetHillDfDsigmaPQsigma(stk,np,p,&h);
	
	// get C df (in CdfQsigma) and df C df (in dfCdf), which needs df/dsig
	ElasticProperties *r = p->ep;
	if(np==THREED_MPM)
	{	p->CdfQsigma.xx = r->C[0][0]*p->dfdsPsigma.xx + r->C[0][1]*p->dfdsPsigma.yy + r->C[0][2]*p->dfdsPsigma.zz;
		p->CdfQsigma.yy = r->C[1][0]*p->dfdsPsigma.xx + r->C[1][1]*p->dfdsPsigma.yy + r->C[1][2]*p->dfdsPsigma.zz;
		p->CdfQsigma.zz = r->C[2][0]*p->dfdsPsigma.xx + r->C[2][1]*p->dfdsPsigma.yy + r->C[2][2]*p->dfdsPsigma.zz;
		p->CdfQsigma.yz = r->C[3][3]*p->dfdsPsigma.yz;
		p->CdfQsigma.xz = r->C[4][4]*p->dfdsPsigma.xz;
		p->CdfQsigma.xy = r->C[5][5]*p->dfdsPsigma.xy;
		p->dfCdf = p->dfdsPsigma.xx*p->CdfQsigma.xx + p->dfdsPsigma.yy*p->CdfQsigma.yy
						+ p->dfdsPsigma.zz*p->CdfQsigma.zz + p->dfdsPsigma.yz*p->CdfQsigma.yz
						+ p->dfdsPsigma.xz*p->CdfQsigma.xz + p->dfdsPsigma.xy*p->CdfQsigma.xy;
	}
	else
	{	p->CdfQsigma.xx = r->C[1][1]*p->dfdsPsigma.xx + r->C[1][2]*p->dfdsPsigma.yy;
		p->CdfQsigma.yy = r->C[1][2]*p->dfdsPsigma.xx + r->C[2][2]*p->dfdsPsigma.yy;
		p->CdfQsigma.xy = r->C[3][3]*p->dfdsPsigma.xy;
		p->CdfQsigma.zz = r->C[4][1]*p->dfdsPsigma.xx + r->C[4][2]*p->dfdsPsigma.yy + r->C[4][4]*p->dfdsPsigma.zz;
		p->dfCdf = p->dfdsPsigma.xx*p->CdfQsigma.xx + p->dfdsPsigma.yy*p->CdfQsigma.yy
							+ p->dfdsPsigma.xy*p->CdfQsigma.xy + p->dfdsPsigma.zz*p->CdfQsigma.zz;
	}
}

#pragma mark AnisoPlasticity::Class Methods

// NEW_HILL properties
// get properties for Hill style 2 (or revert to style 1
// This is a static because it is called by OrthoPlasticSofetning too
void AnisoPlasticity::FillHillStyleProperties(int np,HillProperties &h,ElasticProperties &pr)
{
	if(h.style==SQUARED_TERMS)
	{	// matrix C.P (non-zero elements)
		// Elements of P = 2A are scaled by rho^2 and P by 1/rho so CP is * rho
		if(np==THREED_MPM)
		{	h.CP11 = 2.*( pr.C[0][0]*h.syxx2 - pr.C[0][1]*h.H - pr.C[0][2]*h.G );
			h.CP12 = 2.*( -pr.C[0][0]*h.H + pr.C[0][1]*h.syyy2 - pr.C[0][2]*h.F );
			h.CP13 = 2.*( -pr.C[0][0]*h.G - pr.C[0][1]*h.F + pr.C[0][2]*h.syzz2 );
			h.CP21 = 2.*( pr.C[0][1]*h.syxx2 - pr.C[1][1]*h.H - pr.C[1][2]*h.G );
			h.CP22 = 2.*( -pr.C[0][1]*h.H + pr.C[1][1]*h.syyy2 - pr.C[1][2]*h.F );
			h.CP23 = 2.*( -pr.C[0][1]*h.G - pr.C[1][1]*h.F + pr.C[1][2]*h.syzz2 );
			h.CP31 = 2.*( pr.C[0][2]*h.syxx2 - pr.C[1][2]*h.H - pr.C[2][2]*h.G );
			h.CP32 = 2.*( -pr.C[0][2]*h.H + pr.C[1][2]*h.syyy2 - pr.C[2][2]*h.F );
			h.CP33 = 2.*( -pr.C[0][2]*h.G - pr.C[1][2]*h.F + pr.C[2][2]*h.syzz2 );
			h.CP44 = 2.*pr.C[3][3]*h.L;
			h.CP55 = 2.*pr.C[4][4]*h.M;
			h.CP66 = 2.*pr.C[5][5]*h.N;
		}
		else
		{	h.CP11 = 2.*( pr.C[1][1]*h.syxx2 - pr.C[1][2]*h.H - pr.C[4][1]*h.G );
			h.CP12 = 2.*( -pr.C[1][1]*h.H + pr.C[1][2]*h.syyy2 - pr.C[4][1]*h.F );
			h.CP13 = 2.*( -pr.C[1][1]*h.G - pr.C[1][2]*h.F + pr.C[4][1]*h.syzz2 );
			h.CP21 = 2.*( pr.C[1][2]*h.syxx2 - pr.C[2][2]*h.H - pr.C[4][2]*h.G );
			h.CP22 = 2.*( -pr.C[1][2]*h.H + pr.C[2][2]*h.syyy2 - pr.C[4][2]*h.F );
			h.CP23 = 2.*( -pr.C[1][2]*h.G - pr.C[2][2]*h.F + pr.C[4][2]*h.syzz2 );
			h.CP31 = 2.*( pr.C[4][1]*h.syxx2 - pr.C[4][2]*h.H - pr.C[4][4]*h.G );
			h.CP32 = 2.*( -pr.C[4][1]*h.H + pr.C[4][2]*h.syyy2 - pr.C[4][4]*h.F );
			h.CP33 = 2.*( -pr.C[4][1]*h.G - pr.C[4][2]*h.F + pr.C[4][4]*h.syzz2 );
			h.CP44 = 0.;
			h.CP55 = 0.;
			h.CP66 = 2.*pr.C[3][3]*h.N;
		}
		
		// Q = PZP matrix (scale by rho^2 and symmetric)
		// same 2D and 3D
		h.Q11 = 8.*(h.G*h.G + h.G*h.H + h.H*h.H);
		h.Q12 = 4.*(h.F*(h.G-h.H) - h.H*(h.G+2.*h.H));
		h.Q13 = 4.*(-h.G*(h.F+h.G) + h.F*h.H - h.G*(h.G+h.H));
		h.Q22 = 8.*(h.F*h.F + h.F*h.H + h.H*h.H);
		h.Q23 = 4.*(-2.*h.F*h.F + h.G*h.H - h.F*(h.G+h.H));
		h.Q33 = 8.*(h.F*h.F + h.F*h.G + h.G*h.G);
		h.Q44 = 2.*h.L*h.L;
		h.Q55 = 2.*h.M*h.M;
		h.Q66 = 2.*h.N*h.N;
	}
	else
	{	// revert to valid option is was not set to 1 or 2
		h.style = SQRT_TERMS;
	}
}

// print just yield properties to output window
void AnisoPlasticity::PrintAPYieldProperties(double yld1,double yld2,double yld3,double yld23,double yld13,double yld12)
{
    if(yld1>=0.)
        PrintProperty("syxx",yld1*UnitsController::Scaling(1.e-6),"");
    else
        PrintProperty("syxx= inf",false);
        
    if(yld2>=0.)
        PrintProperty("syyy",yld2*UnitsController::Scaling(1.e-6),"");
    else
        PrintProperty("syyy= inf",false);
        
    if(yld3>=0.)
        PrintProperty("syzz",yld3*UnitsController::Scaling(1.e-6),"");
    else
        PrintProperty("syzz= inf",false);
        
    cout << endl;

    if(yld23>=0.)
        PrintProperty("tyyz",yld23*UnitsController::Scaling(1.e-6),"");
    else
        PrintProperty("tyyz= inf",false);

    if(yld13>=0.)
        PrintProperty("tyxz",yld13*UnitsController::Scaling(1.e-6),"");
    else
        PrintProperty("tyxz= inf",false);
    
    if(yld12>=0.)
        PrintProperty("tyxy",yld12*UnitsController::Scaling(1.e-6),"");
    else
        PrintProperty("tyxy= inf",false);
    
    cout << endl;
}

// NEW_HILL get magnitude (static)
// h->style==SQRT_TERMS: Get sqrt(sAs) where srot is stress in material axis system (negative retruned as zero)
// h->style==SQUARED_TERMS: Get (1/2)sPs where srot is stress in material axis system (negative retruned as zero)
double AnisoPlasticity::GetHillMagnitude(Tensor &srot,const HillProperties *h,int np)
{
    // initialize (with 3D shear)
    double sAs = np==THREED_MPM ? srot.yz*srot.yz*h->L + srot.xz*srot.xz*h->M : 0. ;
    
    // normal and xy shear terms
    double dyz = srot.yy-srot.zz;
    double dxz = srot.xx-srot.zz;
    double dxy = srot.xx-srot.yy;
    sAs += h->F*dyz*dyz + h->G*dxz*dxz + h->H*dxy*dxy + srot.xy*srot.xy*h->N;
    
    // check on negative sAs can happen due to round-off error when stresses near zero
	if(h->style==SQRT_TERMS)
		return sAs>0. ? sqrt(sAs) : 0 ;
	else
		return sAs>0. ? sAs : 0 ;
}

// NEW_HILL get tensors (static)
// h->stylee=SQRT_TERMS: Find dfds = A sigma/sqrt(sigma.Asigma) in material axes
// h->style=SQUARED_TERMS: Find P sigma = 2 A sigma, (2/3)Q sigma and sqrt((2/3)sigma.CdfQsigma)
// load all into plastic properties
void AnisoPlasticity::GetHillDfDsigmaPQsigma(Tensor &stk,int np,AnisoPlasticProperties *p,const HillProperties *h)
{
	if(h->style==SQUARED_TERMS)
	{	// Get P sigma = 2 A sigma
		p->dfdsPsigma.xx = 2.*(h->syxx2*stk.xx - h->H*stk.yy - h->G*stk.zz);
		p->dfdsPsigma.yy = 2.*(-h->H*stk.xx + h->syyy2*stk.yy - h->F*stk.zz);
		p->dfdsPsigma.zz = 2.*(-h->G*stk.xx - h->F*stk.yy + h->syzz2*stk.zz);
		p->dfdsPsigma.xy = 2.*h->N*stk.xy;
		if(np==THREED_MPM)
		{   p->dfdsPsigma.yz = 2.*h->L*stk.yz;
			p->dfdsPsigma.xz = 2.*h->M*stk.xz;
		}
		
		// Get (2/3)Q sigma(lambdak)
		p->CdfQsigma.xx = h->Q11*stk.xx + h->Q12*stk.yy + h->Q13*stk.zz;
		p->CdfQsigma.yy = h->Q12*stk.xx + h->Q22*stk.yy + h->Q23*stk.zz;
		p->CdfQsigma.zz = h->Q13*stk.xx + h->Q23*stk.yy + h->Q33*stk.zz;
		p->CdfQsigma.xy = h->Q66*stk.xy;
		if(np==THREED_MPM)
		{   p->CdfQsigma.yz = h->Q44*stk.yz;
			p->CdfQsigma.xz = h->Q55*stk.xz;
		}
		ScaleTensor(&(p->CdfQsigma),2./3.);
		
		// get sqrt((2/3)sigma.CdfQsigma)
		p->sAQsmag = stk.xx*p->CdfQsigma.xx + stk.yy*p->CdfQsigma.yy + stk.zz*p->CdfQsigma.zz + stk.xy*p->CdfQsigma.xy;
		if(np==THREED_MPM) p->sAQsmag += stk.yz*p->CdfQsigma.yz + stk.xz*p->CdfQsigma.xz;
		p->sAQsmag = sqrt(p->sAQsmag);
	}
	else
	{	// using SQRT_TERMS
		
		// dfds = A.sigma/rootSAS
		if(p->sAQsmag>0.)
		{   p->dfdsPsigma.xx = (h->syxx2*stk.xx - h->H*stk.yy - h->G*stk.zz) / p->sAQsmag;
			p->dfdsPsigma.yy = (-h->H*stk.xx + h->syyy2*stk.yy - h->F*stk.zz) / p->sAQsmag;
			p->dfdsPsigma.zz = (-h->G*stk.xx - h->F*stk.yy + h->syzz2*stk.zz) / p->sAQsmag;
			p->dfdsPsigma.xy = h->N*stk.xy / p->sAQsmag;
			if(np==THREED_MPM)
			{   p->dfdsPsigma.yz = h->L*stk.yz / p->sAQsmag;
				p->dfdsPsigma.xz = h->M*stk.xz / p->sAQsmag;
			}
		}
		
		else
		{   // negative root implies zero stress within roundoff error
			p->dfdsPsigma.xx = p->dfdsPsigma.yy = p->dfdsPsigma.zz = 0.;
			p->dfdsPsigma.yz = p->dfdsPsigma.xz = p->dfdsPsigma.xy = 0.;
		}
	}
}

#pragma mark AnisoPlasticity::Accessors

// store plastic strain in alt strain
int AnisoPlasticity::AltStrainContains(void) const
{	return ENG_BIOT_PLASTIC_STRAIN;
}
