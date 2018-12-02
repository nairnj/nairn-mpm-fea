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
#include "AnisoPlasticity.hpp"
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
	
	return Orthotropic::InputMaterialProperty(xName,input,gScaling);
}

// verify settings and some initial calculations
const char *AnisoPlasticity::VerifyAndLoadProperties(int np)
{
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
	{	double rsxx=0.,rsyy=0.,rszz=0.;
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
    {	syxxred2=rho/syxx;
		syxxred2*=syxxred2;
	}
	else
		syxxred2=0.;		// 1/inf^2
	if(syyy>0.)
    {	syyyred2=rho/syyy;
		syyyred2*=syyyred2;
	}
	else
		syyyred2=0.;		// 1/inf^2
	if(syzz>0.)
	{	syzzred2=rho/syzz;
		syzzred2*=syzzred2;
	}
	else
		syzzred2=0.;		// 1/inf^2
	
	// reciprocals of reduced shear yield stresses
	if(tyxy>0.)
	{	tyxyred2=rho/tyxy;
		tyxyred2*=tyxyred2;
	}
	else
		tyxyred2=0.;		// 1/inf^2
	if(tyxz>0.)
	{	tyxzred2=rho/tyxz;
		tyxzred2*=tyxzred2;
	}
	else
		tyxzred2=0.;		// 1/inf^2
	if(tyyz>0.)
	{	tyyzred2=rho/tyyz;
		tyyzred2*=tyyzred2;
	}
	else
		tyyzred2=0.;		// 1/inf^2
	
	// combination terms
	fTerm = 0.5*(syyyred2 + syzzred2 - syxxred2);
	gTerm = 0.5*(syzzred2 + syxxred2 - syyyred2);
	hTerm = 0.5*(syxxred2 + syyyred2 - syzzred2);
	
	// reference yield stress (it is acutally sqrt(2/3)/sigma(Y,ref)
	sigmaYref = fTerm+gTerm+hTerm;
	if(sigmaYref<=0.) return "(F+G+H) for Hill plastic potential must be positive";
	sigmaYref = sqrt(2.*sigmaYref/3.);
	
	// for convergence problems
	warnNonconvergence=warnings.CreateWarning("anisotropic plastic algorithm failed to converge",-1,3);
	
	// call super class
	return Orthotropic::VerifyAndLoadProperties(np);
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

// print just yield properties to output window
void AnisoPlasticity::PrintYieldProperties(void) const
{
	if(syxx>=0.)
		PrintProperty("yld1",syxx*UnitsController::Scaling(1.e-6),"");
	else
		PrintProperty("yld1= inf",false);
		
	if(syyy>=0.)
		PrintProperty("yld2",syyy*UnitsController::Scaling(1.e-6),"");
	else
		PrintProperty("yld2= inf",false);
		
	if(syzz>=0.)
		PrintProperty("yld3",syzz*UnitsController::Scaling(1.e-6),"");
	else
		PrintProperty("yld3= inf",false);
		
    cout << endl;

	if(tyyz>=0.)
		PrintProperty("yld23",tyyz*UnitsController::Scaling(1.e-6),"");
	else
		PrintProperty("yld23= inf",false);

	if(tyxz>=0.)
		PrintProperty("yld13",tyxz*UnitsController::Scaling(1.e-6),"");
	else
		PrintProperty("yld13= inf",false);
	
	if(tyxy>=0.)
		PrintProperty("yld12",tyxy*UnitsController::Scaling(1.e-6),"");
	else
		PrintProperty("yld12= inf",false);
	
    cout << endl;
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
// Similar, strain increments are rotated back to find work energy (done in global system)
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
	p->sAsmag = GetMagnitudeHill(strial,np);
	double ftrial = p->sAsmag - GetYield(p);
	
	// Step 4: Done if elastic
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual (lazy method is recalculate elastic update)
		Elastic::LRElasticConstitutiveLaw(mptr,de,er,Rtot,dR,Rnm1tot,np,properties,res);
		return;
    }
    
	// Step 5: Solve for plastic strain increment
	Tensor dep = SolveForPlasticIncrement(mptr,np,ftrial,strial,p);
	
	// Step 8: update the particle
	
	// get energies in material axes (Need more plastic energy - c.f. istropic term)
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
		workEnergy = strial.xx*de(0,0) + strial.yy*de(1,1) + strial.zz*de(1,1) + 2.*strial.xy*de(0,1);
		
		// Plastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		dispEnergy = strial.xx*dep.xx + strial.yy*dep.yy + strial.zz*dep.zz + 2.*strial.xy*dep.xy;
		
		// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		resEnergy = strial.xx*er(0,0) + strial.yy*er(1,1) + strial.zz*er(2,2) + 2.*strial.xy*er(0,1);
	}
	
	// add now
	mptr->AddWorkEnergyAndResidualEnergy(workEnergy,resEnergy);
	
	// add dissipated energy to plastic energy to the particle
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

// Solve numerically plastic strain increment by explicit radial return
Tensor AnisoPlasticity::SolveForPlasticIncrement(MPMBase *mptr,int np,double fk,Tensor &stk,AnisoPlasticProperties *p) const
{
	// collect total plastic strain increment from subincrements
	Tensor dep;
	ZeroTensor(&dep);
	// stk is trial stress state assuming elastic increment
	// ftrial>0 is in fk for elastic increment
	// alphak in p->aint
	double cutoff = 1.e-10;
	int nsteps=0;
	//cout << "Start with f = " << fk << ", a = " << p->aint << endl;
	while(nsteps<10)
	{	// find dlambda = f/((df/dsig)C(df/dsig) + g'(alphak))
		GetDfCdf(stk,np,p);
		double dlambda = fk/(p->dfCdf + sigmaYref*GetGPrime(p));
		
		// update variables
		p->aint += sigmaYref*dlambda;
		
		// next subincrement in plastic strain
		Tensor ddep = p->dfds;
		ScaleTensor(&ddep,dlambda);
		
		// update stress
		stk.xx -= dlambda*p->Cdf.xx;
		stk.yy -= dlambda*p->Cdf.yy;
		stk.zz -= dlambda*p->Cdf.zz;
		stk.xy -= dlambda*p->Cdf.xy;
		if(np==THREED_MPM)
		{	stk.xz -= dlambda*p->Cdf.xz;
			stk.yz -= dlambda*p->Cdf.yz;
		}
		
		// total incremental plastic strain accumulated
		AddTensor(&dep,&ddep);
		
		// get new magniture and potential
		p->sAsmag = GetMagnitudeHill(stk,np);
		fk = p->sAsmag - GetYield(p);
		
		// check for convergence
		nsteps++;
		//cout << "   " << nsteps << ": f = " << fk << ", a = " << p->aint << ", dlambda = " << dlambda << endl;;
		if(fabs(fk)<cutoff) break;
		
	}
	
	// check number of steps
	if(nsteps>hstepsMax)
	{	hstepsMax = nsteps;
		cout << "# Took " << nsteps << " steps to converge" << endl;
	}
	
	// return full plastic strain increment
	return dep;
}

// Get sqrt(s As) where srot is stress in material axis system
double AnisoPlasticity::GetMagnitudeHill(Tensor &srot,int np) const
{
	// initialize (with 3D shear)
	double sAs = np==THREED_MPM ? srot.xz*srot.xz*tyxzred2 + srot.yz*srot.yz*tyyzred2 : 0. ;
	
	// normal and xy shear terms
	double dyz = srot.yy-srot.zz;
	double dxz = srot.xx-srot.zz;
	double dxy = srot.xx-srot.yy;
	sAs += fTerm*dyz*dyz + gTerm*dxz*dxz + hTerm*dxy*dxy + srot.xy*srot.xy*tyxyred2;
	
	// check on negative sAs can happen due to round-off error when stresses near zero
	return sAs>0. ? sqrt(sAs) : 0 ;
}

// Find C.df and df.C.df at given stress - store results in plastic property variables active only during the loop
// and only for current material point
void AnisoPlasticity::GetDfCdf(Tensor &stk,int np,AnisoPlasticProperties *p) const
{
	// Get dfds
	GetDfDsigma(stk,np,p);
	
	// get C df and df C df, which needs df/dsig
	ElasticProperties *r = p->ep;
	if(np==THREED_MPM)
	{	p->Cdf.xx = r->C[0][0]*p->dfds.xx + r->C[0][1]*p->dfds.yy + r->C[0][2]*p->dfds.zz;
		p->Cdf.yy = r->C[1][0]*p->dfds.xx + r->C[1][1]*p->dfds.yy + r->C[1][2]*p->dfds.zz;
		p->Cdf.zz = r->C[2][0]*p->dfds.xx + r->C[2][1]*p->dfds.yy + r->C[2][2]*p->dfds.zz;
		p->Cdf.yz = r->C[3][3]*p->dfds.yz;
		p->Cdf.xz = r->C[4][4]*p->dfds.xz;
		p->Cdf.xy = r->C[5][5]*p->dfds.xy;
		p->dfCdf = p->dfds.xx*p->Cdf.xx + p->dfds.yy*p->Cdf.yy + p->dfds.zz*p->Cdf.zz
						+ p->dfds.yz*p->Cdf.yz + p->dfds.xz*p->Cdf.xz + p->dfds.xy*p->Cdf.xy;
	}
	else
	{	p->Cdf.xx = r->C[1][1]*p->dfds.xx + r->C[1][2]*p->dfds.yy;
		p->Cdf.yy = r->C[1][2]*p->dfds.xx + r->C[2][2]*p->dfds.yy;
		p->Cdf.xy = r->C[3][3]*p->dfds.xy;
		p->Cdf.zz = r->C[4][1]*p->dfds.xx + r->C[4][2]*p->dfds.yy + r->C[4][4]*p->dfds.zz;
		p->dfCdf = p->dfds.xx*p->Cdf.xx + p->dfds.yy*p->Cdf.yy + p->dfds.xy*p->Cdf.xy + p->dfds.zz*p->Cdf.zz;
	}
}

// Find dfds = A sigma/sqrt(sigma.Asigma) in material axes
// load it into plastic properties
void AnisoPlasticity::GetDfDsigma(Tensor &stk,int np,AnisoPlasticProperties *p) const
{
	// df = A.sigma/rootSAS
	if(p->sAsmag>0.)
	{	p->dfds.xx = (syxxred2*stk.xx - hTerm*stk.yy - gTerm*stk.zz) / p->sAsmag;
		p->dfds.yy = (-hTerm*stk.xx + syyyred2*stk.yy - fTerm*stk.zz) / p->sAsmag;
		p->dfds.zz = (-gTerm*stk.xx - fTerm*stk.yy + syzzred2*stk.zz) / p->sAsmag;
		p->dfds.xy = tyxyred2*stk.xy / p->sAsmag;
		if(np==THREED_MPM)
		{	p->dfds.yz = tyyzred2*stk.yz / p->sAsmag;
			p->dfds.xz = tyxzred2*stk.xz / p->sAsmag;
		}
	}
	
	else
	{	// negative root implies zero stress within roundoff error
		p->dfds.xx = p->dfds.yy = p->dfds.zz = p->dfds.yz = p->dfds.xz = p->dfds.xy = 0.;
	}
}

#pragma mark AnisoPlasticity::Accessors

// store plastic strain in alt strain
int AnisoPlasticity::AltStrainContains(void) const
{	return ENG_BIOT_PLASTIC_STRAIN;
}
