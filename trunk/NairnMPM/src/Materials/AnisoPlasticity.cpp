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

#include "AnisoPlasticity.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/MPMWarnings.hpp"

#define ITERATIVE_STRESS_UPDATE

// class statics
int AnisoPlasticity::warnNonconvergence;

double maxLambda = 0.;

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
	
	// reciprocals of reduced normal yield stresses
	if(syxx>=0.)
    {	syxxred2=rho/syxx;
		syxxred2*=syxxred2;
	}
	else
		syxxred2=0.;		// 1/inf^2
	if(syyy>=0.)
    {	syyyred2=rho/syyy;
		syyyred2*=syyyred2;
	}
	else
		syyyred2=0.;		// 1/inf^2
	if(syzz>=0.)
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
	
	warnNonconvergence=warnings.CreateWarning("anisotropic plastic algorithm failed to converge",-1,3);
	
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
    return useLargeRotation ? sizeof(LRAnisoPlasticProperties) : sizeof(AnisoPlasticProperties) ;
}

// Get current anisotropic properties (NULL on memory error)
void *AnisoPlasticity::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
{
	// large rotation does not need to rotate stiffness matrix
	if(useLargeRotation)
	{	// fill plastic properties
		LRAnisoPlasticProperties *p = (LRAnisoPlasticProperties *)matBuffer;
		p->ep = (ElasticProperties *)&pr;
		return p;
	}
	
	// full plastic properties
	AnisoPlasticProperties *p = (AnisoPlasticProperties *)matBuffer;
	
	if(np!=THREED_MPM)
	{	double s,c;
		mptr->Get2DSinCos(&s,&c);
		FillElasticProperties2D(&(p->ep),TRUE,s,c,np);
	}
	else
		FillElasticProperties3D(mptr,&(p->ep),np);
	
	return p;
}

// return pointer to elastic properties
ElasticProperties *AnisoPlasticity::GetElasticPropertiesPointer(void *properties) const
{
	if(useLargeRotation)
	{	LRAnisoPlasticProperties *p = (LRAnisoPlasticProperties *)properties;
		return p->ep;
	}
	return Elastic::GetElasticPropertiesPointer(properties);
}

#pragma mark AnisoPlasticity::Methods (Large Rotation)

// Once stress deformation has been decomposed, finish calculations in the material axis system
// When stress are found, they are rotated back to the global axes (using Rtot and dR)
// Similar, strain increments are rotated back to find work energy (done in global system)
void AnisoPlasticity::LRElasticConstitutiveLaw(MPMBase *mptr,Matrix3 de,Matrix3 er,Matrix3 Rtot,Matrix3 dR,
											   Matrix3 *Rnm1tot,int np,void *properties,ResidualStrains *res) const
{
	// Properties
	LRAnisoPlasticProperties *p = (LRAnisoPlasticProperties *)properties;
	ElasticProperties *r = p->ep;
	
	// get stress increment in material axes by first rotating current stress
	// to material and then adding elastic increment
	Tensor *sp = mptr->GetStressTensor();
	
	// Step 1: Get trial stress in material axes assuming increment is elastic
	Matrix3 strial;
	Matrix3 deeff = de-er;							// in material axes
	Tensor *eplast=mptr->GetAltStrainTensor();		// in global axes
	Matrix3 stnm1,stnm1Mat,etn;
	if(np==THREED_MPM)
	{	// rotate current stress (stnm1) to material axes (stnm1Mat) and increment with effective strain increment
		stnm1 = Matrix3(sp->xx,sp->xy,sp->xz,sp->xy,sp->yy,sp->yz,sp->xz,sp->yz,sp->zz);
		stnm1Mat = stnm1.RTMR(*Rnm1tot);
		strial.setIs2D(false);
		strial(0,0) = stnm1Mat(0,0) + r->C[0][0]*deeff(0,0) + r->C[0][1]*deeff(1,1) + r->C[0][2]*deeff(2,2);
		strial(1,1) = stnm1Mat(1,1) + r->C[1][0]*deeff(0,0) + r->C[1][1]*deeff(1,1) + r->C[1][2]*deeff(2,2);
		strial(2,2) = stnm1Mat(2,2) + r->C[2][0]*deeff(0,0) + r->C[2][1]*deeff(1,1) + r->C[2][2]*deeff(2,2);
		
		// deeeff has tensor shear strains
		strial(1,2) = stnm1Mat(1,2) + 2.*r->C[3][3]*deeff(1,2);
		strial(2,1) = strial(1,2);
		strial(0,2) = stnm1Mat(0,2) + 2.*r->C[4][4]*deeff(0,2);
		strial(2,0) = strial(0,2);
		strial(0,1) = stnm1Mat(0,1) + 2.*r->C[5][5]*deeff(0,1);
		strial(1,0) = strial(0,1);
		
		// plastic strain tensor (tensorial shear strains) in current configuration step n-1
		etn = Matrix3(eplast->xx,0.5*eplast->xy,0.5*eplast->xz,0.5*eplast->xy,eplast->yy,0.5*eplast->yz,
					0.5*eplast->xz,0.5*eplast->yz,eplast->zz);
 	}
	else
	{	// rotate current stress (stnm1) to material axes (stnm1Mat) and increment with effective strain increment
		stnm1 = Matrix3(sp->xx,sp->xy,sp->xy,sp->yy,sp->zz);
		stnm1Mat = stnm1.RTMR(*Rnm1tot);
		strial.setIs2D(true);
		strial(0,0) = stnm1Mat(0,0) + r->C[1][1]*deeff(0,0) + r->C[1][2]*deeff(1,1);
		strial(1,1) = stnm1Mat(1,1) + r->C[1][2]*deeff(0,0) + r->C[2][2]*deeff(1,1);
		
		// deeeff has tensor shear strains
		strial(0,1) = stnm1Mat(0,1) + 2.*r->C[3][3]*deeff(0,1);
		strial(1,0) = strial(0,1);
		
		// sigma(zz)
		strial(2,2) = stnm1Mat(2,2) + r->C[4][1]*deeff(0,0) + r->C[4][2]*deeff(1,1) + r->C[4][4]*deeff(2,2);
		if(np==PLANE_STRAIN_MPM)
		{	strial(2,2) += r->C[4][1]*r->alpha[5]*er(2,2) + r->C[4][2]*r->alpha[6]*er(2,2);
		}
		
		// plastic strain tensor (tensorial shear strains) in current configuration step n-1
		etn = Matrix3(eplast->xx,0.5*eplast->xy,0.5*eplast->xy,eplast->yy,eplast->zz);
	}
	
	// Step 2: rotate plastic strain
	// done if elastic, or add new plastic strain later if plastic
	Matrix3 etr = etn.RMRT(dR);
	eplast->xx = etr(0,0);
	eplast->yy = etr(1,1);
	eplast->zz = etr(2,2);
	eplast->xy = 2.*etr(0,1);
	if(np==THREED_MPM)
	{	eplast->xz = 2.*etr(0,2);
		eplast->yz = 2.*etr(1,2);
	}
	
	// Step 3: Rotation for plastic potential (not needed here because using material axes)
	
    // Step 4: Calculate plastic potential f
	UpdateTrialAlpha(mptr,np,&(p->hp));
	double sAstrial = GetMagnitudeHill(strial,np);
	double ftrial = sAstrial - GetYield(&(p->hp));
	
	// Step 5: Done if elastic
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
		Elastic::LRElasticConstitutiveLaw(mptr,de,er,Rtot,dR,Rnm1tot,np,properties,res);
		return;
    }
    
	// Step 6: Solve for lambda
	Matrix3 stk;
	double lambda = LRSolveForLambdaAP(mptr,np,ftrial,strial,stk,p);
	
    // Step 7: Plastic strain (tensorial) increments on particle
	Matrix3 dep;
	double dexyp=0.5*lambda*p->dfds.xy;
	if(np==THREED_MPM)
	{	double deyzp=0.5*lambda*p->dfds.yz;
		double dexzp=0.5*lambda*p->dfds.xz;
		dep = Matrix3(lambda*p->dfds.xx,dexyp,dexzp,dexyp,lambda*p->dfds.yy,deyzp,
					  dexzp,deyzp,lambda*p->dfds.zz);
	}
	else
	{	dep = Matrix3(lambda*p->dfds.xx,dexyp,dexyp,lambda*p->dfds.yy,lambda*p->dfds.zz);
	}
	
	// rotate to current config and add to eplast in current config and make engineering shear strain
	dep = dep.RMRT(Rtot);
	eplast->xx += dep(0,0);
	eplast->yy += dep(1,1);
	eplast->zz += dep(2,2);
	eplast->xy += 2.*dep(0,1);
	if(np==THREED_MPM)
	{	eplast->xz += 2.*dep(0,2);
		eplast->yz += 2.*dep(1,2);
	}
	
	// Step 8: get elastic strain increment
	// Step 9: stress increment on the particle
	// Here stk is final stress in materials axes, the algorithm is
	//		sigma = dR stnm1 dRT + Rtot (stk- stnm1Mat) RtotT = Rtot stk RtotT
	Matrix3 str = stk.RMRT(Rtot);
	sp->xx = str(0,0);
	sp->yy = str(1,1);
	sp->zz = str(2,2);
	sp->xy = str(0,1);
	if(np==THREED_MPM)
	{	sp->xz = str(0,2);
		sp->yz = str(1,2);
	}
	// Step 10: Increment energies
	// rotate strains into global axes
	de = de.RMRT(Rtot);
	er = er.RMRT(Rtot);
	Matrix3 st0 = stnm1.RMRT(dR);
	double workEnergy,dispEnergy,resEnergy;

	if(np==THREED_MPM)
	{	// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		workEnergy = 0.5*((st0(0,0)+sp->xx)*de(0,0) + (st0(1,1)+sp->yy)*de(1,1) + (st0(2,2)+sp->zz)*de(1,1))
								+ (st0(0,1)+sp->xy)*de(0,1) + (st0(0,2)+sp->xz)*de(0,2) + (st0(1,2)+sp->yz)*de(1,2);
		
		// Plastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		dispEnergy = 0.5*((st0(0,0)+sp->xx)*dep(0,0) + (st0(1,1)+sp->yy)*dep(1,1) + (st0(2,2)+sp->zz)*dep(2,2))
								+ (st0(0,1)+sp->xy)*dep(0,1) + (st0(0,2)+sp->xz)*dep(0,2) + (st0(1,2)+sp->yz)*dep(1,2);
		
		// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		resEnergy = 0.5*((st0(0,0)+sp->xx)*er(0,0) + (st0(1,1)+sp->yy)*er(1,1) + (st0(2,2)+sp->zz)*er(2,2))
								+ (st0(0,1)+sp->xy)*er(0,1) + (st0(0,2)+sp->xz)*er(0,2) + (st0(1,2)+sp->yz)*er(1,2);
	}
	else
	{	// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		workEnergy = 0.5*((st0(0,0)+sp->xx)*de(0,0) + (st0(1,1)+sp->yy)*de(1,1) + (st0(2,2)+sp->zz)*de(1,1))
									+ (st0(0,1)+sp->xy)*de(0,1);
		
		// Plastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		dispEnergy = 0.5*((st0(0,0)+sp->xx)*dep(0,0) + (st0(1,1)+sp->yy)*dep(1,1) + (st0(2,2)+sp->zz)*dep(2,2))
								 + (st0(0,1)+sp->xy)*dep(0,1);
		
		// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		resEnergy = 0.5*((st0(0,0)+sp->xx)*er(0,0) + (st0(1,1)+sp->yy)*er(1,1) + (st0(2,2)+sp->zz)*er(2,2))
								+ (st0(0,1)+sp->xy)*er(1,1);
	}
	
	// add now
	mptr->AddWorkEnergyAndResidualEnergy(workEnergy + dispEnergy,resEnergy);
	
	// add dissipated energy to plastic energy to the particle
    mptr->AddPlastEnergy(dispEnergy);
	
    // heat energy is Cv(dT-dTq0) - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adiabatic heating occurs)
    IncrementHeatEnergy(mptr,res->dT,0.,dispEnergy);
	
	// Step 11: update internal variables
	UpdatePlasticInternal(mptr,np,&(p->hp));
}

// Get sqrt(s As) where s is stress in material axis system
double AnisoPlasticity::GetMagnitudeHill(Matrix3 &srot,int np) const
{
	// initialize (with 3D shear)
	double sAs = np==THREED_MPM ? srot(0,2)*srot(0,2)*tyxzred2 + srot(1,2)*srot(1,2)*tyyzred2 : 0. ;
	
	// normal and xy shear terms
	double dyz = srot(1,1)-srot(2,2);
	double dxz = srot(0,0)-srot(2,2);
	double dxy = srot(0,0)-srot(1,1);
	sAs += fTerm*dyz*dyz + gTerm*dxz*dxz + hTerm*dxy*dxy + srot(0,1)*srot(0,1)*tyxyred2;
	
	// check on negative sAs can happen due to round-off error when stresses near zero
	return sAs>0. ? sqrt(sAs) : 0 ;
}

// Find C.df and df.C.df at given stress - store results in plastic property variables active only during the loop
// and only for current material point
void AnisoPlasticity::LRGetDfCdf(Matrix3 &stk,int np,LRAnisoPlasticProperties *p) const
{
	// get C df and df C df, which need df/dsig (Function of yield criteria, and normally only current stress in stk)
	LRGetDfDsigma(stk,np,p);
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

// Find A srot/sqrt(s As) in rotated coordinates (in dfdsijrot) and then rotate the result to
//   get df term such that dep = lamda df in the analysis coordinates. Rotated term
//	 is stored in dfds tensor in plastic properties
// As find R df = -h and store in minush hardning property
void AnisoPlasticity::LRGetDfDsigma(Matrix3 &st0,int np,LRAnisoPlasticProperties *p) const
{
	// clockwise rotation from analysis to material axes
	double rootSAS = GetMagnitudeHill(st0,np);
	
	// df = A.sigma/rootSAS
	if(rootSAS>0.)
	{	p->dfds.xx = (syxxred2*st0(0,0) - hTerm*st0(1,1) - gTerm*st0(2,2)) / rootSAS;
		p->dfds.yy = (-hTerm*st0(0,0) + syyyred2*st0(1,1) - fTerm*st0(2,2)) / rootSAS;
		p->dfds.zz = (-gTerm*st0(0,0) - fTerm*st0(1,1) + syzzred2*st0(2,2)) / rootSAS;
		p->dfds.xy = tyxyred2*st0(0,1) / rootSAS;
		if(np==THREED_MPM)
		{	p->dfds.yz = tyyzred2*st0(1,2) / rootSAS;
			p->dfds.xz = tyxzred2*st0(0,2) / rootSAS;
			
			// for use in alpha upate (it is -h or dalpha/lambda)
			p->hp.minush = p->dfds.xx*p->dfds.xx + p->dfds.yy*p->dfds.yy + p->dfds.zz*p->dfds.zz
							+ 0.5*(p->dfds.yz*p->dfds.yz + p->dfds.xz*p->dfds.xz + p->dfds.xy*p->dfds.xy);
			p->hp.minush = sqrt(p->hp.minush/1.5);
		}
		else
		{	// for use in alpha upate (it is -h or dalpha/lambda)
			p->hp.minush = p->dfds.xx*p->dfds.xx + p->dfds.yy*p->dfds.yy + p->dfds.zz*p->dfds.zz + 0.5*p->dfds.xy*p->dfds.xy;
			p->hp.minush = sqrt(p->hp.minush/1.5);
		}
	}
	
	else
	{	// negative root implies zero stress with roundoff error
		p->dfds.xx = p->dfds.yy = p->dfds.zz = p->dfds.yz = p->dfds.xz = p->dfds.xy = 0.;
		p->hp.minush=0.;
	}
}

// Solve numerically for lambda
// Ouptut is lambdak, final df, final alpha
double AnisoPlasticity::LRSolveForLambdaAP(MPMBase *mptr,int np,double ftrial,Matrix3 &strial,Matrix3 &stk,LRAnisoPlasticProperties *p) const
{
	int step;
	double lambda1,lambda2,f1,f2;
	
	// Step 1: stress in strial, alpha is previous alpha pick a second value
	LRGetDfCdf(strial,np,p);
#ifdef ITERATIVE_STRESS_UPDATE
	p->Cdf0 = p->Cdf;
	p->snorm = strial.DotProduct();
#endif
	lambda2 = ftrial/(p->dfCdf + GetDfAlphaDotH(mptr,np,&(p->hp)));
	double lambdaInitial = lambda2;
	f2 = LRGetFkFromLambdak(mptr,strial,stk,lambda2,np,p);
	
	// Step 2: bracket the solution
	
	/* Numerical Recipies Method */
	
	/*
	// pick second lambda that is lower
	lambda1 = 0.5*lambda2;
	f1 = LRGetFkFromLambdak(mptr,strial,stk,lambda1,np,p);	// pick second lambda that is lower
	
	// bracket the solution
	while(true)
	{	if(f1*f2<0.) break;
		if(fabs(f1)<fabs(f2))
		{	lambda1+=1.6*(lambda1-lambda2);
			if(lambda1>0.)
				f1=LRGetFkFromLambdak(mptr,strial,stk,lambda1,np,p);
			else
			{	f1=ftrial;
				lambda1=0.;
			}
		}
		else
		{	lambda2+=1.6*(lambda2-lambda1);
			f2=LRGetFkFromLambdak(mptr,strial,stk,lambda2,np,p);
		}
		
		// if fails to bracket in 50 tries, return with single step solution
		step++;
		if(step>50)
		{	if(warnings.Issue(warnNonconvergence,1)==GAVE_WARNING)
			{	mptr -> Describe();
				cout << "# Reverting to initial lambda = " << lambdaInitial << endl;
			}
			//GetDfCdf(strial,np,p);
			//UpdateStress(strial,&stk,lambdaInitial,np,p);
			//GetDfCdf(&stk,np,p);
			UpdateTrialAlpha(mptr,np,lambdaInitial,&(p->hp));
			return lambdaInitial;
		}
	}
	*/

	/* Custom method watching for parabolic function, but maybe not more */
	
	if(f2<0)
	{	// solution between 0 (ftrial>0) and lamda2 (f2<0)
		lambda1 = 0.;
		f1 = ftrial;
	}
	else
	{	// pick second lambda that is lower
		lambda1 = 0.5*lambda2;
		f1 = LRGetFkFromLambdak(mptr,strial,stk,lambda1,np,p);
		
		if(f1<0.)
		{	// solution between 0 (ftrial>0) and lamda1 (f1<0)
			lambda2 = 0.;
			f2 = ftrial;
		}
	
		else if(f1>f2)
		{	// solution is between lambda1 and lambda2 or after lambda2
			// increase lambda2 until f2<0 or it starts increasing
			//cout << "# increasing lambda2" << endl;
			double f2min = f2;
			step=1;
			while(true)
			{	// exit if f1 and f2 have different signs so solution between lambda1 and lambda2
				if(f1*f2<0.) break;
				//cout << "#    step " << step << " (" << lambda1 << "," << f1 << ") and (" << lambda2 << "," << f2 << ")" << endl;
				
				// increase lambda2
				lambda2 += 1.6*(lambda2-lambda1);			// increase lambda2
				f2 = LRGetFkFromLambdak(mptr,strial,stk,lambda2,np,p);
				if(f2>f2min) break;
				
				// keep going
				f2min = f2;
				step++;
				if(step>50)
				{	// if fails, use current lambda2, which will be the minimium f found so far
					if(warnings.Issue(warnNonconvergence,1)==GAVE_WARNING)
					{	mptr -> Describe();
						cout << "# Using minimum lambda while increasing = " << lambda2 << " with f = " << f2 << endl;
						LRPrintFk(mptr,strial,stk,lambda2,np,p,ftrial,lambdaInitial);
					}
					return lambda2;
				}
			}
		}
	
		else
		{	// solution is between lambda1 and lambda2 or before lambda1
			// decrease lambda1 until f1<0 or it starts increasing
			//cout << "# increasing lambda2" << endl;
			double f1min = f1;
			step = 1;
			while(true)
			{	// exit if f1 and f2 have different signs so solution between lambda1 and lambda2
				if(f1*f2<0.) break;
				//cout << "#    step " << step << " (" << lambda1 << "," << f1 << ") and (" << lambda2 << "," << f2 << ")" << endl;
				
				// reduce lambda1
				lambda1 += 1.6*(lambda1-lambda2);	
				if(lambda1>0.)
				{	f1=LRGetFkFromLambdak(mptr,strial,stk,lambda1,np,p);
				}
				else
				{	f1=ftrial;							// ftrial > 0
					lambda1=0.;
					break;
				}
				if(f1>f1min) break;
				
				// keep going
				f1min = f1;
				step++;
				if(step>50)
				{	// if fails, use current lambda1, which will be the minimium f found so far
					if(warnings.Issue(warnNonconvergence,1)==GAVE_WARNING)
					{	mptr -> Describe();
						cout << "# Using minimum lambda while decreasing = " << lambda1 << " with f = " << f1 << endl;
						LRPrintFk(mptr,strial,stk,lambda1,np,p,ftrial,lambdaInitial);
					}
					return lambda1;
				}
			}
		}
	
		// If both are positive then appears to be parabolic function and no guarantee that
		// the minimum is below zero, but look for it anyway
		if(f1*f2>0.)
		{	// find solution between lambda1 and lambda2 assume a parabola with local minimum
			// cout << "# Look for local minimum in (" << lambda1 << "," << f1 << ") and (" << lambda2 << "," << f2 << ")" << endl;
			step = 1;
			while(true)
			{	// move 1/4 space from the larger end
				if(f1>f2)
				{	lambda1 = lambda1 + 0.25*(lambda2-lambda1);
					f1=LRGetFkFromLambdak(mptr,strial,stk,lambda1,np,p);
				}
				else
				{	lambda2 = lambda2 - 0.25*(lambda2-lambda1);
					f2=LRGetFkFromLambdak(mptr,strial,stk,lambda2,np,p);
				}
				//cout << "#    step " << step << " (" << lambda1 << "," << f1 << ") and (" << lambda2 << "," << f2 << ")" << endl;
				
				// exit if done
				if(f1*f2<0.) break;
				step++;
				if(step>50)
				{	// Probably solution. Just return the current value which is likely the
					// minimum of the parabola
					if(warnings.Issue(warnNonconvergence,1)==GAVE_WARNING)
					{	mptr -> Describe();
						cout << "# Using minimum lambda in parabola = " << lambda2 << " with f = " << f2 << endl;
						LRPrintFk(mptr,strial,stk,lambda2,np,p,ftrial,lambdaInitial);
					}
					return lambda2;
				}
			}
		}
	}
	
	// The solution is bracketed, now use Newton's method with bracketing
	
	// order solutions so f1<0 and f2>0
	if(f1>0.)
	{	double ftemp=f1;
		double lambdatemp=lambda1;
		f1=f2;
		lambda1=lambda2;
		f2=ftemp;
		lambda2=lambdatemp;
	}
	
	// Step 3: Newton's method with bracketing
	
	// set up
	double lambdak=0.5*(lambda1+lambda2);				// initial guess
	double dxold=fabs(lambda2-lambda1);					// the step size before last
	double dx=dxold;									// the last step sie
	
	// initial guess
	double fk=LRGetFkFromLambdak(mptr,strial,stk,lambdak,np,p);
	double dfkdlam=-(p->dfCdf + GetDfAlphaDotH(mptr,np,&(p->hp)));
	
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
		fk=LRGetFkFromLambdak(mptr,strial,stk,lambdak,np,p);
		dfkdlam=-(p->dfCdf + GetDfAlphaDotH(mptr,np,&(p->hp)));
		
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
	
	if(step>20)
	{	if(warnings.Issue(warnNonconvergence,2)==GAVE_WARNING)
		{	mptr -> Describe();
			cout << "# Using final lambda = " << lambdak << " with f = " << fk << endl;
			LRPrintFk(mptr,strial,stk,lambdak,np,p,ftrial,lambdaInitial);
		}
	}
	
	// output df (from last lambdak), alpha (here with last lambdak), and lambdak (the return value)
	LRUpdateStress(strial,stk,lambdak,np,p);
	UpdateTrialAlpha(mptr,np,lambdak,&(p->hp));
	return lambdak;
}


// Update stress (assumes Cdfij calculated before)
void AnisoPlasticity::LRUpdateStress(Matrix3 &strial,Matrix3 &stk,double lambda,int np,LRAnisoPlasticProperties *p) const
{
#ifdef ITERATIVE_STRESS_UPDATE
	Matrix3 stkprev;
	
	// get stress relative to strial
	stk(0,0) = strial(0,0) - lambda*p->Cdf0.xx;
	stk(1,1) = strial(1,1) - lambda*p->Cdf0.yy;
	stk(2,2) = strial(2,2) - lambda*p->Cdf0.zz;
	stk(0,1) = strial(0,1) - lambda*p->Cdf0.xy;
	stk(1,0) = stk(0,1);
	if(np==THREED_MPM)
	{	stk(1,2) = strial(1,2) - lambda*p->Cdf0.yz;
		stk(2,1) = stk(1,2);
		stk(0,2) = strial(0,2) - lambda*p->Cdf0.xz;
		stk(2,0) = stk(0,2);
	}
	
	// iterate until close or 10 steps
	int nstep = 1;
	double diffnorm;
	while(nstep<10)
	{	// save last
		stkprev = stk;
		
		// get slopes at this first guess stress
		LRGetDfCdf(stkprev,np,p);
		
		// Now get stress using updated slopes
		stk(0,0) = strial(0,0) - lambda*p->Cdf.xx;
		stk(1,1) = strial(1,1) - lambda*p->Cdf.yy;
		stk(2,2) = strial(2,2) - lambda*p->Cdf.zz;
		stk(0,1) = strial(0,1) - lambda*p->Cdf.xy;
		stk(1,0) = stk(0,1);
		if(np==THREED_MPM)
		{	stk(1,2) = strial(1,2) - lambda*p->Cdf.yz;
			stk(2,1) = stk(1,2);
			stk(0,2) = strial(0,2) - lambda*p->Cdf.xz;
			stk(2,0) = stk(0,2);
		}
		
		// find difference normalized to strial norm
		Matrix3 diff = stk-stkprev;
		diffnorm = diff.DotProduct()/p->snorm;
		
		// if sqrt(diff)<1e-12 then return with current stk and slopes in plastic properties
		if(diffnorm < 1.e-24) return;
		
		// do another step
		nstep++;
	}
	
	// if use all steps print warning
	if(warnings.Issue(warnNonconvergence,3)==GAVE_WARNING)
	{	cout << "# Failed to update stress in 10 steps with final Delta = " << sqrt(diffnorm) << endl;
	}
#else
	// Now get stress using current slopes
	stk(0,0) = strial(0,0) - lambda*p->Cdf.xx;
	stk(1,1) = strial(1,1) - lambda*p->Cdf.yy;
	stk(2,2) = strial(2,2) - lambda*p->Cdf.zz;
	stk(0,1) = strial(0,1) - lambda*p->Cdf.xy;
	stk(1,0) = stk(0,1);
	if(np==THREED_MPM)
	{	stk(1,2) = strial(1,2) - lambda*p->Cdf.yz;
		stk(2,1) = stk(1,2);
		stk(0,2) = strial(0,2) - lambda*p->Cdf.xz;
		stk(2,0) = stk(0,2);
	}
#endif
}

// Given stress and lambda, find f (also find the stress and return dfCdf)
double AnisoPlasticity::LRGetFkFromLambdak(MPMBase *mptr,Matrix3 &strial,Matrix3 &stk,
										 double lambda,int np,LRAnisoPlasticProperties *p) const
{
	// Change stress to new value based on strial, lambda, and previous slope
	LRUpdateStress(strial,stk,lambda,np,p);
	
	// update alpha using new lambda and -h from most recent LRGetDfDsigma()
	UpdateTrialAlpha(mptr,np,lambda,&(p->hp));
	
	// update fk using new stress and alpha
	return GetMagnitudeHill(stk,np) - GetYield(&(p->hp));
}

// debugging - print function trying to solve
double AnisoPlasticity::LRPrintFk(MPMBase *mptr,Matrix3 &strial,Matrix3 &stk,
								  double lambda,int np,LRAnisoPlasticProperties *p,double ftrial,double lambdaInitial) const
{
	double lambdamax = log(2.*lambda);
	int npts = 100;
	cout << "0," << ftrial << endl;
	for(int i=1;i<npts;i++)
	{	double lambdak = exp((double)i*lambdamax/(double)npts);
		double fk = LRGetFkFromLambdak(mptr,strial,stk,lambdak,np,p);
		cout << lambdak << "," << fk << endl;
	}
	// resent
	LRUpdateStress(strial,stk,lambda,np,p);
	UpdateTrialAlpha(mptr,np,lambda,&(p->hp));
}

#pragma mark AnisoPlasticity::Methods (Small Rotation)

/* For 2D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, plastic strain, stresses, strain energy, 
		plastic energy, dissipated energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
   For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
   This is general analysis for anisotropic material based on hill criterion. Various
		hardening laws can implement the hardening term.
*/
void AnisoPlasticity::SRConstitutiveLaw2D(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// strain increments
	double dvxx = du(0,0);
	double dvyy = du(1,1);
	double dvxy = du(0,1);
	double dvyx = du(1,0);
	double dvzz = du(2,2);
	
	// properties
	AnisoPlasticProperties *p = (AnisoPlasticProperties *)properties;
	
    // Effective strain by deducting thermal strain
	// (note ep->alpha[1] and ep->beta[1] are reduced in plane strain, but CTE3 and CME3 are not)
	ElasticProperties r = p->ep;
	double erxx = r.alpha[1]*res->dT;
	double eryy = r.alpha[2]*res->dT;
	double erxy = r.alpha[3]*res->dT;
	double erzz = CTE3*res->dT;
	if(DiffusionTask::active)
	{	erxx += r.beta[1]*res->dC;
		eryy += r.beta[2]*res->dC;
		erxy += r.beta[3]*res->dC;
		erzz += CME3*res->dC;
	}
	
	// effective strains
    double dexx = dvxx-erxx;  
    double deyy = dvyy-eryy;
    double dgam = dvxy+dvyx;
    double dgxy = dgam-erxy;
	double dezz = dvzz-erzz;

	// rotational strain
	double dwrotxy=dvyx-dvxy;

    // Step 1: get trial stress by assuming elastic stress increment in current axes
	Tensor *sp = mptr->GetStressTensor();
    Tensor stk = *sp;
    stk.xx += r.C[1][1]*dexx+r.C[1][2]*deyy+r.C[1][3]*dgxy;
    stk.yy += r.C[1][2]*dexx+r.C[2][2]*deyy+r.C[2][3]*dgxy;
    stk.xy += r.C[1][3]*dexx+r.C[2][3]*deyy+r.C[3][3]*dgxy;
	stk.zz += r.C[4][1]*dexx + r.C[4][2]*deyy + r.C[4][3]*dgxy + r.C[4][4]*dezz;
    if(np==PLANE_STRAIN_MPM)
	{	stk.zz += r.C[4][1]*r.alpha[5]*erzz + r.C[4][2]*r.alpha[6]*erzz
					+ r.C[4][3]*r.alpha[7]*erzz;
	}
	
	// Step 2: 1st order rotate plastic strain
	// done if elastic, or add new plastic strain later if plastic
	Tensor *eplast=mptr->GetAltStrainTensor();
	double dwxy2 = dwrotxy*dwrotxy/4.;						// dwxy^2/4
	double shearD = 0.5*dwrotxy*(1.-0.5*(dvxx+dvyy));		// (1/2)*dwxy*(1-(dexx+deyy)/2)
	double diff = eplast->xx-eplast->yy;
	double dnorm = shearD*eplast->xy + dwxy2*diff;
	eplast->xx -= dnorm ;
	eplast->yy += dnorm ;
	double dshear = 2.*(shearD*diff - dwxy2*eplast->xy);
	eplast->xy += dshear ;
	
	// Step 3: Get rotation matrix to rotate stress from current axes to
	// material axes
	double s,c;
	mptr->Get2DSinCos(&s,&c);
	p->rzyx[0][0] = p->rzyx[1][1] = c*c;
	p->rzyx[0][1] = p->rzyx[1][0] = s*s;
	p->rzyx[0][5] = 2.*c*s;
	p->rzyx[1][5] = -p->rzyx[0][5];
	p->rzyx[5][0] = -c*s;
	p->rzyx[5][1] = -p->rzyx[5][0];
	p->rzyx[5][5] = p->rzyx[0][0] - p->rzyx[0][1];

    // Step 4: Calculate plastic potential f
	UpdateTrialAlpha(mptr,np,&(p->hp));
	Tensor srot;
	double sAstrial = GetMagnitudeRotatedHill(&stk,&srot,np,p);
	double ftrial = sAstrial - GetYield(&(p->hp));
	
	// Step 5: Done if elastic
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
		Elastic::SRConstitutiveLaw2D(mptr,du,delTime,np,&r,res);
		return; 
    }
    
	// Step 6: Solve for lambda
	double lambda = SolveForLambdaAP(mptr,np,ftrial,&stk,p);
    
    // Step 7: Plastic strain increments on particle
    double dexxp = lambda*p->dfds.xx;
    double deyyp = lambda*p->dfds.yy;
    double dgxyp = lambda*p->dfds.xy;
	double dezzp = lambda*p->dfds.zz;
    eplast->xx += dexxp;
    eplast->yy += deyyp;
    eplast->xy += dgxyp;
	eplast->zz += dezzp;
    
	// Step 8: Effective elastic strain increment = de - deres
	dexx -= dexxp;
	deyy -= deyyp;
	dgxy -= dgxyp;
	dezz -= dezzp;
	
	// Step 9: increment particle stresses
    Tensor st0=*sp;			// save previous stress for energy updates below
    double c1 = r.C[1][1]*dexx+r.C[1][2]*deyy+r.C[1][3]*dgxy;
    double c2 = r.C[1][2]*dexx+r.C[2][2]*deyy+r.C[2][3]*dgxy;
    double c3 = r.C[1][3]*dexx+r.C[2][3]*deyy+r.C[3][3]*dgxy;
	Hypo2DCalculations(mptr,dwrotxy,dvxx+dvyy,c1,c2,c3);
	
	// out of plane stress
    if(np==PLANE_STRAIN_MPM)
	{	sp->zz += r.C[4][1]*(dexx+r.alpha[5]*erzz) + r.C[4][2]*(deyy+r.alpha[6]*erzz)
					+ r.C[4][3]*(dgxy+r.alpha[7]*erzz) + r.C[4][4]*dezz;
	}

	// Step 10: increment energies
	// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
	double workEnergy = 0.5*((st0.xx+sp->xx)*dvxx
							   + (st0.yy+sp->yy)*dvyy
							   + (st0.xy+sp->xy)*dgam
							   + (st0.zz+sp->zz)*dvzz);

    // Plastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
	double dispEnergy = 0.5*((st0.xx+sp->xx)*dexxp
                        + (st0.yy+sp->yy)*deyyp
                        + (st0.xy+sp->xy)*dgxyp
						+ (st0.zz+sp->zz)*dezzp);

	// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
	double resEnergy = 0.5*((st0.xx+sp->xx)*erxx
							+ (st0.yy+sp->yy)*eryy
							+ (st0.xy+sp->xy)*erxy
							+ (st0.zz+sp->zz)*erzz);
	
	// add now
	mptr->AddWorkEnergyAndResidualEnergy(workEnergy + dispEnergy,resEnergy);
	
	// add dissipated energy to plastic energy to the particle
    mptr->AddPlastEnergy(dispEnergy);
	
    // heat energy is Cv(dT-dTq0) - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adiabatic heating occurs)
    IncrementHeatEnergy(mptr,res->dT,0.,dispEnergy);

	// Step 11: update internal variables
	UpdatePlasticInternal(mptr,np,&(p->hp));
}

/* For 3D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	duij are (gradient rates X time increment) to give deformation gradient change
	Assumes linear elastic, uses hypoelastic correction
*/
void AnisoPlasticity::SRConstitutiveLaw3D(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// strain increments
	double dvxx = du(0,0);
	double dvyy = du(1,1);
	double dvzz = du(2,2);
	
	// engineering shear strain icrements
    double dgamxy = du(0,1)+du(1,0);
    double dgamxz = du(0,2)+du(2,0);
    double dgamyz = du(1,2)+du(2,1);
	
	// rotational strain increments
	double dwrotxy = du(1,0)-du(0,1);
	double dwrotxz = du(2,0)-du(0,2);
	double dwrotyz = du(2,1)-du(1,2);
	
	// anisotropic properties
	AnisoPlasticProperties *p = (AnisoPlasticProperties *)properties;
	
    // Effective strain by deducting thermal strain
	ElasticProperties r = p->ep;
	double erxx=r.alpha[0]*res->dT;
	double eryy=r.alpha[1]*res->dT;
	double erzz=r.alpha[2]*res->dT;
	double eryz=r.alpha[3]*res->dT;
	double erxz=r.alpha[4]*res->dT;
	double erxy=r.alpha[5]*res->dT;
	if(DiffusionTask::active)
	{	erxx+=r.beta[0]*res->dC;
		eryy+=r.beta[1]*res->dC;
		erzz+=r.beta[2]*res->dC;
		eryz+=r.beta[3]*res->dC;
		erxz+=r.beta[4]*res->dC;
		erxy+=r.beta[5]*res->dC;
	}
    double dexx=dvxx-erxx;  
    double deyy=dvyy-eryy; 
	double dezz=dvzz-erzz;
    double dgyz=dgamyz-eryz;
    double dgxz=dgamxz-erxz;
    double dgxy=dgamxy-erxy;
	
    // Step 1: get trial stress by assuming elastic stress increment in current axes
	Tensor *sp=mptr->GetStressTensor();
    Tensor stk=*sp;
    stk.xx += r.C[0][0]*dexx+r.C[0][1]*deyy+r.C[0][2]*dezz+r.C[0][3]*dgyz+r.C[0][4]*dgxz+r.C[0][5]*dgxy;
    stk.yy += r.C[1][0]*dexx+r.C[1][1]*deyy+r.C[1][2]*dezz+r.C[1][3]*dgyz+r.C[1][4]*dgxz+r.C[1][5]*dgxy;
    stk.zz += r.C[2][0]*dexx+r.C[2][1]*deyy+r.C[2][2]*dezz+r.C[2][3]*dgyz+r.C[2][4]*dgxz+r.C[2][5]*dgxy;
    stk.yz += r.C[3][0]*dexx+r.C[3][1]*deyy+r.C[3][2]*dezz+r.C[3][3]*dgyz+r.C[3][4]*dgxz+r.C[3][5]*dgxy;
    stk.xz += r.C[4][0]*dexx+r.C[4][1]*deyy+r.C[4][2]*dezz+r.C[4][3]*dgyz+r.C[4][4]*dgxz+r.C[4][5]*dgxy;
    stk.xy += r.C[5][0]*dexx+r.C[5][1]*deyy+r.C[5][2]*dezz+r.C[5][3]*dgyz+r.C[5][4]*dgxz+r.C[5][5]*dgxy;
	
	// Step 2: 1st order rotate plastic strain
	// done if elastic, or add new plastic strain later if plastic
	Tensor *eplast=mptr->GetAltStrainTensor();
	double depxx = -0.5*(dwrotxy*eplast->xy + dwrotxz*eplast->xz);
	double depyy = 0.5*(dwrotxy*eplast->xy - dwrotyz*eplast->yz);
	double depzz = 0.5*(dwrotxz*eplast->xz + dwrotyz*eplast->yz);
	double depyz = dwrotyz*(eplast->yy-eplast->zz) + 0.5*(dwrotxz*eplast->xz + dwrotxy*eplast->xy);
	double depxz = dwrotxz*(eplast->xx-eplast->zz) + 0.5*(dwrotyz*eplast->yz - dwrotxy*eplast->xy);
	double depxy = dwrotxy*(eplast->xx-eplast->yy) - 0.5*(dwrotyz*eplast->yz + dwrotxz*eplast->xz);
	eplast->xx += depxx ;
	eplast->yy += depyy ;
	eplast->zz += depzz ;
	eplast->yz += depyz ;
	eplast->xz += depxz ;
	eplast->xy += depxy ;

	// Step 3: Get rotation matrix used for plastic potential
	// to rotate stress from current configuration to material axes
	// Gets Rsig^-1 = Re^T
	Matrix3 *Rtot = mptr->GetRtotPtr();
	Matrix3 RtotT = Rtot->Transpose();
	RtotT.GetRStress(p->rzyx);
	
    // Step 4: Calculate plastic potential f
	UpdateTrialAlpha(mptr,np,&(p->hp));
	Tensor srot;
	double sAstrial = GetMagnitudeRotatedHill(&stk,&srot,np,p);
	double ftrial = sAstrial - GetYield(&(p->hp));
	
	// Step 5: Done if elastic
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
		Elastic::SRConstitutiveLaw3D(mptr,du,delTime,np,&r,res);
		return; 
    }
    
	// Step 6: Solve for lambda
	double lambda = SolveForLambdaAP(mptr,np,ftrial,&stk,p);
    
    // Step 7: Plastic strain increments on particle
    double dexxp=lambda*p->dfds.xx;
    double deyyp=lambda*p->dfds.yy;
	double dezzp=lambda*p->dfds.zz;
    double dgyzp=lambda*p->dfds.yz;
    double dgxzp=lambda*p->dfds.xz;
    double dgxyp=lambda*p->dfds.xy;
    eplast->xx+=dexxp;
    eplast->yy+=deyyp;
	eplast->zz+=dezzp;
    eplast->yz+=dgyzp;
    eplast->xz+=dgxzp;
    eplast->xy+=dgxyp;
    
	// Step 8: Effective elastic strain increment = de - deres
	dexx -= dexxp;
	deyy -= deyyp;
	dezz -= dezzp;
	dgyz -= dgyzp;
	dgxz -= dgxzp;
	dgxy -= dgxyp;
	
	// Step 9: increment particle stresses
	Tensor st0=*sp;			// save previous stress for energy updates below
	double dsig[6];
    dsig[XX] = r.C[0][0]*dexx+r.C[0][1]*deyy+r.C[0][2]*dezz+r.C[0][3]*dgyz+r.C[0][4]*dgxz+r.C[0][5]*dgxy;
    dsig[YY] = r.C[1][0]*dexx+r.C[1][1]*deyy+r.C[1][2]*dezz+r.C[1][3]*dgyz+r.C[1][4]*dgxz+r.C[1][5]*dgxy;
    dsig[ZZ] = r.C[2][0]*dexx+r.C[2][1]*deyy+r.C[2][2]*dezz+r.C[2][3]*dgyz+r.C[2][4]*dgxz+r.C[2][5]*dgxy;
    dsig[YZ] = r.C[3][0]*dexx+r.C[3][1]*deyy+r.C[3][2]*dezz+r.C[3][3]*dgyz+r.C[3][4]*dgxz+r.C[3][5]*dgxy;
    dsig[XZ] = r.C[4][0]*dexx+r.C[4][1]*deyy+r.C[4][2]*dezz+r.C[4][3]*dgyz+r.C[4][4]*dgxz+r.C[4][5]*dgxy;
    dsig[XY] = r.C[5][0]*dexx+r.C[5][1]*deyy+r.C[5][2]*dezz+r.C[5][3]*dgyz+r.C[5][4]*dgxz+r.C[5][5]*dgxy;
	Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,dsig);

	// Step 10: Increment energies
    // Elastic work increment per unit mass (dU/(rho0 V0)) (nJ/g)
	double workEnergy = 0.5*((st0.xx+sp->xx)*dvxx
							   + (st0.yy+sp->yy)*dvyy
							   + (st0.zz+sp->zz)*dvzz
							   + (st0.yz+sp->yz)*dgamyz
							   + (st0.xz+sp->xz)*dgamxz
							   + (st0.xy+sp->xy)*dgamxy);
	
    // Plastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
	double dispEnergy=0.5*(0.5*((st0.xx+sp->xx)*dexxp
								+ (st0.yy+sp->yy)*deyyp
								+ (st0.zz+sp->zz)*dezzp
								+ (st0.yz+sp->yz)*dgyzp
								+ (st0.xz+sp->xz)*dgxzp
								+ (st0.xy+sp->xy)*dgxyp));
	
	// Elastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
	double resEnergy = 0.5*((st0.xx+sp->xx)*erxx
							+ (st0.yy+sp->yy)*eryy
							+ (st0.zz+sp->zz)*erzz
							+ (st0.yz+sp->yz)*eryz
							+ (st0.xz+sp->xz)*erxz
							+ (st0.xy+sp->xy)*erxy);
	
	// add now
	mptr->AddWorkEnergyAndResidualEnergy(workEnergy + dispEnergy,resEnergy);
	
	// add dissipated energy to plastic energy to the particle
    mptr->AddPlastEnergy(dispEnergy);
	
    // heat energy is Cv(dT-dTq0) - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adiabatic heating occurs)
    IncrementHeatEnergy(mptr,res->dT,0.,dispEnergy);
	
	// Step 11: update internal variables
	UpdatePlasticInternal(mptr,np,&(p->hp));
}

// Get sqrt(s As) and also return rotated stresses in case caller needs them
double AnisoPlasticity::GetMagnitudeRotatedHill(Tensor *st0,Tensor *srot,int np,AnisoPlasticProperties *p) const
{
	double sAs=0.;
	
	if(np==THREED_MPM)
	{	// rotation from analysis to material axes srot = Rsig^-1 s = Re^T s, where Rsig^-1 is in rzyx array
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
	{	// rotation from analysis to material axes using Rsigma^(-1).st0
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
// Note that engineering plastis strain is 2 f.xy, etc., so need to multiple last three columns of C by 2
void AnisoPlasticity::GetDfCdf(Tensor *stk,int np,AnisoPlasticProperties *p) const
{
	// get C df and df C df, which need df/dsig (Function of yield criteria, and normally only current stress in stk)
	GetDfDsigma(stk,np,p);
	ElasticProperties r = p->ep;
	if(np==THREED_MPM)
	{	p->Cdf.xx = r.C[0][0]*p->dfds.xx+r.C[0][1]*p->dfds.yy+r.C[0][2]*p->dfds.zz
						+r.C[0][3]*p->dfds.yz+r.C[0][4]*p->dfds.xz+r.C[0][5]*p->dfds.xy;
		p->Cdf.yy = r.C[1][0]*p->dfds.xx+r.C[1][1]*p->dfds.yy+r.C[1][2]*p->dfds.zz
						+r.C[1][3]*p->dfds.yz+r.C[1][4]*p->dfds.xz+r.C[1][5]*p->dfds.xy;
		p->Cdf.zz = r.C[2][0]*p->dfds.xx+r.C[2][1]*p->dfds.yy+r.C[2][2]*p->dfds.zz
						+r.C[2][3]*p->dfds.yz+r.C[2][4]*p->dfds.xz+r.C[2][5]*p->dfds.xy;
		p->Cdf.yz = r.C[3][0]*p->dfds.xx+r.C[3][1]*p->dfds.yy+r.C[3][2]*p->dfds.zz
						+r.C[3][3]*p->dfds.yz+r.C[3][4]*p->dfds.xz+r.C[3][5]*p->dfds.xy;
		p->Cdf.xz = r.C[4][0]*p->dfds.xx+r.C[4][1]*p->dfds.yy+r.C[4][2]*p->dfds.zz
						+r.C[4][3]*p->dfds.yz+r.C[4][4]*p->dfds.xz+r.C[4][5]*p->dfds.xy;
		p->Cdf.xy = r.C[5][0]*p->dfds.xx+r.C[5][1]*p->dfds.yy+r.C[5][2]*p->dfds.zz
						+r.C[5][3]*p->dfds.yz+r.C[5][4]*p->dfds.xz+r.C[5][5]*p->dfds.xy;
		p->dfCdf = p->dfds.xx*p->Cdf.xx + p->dfds.yy*p->Cdf.yy + p->dfds.zz*p->Cdf.zz
						+ p->dfds.yz*p->Cdf.yz + p->dfds.xz*p->Cdf.xz + p->dfds.xy*p->Cdf.xy;
	}
	else
	{	p->Cdf.xx = r.C[1][1]*p->dfds.xx + r.C[1][2]*p->dfds.yy + r.C[1][3]*p->dfds.xy;
		p->Cdf.yy = r.C[1][2]*p->dfds.xx + r.C[2][2]*p->dfds.yy + r.C[2][3]*p->dfds.xy;
		p->Cdf.xy = r.C[1][3]*p->dfds.xx + r.C[2][3]*p->dfds.yy + r.C[3][3]*p->dfds.xy;
		p->Cdf.zz = r.C[4][1]*p->dfds.xx + r.C[4][2]*p->dfds.yy + r.C[4][3]*p->dfds.xy + r.C[4][4]*p->dfds.zz;
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
			
			// rotate to analysis coordinates df = Rsig^-1^T dfrot = Re ftot, which is transpose of rzyx array
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
			p->hp.minush = dfdsxxrot*dfdsxxrot + dfdsyyrot*dfdsyyrot + dfdszzrot*dfdszzrot
							+ 0.5*(dfdtyzrot*dfdtyzrot + dfdtxzrot*dfdtxzrot + dfdtxyrot*dfdtxyrot);
			p->hp.minush = sqrt(p->hp.minush/1.5);
		}
		else
		{	// rotate to analysis coordinates df = R^T dfrot
			p->dfds.xx = dfdsxxrot*p->rzyx[0][0] + dfdsyyrot*p->rzyx[1][0] + dfdtxyrot*p->rzyx[5][0];
			p->dfds.yy = dfdsxxrot*p->rzyx[0][1] + dfdsyyrot*p->rzyx[1][1] + dfdtxyrot*p->rzyx[5][1];
			p->dfds.xy = (dfdsxxrot - dfdsyyrot)*p->rzyx[0][5] + dfdtxyrot*p->rzyx[5][5];
			p->dfds.zz = dfdszzrot;
			
			// for use in alpha upate
			p->hp.minush = dfdsxxrot*dfdsxxrot + dfdsyyrot*dfdsyyrot + dfdszzrot*dfdszzrot + 0.5*dfdtxyrot*dfdtxyrot;
			p->hp.minush = sqrt(p->hp.minush/1.5);
		}
	}
	
	else
	{	// negative root implies zero stress with roundoff error
		p->dfds.xx = p->dfds.yy = p->dfds.zz = p->dfds.yz = p->dfds.xz = p->dfds.xy = 0.;
		p->hp.minush=0.;
	}
}

// Solve numerically for lambda
// Ouptut is lambdak, final df, final alpha
double AnisoPlasticity::SolveForLambdaAP(MPMBase *mptr,int np,double ftrial,Tensor *strial,AnisoPlasticProperties *p) const
{
	int step;
	double lambda1,lambda2,f1,f2;
	
	// Step 1: stress in strial, alpha is previous alpha pick a second value
	GetDfCdf(strial,np,p);
#ifdef ITERATIVE_STRESS_UPDATE
	p->Cdf0 = p->Cdf;
	p->snorm = strial->xx*strial->xx+strial->yy*strial->yy+strial->zz*strial->zz+2.*strial->xy*strial->xy;
	if(np==THREED_MPM) p->snorm += 2.*(strial->xz*strial->xz+strial->yz*strial->yz);
#endif
	lambda2 = ftrial/(p->dfCdf + GetDfAlphaDotH(mptr,np,&(p->hp)));
	double lambdaInitial=lambda2;
	Tensor stk;
	f2 = GetFkFromLambdak(mptr,strial,&stk,lambda2,np,p);
	
	// Step 2: bracket the solution
	
	/* Numerical Recipes Method */
	
	/*
	// find df/dsigma and -h at this initial guessed stress rather than at strial as above
	//Tensor stk;
	//UpdateStress(strial,&stk,lambda2,np,p);
	//GetDfCdf(&stk,np,p);
	//UpdateStress(strial,&stk,lambda2,np,p);			// second pass might help
	//GetDfCdf(&stk,np,p);

	// Find f using new slopes at lambda1
	//UpdateTrialAlpha(mptr,np,lambda2,&(p->hp));
	//Tensor srot;
	//f2 = GetMagnitudeRotatedHill(&stk,&srot,np,p) - GetYield(&(p->hp));
	
	// pick second lambda that is lower
	lambda1 = 0.5*lambda2;
	f1 = GetFkFromLambdak(mptr,strial,&stk,lambda1,np,p);
	
	// bracket the solution
	step=1;
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
		{	if(warnings.Issue(warnNonconvergence,1)==GAVE_WARNING)
			{	mptr -> Describe();
				cout << "# Reverting to initial lambda = " << lambdaInitial << endl;
			}
			//GetDfCdf(strial,np,p);
			//UpdateStress(strial,&stk,lambdaInitial,np,p);
			//GetDfCdf(&stk,np,p);
			UpdateTrialAlpha(mptr,np,lambdaInitial,&(p->hp));
			cout << " # Failed to bracket in 50 steps" << endl;
			return lambdaInitial;
		}
	}
	*/
	
	/* Custom Method */
	
	if(f2<0)
	{	// solution between 0 (ftrial>0) and lamda2 (f2<0)
		lambda1 = 0.;
		f1 = ftrial;
	}
	else
	{	// pick second lambda that is lower
		lambda1 = 0.5*lambda2;
		f1 = GetFkFromLambdak(mptr,strial,&stk,lambda1,np,p);
		
		if(f1<0.)
		{	// solution between 0 (ftrial>0) and lamda1 (f1<0)
			lambda2 = 0.;
			f2 = ftrial;
		}
		
		else if(f1>f2)
		{	// solution is between lambda1 and lambda2 or after lambda2
			// increase lambda2 until f2<0 or it starts increasing
			//cout << "# increasing lambda2" << endl;
			double f2min = f2;
			step=1;
			while(true)
			{	// exit if f1 and f2 have different signs so solution between lambda1 and lambda2
				if(f1*f2<0.) break;
				//cout << "#    step " << step << " (" << lambda1 << "," << f1 << ") and (" << lambda2 << "," << f2 << ")" << endl;
				
				// increase lambda2
				lambda2 += 1.6*(lambda2-lambda1);			// increase lambda2
				f2 = GetFkFromLambdak(mptr,strial,&stk,lambda2,np,p);
				if(f2>f2min) break;
				
				// keep going
				f2min = f2;
				step++;
				if(step>50)
				{	// if fails, use current lambda2, which will be the minimium f found so far
					if(warnings.Issue(warnNonconvergence,1)==GAVE_WARNING)
					{	mptr -> Describe();
						cout << "# Using minimum lambda while increasing = " << lambda2 << " with f = " << f2 << endl;
						SRPrintFk(mptr,strial,&stk,lambda2,np,p,ftrial,lambdaInitial);
					}
					return lambda2;
				}
			}
		}
		
		else
		{	// solution is between lambda1 and lambda2 or before lambda1
			// decrease lambda1 until f1<0 or it starts increasing
			//cout << "# increasing lambda2" << endl;
			double f1min = f1;
			step = 1;
			while(true)
			{	// exit if f1 and f2 have different signs so solution between lambda1 and lambda2
				if(f1*f2<0.) break;
				//cout << "#    step " << step << " (" << lambda1 << "," << f1 << ") and (" << lambda2 << "," << f2 << ")" << endl;
				
				// reduce lambda1
				lambda1 += 1.6*(lambda1-lambda2);
				if(lambda1>0.)
				{	f1=GetFkFromLambdak(mptr,strial,&stk,lambda1,np,p);
				}
				else
				{	f1=ftrial;							// ftrial > 0
					lambda1=0.;
					break;
				}
				if(f1>f1min) break;
				
				// keep going
				f1min = f1;
				step++;
				if(step>50)
				{	// if fails, use current lambda1, which will be the minimium f found so far
					if(warnings.Issue(warnNonconvergence,1)==GAVE_WARNING)
					{	mptr -> Describe();
						cout << "# Using minimum lambda while decreasing = " << lambda1 << " with f = " << f1 << endl;
						SRPrintFk(mptr,strial,&stk,lambda1,np,p,ftrial,lambdaInitial);
					}
					return lambda1;
				}
			}
		}
		
		// If both are positive then appears to be parabolic function and no guarantee that
		// the minimum is below zero, but look for it anyway
		if(f1*f2>0.)
		{	// find solution between lambda1 and lambda2 assume a parabola with local minimum
			// cout << "# Look for local minimum in (" << lambda1 << "," << f1 << ") and (" << lambda2 << "," << f2 << ")" << endl;
			step = 1;
			while(true)
			{	// move 1/4 space from the larger end
				if(f1>f2)
				{	lambda1 = lambda1 + 0.25*(lambda2-lambda1);
					f1=GetFkFromLambdak(mptr,strial,&stk,lambda1,np,p);
				}
				else
				{	lambda2 = lambda2 - 0.25*(lambda2-lambda1);
					f2=GetFkFromLambdak(mptr,strial,&stk,lambda2,np,p);
				}
				//cout << "#    step " << step << " (" << lambda1 << "," << f1 << ") and (" << lambda2 << "," << f2 << ")" << endl;
				
				// exit if done
				if(f1*f2<0.) break;
				step++;
				if(step>50)
				{	// Probably solution. Just return the current value which is likely the
					// minimum of the parabola
					if(warnings.Issue(warnNonconvergence,1)==GAVE_WARNING)
					{	mptr -> Describe();
						cout << "# Using minimum lambda in parabola = " << lambda2 << " with f = " << f2 << endl;
						SRPrintFk(mptr,strial,&stk,lambda2,np,p,ftrial,lambdaInitial);
					}
					return lambda2;
				}
			}
		}
	}
	
	// The solution is bracketed, now use Newton's method with bracketing

	// order solutions so f1<0
	if(f1>0.)
	{	double ftemp=f1;
		double lambdatemp=lambda1;
		f1=f2;
		lambda1=lambda2;
		f2=ftemp;
		lambda2=lambdatemp;
	}
	
	// Step 3: Newton's method with bracketing
	
	// set up
	double lambdak=0.5*(lambda1+lambda2);				// initial guess
	double dxold=fabs(lambda2-lambda1);					// the step size before last
	double dx=dxold;									// the last step sie
	
	// initial guess
	double fk=GetFkFromLambdak(mptr,strial,&stk,lambdak,np,p);
	double dfkdlam=-(p->dfCdf + GetDfAlphaDotH(mptr,np,&(p->hp)));
	
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
		dfkdlam=-(p->dfCdf + GetDfAlphaDotH(mptr,np,&(p->hp)));
		
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
	
	if(step>20)
	{	if(warnings.Issue(warnNonconvergence,2)==GAVE_WARNING)
		{	mptr -> Describe();
			cout << "# Using final lambda = " << lambdak << " with f = " << fk << endl;
			SRPrintFk(mptr,strial,&stk,lambdak,np,p,ftrial,lambdaInitial);
		}
	}
	
	// output df (from last stress), alpha (here with latest lambda), and lamda (the return value)
	UpdateStress(strial,&stk,lambdak,np,p);
	UpdateTrialAlpha(mptr,np,lambdak,&(p->hp));
	return lambdak;
}

// Update stress (assumes Cdfij calculated before)
void AnisoPlasticity::UpdateStress(Tensor *strial,Tensor *stk,double lambda,int np,AnisoPlasticProperties *p) const
{
#ifdef ITERATIVE_STRESS_UPDATE
	Tensor stkprev;
	
	// get stress relative to strial
	stk->xx = strial->xx - lambda*p->Cdf0.xx;
	stk->yy = strial->yy - lambda*p->Cdf0.yy;
	stk->xy = strial->xy - lambda*p->Cdf0.xy;
	stk->zz = strial->zz - lambda*p->Cdf0.zz;
	if(np==THREED_MPM)
	{	stk->yz = strial->yz - lambda*p->Cdf0.yz;
		stk->xz = strial->xz - lambda*p->Cdf0.xz;
	}
	
	// iterate until close or 10 steps
	int nstep = 1;
	double diffnorm;
	while(nstep<10)
	{	// save last
		stkprev = *stk;
		
		// get slopes at current stress
		GetDfCdf(stk,np,p);
		
		// Now get stress using updated slopes
		stk->xx = strial->xx - lambda*p->Cdf.xx;
		stk->yy = strial->yy - lambda*p->Cdf.yy;
		stk->xy = strial->xy - lambda*p->Cdf.xy;
		stk->zz = strial->zz - lambda*p->Cdf.zz;
		if(np==THREED_MPM)
		{	stk->yz = strial->yz - lambda*p->Cdf.yz;
			stk->xz = strial->xz - lambda*p->Cdf.xz;
		}
		
		// find difference normalized to strial norm
		diffnorm = (stk->xx-stkprev.xx)*(stk->xx-stkprev.xx)+(stk->yy-stkprev.yy)*(stk->yy-stkprev.yy)+
						+(stk->zz-stkprev.zz)*(stk->zz-stkprev.zz)+2.*(stk->xy-stkprev.xy)*(stk->xy-stkprev.xy);
		if(np==THREED_MPM)
		{	diffnorm += 2.*((stk->xz-stkprev.xz)*(stk->xz-stkprev.xz)+(stk->yz-stkprev.yz)*(stk->yz-stkprev.yz));
		}
		
		// if sqrt(diff)<1e-12 then return wuth current stk and slopes in plastic properties
		if(diffnorm/p->snorm < 1.e-24) return;
		
		// do another step
		nstep++;
	}
	
	// if use all steps print warning
	if(warnings.Issue(warnNonconvergence,3)==GAVE_WARNING)
	{	cout << "# Failed to update stress in 10 steps with final Delta = " << sqrt(diffnorm/p->snorm) << endl;
	}
#else
	// Now get stress using updated slopes
	stk->xx = strial->xx - lambda*p->Cdf.xx;
	stk->yy = strial->yy - lambda*p->Cdf.yy;
	stk->xy = strial->xy - lambda*p->Cdf.xy;
	stk->zz = strial->zz - lambda*p->Cdf.zz;
	if(np==THREED_MPM)
	{	stk->yz = strial->yz - lambda*p->Cdf.yz;
		stk->xz = strial->xz - lambda*p->Cdf.xz;
	}
#endif
}

// Given stress and lambda, find f (also find the stress and return dfCdf)
double AnisoPlasticity::GetFkFromLambdak(MPMBase *mptr,Tensor *strial,Tensor *stk,
										 double lambda,int np,AnisoPlasticProperties *p) const
{
	// Change stress to new value based on strial, lambda, and previous slope
	UpdateStress(strial,stk,lambda,np,p);
	
	// update alpha using new lambda and -h from most recent GetDfDsigma()
	UpdateTrialAlpha(mptr,np,lambda,&(p->hp));
	
	// update fk using new stress and alpha
	Tensor srot;
	return GetMagnitudeRotatedHill(stk,&srot,np,p) - GetYield(&(p->hp));
}

// debugging - print function trying to solve
double AnisoPlasticity::SRPrintFk(MPMBase *mptr,Tensor *strial,Tensor *stk,
								  double lambda,int np,AnisoPlasticProperties *p,double ftrial,double lambdaInitial) const
{
	double lambdamax = log(2.*lambda);
	int npts = 100;
	cout << "0," << ftrial << endl;
	for(int i=1;i<npts;i++)
	{	double lambdak = exp((double)i*lambdamax/(double)npts);
		double fk = GetFkFromLambdak(mptr,strial,stk,lambdak,np,p);
		cout << lambdak << "," << fk << endl;
	}
	UpdateStress(strial,stk,lambda,np,p);
	UpdateTrialAlpha(mptr,np,lambda,&(p->hp));
}

#pragma mark AnisoPlasticity::Accessors

// store plastic strain in alt strain
int AnisoPlasticity::AltStrainContains(void) const
{	return ENG_BIOT_PLASTIC_STRAIN;
}
