/********************************************************************************
	IsoSoftening.cpp
	nairn-mpm-fea

	Created by John Nairn on June 26, 2015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/IsoSoftening.hpp"
#include "Materials/IsotropicMat.hpp"
#include "Materials/Elastic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Materials/FailureSurface.hpp"
#include "Materials/LinearSoftening.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/UnitsController.hpp"
#include "System/ArchiveData.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Boundary_Conditions/InitialCondition.hpp"

// globals (0-0.5 - not damaged, 0.5-1.5 - damage propagation, >1.5 - post failure)
double predamageState = 0.5;
double damageState = 1.5;
extern double mtime;

#pragma mark IsoSoftening::Constructors and Destructors

// Constructor
// throws std::bad_alloc
IsoSoftening::IsoSoftening(char *matName,int matID) : IsotropicMat(matName,matID)
{
	initiationLaw = new FailureSurface(this);
	softeningModeI = new LinearSoftening();
	softeningModeII = new LinearSoftening();
	shearFailureSurface = RECTANGULAR_SURFACE;
	softenCV = 0.;
    softenCVMode = VARY_STRENGTH;
}

#pragma mark IsoSoftening::Initialization

// Read material properties
char *IsoSoftening::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    // look for different damage initiaion law
    if(strcmp(xName,"Initiation")==0)
    {	input = INITIATION_LAW_SELECTION;
        return (char *)this;
    }
	
    // look for different softeningI law
    else if(strcmp(xName,"SofteningI")==0)
    {	input = SOFTI_LAW_SELECTION;
        return (char *)this;
    }
	
    // look for different softeningII law
    else if(strcmp(xName,"SofteningII")==0)
    {	input = SOFTII_LAW_SELECTION;
        return (char *)this;
    }
	
	// type of failure surface in 3D shear
	else if(strcmp(xName,"shearFailureSurface")==0)
	{	input=INT_NUM;
		return((char *)&shearFailureSurface);
	}
	
	// softening coefficient of variation
	else if(strcmp(xName,"coefVariation")==0)
	{	input=DOUBLE_NUM;
		return((char *)&softenCV);
	}
	
    // coefficient of variation for strength (1), toughness (2), or both (3)
    else if(strcmp(xName,"coefVariationMode")==0)
    {	input=INT_NUM;
        return((char *)&softenCVMode);
    }
    
    // check initiation law
    char *ptr = initiationLaw->InputInitationProperty(xName,input,gScaling);
    if(ptr != NULL) return ptr;
	
	// if possible check softening laws
	if(strlen(xName)>3)
	{	if(xName[0]=='I' && xName[1]=='-')
		{	ptr = softeningModeI->InputSofteningProperty(&xName[2],input,gScaling);
			if(ptr != NULL) return ptr;
		}
		else if(xName[0]=='I' && xName[1]=='I' && xName[2]=='-')
		{	ptr = softeningModeII->InputSofteningProperty(&xName[3],input,gScaling);
			if(ptr != NULL) return ptr;
		}
	}
    
    // otherwise get material properties
    return IsotropicMat::InputMaterialProperty(xName,input,gScaling);
}

// Allows selected initiation laws
bool IsoSoftening::AcceptInitiationLaw(FailureSurface *iLaw,int lawID)
{
	// only allows max principle stress
	if(lawID != MAXPRINCIPALSTRESS) return false;
	
	delete initiationLaw;
    initiationLaw = iLaw;
    return true;
}

// Allows mode I and mode II softening laws
bool IsoSoftening::AcceptSofteningLaw(SofteningLaw *sLaw,int lawID,int mode)
{
	if(mode==SOFTI_LAW_SELECTION)
	{	delete softeningModeI;
		softeningModeI = sLaw;
	}
	else if(mode==SOFTII_LAW_SELECTION)
	{	delete softeningModeII;
		softeningModeII = sLaw;
	}
	else
		return false;
    return true;
}

// verify settings and some initial calculations
const char *IsoSoftening::VerifyAndLoadProperties(int np)
{
	if(!useLargeRotation)
		return "IsoSoftening material requires activation of large rotation option";
	
	//if(np==PLANE_STRESS_MPM)
	//	return "IsoSoftening materials cannot yet be used in plane stress calculations";
	
	// call initiation law that is used
    const char *ptr = initiationLaw->VerifyAndLoadProperties(np);
    if(ptr != NULL) return ptr;
	
	// check in superclass (along with its initialization)
	const char *err = IsotropicMat::VerifyAndLoadProperties(np);
	if(err!=NULL) return err;
	
	// constant properties
	// (note both stress and stiffness are divided by rho so strains are absolute)
	if(np==THREED_MPM)
	{	en0 = initiationLaw->sigmaI()/pr.C[0][0];
		gs0 = initiationLaw->sigmaII()/pr.C[3][3];
	}
	else
	{	en0 = initiationLaw->sigmaI()/pr.C[1][1];
		gs0 = initiationLaw->sigmaII()/pr.C[3][3];
	}
    
    if(softenCVMode<VARY_STRENGTH || softenCVMode>VARY_STRENGTH_AND_TOUGHNESS)
        return "Invalid option selected for which softening properties to vary";
	
	if(shearFailureSurface!=RECTANGULAR_SURFACE)
		return "Only the rectangular failure surface is currently implemented";
    
    return NULL;
}

// print mechanical properties to the results
void IsoSoftening::PrintMechanicalProperties(void) const
{
    IsotropicMat::PrintMechanicalProperties();
	cout << "Failure Surface: ";
	initiationLaw->PrintInitiationProperties();
	cout << "Mode I Softening: ";
	softeningModeI->PrintSofteningProperties(rho*initiationLaw->sigmaI());
	cout << "Mode II Softening: ";
	softeningModeII->PrintSofteningProperties(rho*initiationLaw->sigmaII());
    cout << "Coefficient of Variation = " << softenCV << endl;
    if(softenCV>0.)
    {   if(softenCVMode==VARY_STRENGTH_AND_TOUGHNESS)
            cout << "    Vary strength and toughness" << endl;
        else if(softenCVMode==VARY_TOUGHNESS)
            cout << "    Vary toughness" << endl;
        else
            cout << "    Vary strength" << endl;
    }
	cout << "3D Shear Stress Softening: ";
	if(shearFailureSurface == ELLIPTICAL_SURFACE)
		cout << "anisotropic with elliptical failure surface" << endl;
	else
		cout << "anisotropic with rectangular failure surface" << endl;
}

#pragma mark IsoSoftening::History Data Methods

// Create history variables needed for softening behavior
char *IsoSoftening::InitHistoryData(char *pchr,MPMBase *mptr)
{
	// Validate this law (in use) is stable on current grid
	//	This validation uses average values. Large stochastic variations might
	//	cause some particle to be unstable. Enhancement could be to validate
	//	using minimum Gc and maximum strength instead
	// redundant when all particle same size, but general if not
	double delx = mptr->GetMinParticleLength();
	double scale = 1./(rho*initiationLaw->sigmaI()*delx);
	double maxSlope = softeningModeI->GetMaxSlope(scale);
	if (maxSlope > 1. / en0)
	{	cout << "Need to increase GIc, decrease sigmac^2, or decrease particle size "
				<< maxSlope*en0 << " fold for spatial stability." << endl;
		throw CommonException("Normal softening law needs smaller grid cells","IsoSoftening::InitHistoryData");
	}
	scale = 1./(rho*initiationLaw->sigmaII()*delx);
	maxSlope = softeningModeII->GetMaxSlope(scale);
	if(maxSlope > 1./gs0)
	{	cout << "Need to increase GIIc, decrease tauc^2, or decrease particle size "
				<< maxSlope*gs0 << " fold for spatial stability." << endl;
		throw CommonException("Shear softening law needs smaller grid cells","IsoSoftening::InitHistoryData");
	}

	double *p = CreateAndZeroDoubles(pchr,SOFT_NUMBER_HISTORY);
	
	// set damage state to 0.1 instead of zero for help in plotting on a grid
	p[SOFT_DAMAGE_STATE] = 0.1;
	
	// set relative strength
	p[RELATIVE_STRENGTH] = 1.;
	p[RELATIVE_TOUGHNESS] = 1.;
	if(softenCV>0.)
	{	double fract = (double)(rand() % 998 + 1)/1000.;		// .001 to .999, max +/- 3.09 std dev
		if(softenCVMode&VARY_STRENGTH)
			p[RELATIVE_STRENGTH] = fmax(1. + softenCV*NormalCDFInverse(fract),0.1);
		if(softenCVMode&VARY_TOUGHNESS)
			p[RELATIVE_TOUGHNESS] = fmax(1. + softenCV*NormalCDFInverse(fract),0.1);
	}
	
    return (char *)p;
}

// Number of history variables
int IsoSoftening::NumberOfHistoryDoubles(void) const { return SOFT_NUMBER_HISTORY; }

// Initialize damage when requested
void IsoSoftening::SetInitialConditions(InitialCondition *ic,MPMBase *mptr,bool is3D)
{
    double *soft = GetSoftHistoryPtr(mptr);

	if (is3D) {
		Vector dnorm = ic->GetDamageNormal();
		soft[NORMALDIR1] = atan2(dnorm.y,dnorm.x); // Euler angle alpha (should be (y,x))
		soft[NORMALDIR2] = asin(dnorm.z);  //Euler angle beta using asin()Da
		soft[NORMALDIR3] = 0.0; // only 2 angles needed for a normal vector

		// Get failure mode (if not fully damaged) and damage parameters
		double dmode;
		Vector dparams = ic->GetDamageParams(dmode);
		soft[DAMAGENORMAL] = dparams.x;
		soft[DAMAGESHEAR] = dparams.y;
		soft[DAMAGESHEAR2] = dparams.z;

		// Find Ac/Vp
		soft[GCSCALING] = GetAcOverVp(THREED_MPM, mptr, &dnorm) / rho;
		
		// Any damage variable >=1 means failed, otherwise calucate corresponding cracking strains
		// and set provide failure mode
		if(dparams.x>=1. || dparams.y>=1. || dparams.z>=1.)
		{	soft[DAMAGENORMAL] = 1.;
			soft[DAMAGESHEAR] = 1.;
			soft[DAMAGESHEAR2] = 1.;
			soft[SOFT_DAMAGE_STATE] = 4.;
		}
		else
		{	soft[SOFT_DAMAGE_STATE] = dmode;
			double relToughness = soft[RELATIVE_TOUGHNESS];
			double relStrength = soft[RELATIVE_STRENGTH];
			double scale = relToughness*soft[GCSCALING]/(relStrength*initiationLaw->sigmaI());
			double relen0 = relStrength*en0;
			soft[DELTANORMAL] = softeningModeI->GetDeltaFromDamage(dparams.x,scale,relen0);
			scale = soft[GCSCALING]/initiationLaw->sigmaII();
			double relgs0 = relStrength*gs0;
			soft[DELTASHEAR] = softeningModeII->GetDeltaFromDamage(dparams.y,scale,relgs0);
			soft[DELTASHEAR2] = softeningModeII->GetDeltaFromDamage(dparams.z,scale,relgs0);
		}
	}
	else {
		// 2D damage
		Vector dnorm = ic->GetDamageNormal();
		soft[NORMALDIR1] = dnorm.x;
		soft[NORMALDIR2] = dnorm.y;
		dnorm.z = 0.;

		// Get failure mode (if not fully damaged) and damage parameters
		double dmode;
		Vector dparams = ic->GetDamageParams(dmode);
		soft[DAMAGENORMAL] = dparams.x;
		soft[DAMAGESHEAR] = dparams.y;

		// Find Ac/Vp (first parameter anything except THREED_MPM)
		soft[GCSCALING] = GetAcOverVp(PLANE_STRAIN_MPM, mptr, &dnorm) / rho;
		
		// Any damage variable >=1 means failed, otherwise calucate corresponding cracking strains
		// and set provide failure mode
		if(dparams.x>=1. || dparams.y>=1.)
		{	soft[DAMAGENORMAL] = 1.;
			soft[DAMAGESHEAR] = 1.;
			soft[SOFT_DAMAGE_STATE] = 4.;
		}
		else
		{	soft[SOFT_DAMAGE_STATE] = dmode;
			double relToughness = soft[RELATIVE_TOUGHNESS];
			double relStrength = soft[RELATIVE_STRENGTH];
			double scale = relToughness*soft[GCSCALING]/initiationLaw->sigmaI();
			double relen0 = relStrength*en0;
			soft[DELTANORMAL] = softeningModeI->GetDeltaFromDamage(dparams.x,scale,relen0);
			scale = soft[GCSCALING]/initiationLaw->sigmaII();
			double relgs0 = relStrength*gs0;
			soft[DELTASHEAR] = softeningModeII->GetDeltaFromDamage(dparams.y,scale,relgs0);
		}
	}
}

// Return damge normal in global coordinates. It may be used by
// visualization to view normal. The magnitude of the vector
// is set equal to VpAc scaling for visualization use too.
Vector IsoSoftening::GetDamageNormal(MPMBase *mptr,bool threeD) const
{	// zero if no damage
	double *soft = GetSoftHistoryPtr(mptr);

	// Crack rotation
	Matrix3 RToCrack;
	int np = threeD ? THREED_MPM : PLANE_STRAIN_MPM;
	if(!GetRToCrack(&RToCrack,soft,np,0))
		return MakeVector(0.,0.,0.);
	
	// get current rotation matrix
	Matrix3 pF = mptr->GetDeformationGradientMatrix();
    Matrix3 Rtot;
	pF.RightDecompose(&Rtot,NULL);
	Rtot *= RToCrack;
	
	// rotate x axis and scale by Vp/Ac
	double scale = 1./(rho*soft[GCSCALING]);
	return MakeVector(Rtot(0,0)*scale,Rtot(1,0)*scale,Rtot(2,0)*scale);
}

#pragma mark IsoSoftening::Accessors

/* Take increments in strain and calculate new Particle: strains, rotation strain,
		stresses, strain energy,
	dvij are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
*/
void IsoSoftening::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,
									  ResidualStrains *res,int historyOffset,Tensor *gStress) const
{	// set 2D flag
	bool is2D = np == THREED_MPM ? false : true;

    Matrix3 dR,Rnm1,Rtot;
    Matrix3 deT = LRGetStrainIncrement(INITIAL_CONFIGURATION,mptr,du,&dR,NULL,&Rnm1,&Rtot);
    
    // convert matrix to Tensor in initial configuration
	Tensor de = is2D ?
		MakeTensor2D(deT(0,0), deT(1,1), deT(2,2), deT(0,1) + deT(1,0)) :
		MakeTensor(deT(0,0), deT(1,1), deT(2,2), deT(1,2)+deT(2,1), deT(0,2)+deT(2,0), deT(0,1)+deT(1,0));

	// cast pointer to material-specific data
	ElasticProperties *p = GetElasticPropertiesPointer(properties);
	
	// initial stresses
	Tensor *sp = mptr->GetStressTensor();
	
	// residual strain (thermal and moisture)
	// (CTE1 and CME1 are reduced to plane strain, but not CTE3 and CME3)
	double eres = CTE1*res->dT;
	double ezzres = CTE3*res->dT;
	if(DiffusionTask::HasDiffusion())
	{	// Only add for diffusion. Poroelasticity done latter so damage can use solid stress
		eres += CME1*res->dC;
		ezzres += CME3*res->dC;
	}
	// Uses reduced values
	Tensor deres = MakeTensor(eres,eres,eres,0.,0.,0.);
	
#pragma mark ... Code for Previously Undamaged Material
	// get history
	double *soft = GetSoftHistoryPtr(mptr);
	
	// Before cracking, do normal istropic update. If not cracked
	// then done, otherwise initiate the damage
	if(soft[SOFT_DAMAGE_STATE]<predamageState)
	{	// effective strain increment tensor in initial configuration (subtract reduced residual strain)
		Tensor deij = de;
		SubTensor(&deij,&deres);
		
		// Now get stress increment in initial configuration
		// Before rev 2593, this used de, which seems a mistake
		// Using wrong de would affect initiation when there are residual stresses
		Tensor dstr = GetStressIncrement(deij,np,properties);
		
		// Load str with effective stress in initial configuraiton used to look for initiation of failure
		// It is particle stress
		Tensor str = Rnm1.RTVoightR(sp,true,is2D);

		// Add incremental stress
		AddTensor(&str,&dstr);

		// check if has failed
		Vector norm;
		double relStrength = soft[RELATIVE_STRENGTH];
		int failureMode = initiationLaw->ShouldInitiateFailure(&str,&norm,np,relStrength);
		if(failureMode == NO_FAILURE)
		{	// Store str on the particle (rotating by Rtot to updated current config
			AcceptTrialStress(mptr,str,sp,np,&Rtot,properties,de,eres,ezzres);
			return;
		}
		
#pragma mark ...... Undamaged Material Just Initiated Damage
		// initiate failure (code depends on initiation law
		soft[SOFT_DAMAGE_STATE] = 0.01*failureMode ;
		soft[NORMALDIR1] = norm.x;										// cos(theta) (2D) or Euler alpha (3D) to normal
		soft[NORMALDIR2] = norm.y;										// sin(theta) (2D) or Euler beta (3D)  to normal
		soft[NORMALDIR3] = norm.z;										// unused (2D) or Euler gamma (3D) to normal
		
		// get intersection area (for 3D convert angle to normal)
		if(np==THREED_MPM)
		{	Matrix3 RInit;
			GetRToCrack(&RInit, soft, np, 0);
			norm = MakeVector(RInit(0,0),RInit(1,0),RInit(2,0));
		}

		soft[GCSCALING] = GetAcOverVp(np,mptr,&norm)/rho;			// rho because is divided by specific stress
	}
	
#pragma mark ... Code for Damaged Material
	// A crack is present
	
	// Rotate strain increment from initial config to crack axis system
	Matrix3 RToCrack;
	GetRToCrack(&RToCrack, soft, np, 0);
	de = RToCrack.RTVoightR(&de,false,is2D);
	
	// get total rotation to crack axis system
	Rnm1 *= RToCrack;
	Rtot *= RToCrack;
	
	// Get particle stress now rotated to the crack axis system
	Tensor str = Rnm1.RTVoightR(sp,true,is2D);

	// get cracking strain in crack axis system and apply
	// incremental rotation to particle cracking strain
	Tensor *eplast=mptr->GetAltStrainTensor();		// in global axes
	Tensor ecrack = Rnm1.RTVoightR(eplast,false,is2D);
	*eplast = dR.RVoightRT(eplast,false,is2D);
	
	// finish up in separate code
	DamageEvolution(mptr,np,soft,de,str,eres,ezzres,res,dR,Rtot,p,&ecrack,0.);
}

// Finish damage evolution for isotropic materials
// soft is pointer to first sofening history variable
// Tensor de = total elastic strain increment in crack axis system
// eres and ezzres used to find effective strains (from de - eres*I), ezzres for plain strain energy
// str is prior stress in the crack axis system
// res fir temperature and/or concentration/pore pressure change
// dR and Rtot are incemental and total rotation matrices
// p is pointer to isotropic elastic properties
// ecrack is cracking strain in the crack axis system
// dispEnergy is dissipated energy before evolving damage (or zero if none)
Vector IsoSoftening::DamageEvolution(MPMBase *mptr,int np,double *soft,Tensor &de,Tensor &str,double eres,
						double ezzres,ResidualStrains *res,Matrix3 &dR,Matrix3 &Rtot,ElasticProperties *p,
						Tensor *ecrack,double dispEnergy) const
{
	// set 2D flag
	bool is2D = np == THREED_MPM ? false : true;
	
	// some properties
	double relToughness = soft[RELATIVE_TOUGHNESS];
	double relStrength = soft[RELATIVE_STRENGTH];
	double relShearStrength = relStrength;
	
	// get effective strains (effective in plane strain)
	Tensor deres = MakeTensor(eres,eres,eres,0.,0.,0.);
	Tensor deij = de;
	SubTensor(&deij,&deres);
	double nuterm = np==PLANE_STRESS_MPM ? nu : nu/(1.-nu) ;
	double Txy0,dgs,den,Txz0=0.;
	
	if(np==THREED_MPM)
	{	// get dgs, current traction, and sign for traction change
		// sign found when softened
		Txy0 = str.xy;
		Txz0 = str.xz;
	}
	else
	{	// zz effective strain (may get non-zero in above 3D-style code)
		if(np!=AXISYMMETRIC_MPM) deij.zz = 0.;
		
		// get dgs, current traction, and sign for traction change
		Txy0 = fabs(str.xy);
		
		// Thus when tau>0, positive deij.xy may cause damage, but for tau<0 or negative deij.xy can cause damage
		// Finally damage only possible when dgs>0
		dgs = str.xy>0. ? deij.xy : -deij.xy;
	}
	
	// normal strain increments (deij.zz is zero for plane stress and strain)
	den = deij.xx + nuterm*(deij.yy+deij.zz);
	
	// to calculate depending on current state
	double decxx=0.,dgcxy=0.,dgcxz=0.;
	Tensor dsig = MakeTensor(0.,0.,0.,0.,0.,0.);
	double C11 = is2D ? p->C[1][1] : p->C[0][0];
	double Gshear = p->C[3][3];			// C44 in 3D and C66 in 2D are the same
	
#pragma mark ...... Damage Evolution Calculations
										// This section means crack has failed
	if(soft[SOFT_DAMAGE_STATE]>damageState)
	{	// post failure state - return increments in cracking strain (in decxx,dgcxy,dgcxz)
		// and incements in stress (in dsig)
		PostFailureUpdate(decxx,dgcxy,dgcxz,&dsig,&str,ecrack,Rtot,C11,den,deij.xy,deij.xz,is2D);
	}
	
	// The damage is evolving after the crack has initiated
	else
	{	// flag for reached critical strain
		bool criticalStrain = false;
		
		// normal stress and strains
		double sigmaI = relStrength*initiationLaw->sigmaI();
		double scaleI = relToughness*soft[GCSCALING]/sigmaI;
		if(!SoftenAxis(den,soft,DELTANORMAL,DAMAGENORMAL,sigmaI,scaleI,softeningModeI,
					   C11,str.xx,relStrength*en0,decxx,NULL,dispEnergy,criticalStrain))
		{	// Elastic loading but handle contact if den<0
			decxx = soft[DAMAGENORMAL]*den;
			if(den<0.)
			{	// keep ecxx>=0
				if(ecrack->xx+decxx<0.) decxx = -ecrack->xx;
			}
		}

		// shear stress and strains
		double sigmaII = relShearStrength*initiationLaw->sigmaII();
		double scaleII = relToughness*soft[GCSCALING]/sigmaII;
		if(is2D)
		{	// 2D shear strain only for x-y shear
			if(SoftenAxis(dgs,soft,DELTASHEAR,DAMAGESHEAR,sigmaII,scaleII,softeningModeII,
				Gshear,Txy0,relShearStrength*gs0,dgcxy,NULL,dispEnergy,criticalStrain))
			{	// Adjust sign to match shear stress direction
				if(str.xy<0.) dgcxy = -dgcxy;
			}
			else
			{	// elastic loading
				dgcxy = soft[DAMAGESHEAR]*deij.xy;
			}
		}
		else
		{	// 3D shear strains using rectangular failure surfacce
			Txy0 = fabs(str.xy);
			dgs = str.xy>0. ? deij.xy : -deij.xy;

			if(SoftenAxis(dgs,soft,DELTASHEAR,DAMAGESHEAR,sigmaII,scaleII,softeningModeII,
						  Gshear,Txy0,relShearStrength*gs0,dgcxy,NULL,dispEnergy,criticalStrain))
			{	// Adjust sign to match shear stress direction
				if(str.xy<0.) dgcxy = -dgcxy;
			}
			else
			{	// elastic loading
				dgcxy = soft[DAMAGESHEAR]*deij.xy;
			}
			Txz0 = fabs(str.xz);
			dgs = str.xz>0. ? deij.xz : -deij.xz;
			if(SoftenAxis(dgs,soft,DELTASHEAR2,DAMAGESHEAR2,sigmaII,scaleII,softeningModeII,
						  Gshear,Txz0,relShearStrength*gs0,dgcxz,NULL,dispEnergy,criticalStrain))
			{	// Adjust sign to match shear stress direction
				if(str.xz<0.) dgcxz = -dgcxz;
			}
			else
			{	// elastic loading
				dgcxz = soft[DAMAGESHEAR2]*deij.xz;
			}
		}
		
#pragma mark ...... Dissipated Energy and Check for Final Failure
		// Use energy release rate to check on failure
		
		// Get GI/GIc
		double relGI = softeningModeI->GetGoverGc(soft[DELTANORMAL],scaleI);

		// Get (GII/GIIc)'s
		double relGII = softeningModeII->GetGoverGc(soft[DELTASHEAR],scaleII);
		double relGIIxz = softeningModeII->GetGoverGc(soft[DELTASHEAR2],scaleII);
		
		// criterion
		double cutoff = pow(relGI,nmix) + pow(relGII,nmix) + pow(relGIIxz,nmix);
		
		// check limit or energy
		if(criticalStrain || cutoff>=1.)
		{	// report energy released
			double GI = relToughness*relGI*softeningModeI->GetGc()*UnitsController::Scaling(1.e-3);
			double GII = relToughness*(relGII+relGIIxz)*softeningModeII->GetGc()*UnitsController::Scaling(1.e-3);
			double alpha=soft[NORMALDIR1],beta=soft[NORMALDIR2],gamma=soft[NORMALDIR3];
			if(is2D)
			{	alpha = beta>=0. ? acos(alpha) : -acos(alpha) ;
				beta = 0.;
			}
			archiver->Decohesion(mtime,mptr,alpha,beta,gamma,GI,GII,0.0,soft[SOFT_DAMAGE_STATE]);
				
			if(GI!=GI)
			{
#pragma omp critical (output)
				{	cout << "#     enmax=" << soft[DELTANORMAL] << "/" << softeningModeI->GetDeltaMax(scaleI) <<
						", gsmax=" << soft[DELTASHEAR] << "/" << softeningModeII->GetDeltaMax(scaleII) <<
						", dn=" << (soft[DELTANORMAL]/(soft[DELTANORMAL]+relStrength*en0*softeningModeI->GetFFxn(soft[DELTANORMAL],scaleI))) <<
						", ds=" << (soft[DELTASHEAR]/(soft[DELTASHEAR]+relShearStrength*gs0*softeningModeII->GetFFxn(soft[DELTASHEAR],scaleII))) <<
						", scaleI/scaleII=" << scaleI << "/" << scaleII <<
						", relGI/relGII,relGIIxz=" << relGI << "/" << relGII << "/" << relGIIxz <<
						", Ac/Vp=" << soft[GCSCALING] << endl;
				}
			}
			
			// now failed
			// From 2 to 3 will be fraction mode I
			soft[SOFT_DAMAGE_STATE] = 2. + GI/(GI+GII);
			soft[DAMAGENORMAL] = 1.;
			soft[DAMAGESHEAR] = 1.;
			soft[DAMAGESHEAR2] = 1.;
			
			// post failure update
			PostFailureUpdate(decxx,dgcxy,dgcxz,&dsig,&str,ecrack,Rtot,C11,den,deij.xy,deij.xz,is2D);
		}
		else
		{	// crack plane stress increment
			dsig.xx = C11*(den - decxx);
			dsig.xz = Gshear*(deij.xz-dgcxz);
			dsig.xy = Gshear*(deij.xy-dgcxy);
		}
		
		// dissipated energy
		mptr->AddPlastEnergy(dispEnergy);
	}
		
#pragma mark ... Update Cracking Strain and Stresses
	// Have now found increments in cracking strain and stress
	// update cracking strains (in plastic strain)
	// current was rotated by dR (above) now add increment rotated by Rtot
	
	// rotate previous stress by rotation increment
	Tensor *sp = mptr->GetStressTensor();
	str = dR.RVoightRT(sp,true,is2D);

	// get stress incement for other stresses in crack axis system
	if(np==THREED_MPM)
	{	dsig.yy = C11*(deij.yy + nuterm*(deij.xx - decxx + deij.zz));
		dsig.zz = C11*(deij.zz + nuterm*(deij.xx - decxx + deij.yy));
		dsig.yz = Gshear*deij.yz;
	}
	else if(np==AXISYMMETRIC_MPM)
	{	dsig.yy = C11*(deij.yy + nuterm*(deij.xx - decxx - deij.zz));
		dsig.zz = p->C[1][1]*(deij.zz + nuterm*(deij.xx - decxx + deij.yy));
	}
	else
	{	dsig.yy = C11*(deij.yy + nuterm*(deij.xx - decxx));
		if(np==PLANE_STRAIN_MPM)
			dsig.zz = p->C[1][1]*(nuterm*(de.xx - decxx + de.yy - 2.*ezzres) - ezzres);
		else if(np==PLANE_STRESS_MPM)
		{	// zz deformation. zz stress stays at zero for plane stress
			// Here C[4][1] = C[4][2] = -v/(1-v)
			de.zz = p->C[4][1]*(deij.xx-decxx) + p->C[4][2]*deij.yy + eres;
			mptr->IncrementDeformationGradientZZ(de.zz);
		}
	}
	
	// update stress and strain (in which cracking strain rotated to current configuration)
	UpdateCrackingStrainStress(np,mptr->GetAltStrainTensor(),sp,decxx,dgcxy,dgcxz,&dsig,&str,Rtot);
	
#pragma mark ... Work, Residual, and Heat Energies
	// work energy, residual energy, and heat energy
	de = Rtot.RVoightRT(&de,false,is2D);
	if(np==THREED_MPM)
	{	// work and residual energy
		mptr->AddWorkEnergyAndResidualEnergy(DotTensors(sp,&de),(sp->xx + sp->yy + sp->zz)*eres);
		
	}
	else
	{	// work energy, residual energy in current configuration
		double workEnergy = sp->xx*de.xx + sp->yy*de.yy + sp->xy*deij.xy;
		double resEnergy = (sp->xx + sp->yy)*ezzres;
		if(np==PLANE_STRAIN_MPM)
		{	// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			resEnergy += sp->zz*ezzres;
		}
		else if(np==AXISYMMETRIC_MPM)
		{	// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			workEnergy += sp->zz*de.zz;
			resEnergy += sp->zz*eres;
		}
		mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
	}
	
	// Isoentropic temperature rise = -(K 3 alpha T)/(rho Cv) (dV/V) = - gamma0 T (dV/V)
	double delV = de.xx+de.yy+de.zz;
	double dTq0 = -gamma0*mptr->pPreviousTemperature*delV;
	
	// track heat energy
	IncrementHeatEnergy(mptr,dTq0,dispEnergy);
	
	// return cracking strains
	return MakeVector(decxx,dgcxy,dgcxz);
}

// Calculate rotation matrix from crack to initial
bool IsoSoftening::GetRToCrack(Matrix3 *R,double *soft, int np, int Dstyle) const
{	// none if undamaged
	if (soft[SOFT_DAMAGE_STATE] < predamageState) return false;

	// 3D or 2D
	if(np==THREED_MPM)
	{	// get sins an cosines
		double c1 = cos(soft[NORMALDIR1]);
		double s1 = sin(soft[NORMALDIR1]);
		double c2 = cos(soft[NORMALDIR2]);
		double s2 = sin(soft[NORMALDIR2]);
		double c3 = cos(soft[NORMALDIR3]);
		double s3 = sin(soft[NORMALDIR3]);
		R->set(c1*c2*c3-s1*s3,  -c3*s1-c1*c2*s3,  c1*s2,
			c1*s3+c2*c3*s1,   c1*c3-c2*s1*s3,  s1*s2,
			-c3*s2,           s2*s3,        c2    );
	}
	else
	{	// already cos and sin
		R->set(soft[NORMALDIR1], -soft[NORMALDIR2], soft[NORMALDIR2], soft[NORMALDIR1], 1);
	}
	return true;
}

#pragma mark IsoSoftening::Isotropic Elasticity Methods

// For isotropic elastic material, find de(eff). Results same for isotropic material in all coordinate systems
Tensor IsoSoftening::GetStressIncrement(Tensor &deij,int np,void *properties) const
{
	// cast pointer to material-specific data
	ElasticProperties *p = GetElasticPropertiesPointer(properties);

	Tensor dstr;
	ZeroTensor(&dstr);
	
	// stress increments in initial analysis axes
	if(np==THREED_MPM)
	{	// elastic stress update
		dstr.xx = p->C[0][0]*deij.xx + p->C[0][1]*deij.yy + p->C[0][2]*deij.zz;
		dstr.yy = p->C[1][0]*deij.xx + p->C[1][1]*deij.yy + p->C[1][2]*deij.zz;
		dstr.zz = p->C[2][0]*deij.xx + p->C[2][1]*deij.yy + p->C[2][2]*deij.zz;
		dstr.yz = p->C[3][3]*deij.yz;
		dstr.xz = p->C[4][4]*deij.xz;
		dstr.xy = p->C[5][5]*deij.xy;
	}
	
	else
	{	// update stress in initial analysis axes
		if(np==AXISYMMETRIC_MPM)
		{	// hoop stress affect on RR, ZZ, and RZ stresses
			dstr.xx = p->C[1][1]*deij.xx + p->C[1][2]*deij.yy + p->C[4][1]*deij.zz;
			dstr.yy = p->C[1][2]*deij.xx + p->C[2][2]*deij.yy + p->C[4][2]*deij.zz;
		}
		else
		{	dstr.xx += p->C[1][1]*deij.xx + p->C[1][2]*deij.yy;
			dstr.yy += p->C[1][2]*deij.xx + p->C[2][2]*deij.yy;
		}
		dstr.xy += p->C[3][3]*deij.xy;
	}

	// return increment
	return dstr;
}

// For isotropic elastic material, accept new trial stress state
// str, de in initial configuration and need to rotate to current using Rtot
// Need rotate to current axes and update stress and all other terms
void IsoSoftening::AcceptTrialStress(MPMBase *mptr,Tensor &str,Tensor *sp,int np,Matrix3 *Rtot,
									 void *properties,Tensor &de,double eres,double ezzres) const
{
	if(np==THREED_MPM)
	{	// Rotate and store trial stress on the particle
		*sp = Rtot->RVoightRT(&str,true,false);
		
		// work energy increment per unit mass (dU/(rho0 V0)) (OK to use stress and strain in initial configuration)
		mptr->AddWorkEnergyAndResidualEnergy(str.xx*de.xx + str.yy*de.yy + str.zz*de.zz
											 + str.xz*de.yz + str.xz*de.xz + str.xy*de.xy,
											 (str.xx + str.yy + str.zz)*eres);
	}
	
	else
	{	// Rotate to current fconfiguration, but not ready for particle yet
		Tensor stnp1 = Rtot->RVoightRT(&str,true,true);
		sp->xx = stnp1.xx;
		sp->yy = stnp1.yy;
		sp->xy = stnp1.xy;
		
		// cast pointer to material-specific data
		ElasticProperties *p = GetElasticPropertiesPointer(properties);
		
		// work and residual strain energy increments (in initial axes)
		double workEnergy = str.xx*de.xx + str.yy*de.yy + str.xy*de.xy;
		double resEnergy = (str.xx + str.yy)*ezzres;
		if(np==PLANE_STRAIN_MPM)
		{	// need to add back terms to get from reduced cte to actual cte
			sp->zz += p->C[4][1]*(de.xx-ezzres) + p->C[4][2]*(de.yy-ezzres) - p->C[4][4]*ezzres;
			
			// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			resEnergy += sp->zz*ezzres;
		}
		else if(np==PLANE_STRESS_MPM)
		{	// zz deformation. zz stress stays at zero
			// Here C[4][1] = C[4][2] = -v/(1-v)
			de.zz = p->C[4][1]*(de.xx-eres) + p->C[4][2]*(de.yy-eres) + eres;
			mptr->IncrementDeformationGradientZZ(de.zz);
		}
		else
		{	// axisymmetric hoop stress
			sp->zz += p->C[4][1]*(de.xx-eres) + p->C[4][2]*(de.yy-eres) + p->C[4][4]*(de.zz-eres);
			
			// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			workEnergy += sp->zz*de.zz;
			resEnergy += sp->zz*eres;
		}
		mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
	}
	
	// Isoentropic temperature rise = -(K 3 alpha T)/(rho Cv) (dV/V) = - gamma0 T (dV/V)
	double delV = de.xx+de.yy+de.zz;
	double dTq0 = -gamma0*mptr->pPreviousTemperature*delV;
	
	// track heat energy
	IncrementHeatEnergy(mptr,dTq0,0.);
}

#pragma mark IsoSoftening::Accessors

// return material type
const char *IsoSoftening::MaterialType(void) const { return "Isotropic Softening"; }

// store plastic strain in alt strain
int IsoSoftening::AltStrainContains(void) const
{	return ENG_BIOT_PLASTIC_STRAIN;
}

// zero-offset to history data
double *IsoSoftening::GetSoftHistoryPtr(MPMBase *mptr) const
{	return (double *)(mptr->GetHistoryPtr(0));
}

