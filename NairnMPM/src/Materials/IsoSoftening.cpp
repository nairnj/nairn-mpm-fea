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
#include "Exceptions/MPMWarnings.hpp"

#define MAXIMUM_D 0.9999

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
	tractionFailureSurface = CUBOID_SURFACE;
	distributionMode = SOFTDIST_NONE;
	softenStatsMode = VARY_STRENGTH;
	softenCV = 0.;
 	wAlpha = -1.;
	wV0 = -1.;
	frictionCoeff = -1.;				// <0 disables friction
    
    // tolerance of fabs(delta(decexx))/deltaMax when using Ovoid
    // mode with pressure-dependent shear strength
    pdOvoidTolerance = 1.e-9;
    maxOvoidPasses = 10;
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
	
	// type of failure surface (shearFailureSurface for backward compatiility
	//   but better use only documented 0 or 1
	else if(strcmp(xName,"tractionFailureSurface")==0 || strcmp(xName,"shearFailureSurface")==0)
	{	input=INT_NUM;
		return((char *)&tractionFailureSurface);
	}
	
	// softening coefficient of variation
	else if(strcmp(xName,"coefVariation")==0)
	{	input=DOUBLE_NUM;
		distributionMode=SOFTDIST_NORMAL;
		return((char *)&softenCV);
	}
	
    // coefficient of variation for strength (1), toughness (2), or both (3)
    else if(strcmp(xName,"coefVariationMode")==0 || strcmp(xName,"statDistributionMode")==0)
    {	input=INT_NUM;
        return((char *)&softenStatsMode);
    }
	
	// Weibull alpha or scale parameter
	else if(strcmp(xName,"wShape")==0)
	{	input=DOUBLE_NUM;
		distributionMode=SOFTDIST_WEIBULL;
		return((char *)&wAlpha);
	}
	
	// Weibull reference length
	else if(strcmp(xName,"wV0")==0)
	{	input=DOUBLE_NUM;
		distributionMode=SOFTDIST_WEIBULL;
		return((char *)&wV0);
	}
	
	// coefficient of friction (and activates friction)
	else if(strcmp(xName,"coeff")==0)
	{	input=DOUBLE_NUM;
		return (char *)&frictionCoeff;
	}
	
    // pressure dependent loop tolerances
    else if(strcmp(xName,"PDTol")==0)
    {   input=DOUBLE_NUM;
        return (char *)&pdOvoidTolerance;
    }
    
    // pressure dependent loop maximum passes
    else if(strcmp(xName,"PDPasses")==0)
    {   input=INT_NUM;
        return (char *)&maxOvoidPasses;
    }
    
    // check initiation law
	if(initiationLaw!=NULL)
	{	char *ptr = initiationLaw->InputInitationProperty(xName,input,gScaling);
		if(ptr != NULL) return ptr;
	}
	
	// if possible check softening laws
	if(strlen(xName)>2)
	{	if(xName[0]=='I' && xName[1]=='-')
		{	char *ptr = softeningModeI->InputSofteningProperty(&xName[2],input,gScaling);
			if(ptr != NULL) return ptr;
		}
		else if(xName[0]=='I' && xName[1]=='I' && xName[2]=='-' && strlen(xName)>3)
		{	char *ptr = softeningModeII->InputSofteningProperty(&xName[3],input,gScaling);
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
	
	if(initiationLaw!=NULL) delete initiationLaw;
    initiationLaw = iLaw;
    return true;
}

// Allows mode I and mode II softening laws
bool IsoSoftening::AcceptSofteningLaw(SofteningLaw *sLaw,int mode)
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

#ifdef POROELASTICITY
    // Poroelasticity cannot be mixed with pressure dependence
    if(DiffusionTask::HasPoroelasticity())
    {   if(initiationLaw->IsPressureDependent())
            return "Simulations with poroelasticity cannot model pressure-dependent initiation stress";
    }
#endif
	
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
    
    if(softenStatsMode<VARY_STRENGTH || softenStatsMode>VARY_STRENGTH_AND_TOUGHNESS)
        return "Invalid option selected for which softening properties to vary";
	
	// verify distibution parameters
	if(distributionMode==SOFTDIST_NORMAL)
	{	if(softenCV<0.)
			return "Coefficient of variation cannot be negative";
		else if(softenCV==0.)
			distributionMode=SOFTDIST_NONE;
	}
	else if(distributionMode==SOFTDIST_WEIBULL)
	{	if(wAlpha<=0. || wV0<=0.)
			return "Weibull scale and reference volume must be greater than 0";
		wGam1A = gamma_fxn(1.+1./wAlpha);
		softenCV = sqrt(gamma_fxn(1.+2./wAlpha)/(wGam1A*wGam1A) - 1.);
	}
			
	// make sure have a valid setting
	if(tractionFailureSurface!=CUBOID_SURFACE
	   		&& tractionFailureSurface!=CYLINDER_SURFACE
			&& tractionFailureSurface!=OVOID_SURFACE
			&& tractionFailureSurface!=COUPLED_CUBOID_SURFACE)
	{	return "Invalid option selected for the traction failure surface (must be 0 to 3)";
		
	}
	
	// did not implement pressure dependence in the coupled cuboid surface (would not be hard though)
	if(tractionFailureSurface==COUPLED_CUBOID_SURFACE && initiationLaw->IsPressureDependent())
		return "Simulations with coupled cuboid surface cannot model pressure-dependent initiation stress";
	
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
	if(distributionMode==SOFTDIST_WEIBULL)
	{	cout << "Weibull scale = " << wAlpha << " and V0 = " << wV0 <<
                " " << UnitsController::Label(CULENGTH_UNITS) << "^3" << endl;
	}
	cout << "Coefficient of Variation = " << softenCV << endl;
    if(distributionMode!=SOFTDIST_NONE)
    {   if(softenStatsMode==VARY_STRENGTH_AND_TOUGHNESS)
            cout << "    Vary strength and toughness" << endl;
        else if(softenStatsMode==VARY_TOUGHNESS)
            cout << "    Vary toughness" << endl;
        else
            cout << "    Vary strength" << endl;
    }
	
	initiationLaw->SetFailureSurface(tractionFailureSurface);
	cout << "Traction Failure Surface: ";
	if(tractionFailureSurface == CYLINDER_SURFACE)
		cout << "cylindrical" << endl;
	else if(tractionFailureSurface == CUBOID_SURFACE)
		cout << "cuboid" << endl;
	else if(tractionFailureSurface == OVOID_SURFACE)
		cout << "ovoid" << endl;
	else
		cout << "coupled cuboid" << endl;
	
	cout << "Post-failure contact: ";
	if(frictionCoeff>0.)
		cout << "Coulomb coefficient of friction = " << frictionCoeff << endl;
	else
		cout << "frictionless" << endl;
}

#pragma mark IsoSoftening::History Data Methods

// Create history variables needed for softening behavior
// Also check to stability:
// Need cell size < eta p(relG) Gc/(p(rel)^2 ei sigi) where ei is initiation strain and sigi
// is initiation stress and p(rel) is relative strength at specified maxPressure and
// p(relG) is relative toughness. Below finds
//      ratio = (min cell size)(p(rel)^2 ei sigi)/(eta p(relG) Gc)
// If any partivle has ratio>1 for tension or shear, the properties on not stable
// ratio gives factor must decrease (min cell size), p(rel)^2, sig^2 or increase Gc
char *IsoSoftening::InitHistoryData(char *pchr,MPMBase *mptr)
{
	double *p = CreateAndZeroDoubles(pchr,SOFT_NUMBER_HISTORY);
	
	// set damage state to 0.1 instead of zero for help in plotting on a grid
	p[SOFT_DAMAGE_STATE] = 0.1;
	
	// set relative strength and toughness
	// Note that both set to see if used. Future idea would be to use different
	//		stats for strength and toughness
	p[RELATIVE_STRENGTH] = 1.;
	p[RELATIVE_TOUGHNESS] = 1.;
	if(distributionMode==SOFTDIST_NORMAL)
	{	double fract = RandomRange(0.001,0.999);		// .001 to .999, max +/- 3.09 std dev
		double relValue = fmax(1. + softenCV*NormalCDFInverse(fract),0.05);
		if(softenStatsMode&VARY_STRENGTH)
			p[RELATIVE_STRENGTH] = relValue;
		if(softenStatsMode&VARY_TOUGHNESS)
			p[RELATIVE_TOUGHNESS] = relValue;
	}
	else if(distributionMode==SOFTDIST_WEIBULL)
	{	double fract = RandomRange(0.001,0.999);		// .001 to .999 limited range
        double vp = mptr->GetUnscaledVolume();
		double relValue = fmax(pow(-wV0*log(1.-fract)/vp,1./wAlpha)/wGam1A,0.05);
		if(softenStatsMode&VARY_STRENGTH)
			p[RELATIVE_STRENGTH] = relValue;
		if(softenStatsMode&VARY_TOUGHNESS)
			p[RELATIVE_TOUGHNESS] = relValue;
	}

	// Validate this law (in use) is stable on current grid
	// redundant when all particle same size and same relative values, but general if not
    double delx = mptr->GetMinParticleLength();

    // mode I stability
	double sigma0 = rho*initiationLaw->sigmaI();
	double Gc = p[RELATIVE_TOUGHNESS]*softeningModeI->GetGc();
    // get relative strength possibly pressure dependent
	double scale = p[RELATIVE_STRENGTH]*initiationLaw->sigmaIStabilityScale();
	double ratio = delx*en0*sigma0*scale*scale/(softeningModeI->GetEtaStability()*Gc);
	if(ratio>1.)
    {    cout << "Need to increase GIc or decrease particle size " << ratio
                << " fold or decrease sigmac[P] " << sqrt(ratio)
                << " fold for spatial stability." << endl;
        throw CommonException("Normal softening law needs smaller grid cells","IsoSoftening::InitHistoryData");
    }
    
    // mode II
	if(softeningModeII!=NULL)
	{	sigma0 = rho*initiationLaw->sigmaII();
		Gc = p[RELATIVE_TOUGHNESS]*softeningModeII->GetGc();
		// get relative strength possibly pressure dependent
		scale = p[RELATIVE_STRENGTH]*initiationLaw->sigmaIIStabilityScale();
		ratio = delx*gs0*sigma0*scale*scale/(softeningModeII->GetEtaStability()*Gc);
		if(ratio>1.)
		{	cout << "Need to increase GIIc or decrease particle size " << ratio
					<< " fold or decrease tauc[P] " << sqrt(ratio)
					<< " fold for spatial stability." << endl;
			throw CommonException("Shear softening law needs smaller grid cells","IsoSoftening::InitHistoryData");
		}
	}

	return (char *)p;
}

// reset history data
void IsoSoftening::ResetHistoryData(char *pchr,MPMBase *mptr)
{	double *p = (double *)pchr;
	double relStrength = p[RELATIVE_STRENGTH];
	double relToughness = p[RELATIVE_TOUGHNESS];
	ZeroDoubles(pchr,SOFT_NUMBER_HISTORY);
	p[SOFT_DAMAGE_STATE] = 0.1;
	p[RELATIVE_STRENGTH] = relStrength;
	p[RELATIVE_TOUGHNESS] = relToughness;
}

// Number of history variables
int IsoSoftening::NumberOfHistoryDoubles(void) const { return SOFT_NUMBER_HISTORY; }

// Initialize damage when requested
void IsoSoftening::SetInitialConditions(InitialCondition *ic,int ptNum,bool is3D)
{
    MPMBase *mptr = mpm[ptNum-1];
    double *soft = GetSoftHistoryPtr(mptr);

	if(is3D)
    {   Vector dnorm = ic->GetDamageNormal();
		soft[NORMALDIR1] = atan2(dnorm.y,dnorm.x); // Euler angle alpha (should be (y,x))
		soft[NORMALDIR2] = asin(dnorm.z);  //Euler angle beta using asin()Da
		soft[NORMALDIR3] = 0.0; // only 2 angles needed for a normal vector

		// Get failure mode (if not fully damaged) and damage parameters
		double dmode;
		Vector dparams = ic->GetDamageParams(dmode);
		soft[DAMAGENORMAL] = dparams.x;
		soft[DAMAGESHEAR] = dparams.y;
		// otherwise DAMAGESHEAR2=DAMAGEGII is used to store mode II energy
		if(tractionFailureSurface == CUBOID_SURFACE)
			soft[DAMAGESHEAR2] = dparams.z;

		// Find Ac/Vp
		soft[GCSCALING] = GetAcOverVp(THREED_MPM, mptr, &dnorm) / rho;
		
		// Any damage variable >=1 means failed, otherwise calculcate corresponding cracking strains
		// and set provide failure mode
		if(dparams.x>=1. || dparams.y>=1. || dparams.z>=1.)
		{	soft[DAMAGENORMAL] = 1.;
			soft[DAMAGESHEAR] = 1.;
			// otherwise DAMAGESHEAR2 is used to store mode II energy
			if(tractionFailureSurface == CUBOID_SURFACE)
				soft[DAMAGESHEAR2] = 1.;
			soft[SOFT_DAMAGE_STATE] = 4.;
		}
		else
		{	soft[SOFT_DAMAGE_STATE] = dmode;
			double relToughness = soft[RELATIVE_TOUGHNESS];
			double relStrength = soft[RELATIVE_STRENGTH];
			double scale = relToughness*soft[GCSCALING]/(relStrength*initiationLaw->sigmaI());
			double relen0 = relStrength*en0;
			soft[DELTANORMAL] = softeningModeI->GetDeltaFromDamage(dparams.x,scale,relen0,-1.);
			scale = relToughness*soft[GCSCALING]/(relStrength*initiationLaw->sigmaII());
			double relgs0 = relStrength*gs0;
			soft[DELTASHEAR] = softeningModeII->GetDeltaFromDamage(dparams.y,scale,relgs0,-1.);
			// otherwise DELTASHEAR2=DAMAGEGI is used to store mode I energy
			if(tractionFailureSurface == CUBOID_SURFACE)
				soft[DELTASHEAR2] = softeningModeII->GetDeltaFromDamage(dparams.z,scale,relgs0,-1.);
		}
	}
	else
    {   // 2D damage
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
			double scale = relToughness*soft[GCSCALING]/(relStrength*initiationLaw->sigmaI());
			double relen0 = relStrength*en0;
			soft[DELTANORMAL] = softeningModeI->GetDeltaFromDamage(dparams.x,scale,relen0,-1.);
			scale = relToughness*soft[GCSCALING]/(relStrength*initiationLaw->sigmaII());
			double relgs0 = relStrength*gs0;
			soft[DELTASHEAR] = softeningModeII->GetDeltaFromDamage(dparams.y,scale,relgs0,-1.);
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
	if(!GetRToCrack(&RToCrack,soft,!threeD,0))
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

// set relative strength and toughness. Only called by CustomThermalRamp custom tasks
// to set disteibution of strengths and toughness.
void IsoSoftening::SetRelativeStrength(MPMBase *mptr,double rel)
{   double *soft = GetSoftHistoryPtr(mptr);
    soft[RELATIVE_STRENGTH] = rel;
}
void IsoSoftening::SetRelativeToughness(MPMBase *mptr,double rel)
{   double *soft = GetSoftHistoryPtr(mptr);
    soft[RELATIVE_TOUGHNESS] = rel;
}

#pragma mark IsoSoftening::Methods

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
    
    // convert matrix to Tensor in initial configuration - total strain increment (engr shear strains)
	Tensor de = is2D ?
		MakeTensor2D(deT(0,0), deT(1,1), deT(2,2), deT(0,1) + deT(1,0)) :
		MakeTensor(deT(0,0), deT(1,1), deT(2,2), deT(1,2)+deT(2,1), deT(0,2)+deT(2,0), deT(0,1)+deT(1,0));

	// cast pointer to material-specific data
	ElasticProperties *p = GetElasticPropertiesPointer(properties);
	
	// initial stresses
	Tensor *sp = mptr->GetStressTensor();
	
	// residual strain (thermal and moisture)
	// (CTE1 and CME1 are reduced for plane strain, but not CTE3 and CME3)
	double eres = CTE1*res->dT;
	double ezzres = CTE3*res->dT;
	if(DiffusionTask::HasDiffusion())
	{	// Only add for diffusion. Poroelasticity done latter so damage can use solid stress
		eres += CME1*res->dC;
		ezzres += CME3*res->dC;
	}
	// Uses reduced values in plane strain
	Tensor deres = MakeTensor(eres,eres,eres,0.,0.,0.);
	
	// Poroelasticity needs to handle current specific pore pressure
	double rho0 = 1.;			// only defined and used when poroelasticity, calculate once
#ifdef POROELASTICITY
	if(DiffusionTask::HasPoroelasticity())
	{	// get current rho
		int matid = mptr->MatID();
		rho0=theMaterials[matid]->GetRho(mptr);
	}
#endif
	
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
		
		// Load str with effective stress in initial configuration used to look for initiation of failure
		// Need initial, because strains were found in initial configuration
		Tensor str = Rnm1.RTVoightR(sp,true,is2D);

		// Add incremental stress
		AddTensor(&str,&dstr);
		
#ifdef POROELASTICITY
		// Poroelasticity changes to solid stress by adding alphaB*pp*I
		double ppadd = 0.;
		if(DiffusionTask::HasPoroelasticity())
		{	// adding specific pore pressure
			ppadd = alphaPE*fmax(mptr->pDiff[0]->conc,0.)/rho0;
			str.xx += ppadd;
			str.yy += ppadd;
			if(np==THREED_MPM) str.zz += ppadd;
		}
#endif
        
		// check if has failed
		Vector norm;
		double relStrength = soft[RELATIVE_STRENGTH];
        GenADaMVariables alpha;
		int failureMode = initiationLaw->ShouldInitiateFailure(&str,&norm,np,relStrength,&alpha);
		if(failureMode == NO_FAILURE)
		{
#ifdef POROELASTICITY
			// Poroelasticity now subtracts pore pressure added above and pore pressure residual strain term
			if(DiffusionTask::HasPoroelasticity())
			{
				ppadd = ppadd + alphaPE*res->dC/rho0;
				str.xx -= ppadd;
				str.yy -= ppadd;
				if(np==THREED_MPM)
					str.zz -= ppadd;
				else
				{	// increment these to get z terms correct
					eres += CME1*res->dC;
					ezzres += CME3*res->dC;
				}
			}
#endif
			
			// Store str on the particle (rotating by Rtot to updated current config
			AcceptTrialStress(mptr,str,sp,np,&Rtot,properties,de,eres,ezzres);
			return;
		}
		
#pragma mark ...... Undamaged Material Just Initiated Damage
		// initiate failure (code depends on initiation law)
		soft[SOFT_DAMAGE_STATE] = 0.01*failureMode ;
		soft[NORMALDIR1] = norm.x;										// cos(theta) (2D) or Euler alpha (3D) to normal
		soft[NORMALDIR2] = norm.y;										// sin(theta) (2D) or Euler beta (3D)  to normal
		soft[NORMALDIR3] = norm.z;										// unused (2D) or Euler gamma (3D) to normal
		
		// hack to fix the normal
		/*
		soft[NORMALDIR1] = 0.;
		soft[NORMALDIR2] = norm.y>0. ? 1. : 0. ;
		soft[NORMALDIR3] = 0.;
		*/
		
		// get intersection area (for 3D convert angle to normal)
		if(np==THREED_MPM)
		{	Matrix3 RInit;
			GetRToCrack(&RInit, soft, is2D, 0);
			norm = MakeVector(RInit(0,0),RInit(1,0),RInit(2,0));
		}

		soft[GCSCALING] = GetAcOverVp(np,mptr,&norm)/rho;			// rho because is divided by specific stress
    }
	
#pragma mark ... Code for Damaged Material
	// A crack is present
	
	// Rotate strain increment from initial config to crack axis system
	Matrix3 RToCrack;
	GetRToCrack(&RToCrack, soft, is2D, 0);
	de = RToCrack.RTVoightR(&de,false,is2D);
	
	// get total rotation to crack axis system
	Rnm1 *= RToCrack;
	Rtot *= RToCrack;
    
    Tensor str = Rnm1.RTVoightR(sp,true,is2D);

#ifdef POROELASTICITY
	// Poroelasticity bases damage on solid stress by adding alphaB*pp*I
	if(DiffusionTask::HasPoroelasticity())
	{	double ppadd = alphaPE*fmax(mptr->pDiff[0]->conc,0.)/rho0;
		str.xx += ppadd;
		str.yy += ppadd;
		if(np==THREED_MPM) str.zz += ppadd;
	}
#endif

	// get cracking strain in crack axis system (in ecrack)
	// apply incremental rotation to cracking strain on the particle (in eplast)
	Tensor *eplast=mptr->GetAltStrainTensor();		// in global axes
	Tensor ecrack = Rnm1.RTVoightR(eplast,false,is2D);
	*eplast = dR.RVoightRT(eplast,false,is2D);
	
	// finish up in separate code
	DamageEvolution(mptr,np,soft,de,str,eres,ezzres,res,dR,Rtot,p,&ecrack,0.,delTime);
}

// Finish damage evolution for isotropic materials
// soft is pointer to first sofening history variable
// Tensor de = total elastic strain increment in crack axis system
// eres and ezzres used to find effective strains (from de - eres*I), ezzres for plain strain energy
// str is prior stress in the crack axis system (smoothed stress is using that option)
// res for temperature and/or concentration/pore pressure change
// dR and Rtot are incemental and total rotation matrices
// p is pointer to isotropic elastic properties
// ecrack is cracking strain in the crack axis system
// dispEnergy is dissipated energy before evolving damage (or zero if none)
void IsoSoftening::DamageEvolution(MPMBase *mptr,int np,double *soft,Tensor &de,Tensor &str,double eres,
						double ezzres,ResidualStrains *res,Matrix3 &dR,Matrix3 &Rtot,ElasticProperties *p,
						Tensor *ecrack,double dispEnergy,double delTime) const
{
	// set 2D flag
	bool is2D = np == THREED_MPM ? false : true;
	
	// get effective total strain increments (effective in plane strain)
	Tensor deres = MakeTensor(eres,eres,eres,0.,0.,0.);
	Tensor deij = de;
	SubTensor(&deij,&deres);
    double nuterm = np==PLANE_STRESS_MPM ? nu : nu/(1.-nu) ;
    double C11 = is2D ? p->C[1][1] : p->C[0][0];    // C11^r in plane strain, otherwise C11
    double Gshear = p->C[3][3];            // C44 in 3D and C66 in 2D are the same

    // den is increment normal strain in crack axis system
    // alpha.dPe is increment in pressure for elastic material when cracking strain is zero
    GenADaMVariables alpha;
    double den;
    if(np==PLANE_STRAIN_MPM)
    {   den = deij.xx + nuterm*deij.yy;
        alpha.dPe = -C11*( (1+nuterm)*(deij.xx + deij.yy)
                          + nuterm*(de.xx+de.yy-2.*ezzres) - ezzres )/3.;
    }
    else if(np==PLANE_STRESS_MPM)
    {   den = deij.xx + nuterm*deij.yy;
        alpha.dPe = -C11*(1+nuterm)*(deij.xx + deij.yy)/3.;
    }
    else
    {   den = deij.xx + nuterm*(deij.yy+deij.zz);
        alpha.dPe = -C11*(1+2*nuterm)*(deij.xx + deij.yy + deij.zz)/3.;
    }
    
	// to calculate depending on current state
	double decxx=0.,dgcxy=0.,dgcxz=0.;
	Tensor dsig = MakeTensor(0.,0.,0.,0.,0.,0.);
	
#pragma mark ...... Damage Evolution Calculations
	// This section means crack has failed
	if(soft[SOFT_DAMAGE_STATE]>damageState)
	{	// post failure state - return increments in cracking strain (in decxx,dgcxy,dgcxz)
		// and incements in stress (in dsig)
		PostFailureUpdate(decxx,dgcxy,dgcxz,&dsig,&str,ecrack,Rtot,C11,den,Gshear,deij.xy,Gshear,deij.xz,is2D,frictionCoeff);
	}
	
	// The damage is evolving after the crack has initiated
	else
	{	// flag for reached critical strain
		bool criticalStrain = false;
		
		// normal stress and strains
		double relStrength = soft[RELATIVE_STRENGTH];
		double sigmaI = relStrength*initiationLaw->sigmaI();
		double relToughness = soft[RELATIVE_TOUGHNESS];
		double scaleI = relToughness*soft[GCSCALING]/sigmaI;
		
		// increments in tensile and shear dissipation
		double dispEnergyN=0.,dispEnergyXY=0.,dispEnergyXZ=0.;
		
		// independent tensile damage (2D and 3D)
		bool hasDecxx = false;
		if(tractionFailureSurface<OVOID_SURFACE)
		{	hasDecxx = SoftenAxis(mptr,den,soft,DELTANORMAL,DAMAGENORMAL,sigmaI,scaleI,softeningModeI,
								  C11,str.xx,-1.,-1.,decxx,dispEnergyN,criticalStrain);
		}
		if(!hasDecxx)
		{	// Elastic get elastic decxx
			// If tensile not done above, this gets trial decxx for use in pressure effect next
			decxx = soft[DAMAGENORMAL]*den;
			if(den<0.)
			{	// keep ecxx>=0
				if(ecrack->xx+decxx<0.) decxx = -ecrack->xx;
			}
			else if(str.xx<0.)
			{	// start in compression and den>0 - might go into tension
				// reach zero in den1=-str.xx/C11, which is positive
                // zero if den < den1 to reach zero, otherwise subtract part to reach zero
                double den1 = -str.xx/C11;
                decxx = den<den1 ? 0. : soft[DAMAGENORMAL]*(den-den1);
			}
		}
		
		// shear stress and strains (for variable toughness need to call softening law to get new relToughness values)
        double relShearStrength = relStrength;
		double sigmaIIAlpha;
		double sigmaII = initiationLaw->sigmaII(relShearStrength,sigmaIIAlpha,str,&alpha,C11,nuterm,decxx,np);
		double relShearToughness = relToughness;
        
        // scale values for softening law
        // I think this will handle alpha dependent Gc by scaleII *= Gc(alpha)/Gc
        //   and scaleIIAlpha *= Gc(alpha+dalpha)/Gc
		double scaleII = relShearToughness*soft[GCSCALING]/sigmaII;
		double scaleIIAlpha = sigmaIIAlpha>0. ? relShearToughness*soft[GCSCALING]/sigmaIIAlpha : -1.;
		
		if(is2D)
		{	// cuboid and cylinder both use 1D methods for shear update
			// ovoid also uses it if in compression (but decxx was already set above for this case)
            // do not need to iterate pressure effects because shear update does not change decxx
			if(tractionFailureSurface<OVOID_SURFACE)
			{	// 2D shear strain only for x-y shear
				double dgs = str.xy>0. ? deij.xy : -deij.xy;		// * sign(tauxy)
				if(SoftenAxis(mptr,dgs,soft,DELTASHEAR,DAMAGESHEAR,sigmaII,scaleII,softeningModeII,
								Gshear,fabs(str.xy),sigmaIIAlpha,scaleIIAlpha,dgcxy,dispEnergyXY,criticalStrain))
				{	// Adjust sign to match shear stress direction
					if(str.xy<0.) dgcxy = -dgcxy;
				}
				else
				{	// elastic loading
					dgcxy = soft[DAMAGESHEAR]*deij.xy;
				}
			}
			
			// Ovoid 2D coupled shear and tension to a single damage parameter
			else
            {   int pass=1;
                double decxxPrevious = decxx;
                while(true)
                {   DamageState h;
                    h.deltaN = soft[DELTANORMAL];
                    h.deltaS = soft[DELTASHEAR];
                    h.damageN = soft[DAMAGENORMAL];
                    h.damageS = soft[DAMAGESHEAR];
                    h.ecxx = ecrack->xx;
					if(!OvoidSoftening(mptr,is2D,den,deij.xy,deij.xz,&h,str,C11,Gshear,
													  sigmaI,scaleI,sigmaII,scaleII,sigmaIIAlpha,scaleIIAlpha,
													  dispEnergyN,dispEnergyXY,criticalStrain))
					{	// elastic loading, decxx was set above
                        dgcxy = soft[DAMAGESHEAR]*deij.xy;
                        
                        // particle history
                        soft[DELTASHEAR] = h.deltaS;
                    }
                    else
                    {   if(sigmaIIAlpha>0. && !criticalStrain && pass<maxOvoidPasses)
                        {   // pressure dependent and subcritical: is it done?
                            double deltaMax = softeningModeI->GetDeltaMax(scaleI);
                            if(fabs(decxxPrevious-h.decxx)/deltaMax>pdOvoidTolerance)
                            {   // pressure dependence not converged yet
                                decxxPrevious = h.decxx;
                                sigmaII = initiationLaw->sigmaII(relShearStrength,sigmaIIAlpha,str,&alpha,C11,nuterm,h.decxx,np);
                                scaleII = relShearToughness*soft[GCSCALING]/sigmaII;
                                scaleIIAlpha = sigmaIIAlpha>0. ? relShearToughness*soft[GCSCALING]/sigmaIIAlpha : -1.;
                                
                                // repeat calculation
                                pass++;
                                continue;
                            }
                        }
                        
                        // cracking strains
                        decxx = h.decxx;
                        dgcxy = h.dgcxy;
                        
                        // particle history
                        soft[DELTANORMAL] = h.deltaN;
                        soft[DELTASHEAR] = h.deltaS;
                        soft[DAMAGENORMAL] = h.damageN;
                        soft[DAMAGESHEAR] = h.damageS;
                    }
                    
                    // if get here without continue, then done
                    break;
                }
			}
		}
		else if(tractionFailureSurface == CUBOID_SURFACE)
		{	// 3D shear strains using rectangular failure surface
			
			// xy direction
			double dgs = str.xy>0. ? deij.xy : -deij.xy;		// * sign(tauxy)
			if(SoftenAxis(mptr,dgs,soft,DELTASHEAR,DAMAGESHEAR,sigmaII,scaleII,softeningModeII,
						  Gshear,fabs(str.xy),sigmaIIAlpha,scaleIIAlpha,dgcxy,dispEnergyXY,criticalStrain))
			{	// Adjust sign to match shear stress direction
				if(str.xy<0.) dgcxy = -dgcxy;
			}
			else
			{	// elastic loading
				dgcxy = soft[DAMAGESHEAR]*deij.xy;
			}
			
			// xz direction
			dgs = str.xz>0. ? deij.xz : -deij.xz;		// * sign(tauxz)
			if(SoftenAxis(mptr,dgs,soft,DELTASHEAR2,DAMAGESHEAR2,sigmaII,scaleII,softeningModeII,
						  Gshear,fabs(str.xz),sigmaIIAlpha,scaleIIAlpha,dgcxz,dispEnergyXZ,criticalStrain))
			{	// Adjust sign to match shear stress direction
				if(str.xz<0.) dgcxz = -dgcxz;
			}
			else
			{	// elastic loading
				dgcxz = soft[DAMAGESHEAR2]*deij.xz;
			}
		}
		else if(tractionFailureSurface == CYLINDER_SURFACE)
		{	// 3D coupled shear damage
			if(!CoupledShearSoftening(deij.xy,deij.xz,soft,str,Gshear,sigmaII,scaleII,
									  sigmaIIAlpha,scaleIIAlpha,dgcxy,dgcxz,dispEnergyXY,criticalStrain))
			{	// elastic loading
				dgcxy = soft[DAMAGESHEAR]*deij.xy;
				dgcxz = soft[DAMAGESHEAR]*deij.xz;
			}
		}
		else
        {   int pass=1;
            double decxxPrevious = decxx;
            while(true)
            {   DamageState h;
                h.deltaN = soft[DELTANORMAL];
                h.deltaS = soft[DELTASHEAR];
                h.damageN = soft[DAMAGENORMAL];
                h.damageS = soft[DAMAGESHEAR];
                h.ecxx = ecrack->xx;
                if(!OvoidSoftening(mptr,is2D,den,deij.xy,deij.xz,&h,str,C11,Gshear,
                                   sigmaI,scaleI,sigmaII,scaleII,sigmaIIAlpha,scaleIIAlpha,
                                   dispEnergyN,dispEnergyXY,criticalStrain))
                {	// elastic loading, decxx for elastic was set above
                    dgcxy = soft[DAMAGESHEAR]*deij.xy;
                    dgcxz = soft[DAMAGESHEAR]*deij.xz;
                    
                    // particle history
                    soft[DELTASHEAR] = h.deltaS;
                }
                else
                {   if(sigmaIIAlpha>0. && !criticalStrain && pass<maxOvoidPasses)
                    {   // pressure dependent and subscritical: is it done?
                        double deltaMax = softeningModeI->GetDeltaMax(scaleI);
                        if(fabs(decxxPrevious-h.decxx)/deltaMax>pdOvoidTolerance)
                        {   // pressure dependence not converged yet
                            decxxPrevious = h.decxx;
                            sigmaII = initiationLaw->sigmaII(relShearStrength,sigmaIIAlpha,str,&alpha,C11,nuterm,h.decxx,np);
                            scaleII = relShearToughness*soft[GCSCALING]/sigmaII;
                            scaleIIAlpha = sigmaIIAlpha>0. ? relShearToughness*soft[GCSCALING]/sigmaIIAlpha : -1.;
                            
                            // repeat calculation
                            pass++;
                            continue;
                        }
                    }
                    
                    // cracking strains
                    decxx = h.decxx;
                    dgcxy = h.dgcxy;
                    dgcxz = h.dgcxz;
                    
                    // particle history
                    soft[DELTANORMAL] = h.deltaN;
                    soft[DELTASHEAR] = h.deltaS;
                    soft[DAMAGENORMAL] = h.damageN;
                    soft[DAMAGESHEAR] = h.damageS;
                }
                
                // if get here without continue, then done
                break;
            }
		}

#pragma mark ...... Dissipated Energy and Check for Final Failure
		// plastic (if there) and damage dissipated energy
		dispEnergy += dispEnergyN + dispEnergyXY + dispEnergyXZ;
		
        // Use energy release rate to check on failure
		// Get Gbar/Gcbar
		double relGI,relGII,relGIIxz=0.,cutoff=0.;
		if(tractionFailureSurface == CUBOID_SURFACE && !is2D)
		{	// only 3D Cuboid uses softening law (not correct if pressure dependent properties)
			
			// Get GIbar/GIcbar
			relGI = softeningModeI->GetGoverGc(soft[DELTANORMAL],scaleI);
			
			// Get (GIIbar/GIIbarc)
			relGII = softeningModeII->GetGoverGc(soft[DELTASHEAR],scaleII);
			
			// Get (GIIxzbar/GIIbarxzc)
			relGIIxz = softeningModeII->GetGoverGc(soft[DELTASHEAR2],scaleII);
			
			// mixed-mode failure criterion
			cutoff = pow(relGI,nmix) + pow(relGII,nmix) + pow(relGIIxz,nmix);
		}
		else
		{	// for isotropic material, two delta and two d (cuboid (2D only) and cylinder) or one d (ovoid)
			// store GI and GII
			soft[DAMAGEGI] += dispEnergyN/soft[GCSCALING];
			soft[DAMAGEGII] += dispEnergyXY/soft[GCSCALING];
			
			// Get GIbar/GIcbar
			relGI = soft[DAMAGEGI]/(relToughness*softeningModeI->GetGc());
				
			// Get (GIIbar/GIIbarc)
			relGII = soft[DAMAGEGII]/(relShearToughness*softeningModeII->GetGc());
			
			// mixed mode failure criterion
			if(tractionFailureSurface < OVOID_SURFACE)
				cutoff = pow(relGI,nmix) + pow(relGII,nmix);
		}

		// criterion
		
		// check limit or energy
		if(criticalStrain || cutoff>=1.)
		{	// report energy released
            double GI,GII;
            if(tractionFailureSurface == CUBOID_SURFACE && !is2D)
            {   GI = relToughness*relGI*softeningModeI->GetGc()*UnitsController::Scaling(1.e-3);
                GII = relShearToughness*(relGII+relGIIxz)*softeningModeII->GetGc()*UnitsController::Scaling(1.e-3);
            }
            else
            {   GI = soft[DAMAGEGI]*UnitsController::Scaling(1.e-3);
                GII = soft[DAMAGEGII]*UnitsController::Scaling(1.e-3);
            }
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
						", ds=" << (soft[DELTASHEAR]/(soft[DELTASHEAR]+relStrength*gs0*softeningModeII->GetFFxn(soft[DELTASHEAR],scaleII))) <<
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
			// otherwise DAMAGESHEAR2=DAMAGEGII is used to store mode II energy
			if(tractionFailureSurface == CUBOID_SURFACE && !is2D)
				soft[DAMAGESHEAR2] = 1.;

			// post failure update
			PostFailureUpdate(decxx,dgcxy,dgcxz,&dsig,&str,ecrack,Rtot,C11,den,Gshear,deij.xy,Gshear,deij.xz,is2D,frictionCoeff);
		}
		else
		{	// crack plane stress increment
			dsig.xx = C11*(den - decxx);
			dsig.xz = Gshear*(deij.xz-dgcxz);
			dsig.xy = Gshear*(deij.xy-dgcxy);
		}
		
		// dissipated energy from damage (plasticity was added before if it occurred)
        // prior to 10/15/2020, plastic energy was added here too might have double counted
        // 2D tracks damage only in NORMALDIR3 history; subtract it from plastic to get plastic energy
        double dispEnergyDam = dispEnergyN + dispEnergyXY + dispEnergyXZ;
        if(is2D) soft[NORMALDIR3] += UnitsController::Scaling(1.e-9)*mptr->mp*dispEnergyDam;
		mptr->AddPlastEnergy(dispEnergyDam);
	}
		
#pragma mark ... Update Cracking Strain and Stresses
	// Have now found increments in cracking strain and stress
	// update cracking strains (in plastic strain)
	// current was rotated by dR (above) now add increment rotated by Rtot
	
	// rotate previous stress by rotation increment (this switch from smoothed stress, if used, to particle stres now)
	Tensor *sp = mptr->GetStressTensor();
	str = dR.RVoightRT(sp,true,is2D);

#ifdef POROELASTICITY
	// Poroelasticity, needs to add increment due to pp expansion terms
	if(DiffusionTask::HasPoroelasticity())
	{	// Get rho0 for specific Kirchoff stress
		int matid = mptr->MatID();
		double rho0=theMaterials[matid]->GetRho(mptr);
		
		// dsig.xx done above, but now needs pore pressure increment
		dsig.xx += alphaPE*res->dC/rho0;
		
		// substract pore pressure term now for effective increments
		// these changes will get pore pressure added to dsig.yy and dsig.zz below
		double ppadd = CME1*res->dC;
		deij.xx -= ppadd;
		deij.yy -= ppadd;
		eres += ppadd;
		
		// Add increments to eres and deij.zz (only used in 3D or axisymmetric)
		ppadd = CME3*res->dC;
		ezzres += ppadd;
		deij.zz -= ppadd;
	}
#endif
	
	// get stress incement for other stresses in crack axis system
	if(np==THREED_MPM)
	{	dsig.yy = C11*(deij.yy + nuterm*(deij.xx - decxx + deij.zz));
		dsig.zz = C11*(deij.zz + nuterm*(deij.xx - decxx + deij.yy));
		dsig.yz = Gshear*deij.yz;
	}
	else if(np==AXISYMMETRIC_MPM)
	{	dsig.yy = C11*(deij.yy + nuterm*(deij.xx - decxx + deij.zz));
		dsig.zz = C11*(deij.zz + nuterm*(deij.xx - decxx + deij.yy));
	}
	else
	{	dsig.yy = C11*(deij.yy + nuterm*(deij.xx - decxx));
		if(np==PLANE_STRAIN_MPM)
			dsig.zz = C11*(nuterm*(de.xx - decxx + de.yy - 2.*ezzres) - ezzres);
		else if(np==PLANE_STRESS_MPM)
		{	// zz deformation. zz stress stays at zero for plane stress
			// Here C[4][1] = C[4][2] = -v/(1-v)
			de.zz = p->C[4][1]*(deij.xx-decxx) + p->C[4][2]*deij.yy + eres;
			mptr->IncrementDeformationGradientZZ(de.zz);
		}
	}
	
	// update stress and strain (in which cracking strain rotated to current configuration)
	UpdateCrackingStrain(np,mptr->GetAltStrainTensor(),decxx,dgcxy,dgcxz,Rtot,soft);
	UpdateCrackingStress(np,sp,&dsig,&str,Rtot);
	
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
	// But do not inlude increment in cracking strain
	double delV = de.xx-decxx+de.yy+de.zz;
	double dTq0 = -gamma0*mptr->pPreviousTemperature*delV;
	
	// track heat energy with total dissipated energy
	IncrementHeatEnergy(mptr,dTq0,dispEnergy);
	
#ifdef POROELASTICITY
	UndrainedPressIncrement(mptr,delV);
#endif
}

// Calculate rotation matrix from crack to initial
bool IsoSoftening::GetRToCrack(Matrix3 *R,double *soft, bool is2D, int Dstyle) const
{	// none if undamaged
	if (soft[SOFT_DAMAGE_STATE] < predamageState) return false;

	// 3D or 2D
	if(!is2D)
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

// Cylindrical shear softening in an isotropic material
// Note this methods is only called in 3D modeling
bool IsoSoftening::CoupledShearSoftening(double dgxy,double dgxz,double *history,Tensor &str,double Gshear,
										 double sigma,double scale,double sigmaAlpha,double scaleAlpha,
										 double &dgcxy,double &dgcxz,double &dispEnergy,bool &criticalStrain) const
{
	// no damage if negative here and no change in alpha variables
	double Tdotdg = str.xy*dgxy + str.xz*dgxz;
	if(Tdotdg<0. && scaleAlpha<0.) return false;
	
	// get damage parameter, trial stress increment and stress
	double ds = history[DAMAGESHEAR];
	double G1md = Gshear*(1.-ds);
	double dTxyTrial = G1md*dgxy;
	double TxyTrial = str.xy + dTxyTrial;
	double dTxzTrial = G1md*dgxz;
	double TxzTrial = str.xz + dTxzTrial;
	
	// Elastic or damaging?
	double delta = history[DELTASHEAR];
	double Fs = sigma*softeningModeII->GetFFxn(delta,scale);
	
	// look for alpha variables
	if(sigmaAlpha>0.)
	{	// switch to delta+ddeltaElastic, alpha+dAlpha on damage surface
        double ddeltaElastic = softeningModeII->GetDDeltaElastic(delta,sigmaAlpha,scaleAlpha,ds,G1md);
		
		// update strength and delta
		sigma = sigmaAlpha;
		scale = scaleAlpha;
		delta += ddeltaElastic;
		
		// evaluate softening law at new delta and alpha+dalpha
		Fs = sigma*softeningModeII->GetFFxn(delta,scale);
	}
	
	// If no damage, update delta (if needed) and return false
	if(TxyTrial*TxyTrial + TxzTrial*TxzTrial <= Fs*Fs)
	{	history[DELTASHEAR] = delta;			// in case changed
		return false;
	}
	
	// Find scaling of traction to move to failure surface
	// (at delta+ddeltaElast, alpha+dalpha) in general theory
	double Tsmag = sqrt(str.xy*str.xy + str.xz*str.xz);
	double phi = (Fs-Tsmag)/G1md;
	double Thatxy = str.xy/Tsmag;
	double Thatxz = str.xz/Tsmag;
	double dgxy2 = dgxy - Thatxy*phi;
	double dgxz2 = dgxz - Thatxz*phi;
	
	// Find effective damage increment
	double dgs = Thatxy*dgxy2 + Thatxz*dgxz2;
	
	// increment in delta strain (note: returns < 0 if decohesion)
	// note that when sigmaAlpha>0, sigma, scale, and delta are
	//			changed to sigmaAlpha, scaleAlpha, and delta+ddeltaElastic
	double en0 = sigma/Gshear;				 // initiation strain
	double ddelta = softeningModeII->GetDDelta(dgs,en0,delta,ds,scale);
	
	// check if failed
	double deltaPrevious = delta;			// = delta + ddeltaelastic if occurred
	if(ddelta<0.)
	{	// Decohesion: set delta=deltaMax, d=1, and set criticalStrain flag true
		// caller will get cracking strain in post-failure update code
		delta = softeningModeII->GetDeltaMax(scale);
		criticalStrain = true;
		history[DAMAGESHEAR] = 1.;
	}
	else
	{	// damage propagation
		delta += ddelta;
		double fval = softeningModeII->GetFFxn(delta,scale);
		history[DAMAGESHEAR] = delta/(delta+en0*fval);
		
		// Cracking strain increments are sum of elastic part and damage part
		dgcxy = ds*(dgxy - Thatxy*dgs) + Thatxy*ddelta;
		dgcxz = ds*(dgxz - Thatxz*dgs) + Thatxz*ddelta;
	}
	
	// dissipated energy increment per unit mass using dOmega = (1/2) varphi (ddelta-ddeltaElastic)
	// note that deltaPrevious = deltai + ddeltaelastic, so deltaf-deltaPrevious = ddeltaf-ddelti-ddeltaElastic
	double dissipated = 0.5*sigma*softeningModeII->GetPhiFxn(deltaPrevious,scale)*(delta-deltaPrevious);
	
	// This is prior discrete calculation. It is identical to above for linear softening
	//double dissipated = sigma*(softLaw->GetGToDelta(delta,scale) - softLaw->GetGToDelta(deltaPrevious,scale));
	dispEnergy += dissipated;
	
	// new delta paramater
	history[DELTASHEAR] = delta;
	
	return true;
}

// Ovoid damafe evolution surface for isotropic materials
bool IsoSoftening::OvoidSoftening(MPMBase *mptr,bool is2D,double den,double dgxy,double dgxz,
        DamageState *h,Tensor &str,double C11,double Gshear,
        double sigmaI,double scaleI,double sigmaII,double scaleII,double sigmaIIAlpha,double scaleIIAlpha,
        double &dispIEnergy,double &dispIIEnergy,bool &criticalStrain) const
{
	// switch for a different surface
	if(tractionFailureSurface==COUPLED_CUBOID_SURFACE)
	{	return CoupledCuboidSoftening(mptr,is2D,den,dgxy,dgxz,h,str,C11,Gshear,
									  sigmaI,scaleI,sigmaII,scaleII,sigmaIIAlpha,scaleIIAlpha,
									  dispIEnergy,dispIIEnergy,criticalStrain);
		
	}
	
	// get parameter (dn and ds are coupled)
	double d = h->damageN;
	
	// get parameters (deln and dels may differ)
	double deln = h->deltaN;
	double dels = h->deltaS;
	
	// trial stress and increments
	double C111md = C11*(1.-d);
	double dTntrial = C111md*den;
	double G1md = Gshear*(1.-d);
	double dTxytrial = G1md*dgxy;
	double dTxztrial = is2D ? 0. : G1md*dgxz;
	
	// softening laws
	double FI = sigmaI*softeningModeI->GetFFxn(deln,scaleI);
	double FII = sigmaII*softeningModeII->GetFFxn(dels,scaleII);
	
	// look for alpha variables (only allowed in sigmaII for now)
	// Currently shear only
	double FIIstar = FII,ddeltaSElastic=0.;
	if(sigmaIIAlpha>0.)
    {   // switch to delta+ddeltaSElastic, alpha+dAlpha on damage surface
		ddeltaSElastic = softeningModeII->GetDDeltaElastic(dels,sigmaIIAlpha,scaleIIAlpha,d,G1md);
		dels += ddeltaSElastic;
        
		// evaluate softening law at new deltas and alpha+dalpha
		FIIstar = sigmaIIAlpha*softeningModeII->GetFFxn(dels,scaleIIAlpha);
	}
	
	// test for failure
	double stn = str.xx/FI;
	double dstn = dTntrial/FI;
	// sts2 = (|Ts|/FI)^2 and and ststr2 = (|Tstrial|/FIIstar)^2
	double sts2,ststr2;
	if(is2D)
	{	sts2 = str.xy/FIIstar;
		ststr2 = (str.xy+dTxytrial)/FIIstar;
		sts2 *= sts2;
		ststr2 *= ststr2;
	}
	else
	{	double FII2 = FIIstar*FIIstar;
		sts2 = (str.xy*str.xy + str.xz*str.xz)/FII2;
		double Txytrial = str.xy+dTxytrial;
		double Txztrial = str.xz+dTxztrial;
		ststr2 = (Txytrial*Txytrial + Txztrial*Txztrial)/FII2;
	}
	
	// test tension (ststr2->(Tntrial/Fn)^2+(|Tstrial|/FIIstar)^2)
	//			or compression (ststr2 stays (|Tstrial|/FIIstar)^2) for failure
	// If no damage, update delta (if needed) and return false
	double TtrialOverFI = stn+dstn;
	if(TtrialOverFI>0.) ststr2 += TtrialOverFI*TtrialOverFI;
	if(ststr2 <= 1.)
	{	h->deltaS = dels;			// in case changed
		return false;
	}
	
	// Tensile half plane
	double ddeltan,ddeltas,dd;
	if(TtrialOverFI>0.)
	{	// Find phi required to reach the failure envelop at alpha along initiation traction
		double Tmag,Thatn,Thatxy,Thatxz;
		if(is2D)
		{	Tmag = sqrt(str.xx*str.xx + str.xy*str.xy);
			Thatxz = 0.;
		}
		else
		{	Tmag = sqrt(str.xx*str.xx + str.xy*str.xy + str.xz*str.xz);
			Thatxz = str.xz/Tmag;
		}
		Thatn = str.xx/Tmag;
		Thatxy = str.xy/Tmag;
		double Thatn2Diff = 1. - Thatn*Thatn;
		double onePlusPhiTmag = FI*FII/sqrt(FII*FII*Thatn*Thatn + FI*FI*Thatn2Diff);
//#pragma omp critical (output)
//		{
//			cout << "# Tu=(" << Thatxy << "," << Thatxz << "," << Thatn << ") test=" << ststr2
//			<< ", phi=" << (onePlusPhiTmag-Tmag)/Tmag << endl;
//			cout << "#    1-Tn^2=" << Thatn2Diff << ", Txy^2=" << Thatxy*Thatxy
//			<< ", pre=" << str.xx*str.xx/(FI*FI) + str.xy*str.xy/(FIIstar*FIIstar) << endl;
//		}

		// Update traction to be traction just prior to damage step in (Tx,Ty,Tz)
		// Get de2 causing damage
		double Txy,Txz=0.,Tn;
		double den2,dgxy2,dgxz2=0.;
		double phiTmag = onePlusPhiTmag-Tmag;
		if(sigmaIIAlpha>0.)
		{	double FsRatio = onePlusPhiTmag*FIIstar/FII;
			Txy = Thatxy*FsRatio;
			if(!is2D) Txz = Thatxz*FsRatio;
			Tn = Thatn*Tmag;
			
			// damage strain
			FsRatio -= Tmag;
			dgxy2 = dgxy - FsRatio*Thatxy/G1md;
			if(!is2D) dgxz2 = dgxz - FsRatio*Thatxz/G1md;
			den2 = den - phiTmag*Thatn/C111md;
			
			// switch to alpha adjusted properties
			FII = FIIstar;
			sigmaII = sigmaIIAlpha;
			scaleII = scaleIIAlpha;
		}
		else
		{	Txy = onePlusPhiTmag*Thatxy;
			Tn = onePlusPhiTmag*Thatn;
			
			dgxy2 = dgxy - phiTmag*Thatxy/G1md;
			den2 = den - phiTmag*Thatn/C111md;
			
			if(!is2D)
			{	Txz = onePlusPhiTmag*Thatxz;
				dgxz2 = dgxz - phiTmag*Thatxz/G1md;
			}
		}

		// get ovoid traction vectors (magnitude should be one)
		//double Tomag = sqrt((Txy/FII)*(Txy/FII) + (Txz/FII)*(Txz/FII) + (Tn/FI)*(Tn/Fn));
		double Tohatxy = Txy/FII;
		double Tohatxz = 0.;
		double Tohats2 = Tohatxy*Tohatxy;
		double Tdotdgamma = Tohatxy*dgxy2;
		double Tohatn = Tn/FI;
		double Tohatn2 = Tohatn*Tohatn;
		if(!is2D)
		{	Tohatxz = Txz/FII;
			Tohats2 += Tohatxz*Tohatxz;
			Tdotdgamma += Thatxz*dgxz2;
		}
		
		// uniaxial hack
		//Tohatxy=Tohatxz=Tohats2=Tdotdgamma=0.;
		//Tohatn=Tohatn2=1.;
		
		double Fsprime = 1. + sigmaII*softeningModeII->GetFpFxn(dels,scaleII)/Gshear;
		double Fnprime = 1. + sigmaI*softeningModeI->GetFpFxn(deln,scaleI)/C11;
		double Rs = softeningModeII->GetRdFxn(dels,scaleII,sigmaII/Gshear);
		double Rn = softeningModeI->GetRdFxn(deln,scaleI,sigmaI/C11);
//#pragma omp critical (output)
//		{
//			cout << "# To=(" << Tohatxy << "," << Tohatxz << "," << Tohatn << ") mag="
//			<< (Tohats2+Tohatn) << ", Fsp=" << Fsprime << ", Fnp=" << Fnprime
//			<< ", Rs=" << Rs << ", Rn=" << Rn << endl;
//			cout << "#    en0=" << sigmaI/C11 << ", si=" << sigmaI << ", C11=" << C11 <<
//			", G=" << Gshear << endl;
//		}

		// explicit solution
		double Gnorm = Gshear/FII,C11Norm = C11/FI;
		
		// Find which parameters (d, deltan, or deltas) has the smallet fraction increment
		double deltanMax = softeningModeI->GetDeltaMax(scaleI);
		double nRate = deltanMax*Rn;
		double deltasMax = softeningModeII->GetDeltaMax(scaleII);
		double sRate = deltasMax*Rs;
		double newD;
		if(fmax(nRate,sRate)<1.)
		{	// increment d if both 1/(delta rates) are large (this is approaching decohesion)
			dd = ( Gnorm*Tdotdgamma + C11Norm*Tohatn*den2 ) /
					( (Gnorm/Rs)*Tohats2*Fsprime + (C11Norm/Rn)*Tohatn2*Fnprime );
			newD = d+dd;
			if(newD >= MAXIMUM_D)
			{	// failed is handled below
				newD = 1.;
			}
			else
			{	// uppdate delta's and get their increments
				h->deltaN = softeningModeI->GetDeltaFromDamage(newD,scaleI,sigmaI/C11,deln);
                ddeltan = h->deltaN - deln;
				h->deltaS = softeningModeII->GetDeltaFromDamage(newD,scaleII,sigmaII/Gshear,dels);
                ddeltas = h->deltaS - dels;
			}
		}
		else if(nRate>sRate)
		{	// increment deltan
			double RnRs = Rn/Rs;
			ddeltan = ( Gnorm*Tdotdgamma + C11Norm*Tohatn*den2 ) /
						( (RnRs*Gnorm)*Tohats2*Fsprime + C11Norm*Tohatn2*Fnprime );
			double newDeln = deln+ddeltan;
			if(newDeln>=deltanMax)
			{	// failed is handled below
				newD = 1.;
			}
			else
			{	// uppdate delta's and get their increments
				h->deltaN = newDeln;
				double fI = softeningModeI->GetFFxn(newDeln,scaleI);
				newD = newDeln/(newDeln+sigmaI*fI/C11);
				h->deltaS = softeningModeII->GetDeltaFromDamage(newD,scaleII,sigmaII/Gshear,dels);
                ddeltas = h->deltaS - dels;
			}
		}
		else
		{	// increment deltas
			double RsRn = Rs/Rn;
			ddeltas = ( Gnorm*Tdotdgamma + C11Norm*Tohatn*den2 ) /
						( Gnorm*Tohats2*Fsprime + RsRn*C11Norm*Tohatn2*Fnprime );
			double newDels = dels+ddeltas;
			if(newDels>=deltasMax)
			{	// failed is handled below
				newD = 1.;
			}
			else
			{	// uppdate delta's and get their increments
				h->deltaS = newDels;
				double fII = softeningModeII->GetFFxn(newDels,scaleII);
				newD = newDels/(newDels+sigmaII*fII/Gshear);
				h->deltaN = softeningModeI->GetDeltaFromDamage(newD,scaleI,sigmaI/C11,deln);
                ddeltan = h->deltaN - deln;
			}
		}
		// set d parameters to updated values
		h->damageN = newD;
		h->damageS = newD;
		
		// has it failed?
		if(newD >= MAXIMUM_D)
		{	// Decohesion: set delta=deltaMax, d=1, and set criticalStrain flag true
			// caller will get cracking strain in post-failure update code
			h->deltaN = deltanMax;
			h->deltaS = deltasMax;
            ddeltan = deltanMax-deln;
            ddeltas = deltasMax-dels;
			criticalStrain = true;
		}
		else
		{	// Cracking strain increments are sum of elastic part and damage part
			h->dgcxy = d*(dgxy - Tohatxy*Fsprime*ddeltas) + Tohatxy*ddeltas;
			h->decxx = d*(den - Tohatn*Fnprime*ddeltan) + Tohatn*ddeltan;
			if(!is2D)
				h->dgcxz = d*(dgxz - Tohatxz*Fsprime*ddeltas) + Tohatxz*ddeltas;
		}
		

		// dissipated energy increment per unit mass using dOmega = (1/2) varphi (ddelta-ddeltaElastic)
		// note that deltaPrevious = deltai + ddeltaelastic, so deltaf-deltaPrevious = ddeltaf-ddelti-ddeltaElastic
        // no need to add because only one softening method is called
        dispIEnergy = 0.5*Tohatn2*sigmaI*softeningModeI->GetPhiFxn(deln,scaleI)*fmax(ddeltan,0.);
        dispIIEnergy = 0.5*Tohats2*sigmaII*softeningModeII->GetPhiFxn(dels,scaleII)*fmax(ddeltas,0.);
	}
	
	else
	{	// Compression half plane. Reupdate shear damage only, but apply that too normal as well in the coupling
		
		// check for pressure dependent shear
		if(sigmaIIAlpha>0.)
		{	// switch to alpha adjusted properties
			FII = FIIstar;
			sigmaII = sigmaIIAlpha;
			scaleII = scaleIIAlpha;
		}
		
		// Find portion of the step that causes damage
		double dgs2,Thatxy,Thatxz;
		if(is2D)
		{	dgs2 = (str.xy + dTxytrial - FII)/G1md;
		}
		else
		{	double Tsmag = sqrt(str.xy*str.xy + str.xz*str.xz);
			double phi = (FII-Tsmag)/G1md;
			Thatxy = str.xy/Tsmag;
			Thatxz = str.xz/Tsmag;
			double dgxy2 = dgxy - Thatxy*phi;
			double dgxz2 = dgxz - Thatxz*phi;
			
			// Find effective damage increment
			dgs2 = Thatxy*dgxy2 + Thatxz*dgxz2;
		}
		
		// increment in delta strain (note: returns < 0 if decohesion)
		// note that when sigmaAlpha>0, sigma, scale, and delta are
		//			changed to sigmaAlpha, scaleAlpha, and delta+ddeltaElastic
		double gs0 = sigmaII/Gshear;
        // This used softeningModeI unti 9/30/2020
		ddeltas = softeningModeII->GetDDelta(dgs2,gs0,dels,d,scaleII);

		double newD;
		
		// has it failed?
		if(ddeltas<0.)
		{	// Decohesion: set delta=deltaMax, d=1, and set criticalStrain flag true
			// caller will get cracking strain in post-failure update code
			h->deltaN = softeningModeI->GetDeltaMax(scaleI);
            ddeltan = fmax(h->deltaN-deln,0.);
			h->deltaS =  softeningModeII->GetDeltaMax(scaleII);
            ddeltas = fmax(h->deltaS-dels,0.);
			criticalStrain = true;
			newD = 1.;
		}
		else
		{	// get shear direction
			double Rs = softeningModeII->GetRdFxn(dels,scaleII,gs0);
			double deltasMax = softeningModeII->GetDeltaMax(scaleII);
			double sRate = deltasMax*Rs;
			if(sRate<1.)
			{	// increment d if 1/(deltas rate) is large (this is approaching decohesion)
				newD = d + Rs*ddeltas;
				if(newD>MAXIMUM_D)
				{	// Decohesion: set delta=deltaMax, d=1, and set criticalStrain flag true
					// caller will get cracking strain in post-failure update code
					h->deltaN = softeningModeI->GetDeltaMax(scaleI);
                    ddeltan = fmax(h->deltaN-deln,0.);
					h->deltaS =  softeningModeII->GetDeltaMax(scaleII);
                    ddeltas = fmax(h->deltaS-dels,0.);
					criticalStrain = true;
					newD = 1.;
				}
				else
				{	h->deltaS = softeningModeII->GetDeltaFromDamage(newD,scaleII,sigmaII/Gshear,dels);
                    ddeltas = fmax(h->deltaS - dels,0.);
					h->deltaN = softeningModeI->GetDeltaFromDamage(newD,scaleI,sigmaI/C11,deln);
                    ddeltan = fmax(h->deltaN - deln,0.);
				}
			}
			else
			{	// increment deltas (known to be below deltasMax)
				double newDels = dels+ddeltas;
				h->deltaS = newDels;
				double fII = softeningModeII->GetFFxn(newDels,scaleII);
				newD = newDels/(newDels+sigmaII*fII/Gshear);
				h->deltaN = softeningModeI->GetDeltaFromDamage(newD,scaleI,sigmaI/C11,deln);
                ddeltan = fmax(h->deltaN - deln,0.);
			}
			
			// Cracking strain increments are sum of elastic part and damage part
			if(!criticalStrain)
			{	if(is2D)
					h->dgcxy = ddeltas;
				else
				{	h->dgcxy = d*(dgxy - Thatxy*dgs2) + Thatxy*ddeltas;
					h->dgcxz = d*(dgxz - Thatxz*dgs2) + Thatxz*ddeltas;
				}
                h->decxx = -h->ecxx;        // to evolve to zero cracking strain
			}
		}
		
		// finally set updated damage parameters
		h->damageN = newD;
		h->damageS = newD;
		
		// dissipated energy increment per unit mass using dOmega = (1/2) phi (ddelta-ddeltaElastic)
		// note that deltaPrevious = deltai + ddeltaelastic, so deltaf-deltaPrevious = ddeltaf-ddelti-ddeltaElastic
        // no need to add because only one softening method is called
		dispIIEnergy = 0.5*sigmaII*softeningModeII->GetPhiFxn(dels,scaleII)*ddeltas;
	}

    return true;
}

// Coupled cuboid damage evolution surface for isotropic materials
bool IsoSoftening::CoupledCuboidSoftening(MPMBase *mptr,bool is2D,double den,double dgxy,double dgxz,
			  DamageState *h,Tensor &str,double C11,double Gshear,
			  double sigmaI,double scaleI,double sigmaII,double scaleII,double sigmaIIAlpha,double scaleIIAlpha,
			  double &dispIEnergy,double &dispIIEnergy,bool &criticalStrain) const
{
	// get parameter (all coupled)
	double d = h->damageN;
	
	// get parameters (deln and dels may differ)
	double deln = h->deltaN;
	double dels = h->deltaS;
	
	// trial stress and increments
	double C111md = C11*(1.-d);
	double dTntrial = C111md*den;
	double G1md = Gshear*(1.-d);
	double dTxytrial = G1md*dgxy;
	double dTxztrial = is2D ? 0. : G1md*dgxz;

	// softening laws
	double FI = sigmaI*softeningModeI->GetFFxn(deln,scaleI);
	double FII = sigmaII*softeningModeII->GetFFxn(dels,scaleII);
	
	// look for alpha variables (only allowed in sigmaII for now)
	// Pressure dependence not currently allowed in coupled cuboid surface

	// additions for 3D
	if(is2D)
	{   // check for elastic (caller gets cracking strains)
		if((str.xx+dTntrial<=FI) && (fabs(str.xy+dTxytrial)<=FII))
			return false;
	}
	else
	{   // check for elastic (caller gets cracking strains)
		if((str.xx+dTntrial<=FI) && (fabs(str.xy+dTxytrial)<=FII) && (fabs(str.xz+dTxztrial)<=FII))
			return false;
	}
	
	// damage has occurred in at least on direction
	double ddeltaN=0.,ddeltaXY=0.,ddeltaXZ=0.;
	double dDn=0.,dDxy=0.,dDxz=0.;
	double den2,dgxy2,dgxz2;
	double en0 = sigmaI/C11,gs0 = sigmaII/Gshear;
	criticalStrain = false;

	// ratios of current stress to strength (should alays be <= 1)
	double sfN=str.xx/FI,sfXY=str.xy/FII;
	double sfXZ = is2D ? 0. : str.xz/FII;
	
	// normal direction
	if(str.xx+dTntrial>=FI)
	{   // go to elastic surface and then damage (T(trial)-F)/(E(1-D))
		den2 = (str.xx+dTntrial-FI)/C111md;
		ddeltaN = softeningModeI->GetDDelta(den2,en0,deln,d,scaleI);
		if(ddeltaN<0.)
			criticalStrain = true;
		else
		{   // get change in Dn
			double delNew = deln+ddeltaN;
			double Fnew = softeningModeI->GetFFxn(delNew,scaleI);
			dDn = delNew/(delNew + en0*Fnew) - d;
		}
	}
	
	// shear xy direction
	if(!criticalStrain && fabs(str.xy+dTxytrial)>=FII)
	{   // go to elastic surface and then damage (T(trial)-F)/(E(1-D))
		dgxy2 = (fabs(str.xy+dTxytrial)-FII)/Gshear;
		ddeltaXY = softeningModeII->GetDDelta(dgxy2,gs0,dels,d,scaleII);
		if(ddeltaXY<0.)
			criticalStrain = true;
		else
		{   // get change in Dxy
			double delNew = dels+ddeltaXY;
			double Fnew = softeningModeII->GetFFxn(delNew,scaleII);
			dDxy = delNew/(delNew + gs0*Fnew) - d;
		}
	}
	
	// shear xz strain
	if(!criticalStrain && !is2D && fabs(str.xz+dTxztrial)>=FII)
	{   // go to elastic surface and then damage
		dgxz2 = (fabs(str.xz+dTxztrial)-FII)/Gshear;
		ddeltaXZ = softeningModeII->GetDDelta(dgxz2,gs0,dels,d,scaleII);
		if(ddeltaXZ<0.)
			criticalStrain = true;
		else
		{   // get change in Dxz
			double delNew = dels+ddeltaXZ;
			double Fnew = softeningModeII->GetFFxn(delNew,scaleII);
			dDxz = delNew/(delNew + gs0*Fnew) - d;
		}
	}
	
	// see if failed or find direction with most damage
	double dD=0.;
	bool damageN=false,damageXY=false,damageXZ=false;
	if(criticalStrain)
	{   // Decohesion: set delta=deltaMax, D=1
		// caller will get cracking strain in post-failure update code
		h->deltaN = softeningModeI->GetDeltaMax(scaleI);
		ddeltaN = h->deltaN - deln;
		h->deltaS = softeningModeII->GetDeltaMax(scaleII);
		ddeltaXY = h->deltaS - dels;
		ddeltaXZ = ddeltaXY;
		dD = 1.-d;
		damageN = damageXY = damageXZ = true;
		
		// limit to 1, seems rare and not sure how ever, but might help last energy calc
		sfN = fmin(1.,sfN);
		sfXY = fmin(1.,fabs(sfXY));
		if(!is2D) sfXZ = fmin(1.,fabs(sfXZ));
	}
	
	else if(dDn>=fmax(dDxy,dDxz))
	{   // normal direction damages most
		h->deltaN += ddeltaN;
		dD = dDn;
		h->decxx = d*(den-den2) + ddeltaN;
		sfN = 1.;
		damageN = true;
	}
	
	else if(dDxy>=dDxz)
	{   // shear x-y damages most
		h->deltaS += ddeltaXY;
		dD = dDxy;
		h->dgcxy = d*(fabs(dgxy)-dgxy2) + ddeltaXY;
		if(str.xy<0.) h->dgcxy = -h->dgcxy;
		sfXY = 1.;
		damageXY = true;
	}
	
	else if(!is2D)
	{   // shear x-z damages most
		h->deltaS += ddeltaXZ;
		dD = dDxz;
		h->dgcxz = d*(fabs(dgxz)-dgxz2) + ddeltaXZ;
		if(str.xz<0.) h->dgcxz = -h->dgcxz;
		sfXZ = 1.;
		damageXZ = true;
	}
	
	// do non-damaging directions
	if(!damageN)
	{   // one shear direction failed, get deln from the new D
		h->deltaN = softeningModeI->GetDeltaFromDamage(d+dD,scaleI,en0,deln);
		ddeltaN = h->deltaN - deln;
		
		// get cracking strain from initial values and stress at start of damage
		double Fnp = softeningModeI->GetFpFxn(deln,scaleI);
		h->decxx = d*den + sfN*(1.-d*(1+en0*Fnp))*ddeltaN;
		if(h->ecxx+h->decxx<0.) h->decxx = -h->ecxx;
		
		// find shear direction that did not fail (3D only)
		double Fs = softeningModeII->GetFpFxn(dels,scaleII);
		if(!is2D)
		{	if(!damageXY)
			{	ddeltaXY = h->deltaS - dels;
				h->dgcxy = d*dgxy + sfXY*(1.-d*(1+gs0*Fs))*ddeltaXY;
			}
			else
			{	ddeltaXZ = h->deltaS - dels;
				h->dgcxz = d*dgxz + sfXZ*(1.-d*(1+gs0*Fs))*ddeltaXZ;
			}
		}
	}
	else if(!criticalStrain)
	{	// Normal direction damaged, get dels from the new D
		h->deltaS = softeningModeII->GetDeltaFromDamage(d+dD,scaleII,gs0,dels);
		ddeltaXY = h->deltaS - dels;
		
		// get cracking strain from initial values and stress at start of damage
		double Fs = softeningModeII->GetFpFxn(dels,scaleII);
		h->dgcxy = d*dgxy + sfXY*(1.-d*(1+gs0*Fs))*ddeltaXY;
		if(!is2D)
		{	ddeltaXZ = ddeltaXY;
			h->dgcxz = d*dgxz + sfXZ*(1.-d*(1+gs0*Fs))*ddeltaXZ;
		}
	}
	
	// Get energy dissipation
	
	// dissipated energy increment per unit mass using dOmega = (1/2) (S/F)^2 phi ddelta
	if(str.xx>0.)
		dispIEnergy = 0.5*sfN*sfN*sigmaI*softeningModeI->GetPhiFxn(deln,scaleI)*ddeltaN;
	double modeII = sigmaII*softeningModeI->GetPhiFxn(dels,scaleII);
	dispIIEnergy = 0.5*sfXY*sfXY*modeII*ddeltaXY;
	if(!is2D) dispIIEnergy = 0.5*sfXZ*sfXZ*modeII*ddeltaXZ;

	return true;
}

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
	
#ifdef POROELASTICITY
	UndrainedPressIncrement(mptr,delV);
#endif
}

#pragma mark IsoSoftening::Accessors

// return material type
const char *IsoSoftening::MaterialType(void) const { return "Anisotropic Softening of Isotropic Material"; }

// store plastic strain in alt strain
int IsoSoftening::AltStrainContains(void) const
{	return ENG_BIOT_PLASTIC_STRAIN;
}

// return cracking strain and true, or return false if material has no cracking strain
bool IsoSoftening::GetCrackingStrain(MPMBase *mptr,Tensor *ecrack,bool is2D,Matrix3 *Rtot) const
{
	// get history
	double *soft = GetSoftHistoryPtr(mptr);
	
	// Before initiatiation, not damage strain
	if(soft[SOFT_DAMAGE_STATE]<predamageState) return false;
	
	// alternate strain has cracking srain in current configuration
	Tensor *altStrain = mptr->GetAltStrainTensor();
	*ecrack = *altStrain;
	
	// has cracking strain
	return true;
}

// return cracking COD in the crack axis system
// Only used by DeleteDamaged custom task
Vector IsoSoftening::GetCrackingCOD(MPMBase *mptr,bool is2D) const
{
	Vector cod = MakeVector(0.,0.,0.);
	
	// get history
	double *soft = GetSoftHistoryPtr(mptr);
	
	// Before initiatiation, not damage strain
	if(soft[SOFT_DAMAGE_STATE]<predamageState) return cod;
	
	// alternate strain has cracking srain in current configuration
	Tensor *altStrain = mptr->GetAltStrainTensor();
	Tensor ecrack = *altStrain;
	
	// Get R to crack, rotate by Rtot*RToCrack
	Matrix3 F = mptr->GetDeformationGradientMatrix();
	Matrix3 Rtot;
	F.LeftDecompose(&Rtot,NULL);
	Matrix3 RToCrack;
	GetRToCrack(&RToCrack,soft,is2D,0);
	Rtot *= RToCrack;
	ecrack = Rtot.RTVoightR(&ecrack,false,is2D);

	// crack cod in crack axis system (strain scaled by Vp/Ac)
	cod.x = ecrack.xx;
	cod.y = 0.5*ecrack.xy;
	cod.z = 0.5*ecrack.xz;
	ScaleVector(&cod,1./(rho*soft[GCSCALING]));

	// has cracking strain
	return cod;
}

// zero-offset to history data
double *IsoSoftening::GetSoftHistoryPtr(MPMBase *mptr) const
{	return (double *)(mptr->GetHistoryPtr(0));
}

// get traction failure surface (or <0 if not a softening material)
int IsoSoftening::GetTractionFailureSurface(void) const { return tractionFailureSurface; }
