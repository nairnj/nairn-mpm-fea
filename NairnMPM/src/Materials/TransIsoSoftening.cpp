/********************************************************************************
	TransIsoSoftening.cpp
	nairn-mpm-fea

	Created by John Nairn on 4 Aug 2016.
	Copyright (c) 2016 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/TransIsoSoftening.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Materials/TIFailureSurface.hpp"
#include "Materials/LinearSoftening.hpp"
#include "System/UnitsController.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "System/ArchiveData.hpp"

// globals
extern double mtime;

int debugNum=-1;

#pragma mark TransIsoSoftening::Constructors and Destructors

// Constructor
// throws std::bad_alloc
TransIsoSoftening::TransIsoSoftening(char *matName,int matID) : TransIsotropic(matName,matID)
{
	initiationLaw = (FailureSurface *)new TIFailureSurface(this);
	softeningAI = new LinearSoftening();
	softeningI = new LinearSoftening();
	softeningTII = new LinearSoftening();
	softeningAII = new LinearSoftening();
	softeningII = new LinearSoftening();
	tractionFailureSurface = CUBOID_SURFACE;
	distributionMode = SOFTDIST_NONE;
	softenStatsMode = VARY_STRENGTH;
	softenCV = 0.;
	wAlpha = -1.;
	wV0 = -1.;
	frictionCoeff = -1.;					// <=0 disables friction
}

#pragma mark TransIsoSoftening::Initialization

// Read material properties
char *TransIsoSoftening::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    // look for different damage initiaion law
    if(strcmp(xName,"Initiation")==0)
    {	input = INITIATION_LAW_SELECTION;
        return (char *)this;
    }
	
    // look for different softeningAI law
    else if(strcmp(xName,"SofteningAI")==0 || strcmp(xName,"SofteningXX")==0)
    {	input = SOFTAI_LAW_SELECTION;
        return (char *)this;
    }
	
    // look for different softeningI law
    else if(strcmp(xName,"SofteningI")==0 || strcmp(xName,"SofteningYY")==0)
    {	input = SOFTI_LAW_SELECTION;
        return (char *)this;
    }
	
	// look for different softeningTII law
	else if(strcmp(xName,"SofteningTII")==0 || strcmp(xName,"SofteningXYY")==0)
	{	input = SOFTTII_LAW_SELECTION;
		return (char *)this;
	}
	
	// look for different softeningAII law
	else if(strcmp(xName,"SofteningAII")==0 || strcmp(xName,"SofteningXYX")==0)
	{	input = SOFTAII_LAW_SELECTION;
		return (char *)this;
	}
	
	// look for different softeningII law
	else if(strcmp(xName,"SofteningII")==0 || strcmp(xName,"SofteningYZY")==0)
	{	input = SOFTII_LAW_SELECTION;
		return (char *)this;
	}
	
	// type of failure surface in 3D shear
	else if(strcmp(xName,"tractionFailureSurface")==0)
	{	input=INT_NUM;
		return((char *)&tractionFailureSurface);
	}

	// tolerance to converge in 3D shear
	// strength coefficient of variation
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
	
    // check for softening law property
    char *ptr = initiationLaw->InputInitationProperty(xName,input,gScaling);
    if(ptr != NULL) return ptr;
	
	// if possible check softening laws
	if(strlen(xName)>2)
	{	if(xName[0]=='A' && xName[1]=='I' && xName[2]=='-' && strlen(xName)>3)
		{	ptr = softeningAI->InputSofteningProperty(&xName[3],input,gScaling);
			if(ptr != NULL) return ptr;
		}
		else if(xName[0]=='I' && xName[1]=='-')
		{	ptr = softeningI->InputSofteningProperty(&xName[2],input,gScaling);
			if(ptr != NULL) return ptr;
		}
		else if(xName[0]=='I' && xName[1]=='I' && xName[2]=='-' && strlen(xName)>3)
		{	ptr = softeningII->InputSofteningProperty(&xName[3],input,gScaling);
			if(ptr != NULL) return ptr;
		}
		if(strlen(xName)>4)
		{	if(xName[0]=='T' && xName[1]=='I' && xName[2]=='I' && xName[3]=='-')
			{	ptr = softeningTII->InputSofteningProperty(&xName[4],input,gScaling);
				if(ptr != NULL) return ptr;
			}
			else if(xName[0]=='A' && xName[1]=='I' && xName[2]=='I' && xName[3]=='-')
			{	ptr = softeningAII->InputSofteningProperty(&xName[4],input,gScaling);
				if(ptr != NULL) return ptr;
			}
		}
	}
	
    // otherwise get material properties
    return TransIsotropic::InputMaterialProperty(xName,input,gScaling);
}

// Allows selected initiation laws
bool TransIsoSoftening::AcceptInitiationLaw(FailureSurface *iLaw,int lawID)
{
	// only allows max principle stress
	if(lawID != TIFAILURESURFACE) return false;
	
	delete initiationLaw;
    initiationLaw = iLaw;
    return true;
}

// Requies AI, TII, AII, I, or II softening
bool TransIsoSoftening::AcceptSofteningLaw(SofteningLaw *sLaw,int mode)
{
	if(mode==SOFTAI_LAW_SELECTION)
	{	delete softeningAI;
		softeningAI = sLaw;
	}
	else if(mode==SOFTI_LAW_SELECTION)
	{	delete softeningI;
		softeningI = sLaw;
	}
	else if(mode==SOFTTII_LAW_SELECTION)
	{	delete softeningTII;
		softeningTII = sLaw;
	}
	else if(mode==SOFTAII_LAW_SELECTION)
	{	delete softeningAII;
		softeningAII = sLaw;
	}
	else if(mode==SOFTII_LAW_SELECTION)
	{	delete softeningII;
		softeningII = sLaw;
	}
	else
		return false;
    return true;
}

// Swap two softening laws
// do not override
void TransIsoSoftening::SwapLaws(SofteningLaw **sl1,SofteningLaw **sl2)
{   SofteningLaw *tmp = *sl1;
    *sl1 = *sl2;
    *sl2 = tmp;
}

// verify settings and some initial calculations
const char *TransIsoSoftening::VerifyAndLoadProperties(int np)
{
	if(!useLargeRotation)
		return "TransIsoSoftening materials require activation of large rotation option";
	
	if(np==PLANE_STRESS_MPM)
		return "TransIsoSoftening materials cannot yet be used in plane stress calculations";
	
    // 3D requires axial z, 2D can be z or y (3D and axialCode gets done by parent class)
    if(materialID==TRANSISOSOFTENING2)
    {   if(np==THREED_MPM)
            return "3D simulations cannot use TransIsoSoftening 2 material";
        if(swapz>0)
            return "Cannot swap y and z axes of TransIsoSoftening 2 material";
        axialCode = AXIAL_Y;
    }
    
	// call initiation law that is used
    const char *ptr = initiationLaw->VerifyAndLoadProperties(np);
    if(ptr != NULL) return ptr;
	
	// check in superclass (along with its initialization)
	const char *err = TransIsotropic::VerifyAndLoadProperties(np);
	if(err!=NULL) return err;
	
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
			return "Weibull scale and reference volumne must be greater than 0";
		wGam1A = gamma_fxn(1.+1./wAlpha);
		softenCV = sqrt(gamma_fxn(1.+2./wAlpha)/(wGam1A*wGam1A) - 1.);
	}

	// make sure have a valid setting
	if(tractionFailureSurface==CYLINDER_SURFACE)
	{	return "Cylindrical failure surface not allowed for TransIsoSoftening material";
	}
	else if(tractionFailureSurface==COUPLED_CUBOID_SURFACE)
	{	// treat same as OVOID
		tractionFailureSurface = OVOID_SURFACE;
	}
	
#ifdef POROELASTICITY
	// poroelasticity conversions
	if(DiffusionTask::HasPoroelasticity())
	{	return "TransIsoSoftening does not implement poroelasticity";
	}
	else
	{	// diffusion CT is 1
		diffusionCT = 1.;
	}
#else
	// diffusion CT is 1
	//diffusionCT = 1.;
#endif
    
	// finish up in parent class
	return NULL;
}

// print mechanical properties to the results
void TransIsoSoftening::PrintMechanicalProperties(void) const
{
    TransIsotropic::PrintMechanicalProperties();
	
	cout << "Failure Surface: ";
	initiationLaw->PrintInitiationProperties();
	
	cout << "Isotropic Plane Mode I Softening: ";
	softeningI->PrintSofteningProperties(rho*initiationLaw->sigmaI());
	cout << "Isotropic Plane Mode II Softening: ";
	softeningII->PrintSofteningProperties(rho*initiationLaw->sigmaII());
	cout << "Axial Mode I Softening: ";
	softeningAI->PrintSofteningProperties(rho*((TIFailureSurface *)initiationLaw)->sigmaA());
	cout << "Transverse Shear Softening: ";
	softeningTII->PrintSofteningProperties(rho*((TIFailureSurface *)initiationLaw)->tauT());
	cout << "Axial Shear Softening: ";
	softeningAII->PrintSofteningProperties(rho*((TIFailureSurface *)initiationLaw)->tauA());
	
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
	
	// currently only cuboid, but save for future
	initiationLaw->SetFailureSurface(tractionFailureSurface);
	cout << "Traction Failure Surface: ";
	if(tractionFailureSurface == CYLINDER_SURFACE)
		cout << "cylindrical" << endl;
	else if(tractionFailureSurface == CUBOID_SURFACE)
		cout << "decoupled cuboid" << endl;
    else if(tractionFailureSurface == OVOID_SURFACE)
        cout << "coupled cuboid" << endl;
	else
		cout << "unsupported failure surface" << endl;
	
	cout << "Post-failure contact: ";
	if(frictionCoeff>0.)
		cout << "Coulomb coefficient of friction = " << frictionCoeff << endl;
	else
		cout << "frictionless" << endl;
}

#pragma mark TransIsoSoftening::History Data Methods

// Create history variables needed for softening behavior
char *TransIsoSoftening::InitHistoryData(char *pchr,MPMBase *mptr)
{
	// Validate this law (in use) is stable on current grid
	//	This validation uses average values. Large stochastic variations might
	//	cause some particles to be unstable. Enhancement could be to validate
	//	using minimum Gc and maximum strength instead
	// Does for all particles, but this makes sure softening material
	//   is being used and works with variable size particles
	
	// particle size
	double delx = mptr->GetMinParticleLength();
	
	// ET failure (sigmaI/C11)
	double sigMax = initiationLaw->sigmaI();
	double estr0 = fmobj->IsThreeD() ? sigMax/pr.C[0][0] : sigMax/pr.C[1][1] ;
	VerifyStability(sigMax,estr0,delx,softeningI,"I","sigmaT");
	
	// GT failure (sigmaII/C66) not used for AXIAL_Y
	sigMax = initiationLaw->sigmaI();
	estr0 = -1.;
	if(fmobj->IsThreeD())
		estr0 = sigMax/pr.C[5][5];
	else if(AxialDirection()==AXIAL_Z)
		estr0 = sigMax/pr.C[3][3];
	VerifyStability(sigMax,estr0,delx,softeningII,"II","tauRS");
	
	// Axial failure (sigmaA/C33 in 3D or sigmaA/C22 in 2D)
	sigMax = ((TIFailureSurface *)initiationLaw)->sigmaA();
	estr0 = -1.;
	if(fmobj->IsThreeD() ||AxialDirection()==AXIAL_Y)
		estr0 = sigMax/pr.C[2][2];
	VerifyStability(sigMax,estr0,delx,softeningAI,"AI","sigmaA");
	
	// transverse shear failure (tauT/GA - in C[3][3] both 2D and 3D)
	sigMax = ((TIFailureSurface *)initiationLaw)->tauT();
	estr0 = -1.;
	if(fmobj->IsThreeD() || AxialDirection()==AXIAL_Y)
		estr0 = sigMax/pr.C[3][3];
	VerifyStability(sigMax,estr0,delx,softeningTII,"TII","tauT");

	// axial shear failure (tauA/GA - in C[3][3] both 2D and 3D)
	sigMax = ((TIFailureSurface *)initiationLaw)->tauA();
	estr0 = -1.;
	if(fmobj->IsThreeD() || AxialDirection()==AXIAL_Y)
		estr0 = sigMax/pr.C[3][3];
	VerifyStability(sigMax,estr0,delx,softeningAII,"AII","tauA");
	
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
		double relValue = fmax(1. + softenCV*NormalCDFInverse(fract),0.1);
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

    return (char *)p;
}

// reset history data
void TransIsoSoftening::ResetHistoryData(char *pchr,MPMBase *mptr)
{	double *p = (double *)pchr;
	double relStrength = p[RELATIVE_STRENGTH];
	double relToughness = p[RELATIVE_TOUGHNESS];
	ZeroDoubles(pchr,SOFT_NUMBER_HISTORY);
	p[SOFT_DAMAGE_STATE] = 0.1;
	p[RELATIVE_STRENGTH] = relStrength;
	p[RELATIVE_TOUGHNESS] = relToughness;
}

// validate stability for softening law on current particle
void TransIsoSoftening::VerifyStability(double sigMax,double estr0,double delx,SofteningLaw *softLaw,
													const char *softMod,const char *softStrength) const
{	if(estr0<0.) return;
	double maxSize = softLaw->GetEtaStability()*softLaw->GetGc()/(rho*sigMax*estr0);
	if (maxSize < delx)
	{	cout << "Need to increase Gc for " << softMod << ", decrease " << softStrength << "^2, or decrease particle size "
				<< (delx/maxSize) << " fold for spatial stability." << endl;
		throw CommonException("A softening law needs smaller grid cells (see output details)", "TransIsoSoftening::VerifyStability");
	}
}

// Number of history variables
int TransIsoSoftening::NumberOfHistoryDoubles(void) const { return SOFT_NUMBER_HISTORY; }

// Initialize damage when requested
void TransIsoSoftening::SetInitialConditions(InitialCondition *ic,int ptNum,bool is3D)
{
    //MPMBase *mptr = mpm[ptNum-1];
    throw CommonException("Initial conditions not yet implemented for TransIso/OrthoSoftening", "TransIsoSoftening::SetInitialConditions");
}

// Return damage normal of global coordinates. It may be used by
// visualization to view normal. The magnitude of the vector
// is set equal to VpAc scaling for visualization use too.
Vector TransIsoSoftening::GetDamageNormal(MPMBase *mptr,bool threeD) const
{
	double *soft = GetSoftHistoryPtr(mptr);
	
	// Crack rotation (return zero if no damage)
    // This needs to know D form to get R to crack??
    // Need to fix if this material will set initial damage state
	Matrix3 RToCrack;
    if(soft[SOFT_DAMAGE_STATE] < predamageState)
        return MakeVector(0.,0.,0.);
    int DForm = GetDForm(soft[SOFT_DAMAGE_STATE]);
    GetRToCrack(&RToCrack,soft,!threeD,DForm);
	
	// get current rotation matrix
	Matrix3 pF = mptr->GetDeformationGradientMatrix();
    Matrix3 Rtot;
	pF.RightDecompose(&Rtot,NULL);
	
	// initial rotation to material axes and then crack axes
	Matrix3 R0 = mptr->GetInitialRotation();
	Rtot *= R0;
	Rtot *= RToCrack;
	
	// rotate x axis and scale by Vp/Ac
	double scale = 1./(rho*soft[GCSCALING]);
	return MakeVector(Rtot(0,0)*scale,Rtot(1,0)*scale,Rtot(2,0)*scale);
}

// set relative strength and toughness. Only called by CustomThermalRamp custom tasks
// to set disteibution of strengths and toughness.
void TransIsoSoftening::SetRelativeStrength(MPMBase *mptr,double rel)
{   double *soft = GetSoftHistoryPtr(mptr);
    soft[RELATIVE_STRENGTH] = rel;
}
void TransIsoSoftening::SetRelativeToughness(MPMBase *mptr,double rel)
{   double *soft = GetSoftHistoryPtr(mptr);
    soft[RELATIVE_TOUGHNESS] = rel;
}

#pragma mark TransIsoSoftening::Accessors

/* Take increments in strain and calculate new Particle: strains, rotation strain,
		stresses, strain energy,
	dvij are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
*/
void TransIsoSoftening::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,
										   		ResidualStrains *res,int historyOffset,Tensor *gStress) const
{
	// set 2D flag
	bool is2D = np == THREED_MPM ? false : true;

    // get strain increments in material axis system: de = R0T.deT(initial).R0
    // Task 1: Convert input form to effective strains in material configuration
    //-------------------------------------------------------------------------
    Matrix3 dR,Rnm1,Rtot;
    Matrix3 R0 = mptr->GetInitialRotation();
    Matrix3 deT = LRGetStrainIncrement(INITIALMATERIAL,mptr,du,&dR,&R0,&Rnm1,&Rtot);
    
    // get strain increments in material axis system as Tensor type
    Tensor de = is2D ?
        MakeTensor2D(deT(0, 0), deT(1, 1), deT(2, 2), deT(0, 1) + deT(1, 0)) :
        MakeTensor(deT(0,0), deT(1,1), deT(2,2), deT(1,2)+deT(2,1), deT(0,2)+deT(2,0), deT(0,1)+deT(1,0));
    
    // get total rotation (from material axis system to global coordinates)
    // Rnm1 is from material axes to global state n-1
    // Rtot is from material axes to current global axes
    // dR is from state n-1 to n (Rtot = dR*Rnm1)
    Rnm1 *= R0;
    Rtot *= R0;
	
	// cast pointer to material properties in material axes
	ElasticProperties *p = GetElasticPropertiesPointer(properties);
	
	// residual strains (thermal and moisture) in material axes
	double exxr,eyyr,ezzr;
    Tensor er = GetAnisoResStrains(exxr,eyyr,ezzr,p,res,np);
	
	// pointers to stress tensor (in global axes) and history
	Tensor *sp = mptr->GetStressTensor();
    double *soft = GetSoftHistoryPtr(mptr);
	
#pragma mark ... Code for Previously Undamaged Material
	// Before cracking, do provisional transversely istropic update.
	// If not cracked, then done, otherwise initiate the damage
	if(soft[SOFT_DAMAGE_STATE]<predamageState)
	{	// effective strains in mnaterial axis system
		double dvxxeff = de.xx-er.xx;
		double dvyyeff = de.yy-er.yy;
		double dvzzeff = de.zz-er.zz;
		double dgamxy = de.xy,dgamyz=0.,dgamxz=0.;
		
		// stress increments in material analysis axes
		Tensor str;
		if(np==THREED_MPM)
		{	// more effective strains in material axis system
			dgamyz = de.yz;
			dgamxz = de.xz;
			
			// rotate prior stress to material axes
			str = Rnm1.RTVoightR(sp, true, false);
			
			// update stress in material axis system
			str.xx += p->C[0][0]*dvxxeff + p->C[0][1]*dvyyeff + p->C[0][2]*dvzzeff;
			str.yy += p->C[1][0]*dvxxeff + p->C[1][1]*dvyyeff + p->C[1][2]*dvzzeff;
			str.zz += p->C[2][0]*dvxxeff + p->C[2][1]*dvyyeff + p->C[2][2]*dvzzeff;
			str.yz += p->C[3][3]*dgamyz;
			str.xz += p->C[4][4]*dgamxz;
			str.xy += p->C[5][5]*dgamxy;
		}
		
		else
		{	// rotate previous stress to material axes
			str = Rnm1.RTVoightR(sp, true, true);
		
			// update stress in material axis system (no zz stress yet)
			if(np==AXISYMMETRIC_MPM)
			{	// hoop stress affect on RR, ZZ, and RZ stresses
				str.xx += p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff + p->C[4][1]*dvzzeff;
				str.yy += p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff + p->C[4][2]*dvzzeff;
			}
			else
			{   // plain strain
                str.xx += p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff;
				str.yy += p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff;
			}
			str.xy += p->C[3][3]*dgamxy;
		}
		
		// Here: str is stress in material axes after a provisional elastic update

		// check if has failed
		Vector norm;
		double tempRelStrength = soft[RELATIVE_STRENGTH];
		int failureMode = initiationLaw->ShouldInitiateFailure(&str,&norm,np,tempRelStrength,NULL);
		if(failureMode == NO_FAILURE)
		{	// Not failed yet, so finish elastic update
			
			// rotate stress to global coordinates new sigma = Rtot str RtotT
			*sp = Rtot.RVoightRT(&str, true, is2D);
			
			// stresses are in global coordinates so need to rotate strain and residual
			// strain to get work energy increment per unit mass (dU/(rho0 V0))
			de = Rtot.RVoightRT(&de, false, is2D);
			er = Rtot.RVoightRT(&er, false, is2D);
			
			// updates
			if(np==THREED_MPM)
			{	// energy
				mptr->AddWorkEnergyAndResidualEnergy(DotTensors(sp, &de), DotTensors(sp, &er));
			}
			else
			{	// work and residual strain energy increments and sigma or F in z direction
				double workEnergy = sp->xx*de.xx + sp->yy*de.yy + sp->xy*de.xy;
				double resEnergy = sp->xx*er.xx + sp->yy*er.yy + sp->xy*er.xy;
				if(np==PLANE_STRAIN_MPM)
				{	// need to add back terms to get from reduced cte to actual cte
					sp->zz += p->C[4][1]*(dvxxeff+p->alpha[5]*ezzr) + p->C[4][2]*(dvyyeff+p->alpha[6]*ezzr) - p->C[4][4]*ezzr;
					
					// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
					resEnergy += sp->zz*er.zz;
				}
				else
				{	// axisymmetric hoop stress
					sp->zz += p->C[4][1]*dvxxeff + p->C[4][2]*dvyyeff + p->C[4][4]*dvzzeff;
					
					// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
					workEnergy += sp->zz*de.zz;
					resEnergy += sp->zz*er.zz;
				}
				mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
			}
			
			// track heat energy (should get dTq0)
			//IncrementHeatEnergy(mptr,dTq0,0.);
			
			// this elastic update is done
			return;
		}
		
#pragma mark ...... Undamaged Material Just Initiated Damage
		// initiate failure
		// 3D, norm = ZYZ rotation from material axes to crack normal. Currently either (0,pi/2,0)
		//		to indicate normal in axial direction (switches z and x) or (theta,0,0) to keep axial direciton
		//		in z direction, but around that axis
		// 2D, norm = (cos(theta),sin(theta))
		soft[NORMALDIR1] = norm.x;										// cos(theta) (2D) or Euler alpha (3D) to normal
		soft[NORMALDIR2] = norm.y;										// sin(theta) (2D) or Euler beta (3D)  to normal
		soft[NORMALDIR3] = norm.z;										// unused (2D) or Euler gamma (3D) to normal
		
		// initiate failure in material axis system with predamageState=0.5
		// Note that TI, initForm=Di_DAMAGE means axial direction in i direction in crack axis system
		//      but for Ortho subclass, in means material i direction and crack axis x direction
		int initForm = DecodeDamageInitiation(np,&norm,failureMode,soft);
		
		// get intersection area after rotating normal into global axes
		// Rotation from crack to material axes
		Matrix3 RToCrack;
		GetRToCrack(&RToCrack,soft,is2D,initForm);
		
		// initial rotation to material axes and then to the crack
		R0 *= RToCrack;
		Vector globalNorm = MakeVector(R0(0,0),R0(1,0),R0(2,0));
		
		// get intersection with this normal
		soft[GCSCALING] = GetAcOverVp(np,mptr,&globalNorm)/rho;			// divide by rho because divided by specific stress
	}

#pragma mark ... Code for Damaged Material
	// A crack is present - get tensor form and properties
    int DForm = GetDForm(soft[SOFT_DAMAGE_STATE]);
    CrackAxisProperties d;
    LoadCrackAxisProperties(np,&d,DForm,p);

    // Rotate strain increments from material axes to CAS
	Matrix3 RToCrack;
	GetRToCrack(&RToCrack, soft, is2D, DForm);
	de = RToCrack.RTVoightR(&de, false, is2D);
	er = RToCrack.RTVoightR(&er, false, is2D);

	// Rotate from CAS to n-1 or current configuration
	Rnm1 *= RToCrack;
    Rtot *= RToCrack;
    
    // Get current stress in CAS
    Tensor str = Rnm1.RTVoightR(sp, true, is2D);

    // rotate particle cracking strains by increment and get cracking strain in crack axis syste,
    Tensor *eplast=mptr->GetAltStrainTensor();
    Tensor ecrack = Rnm1.RTVoightR(eplast, false, is2D);
    *eplast = dR.RVoightRT(eplast, false, is2D);
    
    // finish evolution
    DamageEvolution(mptr,np,soft,de,str,er,d,dR,Rtot,ecrack,DForm,0.,delTime);
}

// Finish damage evolution for anisotropic materials
// soft is pointer to first sofening history variable
// Tensor de = total elastic strain increment in CAS
// Tensor er = total residual strain increment in CAS
// str is prior stress in the CAS
// dR and Rtot are incemental and total rotation matrices
// ecrack is cracking strain in the CAS
// dispEnergy is dissipated energy before evolving damage (or zero if none)
void TransIsoSoftening::DamageEvolution(MPMBase *mptr,int np,double *soft,Tensor &de,Tensor &str,Tensor &er,CrackAxisProperties &d,
                     Matrix3 &dR,Matrix3 &Rtot,Tensor &ecrack,int DForm,double dispEnergyPlastic,double delTime) const
{
    // set 2D flag
    bool is2D = np == THREED_MPM ? false : true;

    // effective strain
	double dvxxeff = de.xx-er.xx;
	double dvyyeff = de.yy-er.yy;
	double dgamxy = de.xy-er.xy;
	double dvzzeff,dgamxz=0.,dgamyz=0.;
	if(np==THREED_MPM)
	{	// 3D effective strains in crack axis system
		dvzzeff = de.zz-er.zz;
		dgamxz = de.xz-er.xz;
		dgamyz = de.yz-er.yz;
	}
	else
	{	// Effective strains in crack axis system including zz effective strain
		// FINISH UP: Plane stress
		dvxxeff = de.xx-er.xx;
		dvyyeff = de.yy-er.yy;
		dgamxy = de.xy-er.xy;
		dvzzeff = 0.;
		if(np==AXISYMMETRIC_MPM)
			dvzzeff = de.zz-er.zz;
	}
	
	// Up to here:
	// de and er are incremental strain and residual strain in CAS
	// dvijeff and dgamij are effective strain increments in CAS
	// str is previous stress in CAS
	// Rtot rotates between CAS and current configuration
	// Rnm1 rotates between CAS and previous configuration
	// ecrack is cracking strain in CAS
	
	// normal strain increments (dvzzeff=0 in plane strain)
	double den = dvxxeff + d.C12C11*dvyyeff + d.C13C11*dvzzeff;
	
	// to calculate cracking plane changes depending on current state
	double decxx=0.,dgcxy=0.,dgcxz=0.;
	Tensor dsig = MakeTensor(0.,0.,0.,0.,0.,0.);
    
    // increments in tensile and shear dissipation
    double dispEnergyN=0.,dispEnergyXY=0.,dispEnergyXZ=0.,dispEnergyDam=0.;
	
#pragma mark ...... Damage Evolution Calculations
	// This section means crack has failed
	if(soft[SOFT_DAMAGE_STATE]>damageState)
	{	PostFailureUpdate(decxx,dgcxy,dgcxz,&dsig,&str,&ecrack,Rtot,d.C11,den,
                                d.C66,dgamxy,d.C55,dgamxz,np!=THREED_MPM,frictionCoeff);
	}
	
	// The damage is evolving after the crack has initiated
	else
	{	// flag for reached critical strain
		bool criticalStrain = false;
		
		// normal stress and strains
		double relStrength = soft[RELATIVE_STRENGTH];
		double sigmaI = relStrength*d.sigNc;
		double relToughness = soft[RELATIVE_TOUGHNESS];
		double scaleI = relToughness*soft[GCSCALING]/sigmaI;
		// shear terms
		double Txy0 = fabs(str.xy);
		double tau0XY = relStrength*d.tauXYc;
        double scaleXY = relToughness*soft[GCSCALING]/tau0XY;
        double tau0XZ = 1.,scaleXZ = 1.;
		
        if(tractionFailureSurface==CUBOID_SURFACE)
        {	// normal strains
            if(!SoftenAxis(mptr,den,soft,DELTANORMAL,DAMAGENORMAL,sigmaI,scaleI,d.fnLaw,
                                d.C11,str.xx,-1.,-1.,decxx,dispEnergyN,criticalStrain))
            {    // Elastic loading but handle contact if den<0
                decxx = soft[DAMAGENORMAL]*den;
                if(den<0.)
                {    // handle contact to keep keep ecxx>=0
                    if(ecrack.xx+decxx<0.) decxx = -ecrack.xx;
                }
            }

            // 2D and 3D shear strain x-y
            double dgs = str.xy>0. ? dgamxy : -dgamxy;
            if(SoftenAxis(mptr,dgs,soft,DELTASHEAR,DAMAGESHEAR,tau0XY,scaleXY,d.fXYLaw,
                                d.C66,Txy0,-1.,-1.,dgcxy,dispEnergyXY,criticalStrain))
            {	// Adjust sign to match shear stress direction
                if(str.xy<0.) dgcxy = -dgcxy;
            }
            else
            {	// elastic loading
                dgcxy = soft[DAMAGESHEAR]*dgamxy;
            }
            
            // 3D x-z shear strain
            if(!is2D)
            {   double Txz0 = fabs(str.xz);
                tau0XZ = relStrength*d.tauXZc;
                scaleXZ = relToughness*soft[GCSCALING]/tau0XZ;
                dgs = str.xz>0. ? dgamxz : -dgamxz;
                if(SoftenAxis(mptr,dgs,soft,DELTASHEAR2,DAMAGESHEAR2,tau0XZ,scaleXZ,d.fXZLaw,
                                    d.C55,Txz0,-1.,-1.,dgcxz,dispEnergyXZ,criticalStrain))
                {	// Adjust sign to match shear stress direction
                    if(str.xz<0.) dgcxz = -dgcxz;
                }
                else
                {	// elastic loading
                    dgcxz = soft[DAMAGESHEAR2]*dgamxz;
                }
            }
        }
		else if(tractionFailureSurface == OVOID_SURFACE)
        {   // Actually cuboid surface, but directions are coupled like ovoid in isotropic materials
            if(!OvoidSoftening(mptr,is2D,soft,den,dgamxy,dgamxz,str,&d,sigmaI,scaleI,tau0XY,scaleXY,tau0XZ,scaleXZ,
                               dispEnergyN,dispEnergyXY,dispEnergyXZ,criticalStrain,
                               ecrack.xx,decxx,dgcxy,dgcxz))
            {   // elastic, all use same D value
                decxx = soft[DAMAGENORMAL]*den;
                if(den<0.)
                {    // handle contact to keep keep ecxx>=0
                    if(ecrack.xx+decxx<0.) decxx = -ecrack.xx;
                }
                dgcxy = soft[DAMAGENORMAL]*dgamxz;
                if(!is2D) dgcxz = soft[DAMAGENORMAL]*dgamxz;
            }
        }
        else
			throw "Cylindrical surface not available for aniostropic materials";
		
#pragma mark ...... Dissipated Energy and Check for Final Failure
		// Use energy release rate to check on failure
        dispEnergyDam = dispEnergyN+dispEnergyXY+dispEnergyXZ;
		
        // Use energy release rate to check on failure
        // Get Gbar/Gcbar
        double relGI,relGII1,relGII2=0.,cutoff=0.;
        if(tractionFailureSurface == CUBOID_SURFACE)
        {   // Only 3D Cuboid uses softening law to get energy released
            
            // Get GI/GIc
            relGI = d.fnLaw->GetGoverGc(soft[DELTANORMAL],scaleI);
            if(is2D) soft[DAMAGEGI] += dispEnergyN/soft[GCSCALING];
		
            // Get GII/GIIc
            relGII1 = d.fXYLaw->GetGoverGc(soft[DELTASHEAR],scaleXY);
            if(is2D) soft[DAMAGEGII] += dispEnergyXY/soft[GCSCALING];
		
            // Get second GII/GIIc
		    relGII2 = d.fXZLaw!=NULL ? d.fXZLaw->GetGoverGc(soft[DELTASHEAR2],scaleXZ) : 0. ;
		
            // criterion
		    cutoff = pow(relGI,nmix) + pow(relGII1,nmix) + pow(relGII2,nmix);
        }
        else if(is2D)
        {   // for 2D cuboid (above) and 2D ovoid (on cuboid surface), track GI and GII in 3D parameters
            // they are in 3D history parameters
            
            // store GI and GII
            soft[DAMAGEGI] += dispEnergyN/soft[GCSCALING];
            soft[DAMAGEGII] += dispEnergyXY/soft[GCSCALING];
            
            // Get GIbar/GIcbar
            relGI = soft[DAMAGEGI]/(relToughness*d.fnLaw->GetGc());
                
            // Get (GIIbar/GIIbarc)
            relGII1 = soft[DAMAGEGII]/(relToughness*d.fXYLaw->GetGc());
        }
        else
        {   // 3D ovoid (on cuboid surface) combines GXY and GXZ into GII
            // stores on unused damage parameters
            // Future could add new history variable to track GXY and GXZ separately
            
            // store GI and GII
            soft[DAMAGESHEAR] += dispEnergyN/soft[GCSCALING];
            soft[DAMAGESHEAR2] += (dispEnergyXY+dispEnergyXZ)/soft[GCSCALING];
            
            // Get GIbar/GIcbar
            relGI = soft[DAMAGESHEAR]/(relToughness*d.fnLaw->GetGc());
                
            // Get (GIIbar/GIIbarc) (reduces by XY properties, but those removed below)
            relGII1 = soft[DAMAGESHEAR2]/(relToughness*d.fXYLaw->GetGc());
       }

		// check limit or energy
		if(criticalStrain || cutoff>=1.)
		{	// report energy released
			double GI = relToughness*relGI*d.fnLaw->GetGc()*UnitsController::Scaling(1.e-3);
			double GII1 = relToughness*relGII1*d.fXYLaw->GetGc()*UnitsController::Scaling(1.e-3);
			double GII2 = d.fXZLaw!=NULL ? relToughness*relGII2*d.fXZLaw->GetGc()*UnitsController::Scaling(1.e-3) : 0. ;
			double alpha=soft[NORMALDIR1],beta=soft[NORMALDIR2],gamma=soft[NORMALDIR3];
			if(np!=THREED_MPM)
			{	alpha = beta>=0. ? acos(alpha) : -acos(alpha) ;
				beta = 0.;
			}
			archiver->Decohesion(mtime,mptr,alpha,beta,gamma,GI,GII1,GII2,soft[SOFT_DAMAGE_STATE]);
			
			// now failed (this material does not track fraction mode I because initial
			// state sill needed in post failure modeling. Get fraction mode I from
			// decohesion info or add new history variable
			soft[SOFT_DAMAGE_STATE] += 1.0;
			soft[DAMAGENORMAL] = 1.;
            if(tractionFailureSurface == CUBOID_SURFACE)
            {   soft[DAMAGESHEAR] = 1.;
                if(d.fXZLaw!=NULL) soft[DAMAGESHEAR2] = 1.;
            }
			
			// final crack plane stresses should be zero
			dsig.xx = -str.xx;
			decxx = den - dsig.xx/d.C11;
			dsig.xy = -str.xy;
			dgcxy = dgamxy - dsig.xy/d.C66;
			if(d.fXZLaw!=NULL)
			{	dsig.xz = -str.xz;
				dgcxz = dgamxz - dsig.xz/d.C55;
			}
		}
		else
		{	// crack plane stress increment
			dsig.xx = d.C11*(den - decxx);
			dsig.xy = d.C66*(dgamxy-dgcxy);
			if(d.fXZLaw!=NULL) dsig.xz = d.C55*(dgamxz-dgcxz);
		}
		
		// dissipated energy due to damage (if plasticity, it was added before)
        // 2D tracks damage only in NORMALDIR3 history; subtract it from plastic to get plastic energy
        if(is2D) soft[NORMALDIR3] += UnitsController::Scaling(1.e-9)*mptr->mp*dispEnergyDam;
		mptr->AddPlastEnergy(dispEnergyDam);
	}

#pragma mark ... Update Cracking Strain
	// rotate previous stress by rotation increment
    Tensor *sp = mptr->GetStressTensor();
	str = dR.RVoightRT(sp,true,is2D);
	
	// get remaining stress incements in crack axis system
	if(np==THREED_MPM)
	{	// stress increment
		dsig.yy = d.C12*(dvxxeff - decxx) + d.C22*dvyyeff + d.C13*dvzzeff;
		dsig.zz = d.C13*(dvxxeff - decxx) + d.C23*dvyyeff + d.C33*dvzzeff;
		dsig.yz = d.C44*dgamyz;
	}
	else
	{	dsig.yy = d.C12*(dvxxeff - decxx) + d.C22*dvyyeff;
		if(np==AXISYMMETRIC_MPM)
		{	// extra term
			dsig.yy += d.C13*dvzzeff;
			
			// axisymmetric hoop stress
			dsig.zz = d.C13*dvxxeff + d.C23*dvyyeff + d.C33*dvzzeff;
		}
		else if(np==PLANE_STRAIN_MPM)
		{	// need to add back terms to go from reduced cte to actual cte
			// then subtract cracking strain from x and others use current stiffness elements
            dsig.zz = d.C13*(dvxxeff + d.vzxc*er.xx - decxx) + d.C23*(dvyyeff + d.vzyc*er.yy) - d.C33*er.zz;
		}
	}
	
	// update stress and strain (in which cracking strain rotated to current configuration)
	UpdateCrackingStrain(np,mptr->GetAltStrainTensor(),decxx,dgcxy,dgcxz,Rtot,soft);
	UpdateCrackingStress(np,sp,&dsig,&str,Rtot);
	
#pragma mark ... Work, Residual, and Heat Energies
	// work energy, residual energy, and heat energy
	de = Rtot.RVoightRT(&de, false, is2D);		// strain in current config to match stress
	er = Rtot.RVoightRT(&er, false, is2D);		// residual strain to match stress
	if(np==THREED_MPM)
	{	// energy
		mptr->AddWorkEnergyAndResidualEnergy(DotTensors(sp,&de),DotTensors(sp,&er));	
	}
	else
	{	// work energy, residual energy in current configuration
		double workEnergy = sp->xx*de.xx + sp->yy*de.yy + sp->xy*de.xy;
		double resEnergy = sp->xx*er.xx + sp->yy*er.yy + sp->xy*er.xy;
		if(np==PLANE_STRAIN_MPM)
		{	// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
            workEnergy += sp->zz*de.zz;
			resEnergy += sp->zz*er.zz;
		}
		else if(np==AXISYMMETRIC_MPM)
		{	// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			workEnergy += sp->zz*de.zz;
			resEnergy += sp->zz*er.zz;
		}
		mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
		// TIQuery: Plane stress Fzz
	}
	
	// track heat energy (should get dTq0 in second parameter)
    IncrementHeatEnergy(mptr,0.,dispEnergyDam+dispEnergyPlastic);
}

// Ovoid damafe evolution surface for isotropic materials
bool TransIsoSoftening::OvoidSoftening(MPMBase *mptr,bool is2D,double *soft,
        double den,double dgxy,double dgxz,Tensor &str,CrackAxisProperties *d,
        double sigmaN,double scaleN,double sigmaXY,double scaleXY,double sigmaXZ,double scaleXZ,
        double &dispNEnergy,double &dispXYEnergy,double &dispXZEnergy,bool &criticalStrain,
        double ecxx,double &decxx,double &dgcxy,double &dgcxz) const
{
    // get parameter (all coupled)
    double dam = soft[DAMAGENORMAL];
    
    // get parameters
    double deln = soft[DELTANORMAL];
    double delxy = soft[DELTASHEAR];
    
    // trial stress and increments
    double C111md = d->C11*(1.-dam);
    double dTntrial = C111md*den;
    double Gxy1md = d->C66*(1.-dam);
    double dTxytrial = Gxy1md*dgxy;         // may be negative
    
    // softening laws
    double Fn = sigmaN*d->fnLaw->GetFFxn(deln,scaleN);
    double Fxy = sigmaXY*d->fXYLaw->GetFFxn(delxy,scaleXY);

    // additions for 3D
    double delxz,Gxz1md,dTxztrial,Fxz;
    if(is2D)
    {   // check for elastic (caller gets cracking strains)
        if((str.xx+dTntrial<=Fn) && (fabs(str.xy+dTxytrial)<=Fxy))
            return false;
    }
    else
    {   delxz = soft[DELTASHEAR2];
        Gxz1md = d->C55*(1.-dam);
        dTxztrial = Gxz1md*dgxz;            // may be negative
        Fxz = sigmaXZ*d->fXZLaw->GetFFxn(delxz,scaleXZ);
        
        // check for elastic (called gets cracking strain
        if((str.xx+dTntrial<=Fn) && (fabs(str.xy+dTxytrial)<=Fxy) && (fabs(str.xz+dTxztrial)<=Fxz))
            return false;
    }
    
    // damage has occurred in at least on direction
    double ddeltaN=0.,ddeltaXY=0.,ddeltaXZ=0.;
    double dDn=0.,dDxy=0.,dDxz=0.;
    double den2,dgxy2,dgxz2;
    double en0 = sigmaN/d->C11,gxy0 = sigmaXY/d->C66;
    double gxz0 = is2D ? 0. : sigmaXZ/d->C55;
    criticalStrain = false;

    // ratios of current stress to strength (should alays be <= 1)
    double sfN=str.xx/Fn,sfXY=str.xy/Fxy;
    double sfXZ = is2D ? 0. : str.xz/Fxz ;
    
    // normal direction
    if(str.xx+dTntrial>=Fn)
    {   // go to elastic surface and then damage (T(trial)-F)/(E(1-D))
        den2 = (str.xx+dTntrial-Fn)/C111md;
        ddeltaN = d->fnLaw->GetDDelta(den2,en0,deln,dam,scaleN);
        if(ddeltaN<0.)
            criticalStrain = true;
        else
        {   // get change in Dn
            double delNew = deln+ddeltaN;
            double Fnew = d->fnLaw->GetFFxn(delNew,scaleN);
            dDn = delNew/(delNew + en0*Fnew) - dam;
        }
    }
    
    // shear xy direction
    if(!criticalStrain && fabs(str.xy+dTxytrial)>=Fxy)
    {   // go to elastic surface and then damage (T(trial)-F)/(E(1-D))
        dgxy2 = (fabs(str.xy+dTxytrial)-Fxy)/Gxy1md;
        ddeltaXY = d->fXYLaw->GetDDelta(dgxy2,gxy0,delxy,dam,scaleXY);
        if(ddeltaXY<0.)
            criticalStrain = true;
        else
        {   // get change in Dxy
            double delNew = delxy+ddeltaXY;
            double Fnew = d->fXYLaw->GetFFxn(delNew,scaleXY);
            dDxy = delNew/(delNew + gxy0*Fnew) - dam;
        }
    }
    
    // shear xz strain
    if(!criticalStrain && !is2D && fabs(str.xz+dTxztrial)>=Fxz)
    {   // go to elastic surface and then damage
        dgxz2 = (fabs(str.xz+dTxztrial)-Fxz)/Gxz1md;
        ddeltaXZ = d->fXZLaw->GetDDelta(dgxz2,gxz0,delxz,dam,scaleXZ);
        if(ddeltaXZ<0.)
            criticalStrain = true;
        else
        {   // get change in Dxz
            double delNew = delxz+ddeltaXZ;
            double Fnew = d->fXZLaw->GetFFxn(delNew,scaleXZ);
            dDxz = delNew/(delNew + gxz0*Fnew) - dam;
        }
    }
    
    // if failed or find direction with most damage
    double dD=0.;
    bool damageN=false,damageXY=false,damageXZ=false;
    if(criticalStrain)
    {   // Decohesion: set delta=deltaMax, D=1
        // caller will get cracking strain in post-failure update code
        soft[DELTANORMAL] = d->fnLaw->GetDeltaMax(scaleN);
        ddeltaN = soft[DELTANORMAL] - deln;
        soft[DELTASHEAR] = d->fXYLaw->GetDeltaMax(scaleXY);
        ddeltaXY = soft[DELTASHEAR] - delxy;
        if(!is2D)
        {   soft[DELTASHEAR2] = d->fXZLaw->GetDeltaMax(scaleXZ);
            ddeltaXZ = soft[DELTASHEAR2] - delxz;
        }
        dD = 1.-dam;
        damageN = damageXY = damageXZ = true;
        
        // limit to 1, seems rare and not sure how ever, but might help last energy calc
        sfN = fmin(1.,sfN);
        sfXY = fmin(1.,fabs(sfXY));
        if(!is2D) sfXZ = fmin(1.,fabs(sfXZ));
    }
    
    else if(dDn>=fmax(dDxy,dDxz))
    {   // normal direction damages most
        soft[DELTANORMAL] += ddeltaN;
        dD = dDn;
        decxx = dam*(den-den2) + ddeltaN;
        sfN = 1.;
        damageN = true;

    }
    
    else if(dDxy>=dDxz)
    {   // shear x-y damages most
        soft[DELTASHEAR] += ddeltaXY;
        dD = dDxy;
        dgcxy = dam*(fabs(dgxy)-dgxy2) + ddeltaXY;
        if(str.xy<0.) dgcxy = -dgcxy;
        sfXY = 1.;
        damageXY = true;
    }
    
    else if(!is2D)
    {   // shear x-z damages most
        soft[DELTASHEAR2] += ddeltaXZ;
        dD = dDxz;
        dgcxz = dam*(fabs(dgxz)-dgxz2) + ddeltaXZ;
        if(str.xz<0.) dgcxz = -dgcxz;
        sfXZ = 1.;
        damageXZ = true;
    }
    
    // do non-damaging directions
    if(!damageN)
    {   // get delta from new D
        soft[DELTANORMAL] = d->fnLaw->GetDeltaFromDamage(dam+dD,scaleN,en0,deln);
        ddeltaN = soft[DELTANORMAL] - deln;
        
        // get cracking strain from initial values and stress at start of damage
        double Fnp = d->fnLaw->GetFpFxn(deln,scaleN);
        decxx = dam*den + sfN*(1.-dam*(1+en0*Fnp))*ddeltaN;
        if(ecxx+decxx<0.) decxx = -ecxx;
   }
    
    if(!damageXY)
    {   // get delta from new D
        soft[DELTASHEAR] = d->fXYLaw->GetDeltaFromDamage(dam+dD,scaleXY,gxy0,delxy);
        ddeltaXY = soft[DELTASHEAR] - delxy;
        
        // get cracking strain from initial values and stress at start of damage
        double Fxyp = d->fXYLaw->GetFpFxn(delxy,scaleXY);
        dgcxy = dam*dgxy + sfXY*(1.-dam*(1+gxy0*Fxyp))*ddeltaXY;
    }
    
    if(!is2D && !damageXZ)
    {   // get delta from new D
        soft[DELTASHEAR2] = d->fXZLaw->GetDeltaFromDamage(dam+dD,scaleXZ,gxz0,delxz);
        ddeltaXZ = soft[DELTASHEAR2] - delxz;
        
        // get cracking strain from initial values and stress at start of damage
        double Fxzp = d->fXZLaw->GetFpFxn(delxz,scaleXZ);
        dgcxz = dam*dgxz + sfXZ*(1.-dam*(1+gxz0*Fxzp))*ddeltaXZ;
    }

    // common updates
    soft[DAMAGENORMAL] += dD;
    
    // Get energy dissipation
    
    // dissipated energy increment per unit mass using dOmega = (1/2) (S/F)^2 phi ddelta
    if(str.xx>0.)
        dispNEnergy = 0.5*sfN*sfN*sigmaN*d->fnLaw->GetPhiFxn(deln,scaleN)*ddeltaN;
    dispXYEnergy = 0.5*sfXY*sfXY*sigmaXY*d->fXYLaw->GetPhiFxn(delxy,scaleXY)*ddeltaXY;
    if(!is2D) dispXZEnergy = 0.5*sfXZ*sfXZ*sigmaXZ*d->fXZLaw->GetPhiFxn(delxz,scaleXZ)*ddeltaXZ;

    return true;
}

// Convert damage initiation mode and normal to stored setting with information
//      to later get the correct damage tensor
// 3D, norm = ZYZ rotation from material axes to crack normal. For TI softening either (0,pi/2,0)
//		to indicate crack axis normal in axial direction (DX_DAMAGE) (switches z and x) or (theta,0,0) to keep axial direction
//		in z direction (DZ_DAMAGE), but rotate around that axis
// 2D, norm = (cos(theta),sin(theta))
int TransIsoSoftening::DecodeDamageInitiation(int np,Vector *norm,int failureMode,double *soft) const
{
	// initiate failure in material axis system with predamageState=0.5
	int initForm = DX_DAMAGE;
	if(np==THREED_MPM)
	{	if(norm->y > 1.)												// should be pi/2
		{	// Using Dx (axial direction along crack normal)
			if (failureMode == EA_FAILURE)
				soft[SOFT_DAMAGE_STATE] = 0.95;						// axial tension failure
			else
				soft[SOFT_DAMAGE_STATE] = 1.05;						// Transverse shear failure through axial direction
		}
		else
		{	// using Dz (axial direction parallel to crack plane)
			if (failureMode == TENSILE_FAILURE)						// transverse tensile failure
				soft[SOFT_DAMAGE_STATE] = 0.75;
			else if (failureMode == SHEAR_FAILURE)					// rolling shear failure
				soft[SOFT_DAMAGE_STATE] = 0.85;
			else
				soft[SOFT_DAMAGE_STATE] = 0.80;						// GA failure
			initForm = DZ_DAMAGE;
		}
	}
	else
	{	if(AxialDirection() == AXIAL_Z)
		{	// Always using Dz and x-y plane is isotropic
			if (failureMode == TENSILE_FAILURE)
				soft[SOFT_DAMAGE_STATE] = 0.75;						// transverse tensile failure
			else
				soft[SOFT_DAMAGE_STATE] = 0.85;						// Rolling shear failure
			initForm = DZ_DAMAGE;
		}
		else if(norm->y > 0.99)
		{	// Using Dx (norm = (0,1,0)) (axial direction along crack normal)
			if (failureMode == EA_FAILURE)
				soft[SOFT_DAMAGE_STATE] = 0.95;						// axial tension failure
			else
				soft[SOFT_DAMAGE_STATE] = 1.05;						// Transverse shear failure
		}
		else
		{	// Using Dy (norm = (1,0,0)) (axial direction parallel to crack plane)
			if (failureMode == GA_FAILURE)
				soft[SOFT_DAMAGE_STATE] = 1.15;						// axial shear failure
			else
				soft[SOFT_DAMAGE_STATE] = 1.25;						// transverse tension
			initForm = DY_DAMAGE;
		}
	}

	return initForm;
}

// Fill d with properties in the crack axis system depending on the type of damage
// tensor as specified by DForm
void TransIsoSoftening::LoadCrackAxisProperties(int np,CrackAxisProperties *d,int DForm,ElasticProperties *p) const
{
	// some for 3D only are cleared
	d->C44=0.;
	d->C55=0.;
	d->fXZLaw=NULL;
	d->vzxc=0.;
	d->vzyc=0.;
	
	if(DForm==DZ_DAMAGE)
	{	// Use Dz: crack axis system 1=T, 2=T, 3=A
		if(np==THREED_MPM)
		{	d->C11 = p->C[0][0];			// CTT = Cxx
			d->C22 = p->C[1][1];			// CTT = Cyy
			d->C33 = p->C[2][2];			// CAA = Czz
			d->C12 = p->C[0][1];			// CT,T = Cxy
			d->C13 = p->C[0][2];			// CT,A = Cxz
			d->C23 = p->C[1][2];			// CT,A = Cyz
			d->C44 = p->C[3][3];			// GA
		}
		else
		{	d->C11 = p->C[1][1];			// CTT = Cxx (Type 1 only)
			d->C22 = p->C[2][2];			// CTT = Cyy (Type 1 only)
			d->C12 = p->C[1][2];			// CT,T = Cxy (Type 1 only)
			d->C13 = p->C[4][1];			// C13 for sigmaxx in dsig.zz
			d->C23 = p->C[4][2];			// C23 for sigmayy in dsig.zz
			d->C33 = p->C[4][4];			// C33 for ezzr in dsig.zz
			d->vzxc = p->alpha[5];		    // for dsig.zz
			d->vzyc = p->alpha[6];		    // for dsig.zz
		}
		
		// ratios
		d->C12C11 = (nuT + nuA*nuAp)/(1. - nuA*nuAp);
		d->C13C11 = nuA*(1+nuT)/(1. - nuA*nuAp);
		
		// Laws
		d->fnLaw = softeningI;
		d->sigNc = initiationLaw->sigmaI();
		
		d->fXYLaw = softeningII;
		d->tauXYc = initiationLaw->sigmaII();
		d->C66 = (np==THREED_MPM) ? p->C[5][5] : p->C[3][3];
		
		if(np==THREED_MPM)
		{	d->fXZLaw = softeningAII;
			d->tauXZc = ((TIFailureSurface *)initiationLaw)->tauA();
			d->C55 = p->C[3][3];
		}
	}
	else if(DForm==DX_DAMAGE)
	{	// Use Dx: crack axis system 1=A, 2=T, 3=T
		if(np==THREED_MPM)
		{	d->C11 = p->C[2][2];			// CAA = Czz
			d->C22 = p->C[1][1];			// CTT = Cyy
			d->C33 = p->C[0][0];			// CTT = Cxx
			d->C12 = p->C[1][2];			// CA,T = Cyz
			d->C13 = p->C[0][2];			// CA,T = Cxz
			d->C23 = p->C[0][1];			// CT,T = Cxy
			d->C44 = p->C[5][5];			// GT
		}
		else
		{	d->C11 = p->C[2][2];			// CAA = Cyy (Type 2 only)
			d->C22 = p->C[1][1];			// CTT = Cxx (Type 2 only)
			d->C12 = p->C[1][2];			// CA,T = Cxy (Type 2 only)
			d->C13 = p->C[4][2];			// CAT (full) for sigmaxx in dsig.zz
			d->C23 = p->C[4][1];			// CTT (full) for sigmayy in dsig.zz
			d->C33 = p->C[4][4];			// CTT (full) for ezzr in dsig.zz
			d->vzxc = p->alpha[6];		    // for dsig.zz
			d->vzyc = p->alpha[5];		    // for dsig.zz
		}
		
		// ratios
		d->C12C11 = nuAp/(1.-nuT);
		d->C13C11 = d->C12C11;
		
		// laws
		d->fnLaw = softeningAI;
		d->sigNc = ((TIFailureSurface *)initiationLaw)->sigmaA();
		
		d->fXYLaw = softeningTII;
		d->tauXYc = ((TIFailureSurface *)initiationLaw)->tauT();
		d->C66 = p->C[3][3];			// C44 (3D) = C66 (2D)
		
		if(np == THREED_MPM)
		{	d->fXZLaw = softeningTII;
			d->tauXZc = ((TIFailureSurface *)initiationLaw)->tauT();
			d->C55 = p->C[3][3];			// C44 (3D) = C66 (2D)
		}
	}
	else
	{	// Use Dy: crack axis system 1=T, 2=A, 3=T (only 2D and 3 not used)
		d->C11 = p->C[1][1];			// CTT = Cxx (Type 2 only)
		d->C22 = p->C[2][2];			// CAA = Cyy (Type 2 only)
		d->C12 = p->C[1][2];			// CT,A = Cxy (Type 2 only)
		d->C13 = p->C[4][1];			// C13 for sigmaxx in dsig.zz
		d->C23 = p->C[4][2];			// C23 for sigmayy in dsig.zz
		d->C33 = p->C[4][4];			// C33 for ezzr in dsig.zz
		d->vzxc = p->alpha[5];		    // for dsig.zz
		d->vzyc = p->alpha[6];		    // for dsig.zz
		
		// ratios
		d->C12C11 = nuA*(1+nuT)/(1. - nuA*nuAp);
		d->C13C11 = (nuT + nuA*nuAp)/(1. - nuA*nuAp);
		
		// laws
		d->fnLaw = softeningI;
		d->sigNc = initiationLaw->sigmaI();
		d->C66 = p->C[3][3];		// C44 (3D) = C66 (2D)

		d->fXYLaw = softeningAII;
		d->tauXYc = ((TIFailureSurface *)initiationLaw)->tauA();
		d->C55 = d->C66;
	}
}

// Calculate rotation matrix from crack to material axes
bool TransIsoSoftening::GetRToCrack(Matrix3 *R,double *soft, bool is2D, int Dstyle) const
{	// none if undamaged
	if(soft[SOFT_DAMAGE_STATE] < predamageState) return false;
	
	// 3D or 2D
	if(!is2D)
	{	// only DX and DZ are possible
		if(Dstyle==DX_DAMAGE)
		{	// (0,pi/2,0) rotation about y
			R->set(0.,0.,1.,  0.,1., 0.,  -1.,0.,0.);
		}
		else
		{	// (theta,0,0) rotation about z
			double c1 = cos(soft[NORMALDIR1]);
			double s1 = sin(soft[NORMALDIR1]);
			R->set(c1,-s1,0.,  s1,c1,0.,  0.,0.,1.);
		}
	}
	else
	{	// already cos and sin
		R->set(soft[NORMALDIR1], -soft[NORMALDIR2], soft[NORMALDIR2], soft[NORMALDIR1], 1);
	}
	return true;
}

Tensor TransIsoSoftening::GetAnisoResStrains(double &exxr,double &eyyr,double &ezzr,
                                             ElasticProperties *p,ResidualStrains *res,int np) const
{
    if(np==THREED_MPM)
    {   exxr = p->alpha[0]*res->dT;
        eyyr = p->alpha[1]*res->dT;
        ezzr = p->alpha[2]*res->dT;
        if(fmobj->HasFluidTransport())
        {   exxr += p->beta[0]*res->dC;
            eyyr += p->beta[1]*res->dC;
            ezzr += p->beta[2]*res->dC;
        }
    }
    else
    {   // reduced residual strains when in plane strain
        // alpha[1] = alpha(xx)^(r), alpha[2]=alpha(yy)^r, but alpha[5]=alpha(zz)
        exxr = p->alpha[1]*res->dT;
        eyyr = p->alpha[2]*res->dT;
        ezzr = p->alpha[4]*res->dT;
        if(fmobj->HasFluidTransport())
        {   exxr += p->beta[1]*res->dC;
            eyyr += p->beta[2]*res->dC;
            ezzr += p->beta[4]*res->dC;
        }
    }
    Tensor er = MakeTensor2D(exxr, eyyr, ezzr, 0.);
    return er;
}

// A crack is present - get damage tensor form
// IT Softening: Let Di be damage tensor when axial direction of TI material is in the i direction of crack axis system
//         values: .75, .8, .85 using Dz, .95, 1.05 using Dx, 1.15,1.25 using Dy (2D only)
// Ortho softening: Let Di be damage tensor when i material direction is x direction of crack axis system
//         values: .75, .8, .85 using Dz, .95, 1.00, 1.05 using Dx, 1.15, 1.2, 1.25 using Dy (2D only)
int TransIsoSoftening::GetDForm(double damage) const
{
    // dgap = -.25 to -.15 is using Dz, -.05 to .05 is using Dx, .15 to .25 is using Dy
    // damageState is 1.5, predamageState is 0.5, constants
    double dgap = damage>damageState ? damage-damageState-0.5 :  damage-predamageState-0.5;
    int DForm;
    if(dgap<-0.1)
        DForm = DZ_DAMAGE;
    else if(dgap<0.1)
        DForm = DX_DAMAGE;
    else
        DForm = DY_DAMAGE;
    return DForm;
}

#pragma mark TransIsoSoftening::Accessors

// return material type
const char *TransIsoSoftening::MaterialType(void) const
{	if(AxialDirection()==AXIAL_Z)
		return "Tranversely isotropic softening (unrotated A axis along z axis)";
	else
		return "Tranversely isotropic softening (unrotated A axis along y axis)";
}

// store plastic strain in alt strain
int TransIsoSoftening::AltStrainContains(void) const
{	return ENG_BIOT_PLASTIC_STRAIN;
}

// return cracking strain and true, or return false if material has no cracking strain
bool TransIsoSoftening::GetCrackingStrain(MPMBase *mptr,Tensor *ecrack,bool is2D,Matrix3 *Rn) const
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
Vector TransIsoSoftening::GetCrackingCOD(MPMBase *mptr,bool is2D) const
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
double *TransIsoSoftening::GetSoftHistoryPtr(MPMBase *mptr) const
{	return (double *)(mptr->GetHistoryPtr(0));
}

// get traction failure surface (or <0 if not a softening material)
int TransIsoSoftening::GetTractionFailureSurface(void) const { return tractionFailureSurface; }
