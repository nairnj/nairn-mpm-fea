/********************************************************************************
	IsoPlasticSoftening.cpp
	nairn-mpm-fea
 
	Created by John Nairn, Oc. 9, 2017
	Copyright (c) 2017 John A. Nairn, All rights reserved.
  ********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Materials/IsoPlasticSoftening.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Materials/LinearHardening.hpp"
#include "Materials/FailureSurface.hpp"

// for debugging
#include "NairnMPM_Class/NairnMPM.hpp"
int nstepsMax=0;

#pragma mark IsoPlasticity::Constructors and Destructors

// Constructor
// throws std::bad_alloc
IsoPlasticSoftening::IsoPlasticSoftening(char *matName,int matID) : IsoSoftening(matName,matID)
{
	plasticLaw = new LinearHardening(this);
	softHistoryOffset = plasticLaw->HistoryDoublesNeeded();
}

#pragma mark IsoPlasticity::Initialization

// Read material properties
char *IsoPlasticSoftening::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	// look for different plastic law
	if(strcmp(xName,"Hardening")==0)
	{	input = HARDENING_LAW_SELECTION;
		return (char *)this;
	}
	
	// check plastic law
	char *ptr = plasticLaw->InputHardeningProperty(xName,input,gScaling);
	if(ptr != NULL) return ptr;
	
	// otherwise get material properties
	return(IsoSoftening::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *IsoPlasticSoftening::VerifyAndLoadProperties(int np)
{
	// call plastic law that is used
	const char *ptr = plasticLaw->VerifyAndLoadProperties(np);
	if(ptr != NULL) return ptr;
	
	//if(np==PLANE_STRESS_MPM)
	//	return "IsoPlasticSoftening materials cannot yet be used in plane stress calculations";
	
#ifdef POROELASTICITY
	if(DiffusionTask::HasPoroelasticity())
		return "IsoPlasticSoftening materials cannot yet be used when poroelasticity is activated";
#endif
	
	// check in superclass (along with its initialization)
	ptr = IsoSoftening::VerifyAndLoadProperties(np);
	
	// terms used in plasticity - these are undamaged properties
	// Cannot be used with plastic law that changes shear modulus
	Gred = C66/rho;
	Kred = C33/rho - 4.*Gred/3.;							// from C33 = lambda + 2G = K + 4G/3
    
    // only needed for plane stress
	kappa = (1.-2.*nu)/(1.-nu);                             // psRed is materials.tex notes
    psKred = kappa*Kred;
	
	return ptr;
}

// Allows any hardening law
bool IsoPlasticSoftening::AcceptHardeningLaw(HardeningLawBase *pLaw,int lawID)
{
	// don't allow laws that change shear modulus (SCGL and SL)
	if(lawID==4 || lawID==5) return false;
	
	// accept the law
	delete plasticLaw;
	plasticLaw = pLaw;
	softHistoryOffset = plasticLaw->HistoryDoublesNeeded();
	return true;
}

// return plastic law ID (or 0 if none)
HardeningLawBase *IsoPlasticSoftening::GetPlasticLaw(void) const { return plasticLaw; }

// print mechanical properties to the results
void IsoPlasticSoftening::PrintMechanicalProperties(void) const
{
	IsoSoftening::PrintMechanicalProperties();
	plasticLaw->PrintYieldProperties();
}

#pragma mark IsoPlasticity::History Data Methods

// The IsoPlasticity has no history data, but its plasticity law might
// History variables are stores as (plasticity),(Standard Softening),(Cracking Strains)
// Cracking strains history is needed because plastic strain used by plasticity
char *IsoPlasticSoftening::InitHistoryData(char *pchr,MPMBase *mptr)
{	double *p = CreateAndZeroDoubles(pchr,NumberOfHistoryDoubles());
    
    // soften variables (update if parent class updates)
    IsoSoftening::InitHistoryData((char *)p,mptr);
	
	// If has softening history variables, shift other variables by the size
	// of plastic law parameters. Only non-zero ones in parent class needing shift
	// 		are relative strength and toughness
    if(softHistoryOffset>0)
    {	// shift from the end
		p[softHistoryOffset+RELATIVE_TOUGHNESS] = p[RELATIVE_TOUGHNESS];
		p[softHistoryOffset+RELATIVE_STRENGTH] = p[RELATIVE_STRENGTH];
		p[softHistoryOffset+SOFT_DAMAGE_STATE] = p[SOFT_DAMAGE_STATE];
		
		// zero if needed from the beginning
		p[SOFT_DAMAGE_STATE] = 0.;
		p[RELATIVE_STRENGTH] = 0.;
		if(softHistoryOffset>1) p[RELATIVE_TOUGHNESS]=0.;
    }
    
    // do plastic law now, but stored in the beginning
	plasticLaw->InitPlasticHistoryData(p);
	return (char *)p;
}

// reset history data
void IsoPlasticSoftening::ResetHistoryData(char *pchr,MPMBase *mptr)
{	double *p = (double *)pchr;
	double relStrength = p[softHistoryOffset+RELATIVE_STRENGTH];
	double relToughness = p[softHistoryOffset+RELATIVE_TOUGHNESS];
	ZeroDoubles(pchr,NumberOfHistoryDoubles());
	p[softHistoryOffset+RELATIVE_STRENGTH] = relStrength;
	p[softHistoryOffset+RELATIVE_TOUGHNESS] = relToughness;
	p[softHistoryOffset+SOFT_DAMAGE_STATE] = 0.1;
}

// Number of history variables - plastic law plus softening
int IsoPlasticSoftening::NumberOfHistoryDoubles(void) const
{	return softHistoryOffset + SOFT_NUMBER_HISTORY + NUMBER_CRACKING_STRAINS;
}

#pragma mark IsoPlasticSoftening::Methods

// Constitutive Law and it is always in large rotation mode
void IsoPlasticSoftening::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,
											 ResidualStrains *res,int historyOffset,Tensor *gStress) const
{
    // set 2D flag
	bool is2D = np == THREED_MPM ? false : true;

#pragma mark ... Code to get trial stress and check for yielding
	// Task 1: Convert input form to effective strains in initial configuration
	//-------------------------------------------------------------------------
    Matrix3 dR,Rnm1,Rtot;
    Matrix3 deT = LRGetStrainIncrement(INITIAL_CONFIGURATION,mptr,du,&dR,NULL,&Rnm1,&Rtot);
	
	// residual strain (thermal and moisture)
	// Not that ezzres uses unreduced properties while eres is reduced in plain strain
	// For plane stress, 3D, and axisymetric, eres=ezzres
	double eres = CTE1*res->dT;
	double ezzres = CTE3*res->dT;
	if(DiffusionTask::HasFluidTransport())
	{	eres += CME1*res->dC;
		ezzres += CME3*res->dC;
	}
	double dispEnergy = 0.;
	
	// Rotate plastic strain in state n-1 to configuration n
	//-------------------------------------------------------
	Tensor *eplast=mptr->GetAltStrainTensor();
	*eplast = dR.RVoightRT(eplast,false,is2D);
	
	// Task 2: Get trial elastic for deviatoric stress
	//------------------------------------------------
	// cast pointer to material-specific data
	SofteningPlasticProperties *p = (SofteningPlasticProperties *)properties;

    // get history for softening properties (plastic,damage,cracking strains) storage mode
	// Returned value for this material is start of damage properties
	double *soft = GetSoftHistoryPtr(mptr);
	double *ipsoft = soft + SOFT_NUMBER_HISTORY;			// pointer to first cracking strain

#pragma mark ... Plasticity phase input calculations
    // Get de = deT - dec or strain increment in bulk and assuming no damage evolution
	// Here dec is cracking strain stored in history variables starting at ipsoft
	Tensor dec;
	if(soft[SOFT_DAMAGE_STATE]<predamageState)
	{	// will remain in initial configuration
		dec = MakeTensor(0.,0.,0.,0.,0.,0.);
		
		// Rnm1 and Rtot remain rotation from current state to initial configuration
	}
	else
	{	// When damaged, rotate total strain increment from initial axes (found above) to crack axis system
		Matrix3 RToCrack;
		GetRToCrack(&RToCrack, soft, is2D, 0);
		deT = deT.RTMR(RToCrack);
		
		// get cracking strain increments for an elastic update
		// by using ezzres, this works in plane strain (where deT(2,2)=0), 3D, and axisymmetric
		// plane stress is different here, but not yet supported elsewhere
		double den;
		if(np==PLANE_STRESS_MPM)
		{	den = deT(0,0)-ezzres + nu*(deT(1,1)-ezzres) ;
		}
		else
		{	double nuterm = nu/(1.-nu) ;
			den = deT(0,0)-ezzres + nuterm*(deT(1,1)+deT(2,2)-2.*ezzres) ;
		}
		
		// ... and make sure increment does not create negative normal strain
		double decxx = fmax(soft[DAMAGENORMAL]*den,-ipsoft[ECXX_DAMAGE]);
		double dxz = tractionFailureSurface == CUBOID_SURFACE && !is2D ? soft[DAMAGESHEAR2] : soft[DAMAGESHEAR];
		dec = MakeTensor(decxx,0.,0.,0.,dxz*(deT(0,2)+deT(2,0)),soft[DAMAGESHEAR]*(deT(0,1)+deT(1,0)));
		
		// get rotation from current state to crack axis system
		Rnm1 *= RToCrack;
		Rtot *= RToCrack;
	}
	
	// get de = detot-dec (this includes residual strain still)
	Tensor de = is2D ? MakeTensor2D(deT(0,0)-dec.xx, deT(1,1), deT(2,2), deT(0,1)+deT(1,0)-dec.xy) :
			MakeTensor(deT(0,0)-dec.xx, deT(1,1), deT(2,2), deT(1,2)+deT(2,1), deT(0,2)+deT(2,0)-dec.xz, deT(0,1)+deT(1,0)-dec.xy);
	
	// for plasticity phase de(trial)eff = de - ezzres*I
	
	// get deviatoric stress increment (in absence of cracking strain increments and residual strains)
    double delV,delP,dTq0;
    Tensor dels;
    if(np==PLANE_STRESS_MPM)
    {	// see plane stress platicity notes in materials.pdf
		delV = kappa*(de.xx+de.yy-2.*eres);
        double thirdDelV = delV/3.;
        dTq0 = -gamma0*mptr->pPreviousTemperature*(delV+3.*eres);

        // deviatoric stress and pressure increments
		// nets to -(E/(1-nu))(dexxeff+deyyeff)/3
        delP = -Kred*delV;
        dels = MakeTensor2D(2.*Gred*(de.xx-eres-thirdDelV), 2.*Gred*(de.yy-eres-thirdDelV), delP, Gred*de.xy);
    }
    else
    {	// see plane strain platicity notes in materials.pdf
		delV = de.xx+de.yy+de.zz;
        dTq0 = -gamma0*mptr->pPreviousTemperature*delV;
        double thirdDelV = delV/3.;				// with residual strains to find deviatoric stresses
        delV -= 3.*ezzres;						// to get pressure (Needs unreduced here, use ezzres in case plane strain)
        
        // deviatoric stress and pressure increments
        dels = is2D ?
            MakeTensor2D(2.*Gred*(de.xx-thirdDelV), 2.*Gred*(de.yy-thirdDelV), 2.*Gred*(de.zz-thirdDelV), Gred*de.xy) :
            MakeTensor(2.*Gred*(de.xx-thirdDelV), 2.*Gred*(de.yy-thirdDelV), 2.*Gred*(de.zz-thirdDelV),
                                       Gred*de.yz, Gred*de.xz, Gred*de.xy);
        delP = -Kred*delV;
    }
    
	// dels is deviatoric stress increment in crack axis system (which will be initial
	// axes if not damaged). We need to add to deviatoric stress in the same axis system
	
	// rotate previous stress to initial (if undamaged) or crack axis system (if damaged)
	Tensor *sp = mptr->GetStressTensor();
	Tensor st0 = Rnm1.RTVoightR(sp,true,is2D);
	
	// trial deviatoric stress in crack axis system = st0 + P0*I + dels
	double initialP = -(st0.xx+st0.yy+st0.zz)/3.;
	Tensor str = is2D ?
		MakeTensor2D(st0.xx+initialP+dels.xx,st0.yy+initialP+dels.yy,st0.zz+initialP+dels.zz,st0.xy+dels.xy) :
		MakeTensor(st0.xx+initialP+dels.xx,st0.yy+initialP+dels.yy,st0.zz+initialP+dels.zz,
                   		st0.yz+dels.yz,st0.xz+dels.xz,st0.xy+dels.xy);
	double trialP = initialP+delP;
	double trialP0 = trialP;

	// Task 3: Get plastic potential with trial str stress in initial (undamaged) or crack axis (if damaged) system
	//-------------------------------------------------------------------------------------------------------------
	// Calculate plastic potential f = ||s|| - sqrt(2/3)*sy(alpha,rate,...)
	HardeningAlpha alpha;
	plasticLaw->UpdateTrialAlpha(mptr,np,&alpha,0);			// initialize to last value and zero plastic strain rate
	
	// in-line coding of GetPlasticPotential() in IsoPlasticity material
	alpha.strialmag = GetMagnitudeSFromDev(&str,np);
	double fk = alpha.strialmag - SQRT_TWOTHIRDS*plasticLaw->GetYield(mptr,np,delTime,&alpha,p->hardProps);
	
#pragma mark ... Code for elastic or plastic update in plasticity phase
	// Task 4: Check potential and then finish plasticity phase
	//----------------------------------------------------------------------------------------
	Tensor strPlastic;
	if(fk<=0.)
	{	// Task 4a: No plasticity effect on final stress and then get total trial stress
		// -----------------------------------------------------------------------------
		ZeroTensor(&strPlastic);
		str.xx -= trialP;
		str.yy -= trialP;
		str.zz -= trialP;

		// add back cracking strain to get de(elastic) = de(total) = de(trial)+dec (inclusive of residual strains)
		AddTensor(&de,&dec);
		
		// give subclass material chance to update history variables that change in elastic updates
		plasticLaw->ElasticUpdateFinished(mptr,np,delTime,0);
	}
	
	else
	{	// Task 4b: Plasticity - Use Radial Return to find lambdak
		//---------------------------------------------------------
		double lambdak = plasticLaw->SolveForLambdaBracketed(mptr,np,alpha.strialmag,&str,
                                                             Gred,psKred,trialP,delTime,&alpha,p->hardProps,0);

		// Task 4c: plastic strain increment in crack axis system
        //--------------------------------------
        Tensor dep;
        if(np==PLANE_STRESS_MPM)
        {    // Note this plane stress codes assumes J2 plasticity, should generalize when write a subclass
            double d1 = (1. + psKred*lambdak);
            double d2 = (1.+2.*Gred*lambdak);
            double n1 = (str.xx+str.yy-2.*trialP)/d1;
            double n2 = (-str.xx+str.yy)/d2;
            double sxx = (n1-n2)/2.;
            double syy = (n1+n2)/2.;
            double txy = str.xy/d2;
            
			// get final direction (it is not constant in plane stress)
            Tensor dfds;
            dfds.xx = (2.*sxx-syy)/3.;
            dfds.yy = (2.*syy-sxx)/3.;
            dfds.zz = -(dfds.xx+dfds.yy);
            dfds.xy = txy;                // tensorial shear strain
            
            // now the pressure has changed from trialP to new value
            trialP = -n1/3.;                    // final P in plane stress calculations
            
            // -----------------------------------------------
            // dep = lambdak*n
            dep.xx = lambdak*dfds.xx;
            dep.yy = lambdak*dfds.yy;
            dep.zz = lambdak*dfds.zz;
            dep.xy = 2.*lambdak*dfds.xy;     // 2 for engineering plastic shear strain
            dep.xz = 0.;
            dep.yz = 0.;

            // get dezz required to give zero stress in z direction
            // separate out just the dezzp term. Later code will add elastic change
            //  based on de-dep in x and y directions
            mptr->IncrementDeformationGradientZZ(dep.zz);
 
            // Add plastic part to isotropic temperature increment here
            dTq0 -= gamma0*mptr->pPreviousTemperature*dep.zz;

            // Task 4d: Get decrease in deviatoric stress caused by plastic strain
            //    which is equal to ds = strial - s(final) and s(final) = sigma+P
            //-------------------------------------------------------
            strPlastic.xx = str.xx - (sxx + trialP);
            strPlastic.yy = str.yy - (syy + trialP);
            strPlastic.zz = str.zz - trialP;				// = trialP0 - trialP
            strPlastic.xy = str.xy - txy;
            strPlastic.xz = 0.;
            strPlastic.yz = 0.;
        }
        else
        {   // -----------------------------------------------
            // dep = lambdak*(strial/||strial||)
            double pscale = lambdak/alpha.strialmag;
            dep.xx = pscale*str.xx;
            dep.yy = pscale*str.yy;
            dep.zz = pscale*str.zz;
            dep.xy = 2.*pscale*str.xy;     // 2 for engineering plastic shear strain
            if(is2D)
            {   dep.xz = 0.;
                dep.yz = 0.;
            }
            else
            {   dep.xz = 2.*pscale*str.xz;
                dep.yz = 2.*pscale*str.yz;
            }
            
            // Task 4d: Get decrease in stress caused by plastic strain
            //-------------------------------------------------------
            strPlastic.xx = 2.*Gred*dep.xx;
            strPlastic.yy = 2.*Gred*dep.yy;
            strPlastic.zz = 2.*Gred*dep.zz;
            strPlastic.xy = Gred*dep.xy;
            if(is2D)
            {   strPlastic.xz = 0.;
                strPlastic.yz = 0.;
            }
            else
            {   strPlastic.xz = Gred*dep.xz;
                strPlastic.yz = Gred*dep.yz;
            }
        }
		
		// Task 4e: Increment all energies
		//--------------------------------
		// get final deviatoric stress str = str - strPlastic
		// get plastic work from deviatoric stress increment and plastic strains
		SubTensor(&str,&strPlastic);
		double plastEnergy = str.xx*dep.xx + str.yy*dep.yy + str.zz*dep.zz + str.xy*dep.xy;
		if(!is2D) plastEnergy += str.xz*dep.xz + str.yz*dep.yz;
		
		// and subtract q dalpha to get dissipated energy per unit mass
		double qdalphaTerm = lambdak*SQRT_TWOTHIRDS*plasticLaw->GetYieldIncrement(mptr,np,delTime,&alpha,p->hardProps);
		dispEnergy = plastEnergy - qdalphaTerm;
		mptr->AddPlastEnergy(dispEnergy);
		
		// The cumulative dissipated energy is tracked in plastic energy
		
		// Total work energy is sigma.de and remaining code will add only sigma.(de-dep)
		//     To compenstate, we add sigma.dep here. First get total stress
		str.xx -= trialP;
		str.yy -= trialP;
		str.zz -= trialP;                   // will be zero in plane stress

		// correct plastic stress change if in plane stress
		if(np==PLANE_STRESS_MPM)
		{	double deltaP = trialP - trialP0;
			strPlastic.xx += deltaP;
			strPlastic.yy += deltaP;
			strPlastic.zz = 0.;
		}

		double partialWorkEnergy = str.xx*dep.xx + str.yy*dep.yy + str.xy*dep.xy + str.zz*dep.zz;
		if(!is2D) partialWorkEnergy += str.yz*dep.yz + str.xz*dep.xz;
		mptr->AddWorkEnergy(partialWorkEnergy);
		
		// get elastic strain increment (de(elastic) = de(total)-dep = de(trial)+dec-dep) in crack axis system
		AddTensor(&de,&dec);
		SubTensor(&de,&dep);
		
		// rotate plastic strain to current state and add to particle
		dep = Rtot.RVoightRT(&dep,false,is2D);
		AddTensor(eplast,&dep);
		
		// Task 4f: Update Internal Variables
		//----------------------------------
		plasticLaw->UpdatePlasticInternal(mptr,np,&alpha,0);
	}
	
#pragma mark ... Code for Previously Undamaged Material
	// Task 5: Begin damage mechanics phase
	//-------------------------------------
	// These values are in initial (undamaged) or crack axis (damaged) system
	// 1. The trial total stress after plasticity Tensor str
	// 2. The decrease in trial stresses caused by plasticity in strPlastic
	// 3. Total elastic strain increment Tensor de(elastic) = de(total)-dep
    //          (de.zz=0 in plane stress (set below), -dezzp in plane strain, axisymmetric or 3D)

	// Before cracking, do normal istropic update. If not cracked
	// then done, otherwise initiate the damage
	if(soft[SOFT_DAMAGE_STATE]<predamageState)
	{	// check if has failed
		Vector norm;
		double relStrength = soft[RELATIVE_STRENGTH];
		int failureMode = initiationLaw->ShouldInitiateFailure(&str,&norm,np,relStrength,NULL);
		//failureMode = NO_FAILURE;    // Hack to verify plasticity part matches IsoPlasticity material
		if(failureMode == NO_FAILURE)
		{	// Not failed yet, so finish elastic-plastic update

			// rotate to current axes and update stress
			if(np==THREED_MPM)
			{	*sp = Rtot.RVoightRT(&str,true,false);
				
				// work energy increment per unit mass (dU/(rho0 V0)) (in initial axes)
				mptr->AddWorkEnergyAndResidualEnergy(str.xx*de.xx + str.yy*de.yy + str.zz*de.zz
													 + str.xz*de.yz + str.xz*de.xz + str.xy*de.xy,
													 (str.xx + str.yy + str.zz)*eres);
			}
			
			else
			{   // rotate str to current configuration (but keep str in analysis frame)
                Tensor stnp1 = Rtot.RVoightRT(&str,true,true);
				sp->xx = stnp1.xx;
				sp->yy = stnp1.yy;
				sp->xy = stnp1.xy;
				sp->zz = stnp1.zz;
				
				// work and residual strain energy increments (in initial axes)
				double workEnergy = str.xx*de.xx + str.yy*de.yy + str.xy*de.xy;
				double resEnergy = (str.xx + str.yy)*ezzres;
 				if(np==PLANE_STRAIN_MPM)
				{	// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
					resEnergy += sp->zz*ezzres;
				}
				else if(np==PLANE_STRESS_MPM)
                {   // add elastic part in de.zz (because exiting now), plastic part done above
                    de.zz = -(nu/(1.-nu))*(de.xx + de.yy - 2.*eres) + eres;       // elastic ezz part only
                    mptr->IncrementDeformationGradientZZ(de.zz);
                }
                else
				{	// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
					workEnergy += sp->zz*de.zz;
					resEnergy += sp->zz*eres;
				}
				mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
			}
			
			// track heat energy
			IncrementHeatEnergy(mptr,dTq0,dispEnergy);
			
			// this elastic-plastic-no damage update is done
			// plastic strain was rotated above and stored on particle
			return;
		}
        
		// initiate failure - initiate state depends on initiation law
		soft[SOFT_DAMAGE_STATE] = 0.01*failureMode ;
		soft[NORMALDIR1] = norm.x;										// cos(theta) (2D) or Euler alpha (3D) to normal
		soft[NORMALDIR2] = norm.y;										// sin(theta) (2D) or Euler beta (3D)  to normal
		soft[NORMALDIR3] = norm.z;										// unused (2D) or Euler gamma (3D) to normal
		
		// rotation matrix to the new crack axis system
		Matrix3 RToCrack;
		GetRToCrack(&RToCrack, soft, is2D, 0);
		
		// get intersection area (for 3D convert angle to normal)
		if(np==THREED_MPM) norm = MakeVector(RToCrack(0,0),RToCrack(1,0),RToCrack(2,0));
		soft[GCSCALING] = GetAcOverVp(np,mptr,&norm)/rho;			     // rho because is divided by specific stress
		
		// Rotate elastic strain increment from initial to crack axis system
		de = RToCrack.RTVoightR(&de,false,is2D);
		
		// Rotate plastic stress change from initial to crack axis system
		strPlastic = RToCrack.RTVoightR(&strPlastic,false,is2D);
		
		// include rotation in total rotation matrices
		Rnm1 *= RToCrack;
		Rtot *= RToCrack;
	}
	
#pragma mark ... Code for Damaged Material
	// A crack is present
	
	// Elastic strain increment in crack axis system is in de = de(total)-dep
    if(np==PLANE_STRESS_MPM) de.zz=0.;          // Code below assumes plain stress starts at zero

	// Trial elastic-plastic stress in crack axis system after plasticity update is in str
	str = Rnm1.RTVoightR(sp,true,is2D);
    
	if(fk>0.)
	{	// subtract plastic stress from current stress in crack axis system
		SubTensor(&str,&strPlastic);
	
		// subtract plastic stress change from n-1 state currently on the particle
		// DamageEvolution() will finish rotation to current state
		strPlastic = Rnm1.RVoightRT(&strPlastic,true,is2D);
		SubTensor(sp,&strPlastic);
	}
	
	// Get cracking strain in crack axis system
	Tensor ecrack = MakeTensor(ipsoft[ECXX_DAMAGE],0.,0.,0.,ipsoft[GCXZ_DAMAGE],ipsoft[GCXY_DAMAGE]);
    
	// finish up in separate code
    DamageEvolution(mptr,np,soft,de,str,eres,ezzres,res,dR,Rtot,(ElasticProperties *)(p->elasticProps),
                    &ecrack,dispEnergy,delTime);
}

// This material stores cracking strain in history variable
// The cracking strain is stored in the crack axis system
void IsoPlasticSoftening::UpdateCrackingStrain(int np,Tensor *ecrack,double decxx,double dgcxy,double dgcxz,Matrix3 Rtot,double *soft) const
{
	double *ipsoft = soft + SOFT_NUMBER_HISTORY;			// pointer to first cracking strain

	// update cracking strain
	ipsoft[ECXX_DAMAGE] += decxx;
	ipsoft[GCXY_DAMAGE] += dgcxy;
	ipsoft[GCXZ_DAMAGE] += dgcxz;
}

#pragma mark IsoPlasticity::Accessors

// store plastic strain in alt strain
int IsoPlasticSoftening::AltStrainContains(void) const
{	return ENG_BIOT_PLASTIC_STRAIN;
}

// return cracking strain and true, or return false if material has no cracking strain
// (don't change inRtot)
bool IsoPlasticSoftening::GetCrackingStrain(MPMBase *mptr,Tensor *ecrack,bool is2D,Matrix3 *inRtot) const
{
	// get history
	double *soft = GetSoftHistoryPtr(mptr);
	
	// Before initiatiation, no cracking strain
	if(soft[SOFT_DAMAGE_STATE]<predamageState) return false;
	
	// cracking strain in history variable, but in crack axis system
	double *ipsoft = soft + SOFT_NUMBER_HISTORY;
	Tensor eCAS = MakeTensor(ipsoft[ECXX_DAMAGE],0.,0.,0.,ipsoft[GCXZ_DAMAGE],ipsoft[GCXY_DAMAGE]);
	
	// Get R to crack, rotate by Rtot*RToCrack
    Matrix3 Rtot = *inRtot;
	Matrix3 RToCrack;
	GetRToCrack(&RToCrack,soft,is2D,0);
	Rtot *= RToCrack;
	*ecrack = Rtot.RVoightRT(&eCAS,false,is2D);
	
	// has cracking strain
	return true;
}

// return cracking COD in the crack axis system
// Only used by DeleteDamaged custom task
Vector IsoPlasticSoftening::GetCrackingCOD(MPMBase *mptr,bool is2D) const
{
	Vector cod = MakeVector(0.,0.,0.);
	
	// get history
	double *soft = GetSoftHistoryPtr(mptr);
	
	// Before initiatiation, not damage strain
	if(soft[SOFT_DAMAGE_STATE]<predamageState) return cod;
	
	// cracking strain in history variable is already in the crack axis system
	double *ipsoft = soft + SOFT_NUMBER_HISTORY;
	cod.x = ipsoft[ECXX_DAMAGE];
	cod.y = 0.5*ipsoft[GCXY_DAMAGE];
	cod.z = 0.5*ipsoft[GCXZ_DAMAGE];
	
	// scale by Vp/Ac
	ScaleVector(&cod,1./(rho*soft[GCSCALING]));
	
	// has cracking strain
	return cod;
}

// buffer size for mechanical properties
int IsoPlasticSoftening::SizeOfMechanicalProperties(int &altBufferSize) const
{   altBufferSize = plasticLaw->SizeOfHardeningProps();
	return sizeof(SofteningPlasticProperties);
}

// Isotropic material can use read-only initial properties
void *IsoPlasticSoftening::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{
	SofteningPlasticProperties *p = (SofteningPlasticProperties *)matBuffer;
	p->hardProps = plasticLaw->GetCopyOfHardeningProps(mptr,np,altBuffer,offset);
	
	// elastic properties
	p->elasticProps = (void *)&pr;
	
	// should disallow hardening law that changes shear modulus
	//double Gratio = plasticLaw->GetShearRatio(mptr,mptr->GetPressure(),1.,p->hardProps,offset);
	//p->Gred = G0red*Gratio;
	
	return p;
}

// return unique, short name for this material
const char *IsoPlasticSoftening::MaterialType(void) const { return "Isotropic Plastic Softening"; }

// history data for softening after plastic history
double *IsoPlasticSoftening::GetSoftHistoryPtr(MPMBase *mptr) const
{	double *history =  (double *)(mptr->GetHistoryPtr(0));
	return history+softHistoryOffset;
}
