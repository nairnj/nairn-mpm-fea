/********************************************************************************
	IsoDamageMech.cpp
	nairn-mpm-fea

	Created by John Nairn on June 26, 2015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Materials/IsoDamageMech.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Materials/FailureSurface.hpp"
#include "Materials/SofteningLaw.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "System/UnitsController.hpp"
#include "Boundary_Conditions/InitialCondition.hpp"

#include "NairnMPM_Class/NairnMPM.hpp"
#include "System/ArchiveData.hpp"

#pragma mark IsoDamageMech::Constructors and Destructors

// Constructor
// throws std::bad_alloc
IsoDamageMech::IsoDamageMech(char *matName,int matID) : IsoSoftening(matName,matID)
{
	// Initiation law may not always be used, but its mode I will be used
	evolutionStyle = STRESS_ENERGY_METRIC;
	
	// mode I and mode II left intact, but only mode I is used for some
	// Both used in my "improved" approach
	
	// To use my new method, add maximum principle stress law back
	// Enter mode I and mode II softening laws
}

#pragma mark IsoDamageMech::Initialization

// Read material properties - need to pick evolution style
char *IsoDamageMech::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	if(strcmp(xName,"metric")==0)
	{	input=INT_NUM;
		return (char *)&evolutionStyle;
	}
	
	// rest are to handle pressure dependence
	return IsoSoftening::InputMaterialProperty(xName,input,gScaling);
}

const char *IsoDamageMech::VerifyAndLoadProperties(int np)
{
	// verify valid style
	if(evolutionStyle<STRESS_ENERGY_METRIC || evolutionStyle>=LAST_METRIC)
		return "IsoDamageMechanics metric set to an invalid option";

    // ditch mode II if not begin used
	if(evolutionStyle!=MIXED_MODE_METRIC)
	{	delete softeningModeII;
		softeningModeII = NULL;
	}
	else
	{	// if update GetTractionFunction() for ovoid surface style, pick that one here instead
		tractionFailureSurface = CUBOID_SURFACE;
	}
	
	// verify using isosoftening
	const char *msg = IsoSoftening::VerifyAndLoadProperties(np);
	if(msg!=NULL) return msg;
	
#ifdef POROELASTICITY
	// diffusion is OK, but poroelasticity terms not yet implemented
	if(DiffusionTask::HasPoroelasticity())
		return "IsoDamageMechanics cannot be used with poroelasticity";
#endif
	
	Gred = G/rho;
	Kred = E/(3.*(1-2*nu)*rho);
	Ered = E/rho;
	
	return NULL;
}

// print mechanical properties to the results
void IsoDamageMech::PrintMechanicalProperties(void) const
{
    IsotropicMat::PrintMechanicalProperties();
	if(evolutionStyle==MIXED_MODE_METRIC)
	{	// set to use simple principle stress conditions
		initiationLaw->SetFailureSurface(CUBOID_SURFACE);
		cout << "Failure Surface: ";
		initiationLaw->PrintInitiationProperties();
		cout << "Mode I Softening: ";
		softeningModeI->PrintSofteningProperties(rho*initiationLaw->sigmaI());
		cout << "Mode II Softening: ";
		softeningModeII->PrintSofteningProperties(rho*initiationLaw->sigmaII());
	}
	else
	{	if(evolutionStyle==STRESS_ENERGY_METRIC)
			cout << "Failure Surface: Stress metric" << endl;
		else
			cout << "Failure Surface: Tensile principal stress metric" << endl;
		double sigmac = initiationLaw->sigmaI()*rho;
		MaterialBase::PrintProperty("sigc",sigmac*UnitsController::Scaling(1.e-6),"");
		cout << "Damage Softening: ";
		softeningModeI->PrintSofteningProperties(sigmac);
	}
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
}

#pragma mark IsoDamageMech::History Data

// Initialize damage when requested
void IsoDamageMech::SetInitialConditions(InitialCondition *ic,int ptNum,bool is3D)
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

        // Find Ac/Vp
        soft[GCSCALING] = GetAcOverVp(THREED_MPM, mptr, &dnorm) / rho;
        
        // Any damage variable >=1 means failed, otherwise calculcate corresponding cracking strains
        // and set provide failure mode
        double relToughness = soft[RELATIVE_TOUGHNESS];
        double relStrength = soft[RELATIVE_STRENGTH];
        double scaleI = relToughness*soft[GCSCALING]/(relStrength*initiationLaw->sigmaI());
        if(dparams.x>=1.)
        {   soft[DAMAGENORMAL] = 1.;
            soft[SOFT_DAMAGE_STATE] = 4.;
             soft[DELTANORMAL] = softeningModeI->GetDeltaMax(scaleI);
            if(evolutionStyle==MIXED_MODE_METRIC)
            {   double scaleII = relToughness*soft[GCSCALING]/(relStrength*initiationLaw->sigmaII());
                soft[DELTASHEAR] = softeningModeII->GetDeltaMax(scaleII);
            }
        }
        else
        {   soft[SOFT_DAMAGE_STATE] = dmode;
            double relen0 = relStrength*en0;
            soft[DELTANORMAL] = softeningModeI->GetDeltaFromDamage(dparams.x,scaleI,relen0,-1.);
            if(evolutionStyle==MIXED_MODE_METRIC)
            {   double scaleII = relToughness*soft[GCSCALING]/(relStrength*initiationLaw->sigmaII());
                double relgs0 = relStrength*gs0;
                soft[DELTASHEAR] = softeningModeII->GetDeltaFromDamage(dparams.y,scaleII,relgs0,-1.);
            }
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

        // Find Ac/Vp (first parameter anything except THREED_MPM)
        soft[GCSCALING] = GetAcOverVp(PLANE_STRAIN_MPM, mptr, &dnorm) / rho;
        
        // Any damage variable >=1 means failed, otherwise calucate corresponding cracking strains
        // and set provide failure mode
        double relToughness = soft[RELATIVE_TOUGHNESS];
        double relStrength = soft[RELATIVE_STRENGTH];
        double scaleI = relToughness*soft[GCSCALING]/(relStrength*initiationLaw->sigmaI());
        if(dparams.x>=1. || dparams.y>=1.)
        {   soft[DAMAGENORMAL] = 1.;
            soft[SOFT_DAMAGE_STATE] = 4.;
            soft[DELTANORMAL] = softeningModeI->GetDeltaMax(scaleI);
            if(evolutionStyle==MIXED_MODE_METRIC)
            {   double scaleII = relToughness*soft[GCSCALING]/(relStrength*initiationLaw->sigmaII());
                soft[DELTASHEAR] = softeningModeII->GetDeltaMax(scaleII);
            }
        }
        else
        {   soft[SOFT_DAMAGE_STATE] = dmode;
            double relen0 = relStrength*en0;
            soft[DELTANORMAL] = softeningModeI->GetDeltaFromDamage(dparams.x,scaleI,relen0,-1.);
            if(evolutionStyle==MIXED_MODE_METRIC)
            {   double scaleII = relToughness*soft[GCSCALING]/(relStrength*initiationLaw->sigmaII());
                double relgs0 = relStrength*gs0;
                soft[DELTASHEAR] = softeningModeII->GetDeltaFromDamage(dparams.y,scaleII,relgs0,-1.);
            }
        }
    }
}

#pragma mark IsoDamageMech::Methods

/* Take increments in strain and calculate new Particle: strains, rotation strain,
		stresses, strain energy,
	dvij are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
*/
void IsoDamageMech::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,
									  ResidualStrains *res,int historyOffset,Tensor *gStress) const
{	// set 2D flag
	bool is2D = np == THREED_MPM ? false : true;

	// current previous deformation gradient and stretch
	Matrix3 pFnm1 = mptr->GetDeformationGradientMatrix();
	
	// get incremental deformation gradient and decompose it
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
	
	// Update total deformation gradient
	Matrix3 pF = dF*pFnm1;
	mptr->SetDeformationGradientMatrix(pF);
	
	// decompose to get previous Rn and Rn-1
	Matrix3 dR,Rnm1,Rtot;
	
	// two decompositions
	Matrix3 Unm1 = pFnm1.RightDecompose(&Rnm1,NULL);
	pF.LeftDecompose(&Rtot,NULL);
	
	// Rtot = dR*Rnm1 or dR = Rtot*Rnm1^T
	dR = Rtot*Rnm1.Transpose();
	
	// get strain increments in initial configuration Rtot^T(dF-dR)F(n-1)
	Matrix3 dFmdR = dF - dR;
	Matrix3 deT = Rtot.Transpose()*(dFmdR*pFnm1);
    
    // convert matrix to Tensor in initial configuration
	// Note: plane stress will have zero in deT(2,2) to start)
	Tensor de = is2D ?
		MakeTensor2D(deT(0,0), deT(1,1), deT(2,2), deT(0,1) + deT(1,0)) :
		MakeTensor(deT(0,0), deT(1,1), deT(2,2), deT(1,2)+deT(2,1), deT(0,2)+deT(2,0), deT(0,1)+deT(1,0));

	// initial stresses
	Tensor *sp = mptr->GetStressTensor();
	
	// residual strain (thermal and moisture)
	// (CTE1 and CME1 are reduced to plane strain, but not CTE3 and CME3)
	double eres = CTE1*res->dT;
	double ezzres = CTE3*res->dT;
	double dTemp = mptr->pPreviousTemperature-thermal.reference;
	double resStrain = CTE3*dTemp;
	if(DiffusionTask::HasDiffusion())
	{	// Only add for diffusion. Poroelasticity done latter so damage can use solid stress
		eres += CME1*res->dC;
		ezzres += CME3*res->dC;
		double dConc = diffusion->GetDeltaConcentration(mptr);
		resStrain += CME3*dConc;
	}
	
	// Load str with stress in initial configuration used to look for initiation of failure
	// Need initial, because strains were found in initial configuration
	Tensor strn = Rnm1.RTVoightR(sp,true,is2D);

	// previous effective strain in initial configuration (Unm1-(1+resStrain)I)
	deT = Unm1;
	deT(0,0) -= (1+resStrain);
	deT(1,1) -= (1+resStrain);
	deT(2,2) -= (1+resStrain);
	Tensor enm1 = is2D ?
		MakeTensor2D(deT(0,0), deT(1,1), deT(2,2), deT(0,1) + deT(1,0)) :
		MakeTensor(deT(0,0), deT(1,1), deT(2,2), deT(1,2)+deT(2,1), deT(0,2)+deT(2,0), deT(0,1)+deT(1,0));

	// pointer to properties (which are undamaged here)
	ElasticProperties *p = (ElasticProperties *)properties;

#pragma mark ... Code for Previously Undamaged Material
	// get history
	double *soft = GetSoftHistoryPtr(mptr);
	
	// get plane stress zz effective strain and total strain increment
	if(np==PLANE_STRESS_MPM)
	{	// plane stress increment in ezz and previous ezz
		double omDn = 1.-soft[DAMAGENORMAL];
		de.zz = omDn*(p->C[4][1]*(de.xx-ezzres) + p->C[4][2]*(de.yy-ezzres)) + ezzres;
		// previous plane stress ezz was already found, but calculate anyway to be sure
		enm1.zz = omDn*(p->C[4][1]*enm1.xx + p->C[4][2]*enm1.zz);
	}
	
	// Before damage, do normal istropic update. If not damaged yet then done
	// then done, otherwise initiate the damage
	Tensor deij;
	if(soft[SOFT_DAMAGE_STATE]<predamageState)
	{	// Increment with reduced residual strain here
		Tensor deres = MakeTensor(eres,eres,eres,0.,0.,0.);
		deij = de;
		SubTensor(&deij,&deres);

		// Now get stress increment in initial configuration
		Tensor dstr = GetStressIncrement(deij,np,properties);

		// Add incremental stress
		Tensor str = strn;
		AddTensor(&str,&dstr);
		
		// subsequent code needs unreduced 3D strain increments
		deres = MakeTensor(ezzres,ezzres,ezzres,0.,0.,0.);
		deij = de;
		SubTensor(&deij,&deres);

		// check if has failed
		Vector norm;
		double relStrength = soft[RELATIVE_STRENGTH];
		int failureMode;
		if(evolutionStyle==MIXED_MODE_METRIC)
			failureMode = initiationLaw->ShouldInitiateFailure(&str,&norm,np,relStrength,NULL);
		else
		{	// stress metric needs zz stress to work in general
			if(np==PLANE_STRAIN_MPM)
				str.zz += p->C[4][1]*(deij.xx-ezzres) + p->C[4][2]*(deij.yy-ezzres) - p->C[4][4]*ezzres;
			else if(np==AXISYMMETRIC_MPM)
			{	str.zz += p->C[4][1]*(deij.xx-eres) + p->C[4][2]*(deij.yy-eres) +
								p->C[4][4]*(deij.zz-eres);
			}
			failureMode = ShouldInitiateFailure(&str,relStrength,&norm,np);
		}
		if(failureMode == NO_FAILURE)
		{	// Store str on the particle (rotating by Rtot to update current config
			AcceptTrialStress(mptr,str,sp,np,&Rtot,properties,de,eres,ezzres);
			return;
		}
		
#pragma mark ...... Undamaged Material Just Initiated Damage
#pragma omp critical (output)
		// initiate failure (code depends on initiation law)
		soft[SOFT_DAMAGE_STATE] = 0.01*failureMode ;
		soft[NORMALDIR1] = norm.x;										// cos(theta) (2D) or Euler alpha (3D) to normal
		soft[NORMALDIR2] = norm.y;										// sin(theta) (2D) or Euler beta (3D)  to normal
		soft[NORMALDIR3] = norm.z;										// unused (2D) or Euler gamma (3D) to normal
		
		// get intersection area (for 3D convert angle to normal first)
		if(np==THREED_MPM)
		{	Matrix3 RInit;
			GetRToCrack(&RInit, soft, is2D, 0);
			norm = MakeVector(RInit(0,0),RInit(1,0),RInit(2,0));
		}

		soft[GCSCALING] = GetAcOverVp(np,mptr,&norm)/rho;			// rho because is divided by specific stress
	}

	// Now need unreduced 3D strain increments
	Tensor deres = MakeTensor(ezzres,ezzres,ezzres,0.,0.,0.);
	deij = de;
	SubTensor(&deij,&deres);
	
#pragma mark ... Code for Damaged Material
	// Damage is present
	
	// enm1 = initial effective strain (unreduced thermal strains)
	// deij = effective strain increment in initial configuration (unreduced residual strains)
	// plane stress includes ezz, plane strain may have non-zero ezz if residual strains
	
	// strn = initial stress in initial configuration (time step n stress)
	
	// apply incremental rotation to cracking strain on the particle (in eplast)
	Tensor *eplast=mptr->GetAltStrainTensor();		// in global axes
	*eplast = dR.RVoightRT(eplast,false,is2D);
	
	// final cracking strain and strress increments
	Tensor dec;
	Tensor dsig = MakeTensor(0.,0.,0.,0.,0.,0.);
	
	// If damage evolves, find dissipated energy
	double dispEnergy = 0.;
	
#pragma mark ...... Damage Evolution Calculations
	// This section means crack has failed
	if(soft[SOFT_DAMAGE_STATE]>damageState)
	{	// post failure state - increments in cracking strain and stress increment is zero
		dec = deij;
	}
	
	// The damage is evolving after the crack has initiated
	else
	{	// get trial stress and initial D
		double Dn = soft[DAMAGENORMAL];

		// get trial stress increment
		dsig = GetTrialStressIncrement(&deij,Dn,np);
		
		// effective/normal strength properties
		double relStrength = soft[RELATIVE_STRENGTH];
		double relToughness = soft[RELATIVE_TOUGHNESS];
		double sigmaI = relStrength*initiationLaw->sigmaI();
		double scaleI = relToughness*soft[GCSCALING]/sigmaI;
		
		// get shear strength properties (only need for mixed mode methods)
		double sigmaII=0,scaleII=0.;
		if(evolutionStyle==MIXED_MODE_METRIC)
		{	sigmaII = relStrength*initiationLaw->sigmaII();
			scaleII = relToughness*soft[GCSCALING]/sigmaII;
		}

		// trial traction function
		Tensor strial = strn;
		AddTensor(&strial,&dsig);
		double trialPhi = GetTractionFunction(&strial,&enm1,np,Dn,0.,sigmaI,scaleI,sigmaII,scaleII);
		
		// Is it elastic
		if(trialPhi<=0.)
		{	// the update is elastic, we accept above dsig
			// cracking strain in D*de, no energy dissipation
			dec = deij;
			ScaleTensor(&dec,Dn);
		}
		else
		{	// Damage is evolving
			// We need to numerical solve for phi(D+dD)=0
			
			// we know dD=0 has trialPhi>0
			double a = 0,fa = trialPhi;
			
			// Step 1: Small increments in D to find sign change
			//			Failure if equals 1 (or maybe just close to 1)
			double stepSize,b,fb;
			int nsteps;
			if(Dn>0.5)
			{	stepSize = 0.1*(1.-Dn);			// from .1 to .05
				nsteps = 10;
			}
			else if(Dn>0.1)
			{	stepSize = 0.2*(1.-Dn);		   // from .1 to .02
				nsteps = 5;
			}
			else
			{	stepSize = 0.5*(1.-Dn);			// from .05 to 0
				nsteps = 2;
			}
			for(int i=1;i<=nsteps;i++)
			{	b = i*stepSize;
				fb = GetTractionFunction(&strial,&enm1,np,Dn,b,sigmaI,scaleI,sigmaII,scaleII);
				if(fb<0.) break;
				a = b;
				fa = fb;
			}

			// If found bracket, solve for dD, if fb positive or zero, it has failed
			if(fb<0.)
			{	// Step 2: Illinois method to refine the solution between a and b
				//		Sample Code: https://en.wikipedia.org/wiki/Regula_falsi#Example_code
				int side=0,maxSteps=20;
				double c,e = 1e-6;
				for(int n = 0; n < maxSteps; n++)
				{	c = (fb*a - fa*b)/(fb - fa);
					if(fabs(b-a)<e) break;
					double fc = GetTractionFunction(&strial,&enm1,np,Dn,c,sigmaI,scaleI,sigmaII,scaleII);
					if(fc*fb>0.)
					{	// fc and fb have same sign, copy c to b
						b = c; fb = fc;
						if(side==-1) fa /= 2.;
						side = -1;
					}
					else if(fa*fc>0.)
					{	// fc and fa have same sign, copy c to a
						a = c; fa = fc;
						if(side==1) fb /=2.;
						side = 1;
					}
					else
					{	// fc*f_ very small (looks like zero)
						break;
					}
				}
				
				// final answer in c
				double dD = c;
				
				// trial put damage stress increment is total stress increment
				Tensor dstr = GetDamageStressIncrement(&enm1,dD,np);
				AddTensor(&dsig,&dstr);
				
				// Get Ce now for dissipated energy below
				Tensor Ce = GetCe(&enm1,np);
				
				// crack strain increment d(D e) = D*de + enm1*dD
				dec = deij;
				ScaleTensor(&dec,Dn);
				ScaleTensor(&enm1,dD);
				AddTensor(&dec,&enm1);
				
				// Update damage variable
				soft[DAMAGENORMAL] += dD;
				
				// Find delta damage parameters
				double en0 = sigmaI/Ered;
				double deln = softeningModeI->GetDeltaFromDamage(Dn+dD,scaleI,en0,-1.);
				soft[DELTANORMAL] = deln;
				if(evolutionStyle==MIXED_MODE_METRIC)
				{	double gs0 = sigmaII/Gred;
					double delt = softeningModeII->GetDeltaFromDamage(Dn+dD,scaleII,gs0,-1.);
					soft[DELTASHEAR] = delt;
				}
				
				// get dissipated energy (note that enm1 here is enm1*dD)
				dispEnergy = 0.5*DotTensors(&Ce,&enm1);
				soft[DAMAGEGI] += dispEnergy/soft[GCSCALING];
			}
			else
			{	archiver->Decohesion(mtime,mptr,0.,0.,0.,0.,0.,0.,soft[SOFT_DAMAGE_STATE]);

				// Failed, print decohesion
				soft[SOFT_DAMAGE_STATE] += 1.;
				
				// update damage variable
				double dD = 1. - Dn;
				soft[DAMAGENORMAL] = 1.;
				
				// maximum delta values
				soft[DELTANORMAL] = softeningModeI->GetDeltaMax(scaleI);
				if(evolutionStyle==MIXED_MODE_METRIC)
					soft[DELTASHEAR] = softeningModeII->GetDeltaMax(scaleII);

				// Get Ce now for dissipated energy below
				Tensor Ce = GetCe(&enm1,np);
				
				// crack strain increment d(D*e) = D*de + enm1*dD
				dec = deij;
				ScaleTensor(&dec,Dn);
				ScaleTensor(&enm1,dD);
				AddTensor(&dec,&enm1);

				// stress to zero
				dsig = strn;
				ScaleTensor(&dsig,-1.);
				
				// get dissipated energy (note that enm1 here is enm1*dD)
				dispEnergy = 0.5*DotTensors(&Ce,&enm1);
				soft[DAMAGEGI] += dispEnergy/soft[GCSCALING];
			}
		}
	}
	
#pragma mark ... Update Cracking Strain and Stresses
	// Have now found increments in cracking strain and stress
	// update cracking strains (in plastic strain)
	// current was rotated by dR (above) now add increment rotated by Rtot
	
	// rotate previous stress by rotation increment (stress in new current configuration)
	Tensor str = dR.RVoightRT(sp,true,is2D);
	
	// dsig is complete, but is in initial configuration
	// Missing terms is ezz in plane stress
	// get stress incement for other stresses in crack axis system
	if(np==PLANE_STRESS_MPM)
	{	// zz deformation. Here C[4][1] = C[4][2] = -v/(1-v)
		de.zz = p->C[4][1]*(deij.xx-dec.xx) + p->C[4][2]*(deij.yy-dec.yy) + eres;
		mptr->IncrementDeformationGradientZZ(de.zz);
	}
	
	// Update cracking strains
	Tensor *ecrack = mptr->GetAltStrainTensor();
	dec=Rtot.RVoightRT(&dec,false,is2D);
	AddTensor(ecrack, &dec);
	
	// update stress and strain (in which cracking strain rotated to current configuration)
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
	double delV = de.xx-dec.xx+de.yy-dec.yy+de.zz-dec.zz;
	double dTq0 = -gamma0*mptr->pPreviousTemperature*delV;
	
	// track heat energy with total dissipated energy
	IncrementHeatEnergy(mptr,dTq0,dispEnergy);
}

// Get ncrement in stress for constant D for effective strain increment
// Plane stress include deij.zz and can therefore use 3D equation
Tensor IsoDamageMech::GetTrialStressIncrement(Tensor *deij,double Dn,int np) const
{
	// Normal damage term
	double omDnormal = Kred*(1.-Dn);
	double omDshear = Gred*(1.-Dn);
	double dTre = deij->xx+deij->yy+deij->zz;
	
	Tensor dstr;
	ZeroTensor(&dstr);
	
	// normal stress increments in initial analysis axes
	double nTerm = (omDnormal - 2.*omDshear/3.)*dTre;
	dstr.xx = nTerm + 2.*omDshear*deij->xx;
	dstr.yy = nTerm + 2.*omDshear*deij->yy;
	if(np!=PLANE_STRESS_MPM) dstr.zz = nTerm + 2.*omDshear*deij->zz;
	
	// shear stress increments
	dstr.xy = omDshear*deij->xy;
	if(np==THREED_MPM)
	{	// elastic stress update
		dstr.yz = omDshear*deij->yz;
		dstr.xz = omDshear*deij->xz;
	}
	
	// return increment
	return dstr;
}

// Get increment in stress for dD in current state
Tensor IsoDamageMech::GetDamageStressIncrement(Tensor *enm1,double dD,int np) const
{
	Tensor dstr;
	ZeroTensor(&dstr);

	// damage term
	if(dD>0.)
	{	// normal stress with damage increment
		double Tren = enm1->xx+enm1->yy+enm1->zz;
		double nTerm = (Kred - 2.*Gred/3.)*Tren;
		dstr.xx -= (nTerm + 2.*Gred*enm1->xx)*dD;
		dstr.yy -= (nTerm + 2.*Gred*enm1->yy)*dD;
		if(np!=PLANE_STRESS_MPM) dstr.zz -= (nTerm + 2.*Gred*enm1->zz)*dD;
		
		// shear with damage increment
		dstr.xy -= Gred*enm1->xy*dD;
		if(np==THREED_MPM)
		{	// elastic stress update
			dstr.yz -= Gred*enm1->yz*dD;
			dstr.xz -= Gred*enm1->xz*dD;
		}
	}
	
	// return increment
	return dstr;
}

// Get Stiffness times strain for use in energy dissipation
// By leaving reduced, if will result in energy per unit mass
Tensor IsoDamageMech::GetCe(Tensor *enm1,int np) const
{
	Tensor Ce;
	ZeroTensor(&Ce);
	
	// normal stresses
	double Tren = enm1->xx+enm1->yy+enm1->zz;
	Ce.xx = Kred*Tren + 2*Gred*enm1->xx;
	Ce.yy =	Kred*Tren + 2*Gred*enm1->xx;
	if(np!=PLANE_STRESS_MPM) Ce.zz = Kred*Tren + 2.*Gred*enm1->zz;
	
	Ce.xy = Gred*enm1->xy;
	if(np==THREED_MPM)
	{	// elastic stress update
		Ce.yz = Gred*enm1->yz;
		Ce.xz = Gred*enm1->xz;
	}
	
	// return stress
	return Ce;
}


// check stress state in initial axis system to see if failed
// If failed, set normal to principle normal direction, but for
//     3D, find Euler angles in normal
// If not failed, return 0 (NO_FAILURE)
int IsoDamageMech::ShouldInitiateFailure(Tensor *str,double relStrength,Vector *normal,int np) const
{
	// initialize
	int failureMode = NO_FAILURE;
	
	// did it reach initiation value yet?
	double phiCurrent = GetTractionFunction(str,NULL,np,0.,0.,relStrength*initiationLaw->sigmaI(),0.,0.,0.);
	if(phiCurrent<0.) return failureMode;
	
    // failure - take normal in the maximum principal stress direction
    // stress and tensile principal stress metric only
    // only calculated one when initiates, so not computation cost to eigenanalysis
    if(np==THREED_MPM)
    {   // get eigenvalues and vectors
        Matrix3 strmat = TensorToMatrix(str,true);
        Vector s = strmat.Eigenvalues();
        Matrix3 R = strmat.Eigenvectors(s);
        
        // Note the matrix class might reorder eigenvalues
        // Six options - convert to x>y>z
        if(s.y<s.z)
        {   // swap y and z
            R.SwapColumns(1,2);
            double temp = s.y;
            s.y = s.z;
            s.z = temp;
        }
        // Now have x>y>z, y>x>z, or y>z>x
        if(s.x<s.z)
        {   // rotate y->x,z->y,x->z
            R.SwapColumns(1,2);
            R.SwapColumns(0,2);
            Vector t = s;
            s.x = t.y;
            s.y = t.z;
            s.z = t.x;
        }
        else if(s.x<s.y)
        {   // have y>x>z, swap x and y
            R.SwapColumns(0,1);
            double temp = s.y;
            s.y = s.x;
            s.x = temp;
        }

        // flip z if determinant is < 0
        if(R.determinant()<0.)
        {   R(0,2) = -R(0,2);
            R(1,2) = -R(1,2);
            R(2,2) = -R(2,2);
        }
        
        // decode to Euler angles
        // see https://en.wikipedia.org/wiki/Euler_angles - but it has the equation wrong
        // There are two solutions (where second angle is one way or the other) for ZYZ Euler angles
        //   and either can be used for softening material to reconstruct the ZYZ R matrix
        // normal to crack will be first column on reconstructed R matrix
        if(R(2,2)<-1.)
            normal->y = acos(-1.);
        else if(R(2,2)>1.)
            normal->y = acos(1.);
        else
            normal->y = acos(R(2,2));                // beta  = arccos(Z3)
        if(R(1,2)==0. && R(0,2)==0.)
            normal->x = 0.;
        else
            normal->x = atan2(R(1,2),R(0,2));        // alpha = atan2(Z2,Z1)
        if(R(2,1)==0. && R(2,0)==0.)
            normal->z = 0.;
        else
            normal->z = atan2(R(2,1),-R(2,0));        // gamma = atan2(Y3,-X3)
    }
    else
    {   // 2D initiation
        
        // get principal stress and swap if needed to get s.x>s.y
        Vector s = TensorEigenvalues(str,true,true);
        
        // normal in maximum stress direction
        double nx=1.,ny=0.;
        if(s.x<s.y)
        {   nx = 0.;
            ny = 1.;
        }

        // atan2() is ccw angle from max principle direction to crack norm
        // second term is ccw angle from global x axis to max principle
        // sum is ccw angle from global x axis to crack normal
        double theta = atan2(ny,nx) + 0.5*atan2(2.*str->xy,str->xx-str->yy);
        normal->x = cos(theta);
        normal->y = sin(theta);
        normal->z = 0.;
    }
	
	// make as failed
	failureMode = 100;
	return failureMode;
}

// Get effective failure surface
//    Stress Energy: sqrt(E * sigma(D+dD).S.sigma(d+dD)) - sigma_c(D+dD)
//    Principal Stress: tensile stress only
//    Mixed-mode metric:
// Note that enm1 not needed if dD=0 (can be NULL)
double IsoDamageMech::GetTractionFunction(Tensor *strial,Tensor *enm1,int np,double Dn,double dD,
										  double sigmaI,double scaleI,double sigmaII,double scaleII) const
{
	// get strength at D+dD
	// (scale==0 only called for initiation and sigmaI has the initial strength)
	double fnval=sigmaI,ftval=0.;
	if(scaleI>0.)
	{	double en0 = sigmaI/Ered;
		double deln = softeningModeI->GetDeltaFromDamage(Dn+dD,scaleI,en0,-1.);
		fnval = sigmaI*softeningModeI->GetFFxn(deln,scaleI);
		if(evolutionStyle==MIXED_MODE_METRIC)
		{	double gs0 = sigmaII/Gred;
			double delt = softeningModeII->GetDeltaFromDamage(Dn+dD,scaleII,gs0,-1.);
			ftval = sigmaII*softeningModeII->GetFFxn(delt,scaleII);
		}
	}

	// if damage, get new stress
	Tensor str = *strial;
	if(dD>0.)
	{	Tensor dstr = GetDamageStressIncrement(enm1,dD,np);
		AddTensor(&str,&dstr);
	}
	
	if(evolutionStyle==MIXED_MODE_METRIC)
	{	// get principle stress
		Matrix3 strmat = np==THREED_MPM ? TensorToMatrix(&str,true) : TensorToMatrix2D(&str,true);
		Vector s = strmat.Eigenvalues();
		
		// max and min principle stress
		double s1 = fmax(s.x,fmax(s.y,s.z));
		double s3 = fmin(s.x,fmin(s.y,s.z));
		
		double dist;
		
		if(s3<=-s1)
		{	// shear failure (watch for s1=s3,which is along 45 and elastic regime)
			double shearMax = 0.5*(s1-s3);
			dist = shearMax>0.1*ftval ? 1. - ftval/shearMax : -9. ;
		}
		else
		{	double r2 = ftval*ftval/(fnval*fnval);
			if(s3<s1*(1.-2.*r2))
			{	// mixed falure
				double arg = (1.-r2)/((s1-s3)*(s1-s3)+4.*s1*s3*r2);
				dist = (1. - 2.*ftval*sqrt(arg));
			}
			else
			{	// tensile failure
				dist = 1.-fnval/s1;
			}
		}
		
		/* Old code
		if(s1>fnval)
		{	if(s3>fnval)
			{	double dx = s1-fnval;
				double dy = s3-fnval;
				dist = sqrt(dx*dx+dy*dy);
			}
			else if(s3>fnval-2*ftval)
				dist = s1-fnval;
			else if(s3>2.*(fnval-ftval)-s1)
			{	double dx = s1-fnval;
				double dy = s3-fnval+2*ftval;
				dist = sqrt(dx*dx+dy*dy);
			}
			else
				dist = sqrt(2.)*(0.5*(s1-s3)-ftval);
		}
		else if(s1-s3>2.*ftval)
			dist = sqrt(2.)*(0.5*(s1-s3)-ftval);
		else
		{	// negative number inside failure surface
			dist = fmin(sqrt(2.)*(0.5*(s1-s3)-ftval),s1-fnval);
		}
		*/
		
		return dist;
	}
	else if(evolutionStyle==STRESS_ENERGY_METRIC)
	{   // strains here are E S sigma
        double exx = str.xx - nu*(str.yy+str.zz);
		double eyy = str.yy - nu*(str.xx+str.zz);
		double ezz = str.zz - nu*(str.xx+str.yy);
		double arg = str.xx*exx + str.yy*eyy + str.zz*ezz;
		arg += str.xy*str.xy*2.*(1+nu);
		if(np==THREED_MPM)
			arg += (str.xz*str.xz + str.yz*str.yz)*2.*(1+nu);
		return sqrt(arg) - fnval;
	}
	else
    {    // get principle stress
        Matrix3 strmat = np==THREED_MPM ? TensorToMatrix(&str,true) : TensorToMatrix2D(&str,true);
        Vector s = strmat.Eigenvalues();
        
        // strains here are E S sigma and only normal ones in prinple stress space
        double exx = s.x - nu*(s.y+s.z);
        double eyy = s.y - nu*(s.x+s.z);
        double ezz = s.z - nu*(s.x+s.y);
        
        // ignore compressive principle stress
        double arg = 0.;
        if(s.x>0) arg += s.x*exx;
        if(s.y>0) arg += s.y*eyy;
        if(s.z>0) arg += s.z*ezz;
        
		// reture signed metric
        return sqrt(arg) - fnval;
	}
}

#pragma mark IsoDamageMech::Accessors

// return material type
const char *IsoDamageMech::MaterialType(void) const { return "Isotropic Damage Mechanics"; }

