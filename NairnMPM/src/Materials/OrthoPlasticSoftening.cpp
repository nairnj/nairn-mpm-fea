/********************************************************************************
    OrthoPlasticSoftening.hpp
    nairn-mpm-fea
    
    Created by John Nairn, Oct 3, 2020.
    Copyright (c) 2020 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/OrthoPlasticSoftening.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Materials/OrthoFailureSurface.hpp"

int hstepsOPSMax=0;

#pragma mark OrthoPlasticSoftening::Constructors and Destructors

// Thie method should call parent class and then fill in any new initiatizations
// Constructor
OrthoPlasticSoftening::OrthoPlasticSoftening(char *matName,int matID) : OrthoSoftening(matName,matID)
{
    // negative yield stress implies no yielding in that direction
    syxx=-1.;
    syyy=-1.;
    syzz=-1.;
    tyxy=-1.;
    tyxz=-1.;
    tyyz=-1.;
    
    // default value of hardening (elastic plastic)
    Khard=0.;
    nhard=1.;
    exphard=0.;
    alphaMax=-1.;
    hardStyle = AP_LINEAR;
	hMAS.style = SQRT_TERMS;
}

#pragma mark OrthoPlasticSoftening::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *OrthoPlasticSoftening::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"yldxx")==0)
    {   input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&syxx,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"yldyy")==0)
    {   input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&syyy,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"yldzz")==0)
    {   input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&syzz,gScaling,1.e6);
    }

    else if(strcmp(xName,"yldxy")==0)
    {   input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&tyxy,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"yldxz")==0)
    {   input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&tyxz,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"yldyz")==0)
    {   input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&tyyz,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"Khard")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Khard);
    }

    else if(strcmp(xName,"nhard")==0)
    {   input=DOUBLE_NUM;
        hardStyle = AP_UNKNOWN;
        return((char *)&nhard);
    }
    
    // exphard - hardening exponential rate
    else if(strcmp(xName,"exphard")==0)
    {   input=DOUBLE_NUM;
        hardStyle = AP_EXPONENTIAL;
        return((char *)&exphard);
    }

    // exphard - hardening exponential rate
    else if(strcmp(xName,"alphaMax")==0)
    {   input=DOUBLE_NUM;
        return((char *)&alphaMax);
    }

	else if(strcmp(xName,"HillStyle")==0)
	{	input=INT_NUM;
		return (char *)(&hMAS.style);
	}
	
    return OrthoSoftening::InputMaterialProperty(xName,input,gScaling);
}

// Verify input properties do calculations; if problem return string with an error message
// NOTE: This code duplicated in AniosPlasticity. Keep them in sync
const char *OrthoPlasticSoftening::VerifyAndLoadProperties(int np)
{
#ifdef POROELASTICITY
    if(DiffusionTask::HasPoroelasticity())
        return "Anisotropic plasticity materials cannot yet be used when poroelasticity is activated";
#endif
    
    if(hardStyle==AP_UNKNOWN)
    {   hardStyle = nhard<=0. ? AP_NONLINEAR1 : AP_NONLINEAR2;
        nhard = fabs(nhard);
        if(DbleEqual(nhard,1.)) hardStyle=AP_LINEAR;
    }
    else if(hardStyle==AP_EXPONENTIAL)
    {   if(exphard==0.)
            return "The exphard parameter cannot be zero.";
        Kexp = Khard/exphard;
        if(alphaMax>0)
        {   slopeMin = Khard*exp(-exphard*alphaMax);
            hmax = 1+Kexp-slopeMin*(alphaMax+1./exphard);
        }
    }
    
    // requires large rotation mode
    if(useLargeRotation==0) useLargeRotation = 1;
    
    if(swapz==1)
    {   // swap x and z, xx->zz, xy->yz, xz same, yz->xy
        SwapProperties(syxx,syzz,tyyz,tyxy);
    }
    else if(swapz>1)
    {   // swap y and z, yy->zz, xy->xz, xz->xy, yz same
        SwapProperties(syyy,syzz,tyxz,tyxy);
    }

    // check at least some yielding
    if(syzz<0. && syxx<0. && syyy<0. && tyxy<0. && tyxz<0. && tyyz<0.)
        return "No yield stresses were defined";
    
    // check non zeros
    if(syzz==0. || syxx==0. || syyy==0. || tyxy==0. || tyxz==0. || tyyz==0.)
        return "No yield stresses can be zero";

    // check A is positive semi definite (watching for round off error is all the same)
    if(syzz>0. && syxx>0. && syyy>0. && syxx==syyy && syxx==syzz)
    {   // OK if all the same
    }
    else
    {   double rsxx=0.,rsyy=0.,rszz=0.;
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
    {   hMAS.syxx2=rho/syxx;
        hMAS.syxx2*=hMAS.syxx2;
    }
    else
        hMAS.syxx2=0.;        // 1/inf^2
    if(syyy>0.)
    {   hMAS.syyy2=rho/syyy;
        hMAS.syyy2*=hMAS.syyy2;
    }
    else
        hMAS.syyy2=0.;        // 1/inf^2
    if(syzz>0.)
    {   hMAS.syzz2=rho/syzz;
        hMAS.syzz2*=hMAS.syzz2;
    }
    else
        hMAS.syzz2=0.;        // 1/inf^2
    
    // reciprocals of reduced shear yield stresses
    if(tyxy>0.)
    {   hMAS.N=rho/tyxy;
        hMAS.N*=hMAS.N;
    }
    else
        hMAS.N=0.;        // 1/inf^2
    if(tyxz>0.)
    {   hMAS.M=rho/tyxz;
        hMAS.M*=hMAS.M;
    }
    else
        hMAS.M=0.;        // 1/inf^2
    if(tyyz>0.)
    {   hMAS.L=rho/tyyz;
        hMAS.L*=hMAS.L;
    }
    else
        hMAS.L=0.;        // 1/inf^2
    
    // combination terms
    hMAS.F = 0.5*(hMAS.syyy2 + hMAS.syzz2 - hMAS.syxx2);
    hMAS.G = 0.5*(hMAS.syzz2 + hMAS.syxx2 - hMAS.syyy2);
    hMAS.H = 0.5*(hMAS.syxx2 + hMAS.syyy2 - hMAS.syzz2);
    
    // reference yield stress (it is acutally sqrt(2/3)/sigma(Y,ref)
    double sumNormal = hMAS.F+hMAS.G+hMAS.H;
    if(sumNormal<=0.) return "(F+G+H) for Hill plastic potential must be positive";
    sqrt23OversigmaYref = sqrt(2.*sumNormal/3.);
    
    // for convergence problems
    AnisoPlasticity::warnNonconvergence = warnings.CreateWarning("Anisotropic plastic algorithm failed to converge",-1,3);

	// call super class
	const char *errMsg = OrthoSoftening::VerifyAndLoadProperties(np);
	if(errMsg!=NULL) return errMsg;
	
	// set more hill properties
	AnisoPlasticity::FillHillStyleProperties(np,hMAS,pr);

	return NULL;
}

// print mechanical properties to the results
void OrthoPlasticSoftening::PrintMechanicalProperties(void) const
{
	// call superclass
    OrthoSoftening::PrintMechanicalProperties();
	
	// add new properties here
    AnisoPlasticity::PrintAPYieldProperties(syxx,syyy,syzz,tyyz,tyxz,tyxy);
    cout << "Relative yield = ";
    switch(hardStyle)
    {   case AP_LINEAR:
            cout << "1 + " << Khard << "*alpha";
            break;
        case AP_NONLINEAR1:
            cout << "(1 + " << Khard << "*alpha)^" << nhard;
            break;
        case AP_EXPONENTIAL:
            if(alphaMax>0)
            {   cout << "1 + (" << Kexp << ")*(1-exp[-" << exphard << "*alpha]) for alpha<" << alphaMax << endl;
                cout << "                 " << hmax << " + " << slopeMin << "*alpha for alpha>" << alphaMax;
            }
            else if(exphard>0.)
                cout << "(1 + (" << Kexp << ")*(1-exp[-" << exphard << "*alpha])";
            else
            {   // Hill material does not allow negative exphard (I am not sure what it is for?)
                cout << "(1 + (" << Kexp << ")*(1-exp[" << fabs(exphard) << "*alpha])";
            }
            break;
        default:
            cout << "1 + " << Khard << "*alpha^" << nhard;
            break;
    }
    cout << endl;
	
	if(hMAS.style==SQRT_TERMS)
		cout << "Hill style: sqrt(sAs)-relY" << endl;
	else
		cout << "Hill style: sAs - relY^2" << endl;
}


#pragma mark OrthoPlasticSoftening::History Data Methods

// History variables are stores as (1 for Hill strain),(Standard Softening),(Cracking Strains)
// Cracking strains history is needed because plastic strain used by plasticity
char *OrthoPlasticSoftening::InitHistoryData(char *pchr,MPMBase *mptr)
{   // soften variables (update if parent class updates)
    double *p = (double *)OrthoSoftening::InitHistoryData(pchr,mptr);
    
    // If has softening history variables, shift other variables by one
    // Only non-zero ones in parent class needing shift are relative strength and toughness (11 and 12)
    p[1+RELATIVE_TOUGHNESS] = p[RELATIVE_TOUGHNESS];
    p[1+RELATIVE_STRENGTH] = p[RELATIVE_STRENGTH];
    p[RELATIVE_STRENGTH] = 0.;
    
    // add damage state (0)
    p[1+SOFT_DAMAGE_STATE] = p[SOFT_DAMAGE_STATE];
    p[SOFT_DAMAGE_STATE] = 0.;

    return (char *)p;
}

// reset history data
void OrthoPlasticSoftening::ResetHistoryData(char *pchr,MPMBase *mptr)
{	double *p = (double *)pchr;
	double relStrength = p[1+RELATIVE_STRENGTH];
	double relToughness = p[1+RELATIVE_TOUGHNESS];
	ZeroDoubles(pchr,NumberOfHistoryDoubles());
	p[1+RELATIVE_STRENGTH] = relStrength;
	p[1+RELATIVE_TOUGHNESS] = relToughness;
	p[1+SOFT_DAMAGE_STATE] = 0.1;
}

// Number of history variables
// Hill plastic strain, damage variabls, history storage for cracking strains
int OrthoPlasticSoftening::NumberOfHistoryDoubles(void) const
{   return 1 + SOFT_NUMBER_HISTORY + NUMBER_CRACKING_STRAINS;
}

#pragma mark OrthoPlasticSoftening::Methods

// Constitutive Law and it is always in large rotation mode
void OrthoPlasticSoftening::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,
                                             ResidualStrains *res,int historyOffset,Tensor *gStress) const
{
    // set 2D flag
    bool is2D = np == THREED_MPM ? false : true;

#pragma mark ... Code to get trial stress and check for yielding
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
    
    // cast pointer to material-specific data
    AnisoPlasticProperties *ap = (AnisoPlasticProperties *)properties;
    ElasticProperties *ep = GetElasticPropertiesPointer(properties);
	
	// to fill with properties in the crack axis system (or MAS if not damaged)
    CrackAxisProperties d;
    HillProperties h;
    int DForm = -1;
    double dispEnergyPlastic = 0.;
    
    // residual strains (thermal and moisture) in material axes
    double exxr,eyyr,ezzr;
    Tensor er = GetAnisoResStrains(exxr,eyyr,ezzr,ep,res,np);
    
#pragma mark ... Plasticity phase input calculations
    // get history for softening properties (plastic(1),damage,cracking strains) storage mode
    // Returned value for this material is plastic strain, add 1 to get to softening history
    double *soft = GetSoftHistoryPtr(mptr);
    double *ipsoft = soft + SOFT_NUMBER_HISTORY;            // pointer to first cracking strain
    Tensor *sp = mptr->GetStressTensor();

    // Get de = (I-D)deT = deT - dec or strain increment in bulk assuming no damage evolution
    // Here dec is cracking strain stored in history variables starting at ipsoft
    Tensor dec;
    if(soft[SOFT_DAMAGE_STATE]<predamageState)
    {    // will remain in initial configuration
        dec = MakeTensor(0.,0.,0.,0.,0.,0.);
        
        // Rnm1 and Rtot remain rotation from current state to material axes
        // Load material properties into CrackAxisProperties
        LoadCrackAxisProperties(np,&d,DX_DAMAGE,ep);
        LoadCrackHillProperties(np,&h,&d,DX_DAMAGE);
    }
    else
    {   // When damaged, rotate strain increments from material axes (found above) to crack axis system
        DForm = GetDForm(soft[SOFT_DAMAGE_STATE]);
        Matrix3 RToCrack;
        GetRToCrack(&RToCrack, soft, is2D, DForm);
        de = RToCrack.RTVoightR(&de, false, is2D);
        er = RToCrack.RTVoightR(&er, false, is2D);
        
        // Get material properties in crack axis system
        LoadCrackAxisProperties(np,&d,DForm,ep);
        LoadCrackHillProperties(np,&h,&d,DForm);
        
        // get cracking strain increments for an elastic update (de(net) = de-er)
        // plane stress is different here, but not yet supported elsewhere
        double den = de.xx-er.xx +  d.C12C11*(de.yy-er.yy);
        if(np==THREED_MPM) den += d.C13C11*(de.zz-er.zz);
        // ... and make sure increment does not create negative normal strain
        // ... if ecxx + Dn*den < 0, Dn*den = -ecxx
        double decxx = fmax(soft[DAMAGENORMAL]*den,-ipsoft[ECXX_DAMAGE]);
        dec = is2D ? MakeTensor2D(decxx,0.,0.,soft[DAMAGESHEAR]*de.xy) :
                        MakeTensor(decxx,0.,0.,0.,soft[DAMAGESHEAR2]*de.xz,soft[DAMAGESHEAR]*de.xy);
        
        // get rotation from current state to crack axis system
        Rnm1 *= RToCrack;
        Rtot *= RToCrack;
        
        // get de = detot-dec (this includes residual strain still)
        SubTensor(&de,&dec);
    }
    
    // for plasticity phase de(trial)eff = de - er in material/CAS (for undamaged/damaged)
    // Material and Hill properties in material/CAS are in d and h

    // Step 1: Get trial stress (in str) in material/CAS assuming increment is elastic and no cracking strain increment
    
    // effective strains
    double dvxxeff = de.xx-er.xx;
    double dvyyeff = de.yy-er.yy;
    double dvzzeff = de.zz-er.zz;       // not used in plain strain
    // rotate current stress to material axes and then increment
    Tensor stnm1Mat = Rnm1.RTVoightR(sp,true,is2D);
    Tensor strPlastic,str;
    if(np==THREED_MPM)
    {   // rotate current stress (stnm1) to material axes (stnm1Mat) and increment with effective strain increment
        str.xx = stnm1Mat.xx + d.C11*dvxxeff + d.C12*dvyyeff + d.C13*dvzzeff;
        str.yy = stnm1Mat.yy + d.C12*dvxxeff + d.C22*dvyyeff + d.C23*dvzzeff;
        str.zz = stnm1Mat.zz + d.C13*dvxxeff + d.C23*dvyyeff + d.C33*dvzzeff;
        str.yz = stnm1Mat.yz + d.C44*de.yz;
        str.xz = stnm1Mat.xz + d.C55*de.xz;
        str.xy = stnm1Mat.xy + d.C66*de.xy;
    }
    else
    {   // rotate current stress (stnm1) to material/CAS axes (stnm1Mat) and increment with effective strain increment
        str.xx = stnm1Mat.xx + d.C11*dvxxeff + d.C12*dvyyeff;
        str.yy = stnm1Mat.yy + d.C12*dvxxeff + d.C22*dvyyeff;
        str.xy = stnm1Mat.xy + d.C66*de.xy;
        
        // sigma(zz)
        if(np==PLANE_STRAIN_MPM)
            str.zz = stnm1Mat.zz + d.C13*(dvxxeff+d.vzxc*ezzr) + d.C23*(dvyyeff+d.vzyc*ezzr) - d.C33*ezzr;
        else
            str.zz = stnm1Mat.zz + d.C13*dvxxeff + d.C23*dvyyeff + d.C33*dvzzeff;
    }
    
    // Rotate plastic strain in state n-1 to configuration n
    Tensor *eplast=mptr->GetAltStrainTensor();
    *eplast = dR.RVoightRT(eplast,false,is2D);
        
    // Step 3: Calculate plastic potential f
    ap->aint = mptr->GetHistoryDble(0,0);
	// NEW_HILL two surface options
	double halfsPs = AnisoPlasticity::GetHillMagnitude(str,&h,np);
	double yield = GetYield(ap);
	double ftrial;
	if(h.style==SQRT_TERMS)
	{	ap->sAQsmag = AnisoPlasticity::GetHillMagnitude(str,&h,np);
		ftrial = ap->sAQsmag - yield;
	}
	else
		ftrial = halfsPs - yield*yield;
    
    // Step 4: Done if elastic
    Tensor dep;
    if(ftrial<=0.)
    {   // No plasticity effect on final stress
        // add back cracking strain to get de(elastic) = de(total) = de(trial)+dec (inclusive of residual strains)
        AddTensor(&de,&dec);
        dep = MakeTensor(0.,0.,0.,0.,0.,0.);
    }
    else
    {   // Evaluate plastic strains and their change in trial stress
        
        // Solve for plastic strain increment
        strPlastic = str;
        double alpha0 = ap->aint;
        dep = SolveForPlasticIncrement(mptr,np,ftrial,str,ap,&d,&h);
        
        // Get decrease in stress caused by plastic strain
        // str = str0-strPlastic or strPlastic = str0(current strPlastic) - str
        SubTensor(&strPlastic,&str);
        
        // get dissipated energy
        if(np==THREED_MPM)
        {   // Plastic energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
            dispEnergyPlastic = str.xx*dep.xx + str.yy*dep.yy + str.zz*dep.zz
                            + str.xy*dep.xy + str.xz*dep.xz + str.yz*dep.yz;
        }
        else
        {   // Plastic energy increment per unit mass (dU/(rho0 V0)) (Legacy nJ/g)
            dispEnergyPlastic = str.xx*dep.xx + str.yy*dep.yy + str.zz*dep.zz + str.xy*dep.xy;
        }
        
        // add dissipated energy to plastic energy to the particle after subtracting q*dalpha
        dispEnergyPlastic -= (GetYield(ap)-1.)*(ap->aint-alpha0)/sqrt23OversigmaYref;
        mptr->AddPlastEnergy(dispEnergyPlastic);
        
        // Total work energy is sigma.de and softening phase will add only sigma.(de-dep)
        //     To compenstate, we add sigma.dep here (note that using trial elastic-plastic stress
        //     rather than final elastic-plastic-softening stress)
        // FINISH UP: I think correct to add non-zero dep.zz here even in plane strain
        // heat energy incremented here prior to damage or in damage evolution adter damage
        double partialWorkEnergy = str.xx*dep.xx + str.yy*dep.yy + str.xy*dep.xy + str.zz*dep.zz;
        if(!is2D)
            partialWorkEnergy += str.yz*dep.yz + str.xz*dep.xz;
        mptr->AddWorkEnergy(partialWorkEnergy);
        
        // get elastic strain increment (de(elastic) = de(total)-dep = de(trial)+dec-dep) in crack axis system
        AddTensor(&de,&dec);
        SubTensor(&de,&dep);
        
        // rotate plastic strain to current state and add to particle
        Tensor depGlobal = Rtot.RVoightRT(&dep,false,is2D);
        AddTensor(eplast,&depGlobal);
        
        // Update Internal Variable
        mptr->SetHistoryDble(0,ap->aint,0);
    }

#pragma mark ... Code for Previously Undamaged Material
    // Task 5: Begin damage mechanics phase
    //-------------------------------------
    // These values are in material (undamaged) or crack axis (damaged) system
    // 1. The trial stress after plasticity in Tensor str
    // 2. The decrease in trial stresses from elastic trial caused by plasticity in strPlastic
    // 3. Total elastic strain increment Tensor de(elastic) = de(total)-dep

    // Before cracking, do normal anisotropic update. If not cracked
    // then done, otherwise initiate the damage
    if(soft[SOFT_DAMAGE_STATE]<predamageState)
    {   // check if has failed
        Vector norm;
        double tempRelStrength = soft[RELATIVE_STRENGTH];
        int failureMode = initiationLaw->ShouldInitiateFailure(&str,&norm,np,tempRelStrength,NULL);
        if(failureMode == NO_FAILURE)
        {    // Not failed yet, so finish elastic update
            
            // rotate stress to global coordinates new sigma = Rtot str RtotT
            *sp = Rtot.RVoightRT(&str, true, is2D);
            
            // stresses are in global coordinates so need to rotate strain and residual
            // strain to get work energy increment per unit mass (dU/(rho0 V0))
            de = Rtot.RVoightRT(&de, false, is2D);
            er = Rtot.RVoightRT(&er, false, is2D);
            
            // updates
            if(np==THREED_MPM)
            {    // energy
                mptr->AddWorkEnergyAndResidualEnergy(DotTensors(sp, &de), DotTensors(sp, &er));
            }
            else
            {    // work and residual strain energy increments and sigma or F in z direction
                double workEnergy = sp->xx*de.xx + sp->yy*de.yy + sp->xy*de.xy;
                double resEnergy = sp->xx*er.xx + sp->yy*er.yy + sp->xy*er.xy;
                if(np==PLANE_STRAIN_MPM)
                {   // extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
                    workEnergy += sp->zz*de.zz;
                    resEnergy += sp->zz*er.zz;
                }
                else
                {   // extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
                    workEnergy += sp->zz*de.zz;
                    resEnergy += sp->zz*er.zz;
                }
                mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
            }
            
            // track heat energy (should get dTq0), dissipated from plasticity above
            IncrementHeatEnergy(mptr,0.,dispEnergyPlastic);
            
            // this elastic update is done
            return;
        }

#pragma mark ...... Undamaged Material Just Initiated Damage
        // initiate failure
        // 3D, norm = ZYZ rotation from material axes to crack normal. Currently either (0,pi/2,0)
        //        to indicate normal in axial direction (switches z and x) or (theta,0,0) to keep axial direciton
        //        in z direction, but around that axis
        // 2D, norm = (cos(theta),sin(theta))
        soft[NORMALDIR1] = norm.x;                                        // cos(theta) (2D) or Euler alpha (3D) to normal
        soft[NORMALDIR2] = norm.y;                                        // sin(theta) (2D) or Euler beta (3D)  to normal
        soft[NORMALDIR3] = norm.z;                                        // unused (2D) or Euler gamma (3D) to normal
        
        // initiate failure in material axis system with predamageState=0.5
        // Note that DForm=Di_DAMAGE means material i direction in crack axis x direction
        DForm = DecodeDamageInitiation(np,&norm,failureMode,soft);
        LoadCrackAxisProperties(np,&d,DForm,ep);
        
        // get intersection area after rotating normal into global axes
        // Crack rotation to matrix
        Matrix3 RToCrack;
        GetRToCrack(&RToCrack,soft,is2D,DForm);
        
        // initial rotation to material axes and then to the CAS
        R0 *= RToCrack;
        Vector globalNorm = MakeVector(R0(0,0),R0(1,0),R0(2,0));
        
        // get intersection with this normal
        soft[GCSCALING] = GetAcOverVp(np,mptr,&globalNorm)/rho;            // divide by rho because divided by specific stress
        
        // Rotate elastic strain increments from material to CAS
        de = RToCrack.RTVoightR(&de,false,is2D);
        er = RToCrack.RTVoightR(&er,false,is2D);

        // Rotate plastic stress change from material to CAS (should verify)
        if(ftrial>0.)
            strPlastic = RToCrack.RTVoightR(&strPlastic,true,is2D);
        
        // include rotation in total rotation matrices
        Rnm1 *= RToCrack;
        Rtot *= RToCrack;
    }
    
#pragma mark ... Code for Damaged Material
    // A crack is present - get tensor form
    if(DForm<0)
    {   DForm = GetDForm(soft[SOFT_DAMAGE_STATE]);
        LoadCrackAxisProperties(np,&d,DForm,ep);
    }

    // de is elastic strain increment in CAS as de = de(total)-dep
    // er is residual strain incement in CAS
    
    // subtract plastic stress change from n-1 state currently on the particle
    // DamageEvolution() will finish rotation to current state
    if(ftrial>0.)
    {   strPlastic = Rnm1.RVoightRT(&strPlastic,true,is2D);
        SubTensor(sp,&strPlastic);
    }
    
    // Initial elastic stress in CAS after subtracting plastic stress term
    str = Rnm1.RTVoightR(sp,true,is2D);

    // Get cracking strain in CAS
    Tensor ecrack = MakeTensor(ipsoft[ECXX_DAMAGE],0.,0.,0.,ipsoft[GCXZ_DAMAGE],ipsoft[GCXY_DAMAGE]);

    // finish up in separate code
    DamageEvolution(mptr,np,soft,de,str,er,d,dR,Rtot,ecrack,DForm,dispEnergyPlastic,delTime);
}

// This material stores cracking strain in history variable
// The cracking strain is stored in the crack axis system
void OrthoPlasticSoftening::UpdateCrackingStrain(int np,Tensor *ecrack,double decxx,double dgcxy,double dgcxz,Matrix3 Rtot,double *soft) const
{
    double *ipsoft = soft + SOFT_NUMBER_HISTORY;            // pointer to first cracking strain

    // update cracking strain
    ipsoft[ECXX_DAMAGE] += decxx;
    ipsoft[GCXY_DAMAGE] += dgcxy;
    ipsoft[GCXZ_DAMAGE] += dgcxz;
}

// Fill d with properties in the crack axis system depending on the type of damage
// tensor as specified by DForm
void OrthoPlasticSoftening::LoadCrackHillProperties(int np,HillProperties *h,CrackAxisProperties *d,int DForm) const
{
    if(DForm==DZ_DAMAGE)
    {   // crack axis normal in material z direction, swap x and z
        // 1 = Z, 2 = Y, 3 = X (only occurs in 3D)
        h->F = hMAS.H;
        h->G = hMAS.G;
        h->H = hMAS.F;
        h->L = hMAS.N;
        h->M = hMAS.M;
        h->N = hMAS.L;
        h->syxx2 = hMAS.syzz2;
        h->syyy2 = hMAS.syyy2;
        h->syzz2 = hMAS.syxx2;
    }
    else if(DForm==DX_DAMAGE)
    {   // crack axis normal in material x direction, keep original form
        // 1 = X, 2 = Y, 3 = Z
		*h = hMAS;
    }
    else
    {   // crack axis normal in material y direction, swap x and y
        // 1 = Y, 2 = X, 3 = Z
        h->F = hMAS.G;
        h->G = hMAS.F;
        h->H = hMAS.H;
        h->L = hMAS.M;
        h->M = hMAS.L;
        h->N = hMAS.N;
        h->syxx2 = hMAS.syyy2;
        h->syyy2 = hMAS.syxx2;
        h->syzz2 = hMAS.syzz2;
    }
	
	// NEW_HILL properties in CAS
	h->style = hMAS.style;
	
	if(h->style==SQUARED_TERMS && DForm!=DX_DAMAGE)
	{	// matrix C.P (non-zero elements)
		// Elements of P = 2A are scaled by rho^2 and P by 1/rho so CP is * rho
		h->CP11 = 2.*( d->C11*h->syxx2 - d->C12*h->H - d->C13*h->G );
		h->CP12 = 2.*( -d->C11*h->H + d->C12*h->syyy2 - d->C13*h->F );
		h->CP13 = 2.*( -d->C11*h->G - d->C12*h->F + d->C13*h->syzz2 );
		h->CP21 = 2.*( d->C12*h->syxx2 - d->C22*h->H - d->C23*h->G );
		h->CP22 = 2.*( -d->C12*h->H + d->C22*h->syyy2 - d->C23*h->F );
		h->CP23 = 2.*( -d->C12*h->G - d->C22*h->F + d->C23*h->syzz2 );
		h->CP31 = 2.*( d->C13*h->syxx2 - d->C23*h->H - d->C33*h->G );
		h->CP32 = 2.*( -d->C13*h->H + d->C23*h->syyy2 - d->C33*h->F );
		h->CP33 = 2.*( -d->C13*h->G - d->C23*h->F + d->C33*h->syzz2 );
		h->CP66 = 2.*d->C66*h->N;
		if(np==THREED_MPM)
		{	h->CP44 = 2.*d->C44*h->L;
			h->CP55 = 2.*d->C55*h->M;
		}
		
		// Q = PZP matrix (scale by rho^2 and symmetric)
		// same 2D and 3D
		h->Q11 = 8.*(h->G*h->G + h->G*h->H + h->H*h->H);
		h->Q12 = 4.*(h->F*(h->G-h->H) - h->H*(h->G+2.*h->H));
		h->Q13 = 4.*(-h->G*(h->F+h->G) + h->F*h->H - h->G*(h->G+h->H));
		h->Q22 = 8.*(h->F*h->F + h->F*h->H + h->H*h->H);
		h->Q23 = 4.*(-2.*h->F*h->F + h->G*h->H - h->F*(h->G+h->H));
		h->Q33 = 8.*(h->F*h->F + h->F*h->G + h->G*h->G);
		h->Q44 = 2.*h->L*h->L;
		h->Q55 = 2.*h->M*h->M;
		h->Q66 = 2.*h->N*h->N;
	}
}

// Return yield stress for current conditions (p->aint for cum. plastic strain and dalpha/delTime for plastic strain rate)
double OrthoPlasticSoftening::GetYield(AnisoPlasticProperties *p) const
{   switch(hardStyle)
    {   case AP_LINEAR:
            return 1. + Khard*p->aint;
         case AP_NONLINEAR1:
            return pow(1. + Khard*p->aint,nhard);
        case AP_EXPONENTIAL:
            if(alphaMax<0 || p->aint<alphaMax)
                return 1. + Kexp*(1.-exp(-exphard*p->aint));
            else
                return hmax + slopeMin*p->aint;
        default:
            return 1. + Khard*pow(p->aint,nhard);
     }
}

// Get g'(alpha)
double OrthoPlasticSoftening::GetGPrime(AnisoPlasticProperties *p) const
{   switch(hardStyle)
    {   case AP_LINEAR:
            return Khard;
        case AP_NONLINEAR1:
            return Khard*nhard*pow(1. + Khard*p->aint,nhard-1.);
        case AP_EXPONENTIAL:
            if(alphaMax<0 || p->aint<alphaMax)
                return Khard*exp(-exphard*p->aint);
            else
                return slopeMin;
        default:
            if(nhard<1. && p->aint<ALPHA_EPS)
                return Khard*nhard*pow(ALPHA_EPS,nhard-1.);
            else
                return Khard*nhard*pow(p->aint,nhard-1.);
    }
}

// Solve numerically plastic strain increment by explicit radial return
Tensor OrthoPlasticSoftening::SolveForPlasticIncrement(MPMBase *mptr,int np,double fk,Tensor &stk,AnisoPlasticProperties *p,
                                                       CrackAxisProperties *d,HillProperties *h) const
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
	if(h->style==SQUARED_TERMS)
	{	double lambdak = 0;
		Tensor strial = stk;
		double atrial = p->aint;
		
		// P strial = DPhi/dSigma, (2/3)Q strial, and sqrt((2/3)strial Q strial) stored in p->
		AnisoPlasticity::GetHillDfDsigmaPQsigma(strial,np,p,h);
		
		// get (I + lambdak CP)^{-1} for first step when lambdak=0
		Matrix3 IlamCPInv = Matrix3::Identity();
		double lamInv[3][3];
		IlamCPInv.get(lamInv);
		
		while(nsteps<10)
		{	// stk = sigma(lambdak) = (I + lambda CP)^{-1}strial (from above=strial or end of previous step)
			
			// get P stk, (2/3)Q stk, and sqrt((2/3)stk Q stk) (from above using strial or end of previous step)

			// dsigma/dlambda = -(I + lamdak CP)^{-1) CP.sigma(lambdak)
			Tensor CPSigma,dSigdLam;
			CPSigma.xx = h->CP11*stk.xx + h->CP12*stk.yy + h->CP13*stk.zz;
			CPSigma.yy = h->CP21*stk.xx + h->CP22*stk.yy + h->CP23*stk.zz;
			CPSigma.zz = h->CP31*stk.xx + h->CP32*stk.yy + h->CP33*stk.zz;
			CPSigma.xy = h->CP66*stk.xy;
			if(np==THREED_MPM)
			{   CPSigma.yz = h->CP44*stk.yz;
				CPSigma.xz = h->CP55*stk.xz;
			}
			dSigdLam.xx = -lamInv[0][0]*CPSigma.xx - lamInv[0][1]*CPSigma.yy - lamInv[0][2]*CPSigma.zz;
			dSigdLam.yy = -lamInv[1][0]*CPSigma.xx - lamInv[1][1]*CPSigma.yy - lamInv[1][2]*CPSigma.zz;
			dSigdLam.zz = -lamInv[2][0]*CPSigma.xx - lamInv[2][1]*CPSigma.yy - lamInv[2][2]*CPSigma.zz;
			dSigdLam.xy = -CPSigma.xy/(1.+lambdak*h->CP66);
			if(np==THREED_MPM)
			{   dSigdLam.yz = -CPSigma.yz/(1.+lambdak*h->CP44);
				dSigdLam.xz = -CPSigma.xz/(1.+lambdak*h->CP55);
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
			Matrix3 IlamCP = Matrix3(1.+lambdak*h->CP11,lambdak*h->CP12,lambdak*h->CP13,
									lambdak*h->CP21,1.+lambdak*h->CP22,lambdak*h->CP23,
									lambdak*h->CP31,lambdak*h->CP32,1.+lambdak*h->CP33);
			IlamCPInv = IlamCP.Inverse();
			IlamCPInv.get(lamInv);
			stk.xx = lamInv[0][0]*strial.xx + lamInv[0][1]*strial.yy + lamInv[0][2]*strial.zz;
			stk.yy = lamInv[1][0]*strial.xx + lamInv[1][1]*strial.yy + lamInv[1][2]*strial.zz;
			stk.zz = lamInv[2][0]*strial.xx + lamInv[2][1]*strial.yy + lamInv[2][2]*strial.zz;
			stk.xy = strial.xy/(1.+lambdak*h->CP66);
			if(np==THREED_MPM)
			{   stk.yz = strial.yz/(1.+lambdak*h->CP44);
				stk.xz = strial.xz/(1.+lambdak*h->CP55);
			}
			
			// P stk = DPhi/dSigma, (2/3)Q stk, and sqrt((2/3)stk Q stk)
			AnisoPlasticity::GetHillDfDsigmaPQsigma(stk,np,p,h);

			// get new alpha
			p->aint = atrial + lambdak*p->sAQsmag;

			// check convergenve with new values
			double halfsPs = AnisoPlasticity::GetHillMagnitude(stk,h,np);
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
	{	while(nsteps<10)
		{    // find dlambda = f/((df/dsig)C(df/dsig) + g'(alphak))
			GetDfCdf(stk,np,p,d,h);
			double dlambda = fk/(p->dfCdf + sqrt23OversigmaYref*GetGPrime(p));
			
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
			{   stk.xz -= dlambda*p->CdfQsigma.xz;
				stk.yz -= dlambda*p->CdfQsigma.yz;
			}
			
			// total incremental plastic strain accumulated
			AddTensor(&dep,&ddep);
			
			// get new magniture and potential
			p->sAQsmag = AnisoPlasticity::GetHillMagnitude(stk,h,np);
			fk = p->sAQsmag - GetYield(p);
			
			// check for convergence
			nsteps++;
			//cout << "   " << nsteps << ": f = " << fk << ", a = " << p->aint << ", dlambda = " << dlambda << endl;;
			if(fabs(fk)<cutoff) break;
			
		}
	}
    
    // check number of steps
    if(nsteps>hstepsOPSMax)
    {   hstepsOPSMax = nsteps;
#pragma omp critical (output)
        {
            cout << "# OrthoPlasticSoftening took " << nsteps << " steps to converge" << endl;
        }
    }
    
    // return full plastic strain increment
    return dep;
}

// Only called for SQRT_TERMS
// Find C.df and df.C.df at given stress - store results in plastic property variables active only during the loop
// and only for current material point
void OrthoPlasticSoftening::GetDfCdf(Tensor &stk,int np,AnisoPlasticProperties *p,
                                     CrackAxisProperties *r,HillProperties *h) const
{
    // Get dfds (in dfdsPsigma)
    AnisoPlasticity::GetHillDfDsigmaPQsigma(stk,np,p,h);
	
    // get C df (in CdfQsigma) and df C df in (dfCdf), which needs df/dsig
    if(np==THREED_MPM)
    {   p->CdfQsigma.xx = r->C11*p->dfdsPsigma.xx + r->C12*p->dfdsPsigma.yy + r->C13*p->dfdsPsigma.zz;
        p->CdfQsigma.yy = r->C12*p->dfdsPsigma.xx + r->C22*p->dfdsPsigma.yy + r->C23*p->dfdsPsigma.zz;
        p->CdfQsigma.zz = r->C13*p->dfdsPsigma.xx + r->C23*p->dfdsPsigma.yy + r->C33*p->dfdsPsigma.zz;
        p->CdfQsigma.yz = r->C44*p->dfdsPsigma.yz;
        p->CdfQsigma.xz = r->C55*p->dfdsPsigma.xz;
        p->CdfQsigma.xy = r->C66*p->dfdsPsigma.xy;
        p->dfCdf = p->dfdsPsigma.xx*p->CdfQsigma.xx + p->dfdsPsigma.yy*p->CdfQsigma.yy
						+ p->dfdsPsigma.zz*p->CdfQsigma.zz + p->dfdsPsigma.yz*p->CdfQsigma.yz
						+ p->dfdsPsigma.xz*p->CdfQsigma.xz + p->dfdsPsigma.xy*p->CdfQsigma.xy;
    }
    else
    {   p->CdfQsigma.xx = r->C11*p->dfdsPsigma.xx + r->C12*p->dfdsPsigma.yy;
        p->CdfQsigma.yy = r->C12*p->dfdsPsigma.xx + r->C22*p->dfdsPsigma.yy;
        p->CdfQsigma.xy = r->C66*p->dfdsPsigma.xy;
        p->CdfQsigma.zz = r->C13*p->dfdsPsigma.xx + r->C23*p->dfdsPsigma.yy + r->C33*p->dfdsPsigma.zz;
        p->dfCdf = p->dfdsPsigma.xx*p->CdfQsigma.xx + p->dfdsPsigma.yy*p->CdfQsigma.yy
						+ p->dfdsPsigma.xy*p->CdfQsigma.xy + p->dfdsPsigma.zz*p->CdfQsigma.zz;
    }
}

// return pointer to elastic properties
ElasticProperties *OrthoPlasticSoftening::GetElasticPropertiesPointer(void *properties) const
{
    AnisoPlasticProperties *p = (AnisoPlasticProperties *)properties;
    return p->ep;
}

#pragma mark OrthoPlasticSoftening::Accessors

// return unique, short name for this material
const char *OrthoPlasticSoftening::MaterialType(void) const
{    return "Orthotropic plastic softening";
}

// buffer size for mechanical properties
int OrthoPlasticSoftening::SizeOfMechanicalProperties(int &altBufferSize) const
{   altBufferSize = 0;
    return sizeof(AnisoPlasticProperties);
}

// Isotropic material can use read-only initial properties
void *OrthoPlasticSoftening::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{
    AnisoPlasticProperties *p = (AnisoPlasticProperties *)matBuffer;
    p->ep = (ElasticProperties *)&pr;
    
    return p;
}

// history data for softening after one plastic history variable
double *OrthoPlasticSoftening::GetSoftHistoryPtr(MPMBase *mptr) const
{   double *history =  (double *)(mptr->GetHistoryPtr(0));
    return history+1;
}

// store plastic strain in alt strain
int OrthoPlasticSoftening::AltStrainContains(void) const
{    return ENG_BIOT_PLASTIC_STRAIN;
}

// return cracking strain and true, or return false if material has no cracking strain
// (don't change inRtot)
bool OrthoPlasticSoftening::GetCrackingStrain(MPMBase *mptr,Tensor *ecrack,bool is2D,Matrix3 *inRtot) const
{
    // get history
    double *soft = GetSoftHistoryPtr(mptr);
    
    // Before initiatiation, not damage strain
    if(soft[SOFT_DAMAGE_STATE]<predamageState) return false;
    
    // cracking strain in history variable, but in crack axis system
    double *ipsoft = soft + SOFT_NUMBER_HISTORY;
    Tensor eCAS = MakeTensor(ipsoft[ECXX_DAMAGE],0.,0.,0.,ipsoft[GCXZ_DAMAGE],ipsoft[GCXY_DAMAGE]);
    
    // Get R to crack, rotate by Rtot*R0*RToCrack
    Matrix3 Rtot = *inRtot;
    Matrix3 R0 = mptr->GetInitialRotation();
    Rtot *= R0;
    int DForm = GetDForm(soft[SOFT_DAMAGE_STATE]);
    Matrix3 RToCrack;
    GetRToCrack(&RToCrack,soft,is2D,DForm);
    Rtot *= RToCrack;
    *ecrack = Rtot.RVoightRT(&eCAS,false,is2D);
    
    // has cracking strain
    return true;
}

// return cracking COD in the crack axis system
// Only used by DeleteDamaged custom task
Vector OrthoPlasticSoftening::GetCrackingCOD(MPMBase *mptr,bool is2D) const
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

