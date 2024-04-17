/********************************************************************************
    Viscoelastic.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Feb 5 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/Viscoelastic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/MPMWarnings.hpp"

// When uncommented, the constitutive law recalculates stress on each time step
// When commented out, the constitutive law use incremental stress updates
// WARNING: if commented out calculations can not do vertical shifting to implement
//          temperature and moisture dependent properties
#define TOTAL_STRESS_CALC

// class statics
int Viscoelastic::warnExcessiveX = -1;

#pragma mark Viscoelastic::Constructors and Destructors

// Constructor
Viscoelastic::Viscoelastic(char *matName,int matID) : MaterialBase(matName,matID)
{
    ntaus=-1;
    Gk=NULL;
    tauk=NULL;
    G0=0.;
    currentGk=0;
    currentTauk=0;
    aI=40.;
    
    ntausK=0;
    Kk=NULL;
    tauKk=NULL;
    currentKk=0;
    currentTauKk=0;
    
    // MGEOS variables
    pressureLaw = LINEAR_PRESSURE;
    gamma0=1.64;        // dimensionless
    C0=4004000.;        // m/sec
    S1=1.35;            // dimsionless
    S2=0.;                // dimsionless
    S3=0.;                // dimsionless
    Kmax=-1.;            // maxium relative increase allows in K (default no limit)
    
    // WLF parameters
    Tref = -1.;
    C1base10 = 17.44;
    C2 = 51.6;
    
    // Moisture parameters
    mref = -1.;
    Cm1base10 = 10.;
    Cm2base10 = 0.025/0.4;
    
    // vertical shifting (total stress only)
    bTemp.clear();
    bTValue.clear();
    bConc.clear();
    bCValue.clear();
}

#pragma mark Viscoelastic::Initialization

// print mechanical properties to output window
void Viscoelastic::PrintMechanicalProperties(void) const
{
    if(pressureLaw==LINEAR_PRESSURE)
    {   cout << "Pressure law: linear elastic" << endl;
        PrintProperty("K",K*UnitsController::Scaling(1.e-6),"");
        PrintProperty("a",aI,"");
        cout << endl;
    }
    else if(pressureLaw==TIME_DEPENDENT_PRESSURE)
    {   cout << "Pressure law: time-dependent bulk modulus" << endl;
        PrintProperty("K(0)",rho*Kered*UnitsController::Scaling(1.e-6),"");
        PrintProperty("K0",K*UnitsController::Scaling(1.e-6),"");
        PrintProperty("ntausK",(double)ntausK,"");
        cout <<  endl;
        
        for(int i=0;i<ntaus;i++)
        {   PrintProperty("  i",(double)(i+1),"");
            PrintProperty("Kk",Kk[i]*UnitsController::Scaling(1.e-6),"");
            PrintProperty("tauKk",tauKk[i],UnitsController::Label(TIME_UNITS));
            cout << endl;
        }
    }
    else
    {    // core properties
        cout << "Pressure law: MG-EOS" << endl;
        PrintProperty("C0",C0*UnitsController::Scaling(1.e-3),UnitsController::Label(ALTVELOCITY_UNITS));
        PrintProperty("gam0",gamma0,"");
        PrintProperty("K",K*UnitsController::Scaling(1.e-6),"");
        cout << endl;
        
        PrintProperty("S1",S1,"");
        PrintProperty("S2",S2,"");
        PrintProperty("S3",S3,"");
        cout << endl;
        
        // effective volumetric CTE (in ppm/K) alpha = rho0 gamma0 Cv / K
        PrintProperty("aI",1.e6*CTE,"");
        PrintProperty("T0",thermal.reference,"K");
        cout <<  endl;
        
        // Kmax
        if(Kmax>0.)
        {   PrintProperty("Kmax",Kmax," K0");
            PrintProperty("Xmax",Xmax,"");
        }
        else
            cout << "Kmax= no limit";
        cout <<  endl;
    }
    
    PrintProperty("G(0)",rho*Gered*UnitsController::Scaling(1.e-6),"");
    PrintProperty("G0",G0*UnitsController::Scaling(1.e-6),"");
    PrintProperty("ntaus",(double)ntaus,"");
    cout <<  endl;
    
    int i;
    for(i=0;i<ntaus;i++)
    {   PrintProperty("  i",(double)(i+1),"");
        PrintProperty("Gk",Gk[i]*UnitsController::Scaling(1.e-6),"");
        PrintProperty("tauk",tauk[i],UnitsController::Label(TIME_UNITS));
        cout << endl;
    }
    
    // For information, E and nu at time 0
    PrintProperty("E(0)",9.*Kered*Gered*rho/(3.*Kered+Gered)*UnitsController::Scaling(1.e-6),"");
    PrintProperty("nu(0)",(3.*Kered-2.*Gered)/(6.*Kered+2.*Gered),"");
    cout << endl;
    
#ifdef TOTAL_STRESS_CALC
    cout << "Total stress method" << endl;
#else
    cout << "Incremental stress method" << endl;
#endif
    
    // WLF properties
    if(Tref>=0.)
    {   PrintProperty("Tref",Tref,"K");
        PrintProperty("C1",C1base10,"");
        PrintProperty("C2",C2,"");
        cout << endl;
    }
    else if(mref<0.)
        cout << "Isothermal and isosolvent viscoelasticity" << endl;
    else
        cout << "Isothermal viscoelasticity" << endl;

    // WLF moisture properties
    if(mref>=0.)
    {   PrintProperty("mref",mref*concSaturation,"");
        PrintProperty("Cm1",Cm1base10,"");
        PrintProperty("Cm2",Cm2base10,"");
        cout << endl;
    }
    else if(Tref>=0.)
        cout << "Isosolvent viscoelasticity" << endl;
    
#ifdef TOTAL_STRESS_CALC
    if(bTemp.size()==0 || Tref<0)
        cout << "No vertical shifting for temperature";
    else
    {   cout << "Vertical themal shifting:" << endl;
        for(int i=0;i<bTemp.size();i++)
        {   PrintProperty("  T",bTemp[i],"K");
            PrintProperty("  bT",bTValue[i],"");
        }
        cout << endl;
    }
    if(bTemp.size()==0 || mref<0 || diffusion==NULL)
        cout << "No vertical shifting for concentration";
    else
    {   cout << "Vertical concentration shifting:" << endl;
        for(int i=0;i<bConc.size();i++)
        {   PrintProperty("  c",bConc[i]*concSaturation,"K");
            PrintProperty("  bc",bCValue[i],"");
        }
        cout << endl;
    }
#endif
}
    
// Read material properties
// throws std::bad_alloc, SAXException()
char *Viscoelastic::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"G0")==0)
        return UnitsController::ScaledPtr((char *)&G0,gScaling,1.e6);
    
    else if(strcmp(xName,"K")==0)
        return UnitsController::ScaledPtr((char *)&K,gScaling,1.e6);
    
    else if(strcmp(xName,"alpha")==0)
        return((char *)&aI);
    
    else if(strcmp(xName,"ntaus")==0)
    {   input=INT_NUM;
        return((char *)&ntaus);
    }
    
    else if(strcmp(xName,"Gk")==0)
    {   if(Gk==NULL)
        {   if(ntaus<=0)
                ThrowSAXException("Gk found before number of G taus specified.");
            Gk=new double[ntaus];
        }
        currentGk++;
        if(currentGk>ntaus)
            ThrowSAXException("Too many Gk's given.");
        return UnitsController::ScaledPtr((char *)&Gk[currentGk-1],gScaling,1.e6);
    }
    
    else if(strcmp(xName,"tauk")==0)
    {   if(tauk==NULL)
        {   if(ntaus<=0)
                ThrowSAXException("tauk found before number of G taus specified.");
            tauk=new double[ntaus];
        }
        currentTauk++;
        if(currentTauk>ntaus)
            ThrowSAXException("Too many G taus given.");
        return((char *)&tauk[currentTauk-1]);
    }

    else if(strcmp(xName,"ntausK")==0)
    {   input=INT_NUM;
        return((char *)&ntausK);
    }
    
    else if(strcmp(xName,"Kk")==0)
    {   if(Kk==NULL)
        {   if(ntaus<=0)
                ThrowSAXException("Kk found before number of K taus specified.");
            Kk=new double[ntausK];
        }
        currentKk++;
        if(currentKk>ntausK)
            ThrowSAXException("Too many Kk's given.");
        return UnitsController::ScaledPtr((char *)&Kk[currentKk-1],gScaling,1.e6);
    }
    
    else if(strcmp(xName,"tauKk")==0)
    {   if(tauKk==NULL)
        {   if(ntausK<=0)
                ThrowSAXException("tauKk found before number of K taus specified.");
            tauKk=new double[ntausK];
        }
        currentTauKk++;
        if(currentTauKk>ntausK)
            ThrowSAXException("Too many K taus given.");
        return((char *)&tauKk[currentTauKk-1]);
    }
    
    // The follow allow use of MG-EOS law for pressure dependence
    
    else if(strcmp(xName,"pressureLaw")==0)
    {   input=INT_NUM;
        return((char *)&pressureLaw);
    }
    
    else if(strcmp(xName,"gamma0")==0)
    {   input=DOUBLE_NUM;
        return((char *)&gamma0);
    }
    
    else if(strcmp(xName,"C0")==0)
    {   input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&C0,gScaling,1.e3);
    }
    
    else if(strcmp(xName,"S1")==0)
    {   input=DOUBLE_NUM;
        return((char *)&S1);
    }
    
    else if(strcmp(xName,"S2")==0)
    {   input=DOUBLE_NUM;
        return((char *)&S2);
    }
    
    else if(strcmp(xName,"S3")==0)
    {   input=DOUBLE_NUM;
        return((char *)&S3);
    }

    else if(strcmp(xName,"Kmax")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Kmax);
    }
        
    else if(strcmp(xName,"Tref")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Tref);
    }
        
    else if(strcmp(xName,"C1")==0)
    {   input=DOUBLE_NUM;
        return((char *)&C1base10);
    }
        
    else if(strcmp(xName,"C2")==0)
    {   input=DOUBLE_NUM;
        return((char *)&C2);
    }
    
    else if(strcmp(xName,"mref")==0)
    {   input=DOUBLE_NUM;
        return((char *)&mref);
    }
        
    else if(strcmp(xName,"Cm1")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Cm1base10);
    }
        
    else if(strcmp(xName,"Cm2")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Cm2base10);
    }

    else if(strcmp(xName,"bTemp")==0)
    {   bTemp.push_back(0.);
        input=DOUBLE_NUM;
        return (char *)&bTemp[(int)bTemp.size()-1];
    }
    
    else if(strcmp(xName,"bConc")==0)
    {   bConc.push_back(0.);
        input=DOUBLE_NUM;
        return (char *)&bTemp[(int)bConc.size()-1];
    }
    
    else if(strcmp(xName,"bTValue")==0)
    {   bTValue.push_back(0.);
        input=DOUBLE_NUM;
        return (char *)&bTValue[(int)bTemp.size()-1];
    }
    
    else if(strcmp(xName,"bCValue")==0)
    {   bCValue.push_back(0.);
        input=DOUBLE_NUM;
        return (char *)&bCValue[(int)bTemp.size()-1];
    }
    
    return(MaterialBase::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
// throws std::bad_alloc
const char *Viscoelastic::VerifyAndLoadProperties(int np)
{
	if(pressureLaw!=LINEAR_PRESSURE && pressureLaw!=MGEOS_PRESSURE)
		return "Invalid pressure law was selected";
	
	// Shear modulus - must specify ntaus even if it is zero
	if(currentGk<ntaus || currentTauk<ntaus)
		return "Insufficient Gk or tauk for expected number of G taus.";
	if(ntaus<0)
		return "Number of G taus was never entered.";
	
	// bulk modulus. Not entering ntauKs implies time independent bulk modulus
	if(currentKk<ntausK || currentTauKk<ntausK)
		return "Insufficient Kk or tauKk for expected number of K taus.";
	if(ntausK>0)
	{   if(pressureLaw!=LINEAR_PRESSURE)
		return "Time-dependent bulk modulus requires linear pressure law";
		pressureLaw = TIME_DEPENDENT_PRESSURE;
	}
	else
		ntausK = 0;
	
	cout << "Pressure Law: " << pressureLaw << endl;
	
	// Needs non-zero bulk modulus
	if(K<0)
		return "Required bulk modulus must be positive.";
	
	// zero time shear modulus
	Gered = G0;
	TwoGkred = new double[ntaus];
	for(int k=0;k<ntaus;k++)
	{   Gered += Gk[k];
		TwoGkred[k] = 2.*Gk[k]/rho;
	}
	
	// Convert to specific moduli
	Gered /= rho;
	
	// bulk modulus
	if(pressureLaw==LINEAR_PRESSURE)
	{   Kered = K/rho;
	}
	
	else if(pressureLaw==TIME_DEPENDENT_PRESSURE)
	{   // zero time shear modulus
		Kered = K;
		Kkred = new double[ntausK];
		for(int k=0;k<ntausK;k++)
		{   Kered += Kk[k];
			Kkred[k] = Kk[k]/rho;
		}
		
		// Convert to specific moduli
		Kered /= rho;
	}
	
	else if(pressureLaw==MGEOS_PRESSURE)
	{   // Use in place of C0^2. Units are L^2/sec^2 = F/L^2 L^3/mass
		// Equal to reduced bulk modulus
		Kered = C0squared = C0*C0;
		
		// Initial bulk modulus
		K = rho*C0squared;
		
		// expansion coefficients - affect on pressure is handled by eos
		// find for printing (and maybe future large deformation shear)
		double effAlpha = (heatCapacity*gamma0)/C0squared;
		CTE = effAlpha/3.;
		
		// this material not coupled to moisture expansion
		betaI = 0.;
		CME = 0.;
		
		// warning
		if(warnExcessiveX<0)
			warnExcessiveX = warnings.CreateWarning("Compressive strain has exceeded MGEOS law range in Viscoleastic material",-1,0);
		Xmax = GetMGEOSXmax(gamma0,S1,S2,S3,Kmax);
	}
	
	else
		return "Invalid option for the pressure law";
	
	// to absolute CTE and CME
	CTE = 1.e-6*aI;
	CME = betaI*concSaturation;
	
	// for Cp-Cv (units nJ/(g-K^2))
	Ka2sp = 9.*Kered*CTE*CTE;
	
	// WLF coefficients convert to ln aT
	C1 = log(10.)*C1base10;
	
	// Moisture terms input as log ac = -Cm1base(m-mref)/(Cm2base+m-mref)
    // ... convert to use ln ac = - Cm1*(c-cref)/(Cm2+c) where c = m/csat and cref=mref/csat
    // ... Cm1 = Cm1bas*ln(10) and Cm2 = (Cm2base-mref)/csat
    // log ac =
	Cm1 = log(10.)*Cm1base10;
	Cm2 = (Cm2base10-mref)/concSaturation;
	mref /= concSaturation;
	if(mref>=0. && Cm2<=0.)
		return "Cm2 must be greater than mref";
	if(mref>1.)
		return "mref must less than or equal to csat";
	
#ifdef TOTAL_STRESS_CALC
    // vertical shifting temperature same number and sorted
    if(bTemp.size()!=bTValue.size())
        return "Vertical shifting for temperature must have same number of temperature and values";
    for(int i=1;i<bTemp.size();i++)
    {   if(bTemp[i]<=bTemp[i-1])
            return "Vertical shifting temperature points must monotonically increase";
    }

    // vertical shifting temperature same number and sorted (and scaled to csat)
    if(bConc.size()!=bCValue.size())
        return "Vertical shifting for concentration must have same number of concentrations and values";
    if(bConc.size()>0) bConc[0] /= concSaturation;
    for(int i=1;i<bConc.size();i++)
    {   bConc[i] /= concSaturation;
        if(bConc[i]<=bConc[i-1])
            return "Vertical shifting concentration points must monotonically increase";
    }
#endif

    // call super class
    return MaterialBase::VerifyAndLoadProperties(np);
}
    
// plane stress not allows for some viscoelasticity
// throws CommonException()
void Viscoelastic::ValidateForUse(int np) const
{
    if(np==PLANE_STRESS_MPM)
    {   if(pressureLaw==MGEOS_PRESSURE)
        {    throw CommonException("Viscoelastic materials in plane stress require linear pressure model",
                                      "Viscoelastic::ValidateForUse");
        }
        if(artificialViscosity)
        {    throw CommonException("Viscoelastic materials in plane stress do not support artificial viscosity",
                                  "Viscoelastic::ValidateForUse");
        }
    }

    //call super class (why can't call super class?)
    MaterialBase::ValidateForUse(np);
}

#pragma mark Viscoelastic::History Data Methods

// create and return pointer to history variables
// initialize all to zero
// throws std::bad_alloc
char *Viscoelastic::InitHistoryData(char *pchr,MPMBase *mptr)
{
    // count the variables
    numJHistory = pressureLaw==MGEOS_PRESSURE ? 2 : 0;

    if(fmobj->IsThreeD())
        numHistory = numJHistory + 6*ntaus + ntausK;
    else
        numHistory = numJHistory + 4*ntaus + ntausK;
    
    // exit if none
    if(numHistory==0) return NULL;
    
    // all zeros
    double *p = CreateAndZeroDoubles(pchr,numHistory);
    
    // J history starts at 1
    if(numJHistory>0)
    {   p[MGJ_HISTORY] = 1.;
        p[MGJRES_HISTORY] = 1.;
    }
    
    // return pointer
    return (char *)p;
}

// reset history data
void Viscoelastic::ResetHistoryData(char *pchr,MPMBase *mptr)
{    ZeroDoubles(pchr,numHistory);
    // J history starts at 1
    if(numJHistory>0)
    {   double *p = (double *)pchr;
        p[MGJ_HISTORY] = 1.;
        p[MGJRES_HISTORY] = 1.;
    }
}

// Number of history variables - only the plastic law
int Viscoelastic::NumberOfHistoryDoubles(void) const { return numHistory; }

#pragma mark Viscoelastic::Methods

/* Take increments in strain and calculate new Particle: strains, rotation strain,
    stresses, strain energy,
    du are (gradient rates X time increment) to give deformation gradient change
    For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
    This material tracks pressure and stores deviatoric stress only in particle stress tensor
 */
void Viscoelastic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{
    // Note: cannot call generic method because need detdF and Vrot below
    
    // previous deformation gradient and stretch
    Matrix3 pFnm1 = mptr->GetDeformationGradientMatrix();
    
    // get incremental deformation gradient and decompose it
    const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
    
    // Update total deformation gradient (saved on particle at the end)
    Matrix3 pF = dF*pFnm1;

    // two decompositions to get Rn and Rn-1 and U in initial configuration
    Matrix3 Rnm1,Rn;
    pFnm1.RightDecompose(&Rnm1,NULL);
    Matrix3 Umat = pF.RightDecompose(&Rn,NULL);

    // Rn = dR*Rnm1 or dR = Rn*Rnm1^T
    Matrix3 dR = Rn*Rnm1.Transpose();
    
    // get strain and strain increments in initial configuration = Rn^T(dF-dR)F(n-1)
    Matrix3 dFmdR = dF - dR;
    Matrix3 detot = Rn.Transpose()*(dFmdR*pFnm1);
    
#ifdef TOTAL_STRESS_CALC
    // shift of elastic modulus based on mptr->pPreviousTemperature and mptr->pPreviousConcentration
    double bshift = GetVertialShift(mptr,Tref,bTemp,bTValue,mref,bConc,bCValue);
#else
    // incremental method does not support vertical shifting
    double bshift = 1.;
#endif
    
    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
    double eres=CTE*res->dT;
    if(DiffusionTask::HasFluidTransport())
        eres+=CME*res->dC;
    
    // initial terms and get delV
    // this delV is based on dezz=0 for both plane stress and plane strain
    double dispEnergy = 0.,detdF = 1.,J = 1.,dJres = 1.;
    double traceDe = detot.trace();
    double delV = traceDe - 3.*eres;
    
    // history variables start after J history (in MGEOS_PRESSURE)
    double *ak =(double *)(mptr->GetHistoryPtr(0));
    int ai = numJHistory;

    // find dJ and J if needed (only for MGEOS_PRESSURE)
    if(numJHistory>0)
    {    // large strain volume change
        detdF = dF.determinant();
        J = detdF*ak[MGJ_HISTORY];
        ak[MGJ_HISTORY] = J;
        
        // account for residual strains if needed
        double dJres = exp(3.*eres);
        double Jres = dJres*ak[MGJRES_HISTORY];
        ak[MGJRES_HISTORY] = Jres;
    }

    // deviatoric strains increment in de in the initial configuration
    // Actually de is finding 2*d(dev e) to avoid many multiplies by two
    Tensor de;
    double dV = traceDe/3.;
    de.xx = 2.*(detot(0,0) - dV);
    de.yy = 2.*(detot(1,1) - dV);
    de.zz = 2.*(detot(2,2) - dV);
    de.xy = 2.*detot(0,1);
    if(np==THREED_MPM)
    {   de.xz = 2.*detot(0,2);
        de.yz = 2.*detot(1,2);
    }

    // Find initial 2*e(t) (total deviatoric strain) in ed in initial config
    // This has tensorial shear strain from symmetric U
    Tensor ed;
    double thirdV = Umat.trace()/3.;
    ed.xx = 2.*(Umat(0,0)-thirdV);
    ed.yy = 2.*(Umat(1,1)-thirdV);
    ed.zz = 2.*(Umat(2,2)-thirdV);
    ed.xy = 2.*Umat(0,1);                 // = 2 * (1/2)*(Umat(0,1)+Umat(1.0))
    if(np==THREED_MPM)
    {   ed.xz = 2.*Umat(0,2);
        ed.yz = 2.*Umat(1,2);
    }
    // subtract de to get ed as deviatoric strain at start of the time step
    // ednp1 is deviatoric stress at end of the time step
    Tensor ednp1 = ed;
    SubTensor(&ed,&de);

#ifdef TOTAL_STRESS_CALC
    // get new particle deviatoric stresses - elastic part = 2G ednp1
    // note that dsig is actaully total sig in this mode
    Tensor dsig = ednp1;
#else
    // increment particle deviatoric stresses - elastic part = 2G de
    Tensor dsig = de;
#endif
    // strain has factor of 2 so scale by G and not 2G
    ScaleTensor(&dsig,bshift*Gered);

    // get effective time increment
    double delEffTime = GetEffectiveIncrement(mptr,res,delTime,Tref,C1,C2,mref,Cm1,Cm2,0.,1.);
    
    // get internal variable increments, update them, add to incremental stress, and get dissipated energy
    // For plane stress, this gets initial deviatoric stress update only
    Tensor dak;
    int k;
    double omtmp,omtmp2;
    for(k=0;k<ntaus;k++)
    {   GetAlphaArgs(delEffTime,tauk[k],omtmp,omtmp2);
        dak.xx = omtmp*(0.5*ed.xx-ak[ai+XX_HISTORY]) + omtmp2*de.xx;
        dak.yy = omtmp*(0.5*ed.yy-ak[ai+YY_HISTORY]) + omtmp2*de.yy;
        dak.xy = omtmp*(0.5*ed.xy-ak[ai+XY_HISTORY]) + omtmp2*de.xy;
        dak.zz = omtmp*(0.5*ed.zz-ak[ai+ZZ_HISTORY]) + omtmp2*de.zz;
        
        // extra terms for 3D
        if(np==THREED_MPM)
        {    // internal variables
            dak.xz = omtmp*(0.5*ed.xz-ak[ai+XZ_HISTORY]) + omtmp2*de.xz;
            dak.yz = omtmp*(0.5*ed.yz-ak[ai+YZ_HISTORY]) + omtmp2*de.yz;
            
            // update history on particle
            ak[ai+XX_HISTORY] += dak.xx;
            ak[ai+YY_HISTORY] += dak.yy;
            ak[ai+XY_HISTORY] += dak.xy;
            ak[ai+ZZ_HISTORY] += dak.zz;
            ak[ai+XZ_HISTORY] += dak.xz;
            ak[ai+YZ_HISTORY] += dak.yz;
            
            // dissipation (updated deviatoric strain minus updated alpha, remove factor of 2 from strain)
            // except double shear
            dispEnergy += bshift*TwoGkred[k]*(dak.xx*(0.5*ednp1.xx-ak[ai+XX_HISTORY])
                                       + dak.yy*(0.5*ednp1.yy-ak[ai+YY_HISTORY])
                                       + dak.zz*(0.5*ednp1.zz-ak[ai+ZZ_HISTORY])
                                       + dak.xy*(ednp1.xy-2.*ak[ai+XY_HISTORY])
                                       + dak.xz*(ednp1.xz-2.*ak[ai+XZ_HISTORY])
                                       + dak.yz*(ednp1.yz-2.*ak[ai+YZ_HISTORY]));
        }
        else if(np!=PLANE_STRESS_MPM)
        {   // update history on particle
            ak[ai+XX_HISTORY] += dak.xx;
            ak[ai+YY_HISTORY] += dak.yy;
            ak[ai+XY_HISTORY] += dak.xy;
            ak[ai+ZZ_HISTORY] += dak.zz;

            // dissipation  (updated deviatoric strain minus updated alpha, remove factor of 2 from strain)
            // except double shear
            dispEnergy += bshift*TwoGkred[k]*(dak.xx*(0.5*ednp1.xx-ak[ai+XX_HISTORY])
                                       + dak.yy*(0.5*ednp1.yy-ak[ai+YY_HISTORY])
                                       + dak.zz*(0.5*ednp1.zz-ak[ai+ZZ_HISTORY])
                                       + dak.xy*(ednp1.xy-2.*ak[ai+XY_HISTORY]));
        }
        
#ifdef TOTAL_STRESS_CALC
        // add to total stress using updated alphas
        dsig.xx -= bshift*TwoGkred[k]*ak[ai+XX_HISTORY];
        dsig.yy -= bshift*TwoGkred[k]*ak[ai+YY_HISTORY];
        dsig.zz -= bshift*TwoGkred[k]*ak[ai+ZZ_HISTORY];
        dsig.xy -= bshift*TwoGkred[k]*ak[ai+XY_HISTORY];
        if(np==THREED_MPM)
        {   dsig.xz -= bshift*TwoGkred[k]*ak[ai+XZ_HISTORY];
            dsig.yz -= bshift*TwoGkred[k]*ak[ai+YZ_HISTORY];
        }
#else
        // add to stress increments
        dsig.xx -= TwoGkred[k]*dak.xx;
        dsig.yy -= TwoGkred[k]*dak.yy;
        dsig.zz -= TwoGkred[k]*dak.zz;
        dsig.xy -= TwoGkred[k]*dak.xy;
        if(np==THREED_MPM)
        {   dsig.xz -= TwoGkred[k]*dak.xz;
            dsig.yz -= TwoGkred[k]*dak.yz;
        }
#endif
        
        // next history
        ai += (np==THREED_MPM) ? 6 : 4 ;
    }
    
#ifdef TOTAL_STRESS_CALC
    // If pressure is not MGEOS, find elastic pressure (in dP)
#else
    // If pressure is not MGEOS, find elaastic pressure increment (in dP)
#endif
    // Final pressure update, artifical viscosity, and energy done in UpdatePressure() below
    double dP=0.,Vstar=0.;
    if(pressureLaw!=MGEOS_PRESSURE)
    {    // dP and Vstar needed in both these cases
        double dTemp=mptr->pPreviousTemperature-thermal.reference;
        double eresStretch=CTE*dTemp;
        if(DiffusionTask::HasDiffusion())
        {   double dConc = diffusion->GetDeltaConcentration(mptr);
            eresStretch += CME*dConc;
        }
        
        // find current V* (diagonal has 1+eii) from U in updated initial config
        Vstar = Umat(0,0)+Umat(1,1)+Umat(2,2)-3.-3.*eresStretch;
        
#ifdef TOTAL_STRESS_CALC
        // elastic pressure
        dP = -bshift*Kered*Vstar;
#else
        // elastic pressure increment
        dP = -Kered*delV;
#endif
        // Get V* at start of the time step
        Vstar -= delV;
        
        // for time dependent pressure, update alpha and add remaining terms to dP
        if(pressureLaw==TIME_DEPENDENT_PRESSURE)
        {    // viscoelastic history
            for(k=0;k<ntaus;k++)
            {   GetAlphaArgs(delEffTime,tauKk[k],omtmp,omtmp2);
                double dalphaV = omtmp*(Vstar-ak[ai]) + 2.*omtmp2*delV;
                
                if(np!=PLANE_STRESS_MPM)
                {   // update history on particle (plane stress done later)
                    ak[ai] += dalphaV;
                    
                    // dissipation (updated volumetric strain minus updated alpha)
                    dispEnergy += bshift*Kkred[k]*dalphaV*(Vstar+delV-ak[ai]);
                }
                
#ifdef TOTAL_STRESS_CALC
                // add to pressure
                dP += bshift*Kkred[k]*ak[ai];
#else
                // add to pressure increments
                dP += Kkred[k]*dalphaV;
#endif
                
                // next history
                ai++;
            }
        }
    }

    // For plane stress, find dezz and adjust all terms (not allowed for MGEOS_PRESSURE)
    if(np==PLANE_STRESS_MPM)
    {   // equal to Ge*phi in materials.pdf (because GetPSArg() has factor 1/2)
        double phi = Gered;
        for(k=0;k<ntaus;k++) phi -= bshift*TwoGkred[k]*GetPSArg(delEffTime,tauk[k]);

        // equal to Ke*phi in materials.pdf
        double phiK = Kered;
        if(pressureLaw==TIME_DEPENDENT_PRESSURE)
        {   for(k=0;k<ntausK;k++) phiK -= bshift*Kkred[k]*2.*GetPSArg(delEffTime,tauKk[k]);
        }
        
        // get dezz
        double dezz = (dP - dsig.zz)/(phiK + 4.*phi/3.);
        double thirddezz = dezz/3.;

        // adjust deviatoric stress increment
        double ds = 2.*phi*thirddezz;
        dsig.xx -= ds;
        dsig.yy -= ds;
        dsig.zz += 2.*ds;           // should stay zero
        
        // adjust delV for use in pressure update (done below)
        delV += dezz;
        
        // set input strain increment to calculated result (used in work below)
        detot(2,2) = dezz;
        
        // adjust particle deformation gradient (stored below)
        pF(2,2) *= (1.+dezz);
        
        // adjust de for history update and disspated energy (factor of 2 is added)
        de.xx -= 2.*thirddezz;
        de.yy -= 2.*thirddezz;
        de.zz += 4.*thirddezz;
        
        // update history and get dissipation
        ai = numJHistory;
        for(k=0;k<ntaus;k++)
        {   GetAlphaArgs(delEffTime,tauk[k],omtmp,omtmp2);
            dak.xx = omtmp*(0.5*ed.xx-ak[ai+XX_HISTORY]) + omtmp2*de.xx;
            dak.yy = omtmp*(0.5*ed.yy-ak[ai+YY_HISTORY]) + omtmp2*de.yy;
            dak.xy = omtmp*(0.5*ed.xy-ak[ai+XY_HISTORY]) + omtmp2*de.xy;
            dak.zz = omtmp*(0.5*ed.zz-ak[ai+ZZ_HISTORY]) + omtmp2*de.zz;

            // update history on particle
            ak[ai+XX_HISTORY] += dak.xx;
            ak[ai+YY_HISTORY] += dak.yy;
            ak[ai+XY_HISTORY] += dak.xy;
            ak[ai+ZZ_HISTORY] += dak.zz;
            
            // dissipation
            dispEnergy += bshift*TwoGkred[k]*(dak.xx*(0.5*(ed.xx+de.xx)-ak[ai+XX_HISTORY])
                                       + dak.yy*(0.5*(ed.yy+de.yy)-ak[ai+YY_HISTORY])
                                       + dak.zz*(0.5*(ed.zz+de.zz)-ak[ai+ZZ_HISTORY])
                                       + dak.xy*((ed.xy+de.xy)-2.*ak[ai+XY_HISTORY]));
            
            // next history
            ai += 4;
        }
        
        // change pressure
        dP -= phiK*dezz;
        
        if(pressureLaw==TIME_DEPENDENT_PRESSURE)
        {   // update history and get dissipation
            for(k=0;k<ntaus;k++)
            {   GetAlphaArgs(delEffTime,tauKk[k],omtmp,omtmp2);
                double dalphaV = omtmp*(Vstar-ak[ai]) + 2.*omtmp2*delV;
                
                // add to pressure increments
                ak[ai] += dalphaV;
                
                // dissipation (updated volumetric strain minus updated alpha)
                dispEnergy += bshift*Kkred[k]*dalphaV*(Vstar+delV-ak[ai]);
                    
                // next history
                ai++;
            }
        }
    }
    
    // Update particle deviatoric stresses
    bool is2D = np==THREED_MPM ? false : true ;

    // get current particle stress pointer
    Tensor *sp = mptr->GetStressTensor();
    
#ifdef TOTAL_STRESS_CALC
    // get stress in initial configuration (for energy calculations below)
    Tensor sp1 = dsig;
    
    if(numJHistory>0)
    {   // convert sigma(n+1)/rho0 to sigma(n+1)/rho(n+1) using J*V0*sigma(n)/mp
        ScaleTensor(&sp1,J);
    }
#else
    // rotate prior stress to initial configuration (for energy calculations below)
    Tensor sp1 = Rnm1.RTVoightR(sp,true,is2D);

    if(numJHistory>0)
    {   // convert old sigma(n)/rho(n) to new sigma(n)/rho(n+1)
        //  where dF*sigma(n)/rho(n) = (V(n+1)/V(n))*(V(n)*sigma(n)/mp) = sigma(n)/rho(n+1)
        ScaleTensor(&sp1,detdF);
        
        // convert dsig/rho0 to dsig/rho(n+1) using J*V0*dsig/mp
        ScaleTensor(&dsig,J);
    }
    
    // add stress increment in initial configuration
    AddTensor(&sp1,&dsig);
#endif
    
    // incremental work energy = shear energy (dilation and residual energy done in update pressure)
    double shearEnergy = sp1.xx*detot(0,0) + sp1.yy*detot(1,1) + sp1.zz*detot(2,2) + 2.*sp1.xy*de.xy;
    if(np==THREED_MPM)
    {   shearEnergy += 2.*(sp1.xz*de.xz + sp1.yz*de.yz);
    }
    mptr->AddWorkEnergyAndResidualEnergy(shearEnergy,0.);

    // now rotate stress to current configuration
    Tensor spr = Rn.RVoightRT(&sp1,true,is2D);
    *sp = spr;
    
    // finish particle updates
    mptr->SetDeformationGradientMatrix(pF);
    
    // Now update pressure
    double dTq0 = dP;               // pass pressure increment found above, it returns as dTq0
    UpdatePressure(mptr,delV,res,eres,detdF,dJres,delTime,dTq0,dispEnergy);
    
    // dissipated energy per unit mass (dPhi/(rho0 V0)) (Legacy nJ/g)
    mptr->AddPlastEnergy(dispEnergy);
    
    // heat energy
    IncrementHeatEnergy(mptr,dTq0,dispEnergy);
}

// Get parameters to update alpha variables
void Viscoelastic::GetAlphaArgs(double delTime,double tau,double &omtmp,double &omtmp2) const
{   double x = delTime/tau;
    if(x>0.003)
    {   double tmp = exp(-0.5*x);
        omtmp = 1. - tmp*tmp;               // 1-exp(-dt/tau)
        omtmp2 = 0.5*(1. - tmp);            // 0.5*(1-exp(-dt/(2*tau))) because de has factor of 2
    }
    else
    {   // double precision accurate and maybe more precision than 1-exp(-x)
        omtmp = x*(1.-0.5*x*(1-(x/3.)*(1-0.25*x)));         // 1-exp(-dt/tau)
        
        // 0.5*(1-exp(-dt/(2*tau))) because de has factor of 2
        x *= 0.5;
        omtmp2 = 0.5*x*(1.-0.5*x*(1-(x/3.)*(1-0.25*x)));
   }
}

// Get parameter for plane stress calculations
// Get 0.5*(1-exp(-dt/(2*tau))) (0.5 because of 2Gk has factor of 2)
double Viscoelastic::GetPSArg(double delTime,double tau) const
{   double x = 0.5*delTime/tau;
    if(x>0.003)
        return 0.5*(1. - exp(-x));
    
    // double precision accurate and maybe more precision than 1-exp(-x)
    return 0.5*x*(1.-0.5*x*(1.-(x/3.)*(1.-0.25*x)));
}

// This method handles the pressure equation of state.
#ifdef TOTAL_STRESS_CALC
// On input dP is is new pressure except for MGEOS
#else
// On input dP is the new pressure increment except for MGEOS
#endif
// 1. For MGEOS, get new pressure
// 2. Artifical visocist calculations
// 3. Update particle pressure
// 4. Increment the particle energy
void Viscoelastic::UpdatePressure(MPMBase *mptr,double delV,ResidualStrains *res,
                                  double eres,double detdF,double dJres,
                                  double delTime,double &dTq0,double &AVEnergy) const
{
    if(pressureLaw!=MGEOS_PRESSURE)
    {
#ifdef TOTAL_STRESS_CALC
        // new pressure passed in dTq0
#else
        // new pressure increment passsed in dP
#endif
        double dP = dTq0 ;
        
        // artifical viscosity
        // delV is total incremental volumetric strain = total Delta(V)/V
        // not allowed in plane stress because added pressure would change szz
        // may not be valid for time dependent bulk modulus, but is allowed
        if(delV<0. && artificialViscosity)
        {    // Wants K/rho
            double QAVred = GetArtificialViscosity(delV/delTime,sqrt(Kered),mptr);
            AVEnergy += fabs(QAVred*delV);
            dP += QAVred;
        }
        
#ifdef TOTAL_STRESS_CALC
        // get average pressure and set new pressure
        double avgP = 0.5*(mptr->GetPressure()+dP);
        mptr->SetPressure(dP);
#else
        // increment pressure and get average pressure
        mptr->IncrementPressure(dP);
        double avgP = mptr->GetPressure()-0.5*dP;
#endif

        // work energy is dU = -P dV + s.de(total)
        // Here do hydrostatic term
        // Work energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
        mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-3.*avgP*eres);
        
        // Isoentropic temperature rise = -(K 3 alpha T)/(rho Cv) (dV/V) (only spot for this material)
        dTq0 = -3.*Kered*CTE*mptr->pPreviousTemperature*(delV+3*eres)/GetHeatCapacity(mptr);
    }
    else
    {    // J is total volume change - may need to reference to free-swelling volume if that works
        // Note that swelling looks like a problem because the sums of strains needs to be adjusted
        //        to stress-free state
        
        // history pointer
        double *h =(double *)(mptr->GetHistoryPtr(0));
        double J = h[MGJ_HISTORY];
        double Jres = h[MGJRES_HISTORY];
        
        // previous pressure
        double P,P0 = mptr->GetPressure();
        double delVMG = 1. - 1./detdF;            // (Vnp1-Vn)/Vnp1
        double Kred;
        
        // M-G EOS
        // Want specific pressure or pressure over current density (using J = rho0/rho)
        if(J<1.)
        {    // new compression J(k+1) = 1-x(k+1)
            double x = 1.-J;
            
            if(x>Xmax)
            {    // law not valid if gets too high
                if(warnings.Issue(warnExcessiveX,-1)==GAVE_WARNING)
                {
#pragma omp critical (output)
                    {    cout << "# Excessive x = " << x << " causing Kred to increase more than " << Kmax << " fold" << endl;
                        mptr -> Describe();
                    }
                }
                
                // reset
                x = Xmax;
            }
            
            // get reduced bulk modulus
            Kred = C0squared*GetMGEOSKRatio(x,gamma0,S1,S2,S3);
            
            // compression law
            // denominator = 1 - S1*x - S2*x^2 - S3*x^3
            double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));
            
            // current effective and reduced (by rho0) bulk modulus
            double Klawred = C0squared*(1.-0.5*gamma0*x)*denom*denom;
            
            // Pressure from bulk modulus and an energy term
            double e = mptr->GetInternalEnergy();
            P = J*(Klawred*x + gamma0*e);

            // particle isoentropic temperature increment
            dTq0 += -J*gamma0*mptr->pPreviousTemperature*delV;
        }
        else
        {   // In tension hyperelastic law P = - K0(J-1)
            double Jeff = J/Jres;
            Kred = C0squared*Jeff;
            P = -J*C0squared*(Jeff-1.);
            
            // particle isentropic temperature increment
            double Kratio = Jeff;
            dTq0 += -J*Kratio*gamma0*mptr->pPreviousTemperature*delVMG;
        }
        
        // artifical viscosity
        // delVMG is total incremental volumetric strain = total Delta(V)/V
        if(delVMG<0. && artificialViscosity)
        {    double QAVred = GetArtificialViscosity(delVMG/delTime,sqrt(Kred*J),mptr);
            AVEnergy += fabs(QAVred*delVMG);
            P += QAVred;
        }
        
        // set final pressure
        mptr->SetPressure(P);
        
        // work energy is dU = -P dV + s.de(total)
        // Here do hydrostatic terms, deviatoric later
        double avgP = 0.5*(P0+P);
        double delVres = 1. - 1./dJres;
        mptr->AddWorkEnergyAndResidualEnergy(-avgP*delVMG,-avgP*delVres);
    }
}

// convert J to K using isotropic method
Vector Viscoelastic::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{   double nuLS = (3.*Kered-2.*Gered)/(6.*Kered+2.*Gered);
    return IsotropicJToK(d,C,J0,np,nuLS,Gered*rho);
}

// From thermodyanamics Cp-Cv = 9 K a^2 T/rho
// Ka2sp in nJ/(g-K^2) so output in nJ/(g-K)
// Here using K0 and rho0 - could modify if needed
double Viscoelastic::GetCpMinusCv(MPMBase *mptr) const
{   return mptr!=NULL ? Ka2sp*mptr->pPreviousTemperature : Ka2sp*thermal.reference;
}

// Get vertical shift
double Viscoelastic::GetVertialShift(MPMBase *mptr,double T0,vector<double> Txpts,vector<double> Typts,
                                     double m0,vector<double> Cxpts,vector<double> Cypts)
{
    double btot = 1.;
    
    // First check shift to reference temperature
    if(T0>=0.)
    {   btot = PiecewiseInterpolate(mptr->pPreviousTemperature, Txpts, Typts);
    }
    
    // Add concentration shift to reference concentration (requires diffusion task)
    if(m0>=0. && diffusion!=NULL)
    {   btot = PiecewiseInterpolate(mptr->pDiff[0]->prevConc, Cxpts, Cypts);
    }

    return btot;
}

// Get effective time increment for viscoelatic materials
double Viscoelastic::GetEffectiveIncrement(MPMBase *mptr,ResidualStrains *res,double dRealTime,
                                           double T0,double cT1,double cT2,double m0,double cC1,
                                           double cC2,double kms,double amu)
{
    // if reference properties not set, using actual time
    if(T0<0. && m0<0.) return dRealTime;
    
    // we want to find R = aT(T)am(T)/(aT(T+dT)am(c+dc))
    // we start ln R = (ln aT(Told)-ln aT(Tnew)) + (ln am(cold)-ln am(cnew)) = lnR
    double lnR = 0.;

    // First check temperature
	double mlogaold=0.,mlogamold=0.;
    if(T0>=0.)
    {   // previous temperature
        double Tnew = mptr->pPreviousTemperature;
        double Told = Tnew-res->dT;
        
        // freeze viscoleaticity below WLF limit of cT2 degrees below reference temperature
        if(Tnew<=T0-cT2 || Told<=T0-cT2) return 0.;
        
        // only deviates from ln 1=0 when dT changes
        if(!DbleEqual(res->dT,0.))
        {   if(cT1>0)
            {   // WLF equation (using ln) is ln aT = -C1(T-T0)/(C2+T-T0)
                // which leads to this R
                lnR = cT1*cT2*res->dT/((cT2+Tnew-T0)*(cT2+Told-T0));
            }
            else
            {   // Arhenius methods ln aT = -C1(1/T - 1/T0) where C1 = -Delta H/R
                lnR = -cT1*res->dT/(Tnew*Told);
            }
        }
        
        // will need this even if does not change (previous -ln aT)
        if(cT1>0)
        {   // WLF version
            mlogaold = cT1*(Told-T0)/(cT2+Told-T0);
        }
        else
        {   // Arhenius version
            mlogaold= cT1*(1/Told - 1/T0);
        }
    }
    
    // Up to here
    // lnR = ln aT(T)/(aT(T+dT)) and mlogaold = - ln aT(T)

    // now check if moisture effect too (only if simulation has a concentration in pDiff[0])
    // Input properties converted to give ln ac = -cC1*(c-cref)/(cC2+c) where c=m/csat and cref=mref/csat
    if(m0>=0. && diffusion!=NULL)
	{   // concentrations from the grid
		double cnew = mptr->pDiff[0]->prevConc;
		double cold = cnew - res->dC;
		
		// moisture shift (previous -ln ac)
		mlogamold = cC1*(cold-m0)/(cC2+cold);
		
		// only need more when dC changes)
		if(!DbleEqual(res->dC,0.))
		{   // add to temperature changes using (ln ac(cold)-ln ac(cnew))
			// when done, lnR+del = ln aT(T)ac(c)/(aT(T+dT)ac (c+dc))
			double del = cC1*(cC2+m0)*res->dC/((cC2+cnew)*(cC2+cold));

			// add to ln R from temperature
			lnR += del;
		}
		
		// combine temperature and moisture previous shift
		mlogaold += mlogamold;
	}
	
    // When here have
    // lnR = ln aT(T)ac(c)/(aT(T+dT)ac(c+dc)) and mlogaold = -ln aT(T)ac(c)
    
    // return the effective increment = dt*scale/(aT ac)
    return dRealTime*GetDtScale(lnR)*exp(mlogaold);
}

// get scale factor for effective time integration
double Viscoelastic::GetDtScale(double lnR)
{
	double scale = 1.;
	if(lnR!=0.)
	{	double R = exp(lnR);
		if(fabs(R-1.)<0.05)
			scale = (9.+R*(19.-R*(5.-R)))/24.;
		else
			scale = (R-1.)/lnR;
	}
	return scale;
}

#pragma mark Viscoelastic::Accessors

// return material type
const char *Viscoelastic::MaterialType(void) const { return "Viscoelastic"; }

// Calculate wave speed in mm/sec
// Uses sqrt((K +4Ge/3)/rho) which is probably the maximum wave speed possible
double Viscoelastic::WaveSpeed(bool threeD,MPMBase *mptr) const { return sqrt((Kered + 4.*Gered/3.)); }

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatoric stress into full stress
Tensor Viscoelastic::GetStress(Tensor *sp,double pressure,MPMBase *mptr) const
{    return GetStressPandDev(sp,pressure,mptr);
}

// store a new total stress on a particle's stress and pressure variables
void Viscoelastic::SetStress(Tensor *spnew,MPMBase *mptr) const
{    SetStressPandDev(spnew,mptr);
}

// Increment thickness (zz) stress through deviatoric stress and pressure
void Viscoelastic::IncrementThicknessStress(double dszz,MPMBase *mptr) const
{    IncrementThicknessStressPandDev(dszz,mptr);
}

// if a subclass material supports artificial viscosity, override this and return TRUE
bool Viscoelastic::SupportsArtificialViscosity(void) const { return true; }

// Calculate current wave speed in mm/sec. Uses sqrt((K+4G/3)/rho) which is dilational wave speed
// but K/rho = Kred*J and G/rho = Gred*J (in mm^2/sec^2)
double Viscoelastic::CurrentWaveSpeed(bool threeD,MPMBase *mptr,int offset) const
{
    double KcurrRed = Kered;
    
    if(pressureLaw==MGEOS_PRESSURE)
    {   // compressive volumetric strain x = 1-J
        double J = mptr->GetRelativeVolume();
        
        // get K/rho0, but this ignores slope of energy term
        if(J<1.)
        {   double x = 1. - J;
            
            if(x<Xmax)
            {    // compression law
                // denominator = 1 - S1*x - S2*x^2 - S3*x^3
                double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));
            
                // current effective and reduced (by rho0) bulk modulus
                KcurrRed = C0squared*(1.-0.5*gamma0*x)*denom*denom;
            }
            else
            {    // truncate if law seems bad
                KcurrRed = C0squared*Kmax;
            }
        }
        
        //KcurrRed *= J;          // converts to K/rho, but keep K/rho0 for hypoelastic material
    }
    
    // return current save speed
    return sqrt((KcurrRed + 4.*Gered/3.));
}

// Get current relative volume change = J (which this material tracks)
double Viscoelastic::GetCurrentRelativeVolume(MPMBase *mptr,int offset) const
{   if(numJHistory==0) return 1.;
    double *h =(double *)mptr->GetHistoryPtr(offset);
    return h[MGJ_HISTORY];
}

// not supported yet, need to deal with aniostropi properties
bool Viscoelastic::SupportsDiffusion(void) const
{   return DiffusionTask::HasPoroelasticity() ? false : true;
}
