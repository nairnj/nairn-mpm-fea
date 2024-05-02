/********************************************************************************
    TIViscoelastic.cpp
    nairn-mpm-fea
    
    Created by John Nairn, Jan 7, 2021.
    Copyright (c) 2021 John A. Nairn, All rights reserved.
 
	Small, strain linear viscoelastic material
    The materials is transversly isotropic
********************************************************************************/

#include "stdafx.h"
#include "Materials/TIViscoelastic.hpp"
#include "Materials/Viscoelastic.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/CommonException.hpp"
#include "Global_Quantities/ThermalRamp.hpp"

// When uncommented, the constitutive law recalculates stress on each time step
// When commented out, the constitutive law use incremental stress updates
// WARNING: if commented out calculations can not do vertical shifting to implement
//          temperature and moisture dependent properties
#define TOTAL_STRESS_CALC_TI

#pragma mark TIViscoelastic::Constructors and Destructors

// This method should call parent class and then fill in any new initiatizations
// Constructor
TIViscoelastic::TIViscoelastic(char *matName,int matID) : TransIsotropic(matName,matID)
{
    GT0 = -1.;
    GA0 = -1.;
    KT0 = -1.;
    n0 = -1.;
    ell0 = -2.;
    ntaus = -1;
    whichOne = NO_PROP;
    
    // all assume zero terms (but if kept zero, will need positive moduli)
    ntauGT = ntauGA = ntauKT = ntaun = ntauell = 0;
    GTk = tauGTk = NULL;
    GAk = tauGAk = NULL;
    KTk = tauKTk = NULL;
    nk = taunk = NULL;
    ellk = tauellk = NULL;
    
    // WLF
    Tref = -1.;
    C1base10 = 17.44;
    C2 = 51.6;
	C1 = 0.;
    
    // Moisture parameters
    mref = -1.;
	Cm1 = 0.;
	Cm2 = 0.;
    Cm1base10 = 10.;
    Cm2base10 = 0.025/0.04;

	// Visual Studio wants these initiatlized
	GAe = 0.;
	KTe = 0.;
	GTe = 0.;
	elle = 0.;
	ne = 0.;
	betaAconc = 0.;
	betaTconc = 0.;
	currentPk = 0;
	currentTauk = 0;
	numHistory = 0;
    
    // vertical shifting (total stress only)
    bTemp.clear();
    bTValue.clear();
    bConc.clear();
    bCValue.clear();
}

#pragma mark TIViscoelastic::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *TIViscoelastic::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	// read properties for this material
    if(strcmp(xName,"GT0")==0 || strcmp(xName,"G0")==0)
    {   switchSeries(GT_SERIES,ntauGT);
        input = DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&GT0,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"GA0")==0)
    {   switchSeries(GA_SERIES,ntauGA);
        input = DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&GA0,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"KT0")==0)
    {   switchSeries(KT_SERIES,ntauKT);
        input = DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&KT0,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"en0")==0)
    {   switchSeries(N_SERIES,ntaun);
        input = DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&n0,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"ell0")==0)
    {   switchSeries(L_SERIES,ntauell);
        input = DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&ell0,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"ntaus")==0)
    {   // must have one after each "0" property
        input=INT_NUM;
        return((char *)&ntaus);
    }
    
    else if(strcmp(xName,"Pk")==0)
    {   input = DOUBLE_NUM;
        switch(whichOne)
        {   case GT_SERIES:
                return addToPkSeries(&GTk,ntauGT,gScaling);
             case GA_SERIES:
                return addToPkSeries(&GAk,ntauGA,gScaling);
            case KT_SERIES:
                return addToPkSeries(&KTk,ntauKT,gScaling);
            case N_SERIES:
                return addToPkSeries(&nk,ntaun,gScaling);
            case L_SERIES:
                return addToPkSeries(&ellk,ntauell,gScaling);
            default:
                ThrowSAXException("Pk property found when not defining an exponential series");
                break;
        }
    }
    
    else if(strcmp(xName,"tauk")==0)
    {   input = DOUBLE_NUM;
        switch(whichOne)
        {   case GT_SERIES:
                return addToTaukSeries(&tauGTk,ntauGT);
             case GA_SERIES:
                return addToTaukSeries(&tauGAk,ntauGA);
             case KT_SERIES:
                return addToTaukSeries(&tauKTk,ntauKT);
            case N_SERIES:
                return addToTaukSeries(&taunk,ntaun);
             case L_SERIES:
                return addToTaukSeries(&tauellk,ntauell);
             default:
                ThrowSAXException("tauk property found when not defining an exponential series");
                break;
        }
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
        return (char *)&bConc[(int)bConc.size()-1];
    }
    
    else if(strcmp(xName,"bTValue")==0)
    {   bTValue.push_back(0.);
        input=DOUBLE_NUM;
        return (char *)&bTValue[(int)bTValue.size()-1];
    }
    
    else if(strcmp(xName,"bCValue")==0)
    {   bCValue.push_back(0.);
        input=DOUBLE_NUM;
        return (char *)&bCValue[(int)bCValue.size()-1];
    }
    
    return TransIsotropic::InputMaterialProperty(xName,input,gScaling);
}

// Start entering anew series, but verify previous series was completed
void TIViscoelastic::switchSeries(int newSeries,int priorNtaus)
{   // verify current one was done
    switch(whichOne)
    {   case GT_SERIES:
            if(currentPk<ntauGT || currentTauk<ntauGT)
                ThrowSAXException("The GT(t) series did not provide all Pk and tauk values");
            break;
        case GA_SERIES:
            if(currentPk<ntauGA || currentTauk<ntauGA)
                ThrowSAXException("The GA(t) series did not provide all Pk and tauk values");
            break;
        case KT_SERIES:
            if(currentPk<ntauKT || currentTauk<ntauKT)
                ThrowSAXException("The KT(t) series did not provide all Pk and tauk values");
            break;
        case N_SERIES:
            if(currentPk<ntaun || currentTauk<ntaun)
                ThrowSAXException("The n(t) series did not provide all Pk and tauk values");
            break;
        case L_SERIES:
            if(currentPk<ntauell || currentTauk<ntauell)
                ThrowSAXException("The l(t) series did not provide all Pk and tauk values");
            break;
        default:
            // no problem if no previous serties
            break;
    }
    
    // final check passes no property
    if(newSeries==NO_PROP) return;
    
    // each series only allowed once
    if(priorNtaus>0)
        ThrowSAXException("Only one series of exponential allowed for each property");

    // set to new series and reset ntaus
    whichOne = newSeries;
    ntaus = -1;
    currentPk = 0;
    currentTauk = 0;
}

// Start entering anew series, but verify previous series was completed
char *TIViscoelastic::addToPkSeries(double **Pk,int &ntauP,double &gScaling)
{
    if(*Pk==NULL)
    {   if(ntaus<=0)
            ThrowSAXException("Pk found before number of taus specified.");
        *Pk=new double[ntaus];
        ntauP = ntaus;
    }
    currentPk++;
    if(currentPk>ntauP)
        ThrowSAXException("Too many Pk's given.");
    return UnitsController::ScaledPtr((char *)(*Pk+currentPk-1),gScaling,1.e6);
}

// Start entering anew series, but verify previous series was completed
char *TIViscoelastic::addToTaukSeries(double **tauk,int &ntauP)
{
    if(*tauk==NULL)
    {   if(ntaus<=0)
            ThrowSAXException("tauk found before number of taus specified.");
        *tauk=new double[ntaus];
        ntauP = ntaus;
    }
    currentTauk++;
    if(currentTauk>ntauP)
        ThrowSAXException("Too many tauk's given.");
    return (char *)(*tauk+currentTauk-1);
}

// Verify input properties do calculations; if problem return string with an error message
// If OK, MUST pass on to super class. This is called just before PrintMaterial
// (see also ValidateForUse() for checks that depend on MPM calculation mode)
const char *TIViscoelastic::VerifyAndLoadProperties(int np)
{
    // 3D requires axial z, 2D can be z or y
    if(np==THREED_MPM)
    {   // in 3D, require TRANSISO1 and cannot swap
        if(materialID==TIVISCOELASTIC2)
            return "3D simulations cannot use TIViscoelastic 2 material";
        if(swapz>0)
            return "3D simulations cannot swap axial direction, use rotation methods instead";
    }
    else if(materialID==TIVISCOELASTIC2)
    {   if(swapz>0)
            return "Cannot swap axial direction of TIViscoelastic 2 material";
        axialCode = AXIAL_Y;
    }
    else if(swapz>0)
        axialCode = AXIAL_Y;
    
   // Make sure last series was finished (may not like tha error throw)
    switch(whichOne)
    {   case GT_SERIES:
            if(currentPk<ntauGT || currentTauk<ntauGT)
                return "The GT(t) series did not provide all Pk and tauk values";
            break;
        case GA_SERIES:
            if(currentPk<ntauGA || currentTauk<ntauGA)
                return "The GA(t) series did not provide all Pk and tauk values";
            break;
        case KT_SERIES:
            if(currentPk<ntauKT || currentTauk<ntauKT)
                return "The KT(t) series did not provide all Pk and tauk values";
            break;
        case N_SERIES:
            if(currentPk<ntaun || currentTauk<ntaun)
                return "The n(t) series did not provide all Pk and tauk values";
            break;
        case L_SERIES:
            if(currentPk<ntauell || currentTauk<ntauell)
                return "The l(t) series did not provide all Pk and tauk values";
            break;
        default:
            // no problem if no previous serties
            break;
    }
    
    // requires large rotation mode - because no option to rotate properties
    if(useLargeRotation==0) useLargeRotation = 1;
    
    // need five properties
    if(GT0<0.) return "GT0 was not provided (or is negative)";
    if(GA0<0.) return "GA0 was not provided (or is negative)";
    if(KT0<0.) return "KT0 was not provided (or is negative)";
    if(n0<0.) return "n0 was not provided (or is negative)";
    if(ell0<0.) return "ell0 was not provided (or is negative)";

    // zero time moduli and reduced modulus
    int k;
    GTe = GT0;
    for(k=0;k<ntauGT;k++)
    {   GTe += GTk[k];
        GTk[k] /= rho;
    }
    GTe /= rho;
    
    GAe = GA0;
    for(k=0;k<ntauGA;k++)
    {   GAe += GAk[k];
        GAk[k] /= rho;
    }
    GAe /= rho;
    
    KTe = KT0;
    for(k=0;k<ntauKT;k++)
    {   KTe += KTk[k];
        KTk[k] /= rho;
    }
    KTe /= rho;
    
    ne = n0;
    for(k=0;k<ntaun;k++)
    {   ne += nk[k];
        nk[k] /= rho;
    }
    ne /= rho;
    
    elle = ell0;
    for(k=0;k<ntauell;k++)
    {   elle += ellk[k];
        ellk[k] /= rho;
    }
    elle /= rho;
    
    // remaining zero time properties
    nuA = elle/(2.*KTe);
    EA = ne - 4.*KTe*nuA*nuA;
    ET = 1./(1./(4.*KTe) + 1./(4.*GTe) + nuA*nuA/EA);
    nuT = ET/(2.*GTe) - 1.;
    nuAp = ET*nuA/EA;

    // validate Poisson's ratio
    if(nuT<-1. || nuT>1.)
        return "Evaluated value of nuT must be between -1 and 1";
    double nulim = sqrt(EA/ET);
    if(nuA<-nulim || nuA>nulim)
        return "Evaluated value of nuA out of valid range";
    if(ET*nuA*nuA/EA>0.5*(1.-nuT))
        return "Incompatible values of nuA and nuT";

    // convert to analysis units
    aT *= 1.e-6;
    aA *= 1.e-6;
    betaTconc = betaT*concSaturation;
    betaAconc = betaA*concSaturation;

    // diffusion CT is 1
    diffusionCT = 1.;

    // make conductivity (input as (N/(sec-K)) specific (N mm^3/(sec-K-g))
    kCondA /= rho;
    kCondT /= rho;
    
    // always 2 KT (xx and yy or xx and zz), 1 n (yy or zz), 3 ell (xx, yy, zz)
    numHistory = 2*ntauKT + ntaun + 3*ntauell;
    
    // number of shear depends on analysis mode
    if(fmobj->IsThreeD())
    {   // 3 GT (xx, yy, xy) and 2 GA (xz,yz) where z is always axial direction
        numHistory += 3*ntauGT + 2*ntauGA;
    }
    else if(AxialDirection()==AXIAL_Z)
    {   // 2D axial in z direction (3 GT (xx, yy, xy), GA not used where z is axial direction
        numHistory += 3*ntauGT;
    }
    else
    {   // 2D axial in y direction (2 GT (xx, zz) and 1 GA (xy) ) where y is axial direction
        numHistory += 2*ntauGT + ntauGA;
    }

    // WLF coefficients convert to ln aT
    C1 = log(10.)*C1base10;

    // Moisture terms input as log ac = -Cm1base(m-mref)/(Cm2base+m-mref)
    // ... convert to use ln ac = - Cm1*(c-cref)/(Cm2+c) where c = m/csat and cref=mref/csat
    // ... Cm1 = Cm1base*ln(10) and Cm2 = (Cm2base-mref)/csat
    // log ac =
    Cm1 = log(10.)*Cm1base10;
    Cm2 = (Cm2base10-mref)/concSaturation;
    mref /= concSaturation;
    if(mref>=0. && Cm2<=0.)
        return "Cm2 must be greater than mref";
	if(mref>1.)
		return "mref must less than or equal to csat";

#ifdef TOTAL_STRESS_CALC_TI
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

    // superclass call (fills transport property matrices
    return MaterialBase::VerifyAndLoadProperties(np);
}

// Printing the material. The MaterialBase class prints material name and then
// Calls PrintMechanicalProperties(), PrintCommonProperties(), and PrintTransportProperties()
// MaterialBase handles common properties and transport for isotropic materials
void TIViscoelastic::PrintMechanicalProperties(void) const
{
    PrintProperty("EA",rho*EA*UnitsController::Scaling(1.e-6),"");
    PrintProperty("ET",rho*ET*UnitsController::Scaling(1.e-6),"");
    PrintProperty("vA",nuA,"");
    PrintProperty("vT",nuT,"");
    cout << endl;
    
    PrintProperty("GT",rho*GTe*UnitsController::Scaling(1.e-6),"");
    PrintProperty("GA",rho*GAe*UnitsController::Scaling(1.e-6),"");
    PrintProperty("vA'",nuAp,"");
    PrintProperty("KT",rho*KTe*UnitsController::Scaling(1.e-6),"");
    cout << endl;
    
    PrintProperty("n",rho*ne*UnitsController::Scaling(1.e-6),"");
    PrintProperty("ell",rho*elle*UnitsController::Scaling(1.e-6),"");
    double k0 = 1./(1./KTe + (1-2.*nuA)*(1.-2.*nuA)/EA);
    PrintProperty("K",rho*k0*UnitsController::Scaling(1.e-6),"");
    cout << endl;
    
    // print all active series
    PrintPronySeries("GT",GT0,ntauGT,GTk,tauGTk);
    PrintPronySeries("GA",GA0,ntauGA,GAk,tauGAk);
    PrintPronySeries("KT",KT0,ntauKT,KTk,tauKTk);
    PrintPronySeries("n",n0,ntaun,nk,taunk);
    PrintPronySeries("ell",ell0,ntauell,ellk,tauellk);

#ifdef TOTAL_STRESS_CALC_TI
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
        cout << "Isothermal and isosolvent viscoelastic relaxations" << endl;
    else
        cout << "Isothermal viscoelastic relaxations" << endl;

    // WLF moisture properties
    if(mref>=0.)
    {   PrintProperty("mref",mref*concSaturation,"");
        PrintProperty("Cm1",Cm1base10,"");
        PrintProperty("Cm2",Cm2base10,"");
        cout << endl;
    }
    else if(Tref>=0.)
        cout << "Isosolvent viscoelastic relaxations" << endl;

#ifdef TOTAL_STRESS_CALC_TI
    if(bTemp.size()==0)
        cout << "No vertical shifting for temperature" << endl;
    else
    {   cout << "Vertical themal shifting:" << endl;
        for(int i=0;i<bTemp.size();i++)
        {   PrintProperty("  T",bTemp[i],"K");
            PrintProperty("  bT",bTValue[i],"");
        }
        cout << endl;
    }
    if(bConc.size()==0 || diffusion==NULL)
        cout << "No vertical shifting for concentration" << endl;
    else
    {   cout << "Vertical concentration shifting:" << endl;
        for(int i=0;i<bConc.size();i++)
        {   PrintProperty("  c",bConc[i]*concSaturation,"K");
            PrintProperty("  bc",bCValue[i],"");
			cout << endl;
        }
    }
#endif
    PrintProperty("aA",aA*1.e6,"");
    PrintProperty("aT",aT*1.e6,"");
    cout << endl;
}

// Print on prony series
void TIViscoelastic::PrintPronySeries(const char *base,double zeroVal,double ntau,double *Pk,double *tk) const
{   // don't skip zero term series
    if(ntau==0) return;
    
    int i;
    char pname[10];
    
    // Long-term value
    strcpy(pname,base);
    strcat(pname,"0");
    PrintProperty(pname,zeroVal*UnitsController::Scaling(1.e-6),"");
    PrintProperty("ntaus",(double)ntau,"");
    cout <<  endl;
    
    // each term in serties
    strcpy(pname,base);
    strcat(pname,"k");
    for(i=0;i<ntau;i++)
	{	double pnum = (double)i + 1.;
		PrintProperty("  i",pnum,"");
        PrintProperty(pname,rho*Pk[i]*UnitsController::Scaling(1.e-6),"");
        PrintProperty("tauk",tk[i],UnitsController::Label(TIME_UNITS));
        cout << endl;
    }
}

#pragma mark Viscoelastic::History Data Methods

// create and return pointer to history variables
// initialize all to zero
// throws std::bad_alloc
char *TIViscoelastic::InitHistoryData(char *pchr,MPMBase *mptr)
{   // done if none
    if(numHistory==0) return NULL;
    
    // all zeros
    double *p = CreateAndZeroDoubles(pchr,numHistory);
    return (char *)p;
}

// Number of history variables - only the plastic law
int TIViscoelastic::NumberOfHistoryDoubles(void) const { return numHistory; }

#pragma mark TIViscoelastic:Step Methods

// Apply Constitutive law, check np to know what type
void TIViscoelastic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{
    // always in large rotation mode
    // Note: cannot call generic method because need Umat from decomposition
    
    // get current deformation gradient at the start of this time step
    Matrix3 pFnm1 = mptr->GetDeformationGradientMatrix();
    
    // get incremental deformation gradient and decompose it
    const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
    
    // Update total deformation gradient
    // Both plane strain and plane stress will have dezz=0 here, so no change in Fzz
    // ... plane strews adds calculated dezz later
    Matrix3 pF = dF*pFnm1;
    //mptr->SetDeformationGradientMatrix(pF);       // done at the end instead of hear
    
    // two decompositions to get rotation to initial config: Rn-1 and Rn
    // Umat is strength at start of the time step in initial config
    Matrix3 Rnm1,Rn;
    Matrix3 Umat = pFnm1.RightDecompose(&Rnm1,NULL);
    pF.LeftDecompose(&Rn,NULL);

    // Rn = dR*Rnm1 or dR = Rn*Rnm1^T
    Matrix3 dR = Rn*Rnm1.Transpose();
    
    // Material configuration is rotation R0 from initial configuration
    // get strain and strain increments in material configuration using R0^T (initial) R0
    //    where (initial) de = Rn^T(dF-dR)F(n-1) and Umat above
    Matrix3 dFmdR = dF - dR;
    Matrix3 de = Rn.Transpose()*(dFmdR*pFnm1);
    
    // Rotate initial de and Umat to material configuration
    Matrix3 R0 = mptr->GetInitialRotation();
    de = de.RTMR(R0);
    Umat = Umat.RTMR(R0);
    
    // Full rotation from material to current configuration
    Matrix3 Rtot = Rn*R0;
    
#ifdef TOTAL_STRESS_CALC_TI
    // shift of elastic modulus based on mptr->pPreviousTemperature and mptr->pPreviousConcentration
    double bshift = Viscoelastic::GetVertialShift(mptr,bTemp,bTValue,bConc,bCValue);
#else
    // incremental method does not support elastic modulus variations (or vertical shifting)
    double bshift=1.;
#endif
    
    // residual strains (thermal and moisture) in material axes
    double eTr,eAr,eTstretch,eAstretch;
    eTr = aT*res->dT;
    eAr = aA*res->dT;
    double dTemp=mptr->pPreviousTemperature-thermal.reference;
    eTstretch = aT*dTemp;
    eAstretch = aA*dTemp;
    if(fmobj->HasDiffusion())
    {   eTr += betaTconc*res->dC;
        eAr += betaAconc*res->dC;
		double dConc = diffusion->GetDeltaConcentration(mptr);
        eTstretch += betaTconc*dConc;
        eAstretch += betaAconc*dConc;
    }
    
    // history pointer and address
    double *ak = (double *)mptr->GetHistoryPtr(0);
    
    // plane stress needs to remember some dazz
    double *dazzHold = NULL;
    if(np==PLANE_STRESS_MPM)
    {   int nsave = AxialDirection()==AXIAL_Z ? ntaun+3*ntauell : 2*ntauGT+2*ntauKT+3*ntauell ;
        if(nsave>0) dazzHold = new double[nsave];
    }

    // terms for different axes and other variables
    double arg,omtmp,omtmp2,*sT,*sA,dispEnergy=0.,dTq0=0.;
    int k,eT,eA;
    // For TOTAL_STRESS_CALC_TI, dsig will be total stress, otherwise it is stress increment
    Tensor dsig;
    if(AxialDirection()==AXIAL_Z)
    {   eT = 1;         // y direction in MAS
        sT = &dsig.yy;
        eA = 2;         // z direction in MAS
        sA = &dsig.zz;
    }
    else
    {   eT = 2;         // z direction in MAS
        sT = &dsig.zz;
        eA = 1;         // y direction in MAS
        sA = &dsig.yy;
    }

    // effective strain increment
    de(0,0) -= eTr;
    de(eT,eT) -= eTr;
    de(eA,eA) -= eAr;
    
    // Umat is strain at start of the time step
    Umat(0,0) -= (1.+eTstretch);
    Umat(eT,eT) -= (1.+eTstretch);
    Umat(eA,eA) -= (1.+eAstretch);

#ifdef TOTAL_STRESS_CALC_TI
    // elastic terms in each stress component (scaled by bshift when done)
    Matrix3 etot = Umat+de;
    dsig.xx = (KTe+GTe)*etot(0,0) + (KTe-GTe)*etot(eT,eT) + elle*etot(eA,eA);
    *sT =     (KTe-GTe)*etot(0,0) + (KTe+GTe)*etot(eT,eT) + elle*etot(eA,eA);
    *sA =       elle*etot(0,0)    +     elle*etot(eT,eT)  +  ne*etot(eA,eA);
    if(np==THREED_MPM)
    {   // axial direction always in z direction
        dsig.xy = 2.*GTe*etot(0,1);
        dsig.xz = 2.*GAe*etot(0,2);
        dsig.yz = 2.*GAe*etot(1,2);
    }
    else
    {   // 2D
        dsig.xy = AxialDirection()==AXIAL_Z ? 2.*GTe*etot(0,1) : 2.*GAe*etot(0,1) ;
        dsig.xz = dsig.yz = 0.;
    }
    ScaleTensor(&dsig,bshift);
#else
    // increment with zero-time modulis
    dsig.xx = (KTe+GTe)*de(0,0) + (KTe-GTe)*de(eT,eT) + elle*de(eA,eA);
    *sT =     (KTe-GTe)*de(0,0) + (KTe+GTe)*de(eT,eT) + elle*de(eA,eA);
    *sA =       elle*de(0,0)    +     elle*de(eT,eT)  +  ne*de(eA,eA);
    if(np==THREED_MPM)
    {   // axial direction always in z direction
        dsig.xy = 2.*GTe*de(0,1);
        dsig.xz = 2.*GAe*de(0,2);
        dsig.yz = 2.*GAe*de(1,2);
    }
    else
    {   // 2D
        dsig.xy = AxialDirection()==AXIAL_Z ? 2.*GTe*de(0,1) : 2.*GAe*de(0,1) ;
        dsig.xz = dsig.yz = 0.;
    }
#endif
    
    // get effective time increment
    double delEffTime = Viscoelastic::GetEffectiveIncrement(mptr,res,delTime,Tref,C1,C2,mref,Cm1,Cm2,0.,1.);
    
    // subtract each serires
    int ai=0,asave=0;
    
#pragma mark GT Series
    // 3 GT (xT, TT, xy) for 3D or for AXIAL_Z, (xx and TT) for AXIAL_Y (in 2D)
    for(k=0;k<ntauGT;k++)
    {   GetAlphaArgs(delEffTime,tauGTk[k],omtmp,omtmp2);
        // xx, yy(T), xy(T) (2X for shear), but xx, zz(T) for AXIAL_Y (always 2D)
        double daxx = omtmp*(Umat(0,0)-ak[ai]) + omtmp2*de(0,0);
        double daT = omtmp*(Umat(eT,eT)-ak[ai+1]) + omtmp2*de(eT,eT);
        
#ifdef TOTAL_STRESS_CALC_TI
        // add to stress
        arg = bshift*GTk[k]*(ak[ai]+daxx-ak[ai+1]-daT);
#else
        // add to stress increments
        arg = GTk[k]*(daxx-daT);
#endif
        dsig.xx -= arg;
        *sT += arg;
        
        // 2X for shear terms
        double daxy = 0.;
        if(AxialDirection()==AXIAL_Z)
        {   // includes 3D, which is always axial z
            daxy = omtmp*(2.*Umat(0,1)-ak[ai+2]) + omtmp2*2.*de(0,1);
#ifdef TOTAL_STRESS_CALC_TI
            // add to stress
            dsig.xy -= bshift*GTk[k]*(ak[ai+2]+daxy);
#else
            // add to stress increments
            dsig.xy -= GTk[k]*daxy;
#endif
        }
         
        // update history
        ak[ai] += daxx;
        ak[ai+1] += daT;

        // dissipated energy
        if(np!=PLANE_STRESS_MPM || AxialDirection()==AXIAL_Z)
        {   // (exx-axx) - (eT-aT) where T is always y here
            arg = Umat(0,0)+de(0,0)-Umat(eT,eT)-de(eT,eT)-ak[ai]+ak[ai+1];
            dispEnergy += bshift*GTk[k]*( arg*(daxx-daT) );
            if(AxialDirection()==AXIAL_Z)
            {   // includes 3D, which is always axial z and eT=1 here
                ak[ai+2] += daxy;
                dispEnergy += bshift*GTk[k]*(2.*(Umat(0,eT)+de(0,eT))-ak[ai+2])*daxy;
                ai++;
            }
        }
        else
        {   // only save for PLANE_STRESS_MPM and AXIAL_Y
            dazzHold[asave++] = daxx;
            dazzHold[asave++] = daT;
        }
         
        // next history variable (shear one added above when used)
        ai += 2 ;
    }

#pragma mark GA Series
    //  2 GA (xA=xz,TA=yz) for 3D, none for 2D/AXIAL_Z, and 1 GA (xy) for 2D/AXIAL_Y
    if(np==THREED_MPM)
    {   for(k=0;k<ntauGA;k++)
        {   GetAlphaArgs(delEffTime,tauGAk[k],omtmp,omtmp2);
            // xz amd yz (2X for shear)
            double daxz = omtmp*(2.*Umat(0,2)-ak[ai]) + 2.*omtmp2*de(0,2);
            double dayz = omtmp*(2.*Umat(1,2)-ak[ai+1]) + 2.*omtmp2*de(1,2);
            
			// update history
			ak[ai] += daxz;
			ak[ai+1] += dayz;
			
#ifdef TOTAL_STRESS_CALC_TI
            // add to stress
            dsig.xz -= bshift*GAk[k]*ak[ai];
            dsig.yz -= bshift*GAk[k]*ak[ai+1];
#else
            // add to stress increments
            dsig.xz -= GAk[k]*daxz;
            dsig.yz -= GAk[k]*dayz;
#endif
                
            // dissipated energy
            dispEnergy += bshift*GAk[k]*( (2.*(Umat(0,2)+de(0,2))-ak[ai])*daxz
                                  + (2.*(Umat(1,2)+de(1,2))-ak[ai+1])*dayz );
            
            // next history variable
            ai += 2;
        }
    }
    else if(AxialDirection()==AXIAL_Y)
    {   // 1 GA (xy) for AXIAL_Y and always 2D
        for(k=0;k<ntauGA;k++)
        {   GetAlphaArgs(delEffTime,tauGAk[k],omtmp,omtmp2);
            // xy (2X for shear)
            double daxy = omtmp*(2.*Umat(0,1)-ak[ai]) + omtmp2*2.*de(0,1);
            
			// update history
			ak[ai] += daxy;
			
#ifdef TOTAL_STRESS_CALC_TI
            // add to stress
            dsig.xy -= bshift*GAk[k]*ak[ai];
#else
            // add to stress increments
            dsig.xy -= GAk[k]*daxy;
#endif
                     
            // dissipated energy - here eA=1 is y direction
            dispEnergy += bshift*GAk[k]*(2.*(Umat(0,eA)+de(0,eA))-ak[ai])*daxy;
            
            // next history variable
            ai++;
        }
    }
   
#pragma mark KT Series
    // 2 for KT series
    for(k=0;k<ntauKT;k++)
    {   // x and T
        GetAlphaArgs(delEffTime,tauKTk[k],omtmp,omtmp2);
        double daxx = omtmp*(Umat(0,0)-ak[ai]) + omtmp2*de(0,0);
        double daT = omtmp*(Umat(eT,eT)-ak[ai+1]) + omtmp2*de(eT,eT);

#ifdef TOTAL_STRESS_CALC_TI
        // add to stress
        arg = bshift*KTk[k]*(ak[ai]+daxx+ak[ai+1]+daT);
#else
        // add to stress increments
        arg = KTk[k]*(daxx+daT);
#endif
        dsig.xx -= arg;
        *sT -= arg;
        
        // update history
        ak[ai] += daxx;
        ak[ai+1] += daT;

        // dissipated energy
        if(np!=PLANE_STRESS_MPM || AxialDirection()==AXIAL_Z)
        {	// Here eT=1 is y direction or eT=2 for 2D/AXIAL_Y/Plane strain
			arg = Umat(0,0)+de(0,0)+Umat(eT,eT)+de(eT,eT)-ak[ai]-ak[ai+1];
            dispEnergy += bshift*KTk[k]*arg*(daxx+daT);
        }
        else
        {   // only save for PLANE_STRESS_MPM and AXIAL_Y
            dazzHold[asave++] = daxx;
            dazzHold[asave++] = daT;
        }

        // next history variable
        ai += 2;
    }

#pragma mark n Series
    // 1 for n series (axial direction only)
    for(k=0;k<ntaun;k++)
    {   GetAlphaArgs(delEffTime,taunk[k],omtmp,omtmp2);
        double daA = omtmp*(Umat(eA,eA)-ak[ai]) + omtmp2*de(eA,eA);
        
		// update history
		ak[ai] += daA;
		
#ifdef TOTAL_STRESS_CALC_TI
        // add to stress
        *sA -= bshift*nk[k]*ak[ai];
#else
        // add to stress increments
        *sA -= nk[k]*daA;
#endif
             
        // dissipated energy
        if(np!=PLANE_STRESS_MPM || AxialDirection()==AXIAL_Y)
		{	// Here eA=2 unless axial y then eA=1
			dispEnergy += bshift*nk[k]*(Umat(eA,eA)+de(eA,eA)-ak[ai])*daA;
		}
        else
        {   // only save for PLANE_STRESS_MPM and AXIAL_Z
            dazzHold[asave++] = daA;
        }
        
        // next history variable
        ai++;
    }
        
#pragma mark ell Series
    // 3 for ell series (xx, yy, zz)
    for(k=0;k<ntauell;k++)
    {   GetAlphaArgs(delEffTime,tauellk[k],omtmp,omtmp2);
        // xx, yy, and zz
        double daxx = omtmp*(Umat(0,0)-ak[ai]) + omtmp2*de(0,0);
        double daT = omtmp*(Umat(eT,eT)-ak[ai+1]) + omtmp2*de(eT,eT);
        double daA = omtmp*(Umat(eA,eA)-ak[ai+2]) + omtmp2*de(eA,eA);
        
		// update history
		ak[ai] += daxx;
		ak[ai+1] += daT;
		ak[ai+2] += daA;
		
#ifdef TOTAL_STRESS_CALC_TI
        // add to stress increments
        dsig.xx -= bshift*ellk[k]*ak[ai+2];
        *sT -= bshift*ellk[k]*ak[ai+2];
        *sA -= bshift*ellk[k]*(ak[ai]+ak[ai+1]);
#else
        // add to stress increments
        dsig.xx -= ellk[k]*daA;
        *sT -= ellk[k]*daA;
        *sA -= ellk[k]*(daxx+daT);
#endif

        // dissipated energy
        if(np!=PLANE_STRESS_MPM)
        {   dispEnergy += bshift*ellk[k]*( (Umat(eA,eA)+de(eA,eA)-ak[ai+2])*(daxx+daT)
                            + (Umat(0,0)+de(0,0)+Umat(eT,eT)+de(eT,eT)-ak[ai]-ak[ai+1])*daA );
        }
        else
        {   // only save for PLANE_STRESS_MPM
            dazzHold[asave++] = daxx;
            dazzHold[asave++] = daT;
            dazzHold[asave++] = daA;
        }
        
        // next history variable
        ai += 3;
    }
    
#pragma mark Plane stress corrections
    //  plane stress corrections
    if(np==PLANE_STRESS_MPM)
    {   double dezz;
        
        if(AxialDirection()==AXIAL_Z)
        {   double phi = bshift*ne;
            for(k=0;k<ntaun;k++) phi -= bshift*nk[k]*GetPSArg(delEffTime,taunk[k]);
            dezz = -dsig.zz/phi;
            pF(2,2) *= (1.+dezz);
            de(2,2) += dezz;
            dsig.zz = 0.;

            // skip GT and KT (no GA in this 2D with A in z direction)
            // For AXIAL_Z, will have saves 1 for n and 3 for ell
            ai = 3*ntauGT + 2*ntauKT ;
            asave = 0;
                
            // n series
            for(k=0;k<ntaun;k++)
            {   // alpha zz has changed
                double daA = dazzHold[asave++];
                double dazz = GetPSArg(delEffTime,taunk[k])*dezz;
                
                // update history
				daA += dazz;
                ak[ai] += dazz;
                     
                // dissipated energy with final daA not added above
                dispEnergy += bshift*nk[k]*(Umat(eA,eA)+de(eA,eA)-ak[ai])*daA;
                
                // next history
                ai++;
            }
                
            // ell series
            dsig.xx += bshift*elle*dezz;
            dsig.yy += bshift*elle*dezz;
            for(k=0;k<ntauell;k++)
            {   // alpha xx, yy, and zz
                double daxx = dazzHold[asave++];
                double daT = dazzHold[asave++];
                double daA = dazzHold[asave++];
                double dazz = GetPSArg(delEffTime,tauellk[k])*dezz;
                
                // stress increments
                dsig.xx -= bshift*ellk[k]*dazz;
                dsig.yy -= bshift*ellk[k]*dazz;
                
                // update history
				daA += dazz;
				ak[ai+2] += dazz;

                // dissipated energy with final daA not added above
                dispEnergy += bshift*ellk[k]*( (Umat(eA,eA)+de(eA,eA)-ak[ai+2])*(daxx+daT)
                                    + (Umat(0,0)+de(0,0)+Umat(eT,eT)+de(eT,eT)-ak[ai]-ak[ai+1])*daA );
                
                // next history
                ai += 3;
            }
        }
        else
        {   // axial in y direction, z is a transverse direction
			// In MAS, x is x, A is y, and T is z
            double phiG = bshift*GTe,phiK = bshift*KTe;
            for(k=0;k<ntauGT;k++) phiG -= bshift*GTk[k]*GetPSArg(delEffTime,tauGTk[k]);
            for(k=0;k<ntauKT;k++) phiK -= bshift*KTk[k]*GetPSArg(delEffTime,tauKTk[k]);
            dezz = -dsig.zz/(phiG+phiK);
            pF(2,2) *= (1.+dezz);
            de(2,2) += dezz;
            
            // stress increments
            dsig.xx += (phiK-phiG)*dezz;
            dsig.zz = 0.;
            
            // history indexing
            ai = asave = 0;
            
            // GT series (xx and zz) for AXIAL_Y
            // Here T is z in MAS (a transverse direction). Update it with dezz
            for(k=0;k<ntauGT;k++)
            {   double daxx = dazzHold[asave++];
                double daT = dazzHold[asave++];
                double dazz = GetPSArg(delEffTime,tauGTk[k])*dezz;
                 
                // update history for daT
				daT += dazz;
                ak[ai+1] += dazz;
                    
                // dissipated energy, eT=2 is z in MAS
                arg = Umat(0,0)+de(0,0)-Umat(eT,eT)-de(eT,eT)-ak[ai]+ak[ai+1];
                dispEnergy += bshift*GTk[k]*( arg*(daxx-daT) );
                  
                // next history variable
                ai += 2 ;
            }
            
            // skip GA
            ai += ntauGA;
            
            // 2 for KT series
			// Here T is z in MAS (a transverse direction). Update it with dezz
            for(k=0;k<ntauKT;k++)
            {   // xx and T
                double daxx = dazzHold[asave++];
                double daT = dazzHold[asave++];
                double dazz = GetPSArg(delEffTime,tauKTk[k])*dezz;

                // update history
				daT += dazz;
                ak[ai+1] += dazz;

                // dissipated energy, eT=2 is z in MAS
				arg = Umat(0,0)+de(0,0)+Umat(eT,eT)+de(eT,eT)-ak[ai]-ak[ai+1];
                dispEnergy += bshift*KTk[k]*arg*(daxx+daT);

                // next history variable
                ai += 2;
            }
            
            // skip n
            ai += ntaun;
            
            // ell series
			// Here T is z in MAS (a transverse direction). Update it with dezz
            dsig.yy += bshift*elle*dezz;
            for(k=0;k<ntauell;k++)
            {   // xx, yy, and zz
                double daxx = dazzHold[asave++];
                double daT = dazzHold[asave++];
                double daA = dazzHold[asave++];
                double dazz = GetPSArg(delEffTime,tauellk[k])*dezz;
                
                // stresses
                dsig.yy -= bshift*ellk[k]*dazz;
                
                // update history
				daT += dazz;
                ak[ai+1] += dazz;

                // dissipated energy
                dispEnergy += bshift*ellk[k]*( (Umat(eA,eA)+de(eA,eA)-ak[ai+2])*(daxx+daT)
                                    + (Umat(0,0)+de(0,0)+Umat(eT,eT)+de(eT,eT)-ak[ai]-ak[ai+1])*daA );
                
                // next history
                ai += 3;
            }
        }
        
        // plane stress
        if(dazzHold!=NULL) delete [] dazzHold;
    }

    // Update particle stresses
    Tensor *sp=mptr->GetStressTensor();
    
    // update stresses
    double workEnergy,resEnergy;
#ifdef TOTAL_STRESS_CALC_TI
    if(np==THREED_MPM)
    {   // energy uses total strain, so add back residual strains
        workEnergy = dsig.xx*(de(0,0)+eTr) + dsig.yy*(de(1,1)+eTr) + dsig.zz*(de(2,2)+eAr)
               + 2.*(dsig.xz*de(0,2) + dsig.yz*de(1,2) + dsig.xy*de(0,1));
        resEnergy = dsig.xx*eTr + dsig.yy*eTr + dsig.zz*eAr;
		
		// now get stress in current configuration
		*sp = Rtot.RVoightRT(&dsig,true,false);
    }
    else
    {   // energy uses total strain, so add back residual strains
        workEnergy = dsig.xx*(de(0,0)+eTr) + 2.*dsig.xy*de(0,1);
        resEnergy = dsig.xx*eTr;
        if(AxialDirection()==AXIAL_Z)
        {   workEnergy += dsig.yy*(de(1,1)+eTr) + dsig.zz*(de(2,2)+eAr);
            resEnergy += dsig.yy*eTr + dsig.zz*eAr;
        }
        else
        {   workEnergy += dsig.yy*(de(1,1)+eAr) + dsig.zz*(de(2,2)+eTr);
            resEnergy += dsig.yy*eAr + dsig.zz*eTr;
        }
		
		// now get stress in current configuration
		*sp = Rtot.RVoightRT(&dsig,true,true);
    }
    
#else
    if(np==THREED_MPM)
    {   // rotate increment to current configuration and prior stress by dR
        //... then add to get new stress in current configuration
        Tensor dsigma = Rtot.RVoightRT(&dsig,true,false);
        *sp = dR.RVoightRT(sp, true, false);
        AddTensor(sp, &dsigma);
        
        // rotate back to material axes (where strain is known) for energy calculations)
        // Add residual strain increment back to get total strain increment
        Tensor matSig = Rtot.RTVoightR(sp,true,false);
        workEnergy = matSig.xx*(de(0,0)+eTr) + matSig.yy*(de(1,1)+eTr) + matSig.zz*(de(2,2)+eAr)
               + 2.*(matSig.xz*de(0,2) + matSig.yz*de(1,2) + matSig.xy*de(0,1));
        resEnergy = matSig.xx*eTr + matSig.yy*eTr + matSig.zz*eAr;
    }
    else
    {   // rotate increment to current configuration and prior stress by dR
        //... then add to get new stress in current configuration
        Tensor dsigma = Rtot.RVoightRT(&dsig,true,true);
        *sp = dR.RVoightRT(sp, true, true);
        AddTensor(sp, &dsigma);
        
        // rotate back to material axes (where strain is known) for energy calculations)
        // Add residual strain increment back to get total strain increment
        Tensor matSig = Rtot.RTVoightR(sp,true,true);
        workEnergy = matSig.xx*(de(0,0)+eTr) + 2.*matSig.xy*de(0,1);
        resEnergy = matSig.xx*eTr;
        if(AxialDirection()==AXIAL_Z)
        {   workEnergy += matSig.yy*(de(1,1)+eTr) + matSig.zz*(de(2,2)+eAr);
            resEnergy += matSig.yy*eTr + matSig.zz*eAr;
        }
        else
        {   workEnergy += matSig.yy*(de(1,1)+eAr) + matSig.zz*(de(2,2)+eTr);
            resEnergy += matSig.yy*eAr + matSig.zz*eTr;
        }
    }
#endif
    
    // incremental work energy
    // Whan it residual energy?
    mptr->AddWorkEnergyAndResidualEnergy(workEnergy,resEnergy);
    
    // finish particle updates
    mptr->SetDeformationGradientMatrix(pF);
    
    // dissipated energy per unit mass (dPhi/(rho0 V0)) (Legacy nJ/g)
    mptr->AddPlastEnergy(dispEnergy);
    
    // heat energy (should get dTq0)
    IncrementHeatEnergy(mptr,dTq0,dispEnergy);
}

// Get parameters to update alpha variables
void TIViscoelastic::GetAlphaArgs(double delTime,double tau,double &omtmp,double &omtmp2) const
{   double x = delTime/tau;
    if(x>0.003)
    {   double tmp = exp(-0.5*x);
        omtmp = 1. - tmp*tmp;           // 1-exp(-dt/tau)
        omtmp2 = 1. - tmp;              // 1-exp(-dt/(2*tau))
    }
    else
    {   // double precision accurate and maybe more precision than 1-exp(-x)
        omtmp = x*(1.-0.5*x*(1-(x/3.)*(1-0.25*x)));         // 1-exp(-dt/tau)
        x *= 0.5;
        omtmp2 = x*(1.-0.5*x*(1-(x/3.)*(1-0.25*x)));        // 1-exp(-dt/(2*tau))
   }
}

// Get parameter for plane stress calculations
// Get 1-exp(-dt/(2*tau))
double TIViscoelastic::GetPSArg(double delTime,double tau) const
{   double x = 0.5*delTime/tau;
    if(x>0.003)
        return 1. - exp(-x);
    
    // double precision accurate and maybe more precision than 1-exp(-x)
    return x*(1.-0.5*x*(1-(x/3.)*(1-0.25*x)));
}

#pragma mark TIViscoelastic::Accessors

// return material type
const char *TIViscoelastic::MaterialType(void) const
{    if(AxialDirection()==AXIAL_Z)
        return "Tranversely isotropic viscoelastic (unrotated A axis along z axis)";
    else
        return "Tranversely isotropic viscoelastic (unrotated A axis along y axis)";
}

/* Calculate maximum wave speed in mm/sec (moduli in MPa, rho in g/mm^3)
    AXIAL_Z and in 2D
        wave speeds are GT/rho and (KT+GT)/rho - return larger one
    All other cases
        wave speeds are GT/rho, GA/rho, (KT+GT)/rho, and (EA + 4KT nuA^2)/rho
        return largest (assumes shear ones are not largest)
*/
#ifdef TOTAL_STRESS_CALC_TI
// Should get shift in KTe and GTe and GAe and ne to do this right
#endif
double TIViscoelastic::WaveSpeed(bool threeD,MPMBase *mptr) const
{
    if(AxialDirection()==AXIAL_Z && !threeD)
        return sqrt(KTe+GTe);
    else
        return sqrt(fmax(GAe,fmax(KTe+GTe,ne)));
}


// not supported yet, need to deal with aniostropi properties
bool TIViscoelastic::SupportsDiffusion(void) const
{   return DiffusionTask::HasPoroelasticity() ? false : true;
}
