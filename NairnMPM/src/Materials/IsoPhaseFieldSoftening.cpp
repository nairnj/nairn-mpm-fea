/********************************************************************************
	IsoPhaseFieldSoftening.hpp
	nairn-mpm-fea
	
	Created by John Nairn, Sep 27, 2021.
	Copyright (c) 202` John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/IsoPhaseFieldSoftening.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "System/UnitsController.hpp"
#include "Custom_Tasks/PhaseFieldDiffusion.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Boundary_Conditions/InitialCondition.hpp"

extern double mtime;

#pragma mark IsoPhaseFieldSoftening::Constructors and Destructors

// Constructor
// throws std::bad_alloc
IsoPhaseFieldSoftening::IsoPhaseFieldSoftening(char *matName,int matID) : IsotropicMat(matName,matID)
{
	viscosity = -1.;
	ell = -1.;
	JIc = -1.;
	gdMode = GEN_QUADRATIC;
	kStability = 0.;
	hasStability = false;
	garg = 1;
	gargPicked = false;
	partition = PRESSURE_SHEAR;
	phaseTask = NULL;
	taskNum = -1;
}

#pragma mark IsoPhaseFieldSoftening::Initialization

// Read material properties
char *IsoPhaseFieldSoftening::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	if(strcmp(xName,"viscosity")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&viscosity,gScaling,1.e-3);
	}
	
	else if(strcmp(xName,"ell")==0)
	{	input=DOUBLE_NUM;
		return (char *)&ell;
	}

	else if(strcmp(xName,"stability")==0)
	{	input=DOUBLE_NUM;
		return (char *)&kStability;
	}
	
	else if(strcmp(xName,"gd")==0)
	{	input=INT_NUM;
		return (char *)&gdMode;
	}
	
	else if(strcmp(xName,"garg")==0)
	{	input=DOUBLE_NUM;
		gargPicked = true;
		return (char *)&garg;
	}
	
	else if(strcmp(xName,"partition")==0)
	{	input=INT_NUM;
		return (char *)&partition;
	}

	return IsotropicMat::InputMaterialProperty(xName,input,gScaling);
}

// calculate properties used in analyses
const char *IsoPhaseFieldSoftening::VerifyAndLoadProperties(int np)
{
	// required parameters
	if(viscosity<=0. || ell<=0. || JIc<=0)
		return "Phase field material must specify JIc, ell, and viscosity";
	
	if(useLargeRotation!=1)
		return "Phase field material requires large rotation mode";
	
	// stability factor
	if(kStability<0.)
		return "The stability factor must be non-negative";
	if(np==PLANE_STRESS_MPM)
	{	// I am not sure if plane stresses needs stability factor?
		//if(kStability<=0.)
		//	return "Plane stress analysis requires a positive stability factor";
	}
	if(kStability>0.) hasStability = true;
	
	// only two modes allowed
	if(gdMode<GEN_QUADRATIC || gdMode>PF_LINEAR_SOFTENING)
		return "'gd' parameter must be 0, 1, or 2";
	
	// only two partition modes
	if(partition<TENSILE_EIGENVALUES || partition>PRESSURE_SHEAR)
		return "'partition' parameter must be 0 or 1";

    // requires phase field transport task
    phaseTask = (PhaseFieldDiffusion *)DiffusionTask::FindDiffusionTask(FRACTURE_PHASE_FIELD);
    if(phaseTask!=NULL)
        taskNum = phaseTask->GetNumber();
    else
        return "Phase field material requires phase field diffusion task to be active";
    
	// basic isotropic stuff
	const char *msg = IsotropicMat::VerifyAndLoadProperties(np);
	if(msg!=NULL) return msg;
	
	// Gc*ell is diffusion term (put in unrotated properties too)
	ZeroTensor(&pfDTensor);
	pfDTensor.xx = JIc*ell;
	pfDTensor.yy = pfDTensor.xx;
	pfDTensor.zz = pfDTensor.xx;
	
	// pick default value for each gdMode (if was not picked)
	Psii = 0.;
	if(gdMode==GEN_QUADRATIC)
	{	if(!gargPicked)
			garg = 1;
		else if(garg<0 || garg>1)
			return "For quadratic g(d) (gd=0), must set 0 <= garg <= 1";
		else
		{	// input is ductility parameter from 0 to 1
			//	code was written using a = 2*xi-1 from -1 to 1
			// if needed xi = 0.5*(1+a)
			// g(d) = (1-d)(1-a*d)
			// g'(d) = -(1+a-2*a*d)
			garg = 2.*garg - 1.;
		}
	}
	else if(gdMode==FINITE_EXPONENTIAL)
	{	if(!gargPicked)
			garg = 3.;
		else if(garg<=0.)
			return "For exponential g(d) (gd=1), must set garg > 0";
	
		// terms for exponential g(d) function
		double arg = exp(-garg);
		// g(d) = (exp(-k*d)-exp(-k)/(1-exp(-k))
		//		= (ec1*exp(-k*d)+ec2)
		ec1 = 1./(1.-arg);
		ec2 = -arg/(1.-arg);
		// g'(d) = -k exp(-k*d)/(1-exp(-k)) = -ec3*exp(-k*d)
		ec3 = garg*ec1;
	}
	else
	{   // Note change garg (=chi in notes) in revision 3490 to be twice one used in earlier version
        // code is using earlier version. Input garg is ratio of area after peak ot area before the peak
        // in vode garg is half that value (see JANOSE-026-5 and following)
        if(!gargPicked)
            garg = 4;		// give re = 1/5
		else if(garg<=0.)
			return "For linear softening g(d) (gd=2), must set garg > 0";
		
		// initiation energy
        // change to half the area ratio
        garg *= 0.5;
		Psii = JIc/(4.*garg*ell);
	}

	return NULL;
}

// print mechanical properties to the results
void IsoPhaseFieldSoftening::PrintMechanicalProperties(void) const
{
	IsotropicMat::PrintMechanicalProperties();
    
    // basic properties
    PrintProperty("Gc",JIc*UnitsController::Scaling(0.001),UnitsController::Label(ERR_UNITS));
    PrintProperty("ell",ell,UnitsController::Label(CULENGTH_UNITS));
    PrintProperty("eta",viscosity*UnitsController::Scaling(1.e3),UnitsController::Label(VISCOSITY_UNITS));
    cout << "(task #" << taskNum << ")";
    cout << endl;

    // custom softening law (literature is (1-d)^2 or GEN_QUADRATIC with garg=xi=1
	cout << "g(d): ";
	if(gdMode==GEN_QUADRATIC)
		cout << "(1-d)*(1-(2*xi-1)*d) with ductility factor xi = " << (1.+garg)/2. << endl;
	else if (gdMode==FINITE_EXPONENTIAL)
		cout << "(exp(-alpha*d)-exp(-alpha))/(1-exp(-alpha)) with alpha = " << garg << endl;
	else
	{   // See JANOSU-025-5 and following for change in chi to chi' (output here is chi')
        cout << "1+chi*d^2/2-d*sqrt(1+chi+chi^2*d^2/4) with chi = " << (2.*garg) << endl;
        double Ered = fmobj->np==PLANE_STRAIN_MPM ? E/(1-nu*nu) : E ;
		double exi = sqrt(2*Psii/Ered);
		cout << "   re = " << 1/(1+2*garg) << ", ei = " << 100*exi << "%" << endl;
		double sigmax = Ered*exi*UnitsController::Scaling(1.e-6);
		cout << "   sigma(max) = " << sigmax << ", tau(max) = " << sigmax/sqrt(2.*(1+nu)) << endl;
	}
    
    // stability factor
	cout << "Stability factor for (1-k)*g(d)+k with k = " << kStability << endl;
    
    // partitioning method
	cout << "Strain partitioning: ";
	if(partition==TENSILE_EIGENVALUES)
		cout << "tensile eigenstrains" << endl;
	else
		cout << "volumetric/deviatoric strains" << endl;
}

#pragma mark IsoPhaseFieldSoftening::History Data Methods

// store phase field - with diffusion, it is both on pDiff[] and in history variables
// 3 (phase) and 4 (change in phase field)
char *IsoPhaseFieldSoftening::InitHistoryData(char *pchr,MPMBase *mptr)
{	// Check for using diffusion or direct linear solver
	if(taskNum>=0)
	{	mptr->pDiff[taskNum]->conc = 0.;
		mptr->pDiff[taskNum]->prevConc = mptr->pDiff[taskNum]->conc;
	}
	
	// all zeros
	double *p = CreateAndZeroDoubles(pchr,NumberOfHistoryDoubles());
	return (char *)p;
}

// Number of history variables
// 1: Maximum psiPlus value (i.e. tensile energy) (HISTORY_VALUE)
// 2: 0 if not failed, 1 when failed (PHASE_FAILURE_STATE)
// 3: Current phase field value (PHASE_SOLVED)
// 4: Change in phase field on the last time step (DELTA_PHASE_SOLVED)
int IsoPhaseFieldSoftening::NumberOfHistoryDoubles(void) const { return 4; }

// Initialize damage when requested
void IsoPhaseFieldSoftening::SetInitialConditions(InitialCondition *ic,int ptNum,bool is3D)
{
    MPMBase *mptr = mpm[ptNum-1];
	double initPhase = ic->GetPhaseField();
    if(initPhase<0.)
        initPhase = 0.;
    else if(initPhase>1.)
        initPhase = 1.;
	if(taskNum>=0)
	{   mptr->pDiff[taskNum]->conc = initPhase;
		mptr->pDiff[taskNum]->prevConc = initPhase;
	}
	
	// on history variables too
	mptr->SetHistoryDble(PHASE_SOLVED,initPhase,0);
	if(initPhase>=1.) mptr->SetHistoryDble(PHASE_FAILURE_STATE,1.,0);
}

#pragma mark IsoPhaseFieldSoftening::Methods

// Revised code for isotropic material with phase field
// Only plane strain option is implemented
// axisymmetric, 3D possible, or plane stress not updated
// Energy dissipation not calculated = 2(1-phi)C(e-ers)-(e-eres)dPhi
// dV/V should adjust for opening not due to cracks (used for compression heating)
void IsoPhaseFieldSoftening::LRConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// current previous deformation gradient and stretch
	Matrix3 pFnm1 = mptr->GetDeformationGradientMatrix();
	
	// get incremental deformation gradient and decompose it
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
	
	// Update total deformation gradient (saved on particle at the end)
	Matrix3 pF = dF*pFnm1;
	mptr->SetDeformationGradientMatrix(pF);
	
	// decompose to get previous Rn and Rn-1 and current V = I + Biot Strain
	Matrix3 Rnm1,Rn;
	pFnm1.LeftDecompose(&Rnm1,NULL);
	Matrix3 Vn = pF.LeftDecompose(&Rn,NULL);
	Matrix3 dR = Rn*Rnm1.Transpose();
	
	// get strain increments in current configuration (dF-dR)F(n-1)Rn^T
	Matrix3 dFmdR = dF - dR;
	Matrix3 de = dFmdR*(pFnm1*Rn.Transpose());

	// cast pointer to material-specific data
	ElasticProperties *p = (ElasticProperties *)properties;
	double lambda,twoG,C11;
	if(np==THREED_MPM)
	{	twoG = p->C[0][0]-p->C[0][1];
		lambda = p->C[0][1];
		C11 = p->C[0][0];
	}
	else if(np==PLANE_STRESS_MPM)
	{	twoG = 2.*p->C[3][3];
		C11 = p->C[4][4];
		lambda = C11 - twoG;
	}
	else
	{	twoG = p->C[1][1]-p->C[1][2];
		lambda = p->C[1][2];
		C11 = p->C[1][1];
	}
	double Ksp = lambda + twoG/3.;
	
	// save initial stresses
	Tensor *sp = mptr->GetStressTensor();
	
    // residual strains (thermal)
    double eres = CTE1*res->dT;
    double ezzres = CTE3*res->dT;
    
    // initial phase field effective strain  = en - de - alpha(T-Tref) (get last term here)
    double resStrain = CTE3*(mptr->pPreviousTemperature-thermal.reference);
    
    // residual strains (moisture always in pDiff[0])
    if(DiffusionTask::HasFluidTransport())
    {   eres += CME1*res->dC;
        ezzres += CME3*res->dC;
        resStrain += CME3*(mptr->pDiff[0]->prevConc-diffusion->reference);
    }
	
	// Trace terms
	double trenm1,trenm1p,trenm1m,tren,trenp,trenm;
	
	// previous damage strate
	double gdp;
	double gd = GetDegradation(mptr, gdp);
	double dPhase = mptr->GetHistoryDble(DELTA_PHASE_SOLVED,0);
	
	// Effective 3D strains previous step (plane strain enm1(2,2) = -resStrain)
	// Plane stress correct because de(2,2)=0 here and Vn(2,2) is previous z deformation
	Matrix3 enm1 = Vn-de;
	enm1(0,0) -= (1. + resStrain);
	enm1(1,1) -= (1. + resStrain);
	enm1(2,2) -= (1. + resStrain);
	trenm1 = enm1.trace();					// using effective trains
	trenm1p = fmax(trenm1,0.);				// +/- versions
	trenm1m = trenm1 - trenm1p;

	// get +/- initial strains
	// These correct all modes
	Matrix3 enm1p;
	if(partition==TENSILE_EIGENVALUES)
	{	Vector lam = enm1.Eigenvalues();
		Matrix3 Qcol = enm1.Eigenvectors(lam);
		Matrix3 QcolT = Qcol.Transpose();
		Matrix3 Lam(fmax(lam.x,0.),0.,0.,fmax(lam.y,0.),fmax(lam.z,0.));
		Matrix3 LamQcolT = Lam*QcolT;
		enm1p = Qcol*LamQcolT;
	}

	// current phase field effective strain = en - alpha*(T+dT-Tref) = en - eres - resStrain
	// 2D: en(2,2) = ezz - resStrain - eres
	// Plane strain: ezz=0; Plane stress: ezz is previous ezz; needs to find dezz and add to en(2,2)
	Matrix3 en = Vn;
	en(0,0) -= (1. + eres + resStrain);
	en(1,1) -= (1. + eres + resStrain);
	en(2,2) -= (1. + eres + resStrain);
	
	// current trace
	// 2D has en(2,2) = -eres-resStrain (plane strain), en(2,2)=ezz(n-1)-eres-resStrain, and de(2,2) = 0
	// Plane stress needs to find total dezz and add to tren and dVoverV
	tren = en.trace();
	double dVoverV = de.trace();
		
	// Plane stress needs to make dszz = 0
	// for TENSILE_EIGENVALUES
	//   dszz = gd*(lambda*dTrep+twoG*dvzzeffp) + (lambda*dTrem+twoG*dvzzeffm)
	//					+ (gdp*(lambda*trenm1p+twoG*enm1p(2,2))*dPhase)
	// for PRESSURE_SHEAR
	//   dszz = gd(Ksp*dTrep + 2G*dezzdev) + (Ksp*dTrem) + gdp*(Ksp*trenm1p+twoG*ezzdevnm1)
	if(np==PLANE_STRESS_MPM)
	{	double dezz;
		if(partition==TENSILE_EIGENVALUES)
		{	// constant terms depending only on previous strain
			double prevTerm = gdp*(lambda*trenm1p + twoG*enm1p(2,2))*dPhase;
			double ezm1p = enm1p(2,2);
			double ezm1m = enm1(2,2)-ezm1p;
		
			// first assume tren+dezz>0 and en(2,2)+dezz>0
			// then dTrepTest=tren+dezz-trenm1p, dTrem=-trenm1m, dvzzeffp=en(2,2)+dezz-ezm1p, dvzzeffm=-ezm1m
			dezz = (-gd*(lambda*(tren-trenm1p)+twoG*(en(2,2)-ezm1p)) - (-lambda*trenm1m-twoG*ezm1m) - prevTerm)
										/(gd*C11);
			if(tren+dezz<0. || en(2,2)+dezz<0.)
			{	// first try failed, now assume tren+dezz<0 and en(2,2)+dezz<0
				// then dTrep=-trenm1p, dTrem=tren+dezz-trenm1m, dvzzeffp=-enm1p, dvzzeffm=en(2,2)+dezz-ezm1m
				dezz = (-gd*(-lambda*trenm1p-twoG*ezm1p) - (lambda*(tren-trenm1m)+twoG*(en(2,2)-ezm1m)) - prevTerm)
										/C11;
				if(tren+dezz>0. || en(2,2)+dezz>0.)
				{	// second try failed, now assume tren+dezz>0 and en(2,2)+dezz<0
					// then dTrep=tren+dezz-trenm1p, dTrem=-trenm1m, dvzzeffp=-enm1p, dvzzeffm=en(2,2)+dezz-ezm1m
					dezz = (-gd*(lambda*(tren-trenm1p)-twoG*ezm1p) - (-lambda*trenm1m+twoG*(en(2,2)-ezm1m)) - prevTerm)
										/(gd*lambda+twoG);
					if(tren+dezz<0. || en(2,2)+dezz>0.)
					{	// third try failed, now better be tren+dezz<0 and en(2,2)+dezz>0
						// then dTrep=-trenm1p, dTrem=tren+dezz-trenm1m, dvzzeffp=en(2,2)+dezz-ezm1p, dvzzeffm=-ezm1m
						dezz = (-gd*(-lambda*trenm1p+twoG*(en(2,2)-ezm1p)) - (lambda*(tren-trenm1m)-twoG*ezm1m) - prevTerm)
										/(lambda+gd*twoG);
					}
				}
			}
		}
		else
		{	// constant terms depending only on previous strain
			double ezzdevnm1 = enm1(2,2)-trenm1/3.;
			double prevTerm = gdp*(Ksp*trenm1p + twoG*ezzdevnm1)*dPhase;
			double dezzdev = en(2,2) - tren/3. - ezzdevnm1;
			
			// first assume tren+dezz>0
			// then dTrep=tren+dezz-trenm1p, dTrem=-trenm1m
			dezz = (-gd*(Ksp*(tren-trenm1p)+twoG*dezzdev) + Ksp*trenm1m - prevTerm)
								/(gd*C11);
			if(tren+dezz<0)
			{	// first try failed, we must have dtren+dezz<0
				// then dTrep=-trenm1p, dTrem=tren+dezz-trenm1m
				dezz = (-gd*(-Ksp*trenm1p+twoG*dezzdev) - Ksp*(tren-trenm1m) - prevTerm)
							/(Ksp+2.*gd*twoG/3.);
			}
		}
		
		// Update terms that need dezz
		tren += dezz;
		en(2,2) += dezz;
		dVoverV += dezz;
		de(2,2) = dezz;
		Vn(2,2) *= (1.+dezz);
		mptr->IncrementDeformationGradientZZ(dezz);
	}
	
	// see if failed
	if(gd<=kStability)
	{	// mark and display as failed
		if(mptr->GetHistoryDble(PHASE_FAILURE_STATE,0)<0.5)
		{	mptr->SetHistoryDble(PHASE_FAILURE_STATE,1.0,0);
#pragma omp critical (output)
			{	// assumes using Legacy units and in format of ADaM output
				cout << "# Decohesion: t=" << mtime*1e3;
				cout << " x=(" << mptr->pos.x << "," << mptr->pos.y << "," << mptr->pos.z << ")";
				cout << "  (GI,GII1,GII2,Gtot)=(0.0,0.0,0.0," <<
						1.e-9*mptr->mp*mptr->GetPlastEnergy() << ")" << endl;
			}
		}
	}
	
	// signed traces
	trenp = fmax(tren,0.);
	trenm = tren - trenp;
	
	// trace increments
	double dTrep = trenp - trenm1p;
	double dTrem = trenm - trenm1m;

	// Terms to find in next sections
	Tensor delsp;				// stress increment
	double psiPlus;				// dissipation energy term
	
	if(partition==TENSILE_EIGENVALUES)
	{	// get current eigenvalues
		// Plane strain: lam.z = -eres-resStrain
		// Plane stress: lam.z = ezz-resStrain+dezz-eres
		Vector lam = en.Eigenvalues();
		Matrix3 Qcol = en.Eigenvectors(lam);
		Matrix3 QcolT = Qcol.Transpose();
		Matrix3 Lam(fmax(lam.x,0.),0.,0.,fmax(lam.y,0.),fmax(lam.z,0.));
		Matrix3 LamQcolT = Lam*QcolT;
		Matrix3 enp = Qcol*LamQcolT;
		
		// effective + and - strains and increments (as needed)
		double dvxxeff = de(0,0)-eres;
		double dvyyeff = de(1,1)-eres;
		double dvzzeff = de(2,2)-eres;
		double dgamxy = de(0,1)+de(1,0);
		double dvxxeffp = enp(0,0) - enm1p(0,0),dvxxeffm = dvxxeff - dvxxeffp;
		double dvyyeffp = enp(1,1) - enm1p(1,1),dvyyeffm = dvyyeff - dvyyeffp;
		double dvzzeffp = enp(2,2) - enm1p(2,2),dvzzeffm = dvzzeff - dvzzeffp;
		double dgamxyp = enp(0,1) + enp(1,0) - enm1p(0,1) - enm1p(1,0),dgamxym = dgamxy - dgamxyp;
		double dgamxz=0.,dgamyz=0.,dgamxzp=0.,dgamxzm=0.,dgamyzp=0.,dgamyzm=0.;
		if(np==THREED_MPM)
		{	dgamxz = de(0,2)+de(2,0);
			dgamxzp = enp(0,2) + enp(2,0) - enm1p(0,2) - enm1p(2,0);
			dgamxzm = dgamxz - dgamxzp;
			dgamyz = de(1,2)+de(2,1);
			dgamyzp = enp(1,2) + enp(2,1) - enm1p(1,2) - enm1p(2,1);
			dgamyzm = dgamyz - dgamyzp;
		}

		// stress increment
		double gdpdd = gdp*dPhase;
		double dsigCon = gd*lambda*dTrep + lambda*dTrem + gdpdd*lambda*trenm1p;
		delsp.xx = dsigCon + twoG*(gd*dvxxeffp + dvxxeffm + gdpdd*enm1p(0,0));
		delsp.yy = dsigCon + twoG*(gd*dvyyeffp + dvyyeffm + gdpdd*enm1p(1,1));
		delsp.xy = p->C[3][3]*(gd*dgamxyp + dgamxym + gdpdd*(enm1p(0,1)+enm1p(1,0)));
		
		// positive energy as well as disspation energy
		psiPlus = 0.5*lambda*trenp*trenp + 0.5*twoG*(enp(0,0)*enp(0,0) + enp(1,1)*enp(1,1)
										+ enp(2,2)*enp(2,2) + 2.*enp(0,1)*enp(0,1));

		if(np==PLANE_STRAIN_MPM || np==AXISYMMETRIC_MPM)
		{	delsp.zz = dsigCon + twoG*(gd*dvzzeffp + dvzzeffm + gdpdd*enm1p(2,2));
		}
		else if(np==THREED_MPM)
		{	delsp.zz = dsigCon + twoG*(gd*dvzzeffp + dvzzeffm + gdpdd*enm1p(2,2));
			delsp.xz = p->C[3][3]*(gd*dgamxzp + dgamxzm + gdpdd*(enm1p(0,2)+enm1p(2,0)));
			delsp.yz = p->C[3][3]*(gd*dgamyzp + dgamyzm + gdpdd*(enm1p(1,2)+enm1p(2,1)));
			psiPlus += twoG*(enp(0,2)*enp(0,2) + enp(1,2)*enp(1,2));
		}
		else
		{	// plane stress
			delsp.zz = 0.;
		}
	}
	else
	{	// effective strains and increments (as needed)
		double devThird = (tren-trenm1)/3.;
		double evnm1Third = trenm1/3.;
		double dexxdev = de(0,0)-eres - devThird;
		double deyydev = de(1,1)-eres - devThird;
		double dezzdev = de(2,2)-eres - devThird;
		double dgamxy = de(0,1)+de(1,0);

		// stress increment
		double gdpdd = gdp*dPhase;
		double dsigNorm = gd*Ksp*dTrep + Ksp*dTrem + gdpdd*(Ksp*trenm1p-twoG*evnm1Third);
		delsp.xx = dsigNorm + twoG*(gd*dexxdev + gdpdd*enm1(0,0));
		delsp.yy = dsigNorm + twoG*(gd*deyydev + gdpdd*enm1(1,1));
		delsp.xy = p->C[3][3]*(gd*dgamxy + gdpdd*(enm1(0,1)+enm1(1,0)));
		
		// positive energy as well as disspation energy
		double enxxdev = en(0,0)-tren/3.;
		double enyydev = en(1,1)-tren/3.;
		double enzzdev = en(2,2)-tren/3.;
		psiPlus = 0.5*Ksp*trenp*trenp + 0.5*twoG*(enxxdev*enxxdev + enyydev*enyydev
										+ enzzdev*enzzdev + 2.*en(0,1)*en(0,1));

		if(np==PLANE_STRAIN_MPM || np==AXISYMMETRIC_MPM)
		{	delsp.zz = dsigNorm + twoG*(gd*dezzdev + gdpdd*enm1(2,2));
		}
		else if(np==THREED_MPM)
		{	delsp.zz = dsigNorm + twoG*(gd*dezzdev + gdpdd*enm1(2,2));
			double dgamxz = de(0,2)+de(2,0);
			double dgamyz = de(1,2)+de(2,1);
			delsp.xz = p->C[3][3]*(gd*dgamxz + gdpdd*(enm1(0,2)+enm1(2,0)));
			delsp.yz = p->C[3][3]*(gd*dgamyz + gdpdd*(enm1(1,2)+enm1(2,1)));
			psiPlus += twoG*(en(0,2)*en(0,2) + en(1,2)*en(1,2));
		}
		else
		{	// plane stress
			delsp.zz = 0.;
		}
	}

	// increment stress: rotate previous stress and add in-plane components
	*sp = dR.RVoightRT(sp,true,true);
	sp->xx += delsp.xx;
	sp->yy += delsp.yy;
	sp->xy += delsp.xy;
	sp->zz += delsp.zz;
	if(np==THREED_MPM)
	{	sp->xz += delsp.xz;
		sp->yz += delsp.yz;
	}
	
	// work and resdidual strain energy increments (2D and axisymetric)
	double workEnergy = sp->xx*de(0,0) + sp->yy*de(1,1) + sp->xy*(de(0,1)+de(1,0));
	if(np==THREED_MPM) workEnergy += sp->xz*(de(0,2)+de(2,0)) + sp->yz*(de(1,2)+de(2,1));
	double resEnergy = (sp->xx + sp->yy + sp->zz)*eres;
	mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
	
	// dissipated energy per unit mass
	// stored in plastic energy (viz tools will convert to Joules in Legacy or energy units for consistent)
	double dispEnergy = -gdp*psiPlus*dPhase;
	mptr->AddPlastEnergy(dispEnergy);
	
	// Isoentropic temperature rise = -(K 3 alpha T)/(rho Cv) (dV/V) = - gamma0 T (dV/V)
	// Used reduced bulk modulus in tension
	double dTq0 = -gamma0*mptr->pPreviousTemperature*dVoverV;
	if(dVoverV>0.) dTq0 *= gd;
	
	// track heat energy
	IncrementHeatEnergy(mptr,dTq0,dispEnergy);

	// update phase energy history variable not using reduced properties
	psiPlus *= rho;
	double psiP = mptr->GetHistoryDble(HISTORY_VALUE,0);
	if(psiPlus-Psii>psiP) mptr->SetHistoryDble(HISTORY_VALUE,psiPlus-Psii,0);
	
	// add particle source term to phase field task
	if(taskNum>=0)
    {   double psource;
        if(GetParticleDiffusionSource(phaseTask,mptr,delTime,&psource))
            mptr->AddParticleDiffusionSource(taskNum,psource);
    }
}

#pragma mark IsoPhaseFieldSoftening::Accessors

// return material type
const char *IsoPhaseFieldSoftening::MaterialType(void) const
{	return "Isotropic Phase Field Softening";
}

// Get g(d) and g'(d) using current phase field value
double IsoPhaseFieldSoftening::GetDegradation(MPMBase *mptr,double &gdp) const
{
    double phi = mptr->pDiff[taskNum]->prevConc;
    double gd = 0.;
    if(phi<1.)
    {   switch(gdMode)
        {   case GEN_QUADRATIC:
            {   gd = (1.-phi)*(1.-garg*phi);
                gdp = -(1.+garg-2*garg*phi);
                break;
            }
            case FINITE_EXPONENTIAL:
            {   double arg = exp(-garg*phi);
                gd = ec1*arg+ec2;
                gdp = -ec3*arg;
                break;
            }
            case PF_LINEAR_SOFTENING:
            {   double arg = sqrt(1.+garg*(2.+garg*phi*phi));
                double numer = garg*phi-arg;
                gd = (1 + garg*phi*phi - phi*arg);
                gdp = -numer*numer/arg;
                break;
            }
            default:
                gdp = 0.;
                break;
        }
    }
    else
        gdp = 0.;
    
    // add stability factor if used
    if(hasStability)
    {   gdp *= (1-kStability);
        gd = (1.-kStability)*gd + kStability;
    }
    
    // return g(d)
    return gd;
}

#pragma mark IsoPhaseFieldSoftening::Phase Field Diffusion Accessors

// This material support fracture phase field
// None of the subsequent methods are called unless this is true
bool IsoPhaseFieldSoftening::SupportsPhaseField(int phaseStyle) const
{   return phaseStyle==FRACTURE_PHASE_FIELD;
}

// Get maximum phase field diffusivity for time step calculations JIc*ell/viscocity
double IsoPhaseFieldSoftening::MaximumPhaseDiffusivity(int phaseStyle) const
{   return pfDTensor.xx/viscosity ;
}

// Viscosity for phase field diffusion (do not call unless doing fracture phase field)
double IsoPhaseFieldSoftening::GetPhaseFieldViscosity(int phaseStyle) const
{   return viscosity;
}

// Diffusion constant for phase field (do not call unless doing fracture phase field)
Tensor IsoPhaseFieldSoftening::GetPhaseFieldDTensor(int phaseStyle,MPMBase *mptr) const
{    return pfDTensor;
}

// For JIc/ell and g'(phi) for phase field diffusion (g'(phi) in arg)
// Note this is only called by phase field diffusion task and only with mechanics activate
// Because part of diffusion equation, needs particle values (not smoothed grid values)
bool IsoPhaseFieldSoftening::GetParticleDiffusionSource(DiffusionTask *task,MPMBase *mptr,
                                                        double delTime,double *source) const
{
    // part or all as particle source (need phi and derivative
    
    // value on the particle
    // FLIP might might but not work well. They are same with FMPM(k) (or close if periodic)
    //double phi = task->GetParticleValue(mptr);
    double phi = task->GetPrevParticleValue(mptr);

    // the derivative
    double gdp = 0.;
    if(phi<1.)
    {    switch(gdMode)
        {    case GEN_QUADRATIC:
                gdp = -(1.+garg-2*garg*phi);
                break;
            case FINITE_EXPONENTIAL:
                gdp = -ec3*exp(-garg*phi);
                break;
            case PF_LINEAR_SOFTENING:
            {   double arg = sqrt(1.+garg*(2.+garg*phi*phi));
                double numer = garg*phi-arg;
                gdp = -numer*numer/arg;
                break;
            }
            default:
                gdp = 0.;
                break;
        }
        if(hasStability) gdp *= (1.-kStability);
    }
    
    // all source terms
	double psource = -JIc*phi/ell - gdp*mptr->GetHistoryDble(HISTORY_VALUE,0);
    
    *source = psource*delTime/viscosity;
    return true;
}

// For JIc/ell and g'(phi) for phase field diffusion
// Only called by phase field diffusion and after verifed OK to call
void IsoPhaseFieldSoftening::AdjustPhaseFieldValue(DiffusionTask *phaseTask,double pConc,double pPrevConc,
                                        double &value,double &rate,double &lumpedValue,double deltime) const
{ 
    if(pConc>=1.)
    {   // already failed, keep at 1
        value = 1.;
        rate = 0.;
		lumpedValue = 1.;
    }
    else
    {   // extrapolated fracture phase field must increase and cannot exceed 1
		// 10/7/2022 change to check pPrevCon here instead of pConc
		value = fmin(fmax(pPrevConc,value),1.);
		
		// do lumped value to in case blended FLIP/FMPM(k>1)
		lumpedValue = fmin(fmax(pPrevConc,lumpedValue),1.);

        // make sure rate positive and does not cause value to exceed 1
		// Only needed for FLIP OR blended FLIP/FMPM
		double fractFMPM;
		bool usingFMPM = phaseTask->IsUsingTransportXPIC(fractFMPM);
		if(!usingFMPM || fractFMPM<1.)
		{	if(rate<=0.)
				rate = 0.;
			else
			{	// FLIP rate
				rate = fmin(rate,(1.-pConc)/deltime);
			}
		}
    }
}

// Store phase field value in history (only call after verfied supports phase diffusion task)
// This material only supports one type of phase field
void IsoPhaseFieldSoftening::StorePhaseFieldValue(int phaseStyle,MPMBase *mptr,double newValue) const
{   mptr->SetHistoryDble(PHASE_SOLVED,newValue,0);
}

// Store delta phase field value in history (only call after verfied supports phase diffusion task)
// This material only supports one type of phase field
void IsoPhaseFieldSoftening::StorePhaseFieldDelta(int phaseStyle,MPMBase *mptr,double newDelta) const
{	// for USAVG, need to cut delta in half because it will be used twice
	// Note that this approach does not support fractureUSF if that is ever an option that
	//		deviated from hard coded 0.5
	if(fmobj->mpmApproach==USAVG_METHOD) newDelta *= 0.5;
	
	// exponential softening test
	// delta = alpha*newDelta + (1-alpha)*prevDelta = newDelta + (1-alpha)*(prevDelta-newDelta)
	//double alpha = 0.1;
	//newDelta += (1.-alpha)*(mptr->GetHistoryDble(DELTA_PHASE_SOLVED,0) - newDelta);
	mptr->SetHistoryDble(DELTA_PHASE_SOLVED,newDelta,0);
}

