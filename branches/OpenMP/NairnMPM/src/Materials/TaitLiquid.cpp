/********************************************************************************
    TaitLiquid.cpp
    nairn-mpm-fea

    Created by John Nairn, Dec 4, 2013.
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    See Create_MPM_Material for details in google code web site wiki
 ********************************************************************************/

#include "Materials/TaitLiquid.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark TaitLiquid::Constructors and Destructors

// Constructors
TaitLiquid::TaitLiquid() {}

// The default contructor should call a parent class constructor and
// then fill in any new initialization.
TaitLiquid::TaitLiquid(char *matName) : HyperElastic(matName)
{
    viscosity = -1.;
}

#pragma mark TaitLiquid::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *TaitLiquid::InputMat(char *xName,int &input)
{
	// read properties for this material
    if(strcmp(xName,"viscosity")==0)
    {	input=DOUBLE_NUM;
        return((char *)&viscosity);
    }
	
    return HyperElastic::InputMat(xName,input);
}

// Verify input properties do calculations; if problem return string with an error message
// If OK, MUST pass on to super class. This is called just before PrintMaterial
// (see also ValidateForUse() for checks that depend on MPM calculation mode)
const char *TaitLiquid::VerifyAndLoadProperties(int np)
{
	// make sure all were set
    if(Kbulk <= 0. || viscosity<0.)
		return "TaitLiquid material model needs positive K and zero or positive viscosity";
    
	// Viscosity in Specific units using initial rho (times 2)
	// for MPM (units N sec/m^2 cm^3/g) (1 cP = 0.001 N sec/m^2)
	Etasp = 0.002*viscosity/rho;
    
    // heating gamma0
    double alphaV = 3.e-6*aI;
    gamma0 = 1000.*Kbulk*alphaV/(rho*heatCapacity);
    
	// must call super class
	return HyperElastic::VerifyAndLoadProperties(np);
}

// print mechanical properties to the results
void TaitLiquid::PrintMechanicalProperties(void) const
{
	// print properties
	PrintProperty("K",Kbulk,"");
	PrintProperty("eta",viscosity,"cP");
	PrintProperty("a",aI,"");
    cout << endl;
    PrintProperty("gam0",gamma0,"");
}

// if analysis not allowed, throw a CommonException
void TaitLiquid::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
    {	throw CommonException("TaitLiquid material cannot do 2D plane stress analysis",
                          "TaitLiquid::ValidateForUse");
    }

    return HyperElastic::ValidateForUse(np);
}

#pragma mark TaitLiquid:HistoryVariables

// Particle J
char *TaitLiquid::InitHistoryData(void)
{	J_history = 0;
	double *p = CreateAndZeroDoubles(J_history+1);
	p[J_history]=1.;					// J
	return (char *)p;
}

// this material has one
double TaitLiquid::GetHistory(int num,char *historyPtr) const
{
    double history=0.;
	if(num>0 && num<=J_history+1)
	{	double *p=(double *)historyPtr;
		history=p[num-1];
	}
    return history;
}

#pragma mark TaitLiquid:Step Methods

// Apply Constitutive law, check np to know what type
void TaitLiquid::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// Update total deformation gradient, and calculate trial B
	double detdF = IncrementDeformation(mptr,du,NULL,np);
    
	// Deformation gradients and Cauchy tensor differ in plane stress and plane strain
    // This code handles plane strain, axisymmetric, and 3D - Plane stress is blocked
	double Jtot = detdF * mptr->GetHistoryDble(J_history);
    
    // save new J
    mptr->SetHistoryDble(J_history,Jtot);  // Stocking J
	
    // account for residual stresses
	double dresStretch,resStretch = GetResidualStretch(mptr,dresStretch,res);
	double Jres = resStretch*resStretch*resStretch;
    double detdFres = dresStretch*dresStretch*dresStretch;
    double J = Jtot/Jres;
    detdF /= detdFres;

    // new pressure from Tait equation
    double pressure = TAIT_C*Ksp*(exp((1.-J)/TAIT_C)-1.);
    //double pressure = Ksp*(pow(Jtot,-7.)-1.);
    mptr->SetPressure(pressure);
    
    // viscosity term = 2 eta (0.5(grad v) + 0.5*(grad V)^T - (1/3) tr(grad v) I)
    Matrix3 shear;
    double c[3][3];
    c[0][0] = (2.*du(0,0)-du(1,1)-du(2,2))/3.;
    c[1][1] = (2.*du(1,1)-du(0,0)-du(2,2))/3.;
    c[2][2] = (2.*du(2,2)-du(0,0)-du(1,1))/3.;
    c[0][1] = 0.5*(du(0,1)+du(1,0));
    c[1][0] = c[0][1];
    if(np==THREED_MPM)
    {   c[0][2] = 0.5*(du(0,2)+du(2,0));
        c[2][0] = c[0][2];
        c[1][2] = 0.5*(du(1,2)+du(2,1));
        c[2][1] = c[1][2];
        shear.set(c);
    }
    else
        shear.set(c[0][0],c[0][1],c[1][0],c[1][1],c[2][2]);
    shear.Scale(Etasp/delTime);
    
    // update deviatoric stress
	Tensor *sp=mptr->GetStressTensor();
    sp->xx = shear(0,0);
    sp->yy = shear(1,1);
    sp->zz = shear(2,2);
    sp->xy = shear(0,1);
    if(np==THREED_MPM)
    {   sp->xz = shear(0,2);
        sp->yz = shear(1,2);
    }
    
    // particle isentropic temperature increment dT/T = - gamma0 (del(J)/J)
    // dJ/J = 1. - 1/detdF (relative to stress free change)
    double delV = 1. - 1./detdF;
    double dTq0 = -gamma0*mptr->pPreviousTemperature*delV;
    
    // heat energy is Cv (dT - dTq0) - dPhi
	// Here do Cv (dT - dTq0)
    // dPhi is lost due to shear term (but should it be added)
    IncrementHeatEnergy(mptr,res->dT,dTq0,0.);

}

#pragma mark TaitLiquid::Custom Methods

#pragma mark TaitLiquid::Accessors

// return unique, short name for this material
const char *TaitLiquid::MaterialType(void) const { return "Tait Liquid"; }

// Return the material tag
int TaitLiquid::MaterialTag(void) const { return TAITLIQUID; }

// Calculate wave speed for material in mm/sec.
double TaitLiquid::WaveSpeed(bool threeD,MPMBase *mptr) const
{
    return sqrt(1.e9*(Kbulk)/rho);
}

// Copy stress to a read-only tensor variable after combininng deviatoric and pressure
Tensor TaitLiquid::GetStress(Tensor *sp,double pressure) const
{   Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}
