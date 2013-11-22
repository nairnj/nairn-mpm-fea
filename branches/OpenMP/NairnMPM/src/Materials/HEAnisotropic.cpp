/********************************************************************************
    HEAnisotropic.cpp
    nairn-mpm-fea
    
    Created by John Nairn, Sept 27, 2011.
    Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		HyperElastic.hpp (MaterialBase.hpp)
********************************************************************************/

#include "HEAnisotropic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark HEAnisotropic::Constructors and Destructors

// Constructors
HEAnisotropic::HEAnisotropic() {}

/* The default contructor should call a parent class constructor and
	then fill in any new initialization.
	*/
// Constructors
HEAnisotropic::HEAnisotropic(char *matName) : HyperElastic(matName)
{
	// inintialize unique properties
    THETA0= 0.;			// Initial orientation - Direction 1
	Eyarn = 0.;			// Equivalent Elastic Modulus of the fabric
    nu = 0.;			// Equivalent Elastic Modulus of the fabric
    T10= 0.;			// Traction Interpolation coeff in direction 1
	T11= 0.;
    T12= 0.;
    T13= 0.;
    T14= 0.;
    T20= 0.;			// Traction Interpolation coeff in direction 2
	T21= 0.;
    T22= 0.;
    T23= 0.;
    T24= 0.;
    SH0= 0.;			// Shearing Interpolation coeff in XY plane
	SH1= 0.;
    SH2= 0.;
    SH3= 0.;
    SH4= 0.;
    balanced = TRUE;
}

#pragma mark HEAnisotropic::Initialization

// Read material properties
char *HEAnisotropic::InputMat(char *xName,int &input)
{
	input=DOUBLE_NUM;
    
    if(strcmp(xName,"THETA0")==0)
       return((char *)&THETA0);
    else if(strcmp(xName,"Eyarn")==0)
        return((char *)&Eyarn);
    else if(strcmp(xName,"nu")==0)
        return((char *)&nu);
    else if(strcmp(xName,"T10")==0)
        return((char *)&T10);
    else if(strcmp(xName,"T11")==0)
        return((char *)&T11);
    else if(strcmp(xName,"T12")==0)
        return((char *)&T12);
    else if(strcmp(xName,"T13")==0)
        return((char *)&T13);
    else if(strcmp(xName,"T14")==0)
        return((char *)&T14);
    else if(strcmp(xName,"T20")==0)
    {   balanced = FALSE;
        return((char *)&T20);
    }
    else if(strcmp(xName,"T21")==0)
    {   balanced = FALSE;
        return((char *)&T21);
    }
    else if(strcmp(xName,"T22")==0)
    {   balanced = FALSE;
        return((char *)&T22);
    }
    else if(strcmp(xName,"T23")==0)
    {   balanced = FALSE;
        return((char *)&T23);
    }
    else if(strcmp(xName,"T24")==0)
    {   balanced = FALSE;
        return((char *)&T24);
    }
    else if(strcmp(xName,"SH0")==0)
        return((char *)&SH0);
    else if(strcmp(xName,"SH1")==0)
        return((char *)&SH1);
    else if(strcmp(xName,"SH2")==0)
        return((char *)&SH2);
    else if(strcmp(xName,"SH3")==0)
        return((char *)&SH3);
    else if(strcmp(xName,"SH4")==0)
        return((char *)&SH4);
    
    // read properties for this material
    //if(strcmp(xName,"mypropertyname")==0)
    //{	input=DOUBLE_NUM;
    //    return((char *)&newproperty);
    //}
	
    return(HyperElastic::InputMat(xName,input));
}

// verify settings and maybe some initial calculations
const char *HEAnisotropic::VerifyAndLoadProperties(int np)
{
	// check properties

	// must call super class
	return HyperElastic::VerifyAndLoadProperties(np);
}

// plane stress not allowed in viscoelasticity
void HEAnisotropic::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)// || np==THREED_MPM)
    {	throw CommonException("HEAnisotropic materials can only be used in 3D calculation",
                          "HEAnisotropic::ValidateForUse");
    }
	
	//call super class (why can't call super class?)
	return HyperElastic::ValidateForUse(np);
}

// print mechanical properties to the results
void HEAnisotropic::PrintMechanicalProperties(void) const
{	
	// call superclass here if it is not Material base
	HyperElastic::PrintMechanicalProperties();
	
	// add new properties here
}

// First ones for hardening law. Particle J appended at the end
char *HEAnisotropic::InitHistoryData(void)
{   double *p = CreateAndZeroDoubles(NUM_HISTORY);
	p[J_HISTORY]=1.;					// J
    p[THETA12_HISTORY]=90.;
	return (char *)p;
}

#pragma mark HEAnisotropic:Methods

// buffer size for mechanical properties
//int HEAnisotropic::SizeOfMechanicalProperties(int &altBufferSize) const
//{   
//}

// Get elastic and plastic properties, return null on error
//void *HEAnisotropic::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
//{
//}

/* Take increments in strain and calculate new Particle: strains, rotation strain,
 stresses, strain energy,
 dvij are (gradient rates X time increment) to give deformation gradient change
 For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
 This material tracks pressure and stores deviatoric stress only in particle stress
 tensor
 */
void HEAnisotropic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
   
    // Getting Deformation Gradient in Global Basis GxG
    Matrix3 pF = mptr->GetDeformationGradientMatrix();
    
    double F33=1./(pF(0,0)*pF(1,1)-pF(0,1)*pF(1,0));
    pF(2,2) = pow(F33, nu/(1-nu));
    //pF.set(2,2,F33_p);

    
    Tensor Btrial;
    double detdF = IncrementDeformation(mptr,du,&Btrial,np);
    
    // Deformation gradients and Cauchy tensor differ in plane stress and plane strain
    // This code handles plane strain, axisymmetric, and 3D - Plane stress is blocked
	// Increment J and save it in history data
    double J = detdF*mptr->GetHistoryDble(J_HISTORY);
    mptr->SetHistoryDble(J_HISTORY,J);
    
	
    // Getting Initial Orientation of the fabrib;
    //double THETA0 = mptr->GetAnglez0InRadians();   ACTIVATE AFTER ***********

    double pi = 3.14159265358979;
    double THETA0_Rad= THETA0*pi/180.;
    
    // Initial Orientation of the Fibers in the Global Basis
    Matrix3 O(cos(THETA0_Rad),-sin(THETA0_Rad),0,sin(THETA0_Rad),cos(THETA0_Rad),0,0,0,1);
    
    
    // Right Cauchy-Green in the Local Orientation of the Fabric LxL
    Matrix3 F_GxL = pF*O;
    Matrix3 FT_LxG = F_GxL.Transpose();
    Matrix3 C_LxL = FT_LxG*F_GxL;
    
    
    // Actual Elongations of the yarns    
    double lam1=sqrt(C_LxL(0,0));
    double lam2=sqrt(C_LxL(1,1));
    
    double lam1_2=pow(lam1, 2.);
    double lam1_3=pow(lam1, 3.);
    
    double lam2_2=pow(lam2, 2.);
    double lam2_3=pow(lam2, 3.);
    
    
    // Actual relative Rotation of the yarns
    
    double CTHETA12=C_LxL(0,1)/(lam1*lam2);
    double THETA12=acos(CTHETA12);
    double THETA12_DEG=THETA12*180./pi;
        
    mptr->SetHistoryDble(THETA12_HISTORY,THETA12_DEG);

    
    // Tension Material properties   -  Balanced Fabric
    double T13sp=(T13*1.0e+06)/rho;
    double T12sp=(T12*1.0e+06)/rho;
    double T11sp=(T11*1.0e+06)/rho;
    double T10sp=(T10*1.0e+06)/rho;
    
    // Tension Material properties   -  for Unbalanced Fabric
    double T23sp=(T23*1.0e+06)/rho;
    double T22sp=(T22*1.0e+06)/rho;
    double T21sp=(T21*1.0e+06)/rho;
    double T20sp=(T20*1.0e+06)/rho;
    
    // Tensile Components of PK2 in Local Orientation LxL of the Fabric 
    double S11=T13sp*lam1_3+T12sp*lam1_2+T11sp*lam1+T10sp;
    double S22=T23sp*lam2_3+T22sp*lam2_2+T21sp*lam2+T20sp;
    
    // Shear Material properties
    double SH2sp=(SH2*1.0e+06)/rho;
    double SH1sp=(SH1*1.0e+06)/rho;
    double SH0sp=(SH0*1.0e+06)/rho;    
    
    // Shear Components of PK2 in Local Orientation LxL of the Fabric
    double STHETA12=sin((180-THETA12_DEG)*pi/180);
    double STHETA12_2=pow(STHETA12, 2.);
    
    double S12=SH2sp*STHETA12_2+SH1sp*STHETA12+SH0sp;
    
    // Component of PK2 in Local Orientation LxL of Fabric
    Matrix3 PK2(S11,S12,0,S12,S22,0,0,0,0);

    if (THETA12_DEG>90.)
    {
        PK2(1,2)=+PK2(1,2);// MAYBE NEEDED TO PUT -PK2(1,2) AFTER REMENDER**************;
    }
        if (THETA12_DEG==90.)
        {  PK2(1,2)=0.;}
        
        PK2(2,1)=PK2(1,2);
    
    if (THETA12_DEG>90.)
    {
       THETA12_DEG=180.-THETA12_DEG;//  0 <= THETA12<= 90 ;
    }
    
    // PK2 in the Global Basis GxG
    Matrix3 PK2_GxG=O*PK2*O.Transpose();
    
    // Specific Cauchy Stress in the Global Basis GxG
    Matrix3 sp_GxG=pF*PK2_GxG*pF.Transpose();

        
    // store initial stress
    Tensor *sp = mptr->GetStressTensor();
    
    // save new J
    //mptr->SetHistoryDble(J_history,J2);  // Stocking J
    
    
    sp->xx = sp_GxG(0,0);
    sp->yy = sp_GxG(1,1);
    sp->zz = 0.;
    sp->xy = sp_GxG(1,0);
    sp->xz = 0.;
    sp->yz = 0.;
    
}

#pragma mark HEAnisotropic::Custom Methods

#pragma mark HEAnisotropic::Accessors

// this material has two history variables
double HEAnisotropic::GetHistory(int num,char *historyPtr) const
{
    double history=0.;
	if(num>0 && num<=NUM_HISTORY)
	{	double *p=(double *)historyPtr;
		history=p[num-1];
	}
    return history;
}

// Return the material tag
int HEAnisotropic::MaterialTag(void) const { return HEANISOTROPIC; }

// return unique, short name for this material
const char *HEAnisotropic::MaterialType(void) const { return "Hyperelastic Anisotropic"; }

// calculate wave speed in mm/sec - needs to be implemented
// calculate wave speed in mm/sec - needs to be implemented
double HEAnisotropic::WaveSpeed(bool threeD,MPMBase *mptr) const
{
	return sqrt(1.e9*Eyarn/rho);
    
}

