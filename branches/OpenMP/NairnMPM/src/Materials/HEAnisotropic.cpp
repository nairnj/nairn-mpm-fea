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
}

#pragma mark HEAnisotropic::Initialization

// Read material properties
char *HEAnisotropic::InputMat(char *xName,int &input)
{
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
{	if(np==PLANE_STRESS_MPM || np==THREED_MPM)
    {	throw CommonException("HEAnisotropic materials can only be used in plane strain",
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
{   double *p = CreateAndZeroDoubles(2);
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
    double THETA0 = 45.;
    
    Matrix3 pF = mptr->GetDeformationGradientMatrix();
    Tensor Btrial;
    double detdF = IncrementDeformation(mptr,du,&Btrial,np);
    
    // Deformation gradients and Cauchy tensor differ in plane stress and plane strain
    // This code handles plane strain, axisymmetric, and 3D - Plane stress is blocked
	// Increment J and save it in history data
    double J = detdF*mptr->GetHistoryDble(J_HISTORY);
    mptr->SetHistoryDble(J_HISTORY,J);
	
	//double JJ = pF.determinant();//detdF * mptr->GetHistoryDble(J_history);
    //cout << "   Det J =    " << J <<endl;
    
    //cout << "   Deformation Gradient in Basis GxG  =    " << pF <<endl;
    //cout << " Verif in HEAniso J =  " << JJ << endl;
    
    
    
    // cout << " Verif in HEAniso pF =  " << pF << endl;
    
    //Matrix3 pFT = pF.Transpose();
    //Matrix3 C = pFT*pF;
    
    // Right Stretch Tensor - Calculation with CG
    //Matrix3 Eigenvals = C.Eigenvalues();
    //Matrix3 U = C.RightStretch(Eigenvals);
    
    // Rotation Tensor - Calculation with FG
    //Matrix3 R = pF.Rotation(Eigenvals,U);
    
    //  Verification of the Polar Decomposition
    //Matrix3 F_Verif = R*U;
    
    
    //   Initial Orientation
    double pi = 3.14159265358979;
    double THETA0_Rad= THETA0*pi/180.;
    
    // Initial Orientation of the Fibers in the Global Basis
    Matrix3 O(cos(THETA0_Rad),-sin(THETA0_Rad),0,sin(THETA0_Rad),cos(THETA0_Rad),0,0,0,1);
    
    
    
    // Right Cauchy-Green in the Local Orientation of the fibers
    
    Matrix3 F_GxL = pF*O;
    //cout << "   Deformation Gradient Tensor in mixed GxL =    " << F_GxL <<endl;
    
    Matrix3 FT_LxG = F_GxL.Transpose();
    //cout << "   Deformation Gradient Tensor Transpose, so in LxG =    " << FT_LxG <<endl;
    
    Matrix3 C_LxL = FT_LxG*F_GxL;
    
    //Matrix3 C_LxL = O.Transpose()*C*O;  For verification
    //cout << "   Right Cauchy-Green Tensor in Local Basis 2 =    " << C_LxL <<endl;
    
    // Actual Elongations in the Meches Directions
    
    double lam1=sqrt(C_LxL(0,0));
    double lam2=sqrt(C_LxL(1,1));
    
    double lam1_2=pow(lam1, 2.);
    double lam1_3=pow(lam1, 3.);
    
    double lam2_2=pow(lam2, 2.);
    double lam2_3=pow(lam2, 3.);
    
    //cout << "   Elongation 1 =    " << lam1 << "   Elongation 2 =    " << lam2 << endl;
    //========================================
    
    
    // Actual Rotation of the Meches
    
    double CTHETA12=C_LxL(0,1)/(lam1*lam2);
    double THETA12=acos(CTHETA12);
    double THETA12_DEG=THETA12*180./pi;
    mptr->SetHistoryDble(THETA12_HISTORY,THETA12_DEG);
    //cout << "   THETA12_DEG =    " << THETA12_DEG << endl;
    //cout << "   THETA0 =    " << THETA0 << endl;
    
    // PK2 in the LxL
    
    //Matrix3 PK2(E1sp*(lam1*lam1-1),(G1sp*CTHETA12+G2sp*CTHETA12*CTHETA12+G3sp*CTHETA12*CTHETA12*CTHETA12)/(lam1*lam2),0,
    //            (G1sp*CTHETA12+G2sp*CTHETA12*CTHETA12+G3sp*CTHETA12*CTHETA12*CTHETA12)/(lam1*lam2), E2sp*(lam2*lam2-1),0,0,0,0);
    
    // Tension Material properties   -  Balanced Fabric
    double T3sp=112.55*1.0e+06/rho;
    double T2sp=-670.34*1.0e+06/rho;
    double T1sp=1472.89*1.0e+06/rho;
    double T0sp=-915.1*1.0e+06/rho;
    
    // Diagonal Components of PK2 in Local Orientation of Fabric 
    double S11=T3sp*lam1_3+T2sp*lam1_2+T1sp*lam1+T0sp;
    double S22=T3sp*lam2_3+T2sp*lam2_2+T1sp*lam2+T0sp;
    
    // Shear Material properties
    double C2sp=-646.85*1.0e+06/rho;
    double C1sp=-823.49*1.0e+06/rho;
    double C0sp=1470.34*1.0e+06/rho;
    //cout << "   rho =    " << rho << endl;
    
    
    // Shear Component of PK2 in Local Orientation of Fabric
    double STHETA12=sin((180-THETA12_DEG)*pi/180);
    double STHETA12_2=pow(STHETA12, 2.);
    
    double S12=C2sp*STHETA12_2+C1sp*STHETA12+C0sp;
    
    
    // PK2 in Local Orientation LxL
    Matrix3 PK2(S11,S12,0,S12,S22,0,0,0,0);
    //cout << "   PK2 Tensor in Local Basis =    " << PK2 <<endl;
    
    
    
    if (THETA12_DEG>90.)
    {
        PK2(1,2)=PK2(1,2);// Don't forget to put -PK2(1,2) after;
    }
    
    if (THETA12_DEG==90.)
    {
        PK2(1,2)=0.;
    }
    
    PK2(2,1)=PK2(1,2);
    
    
    // PK2 in GxG
    Matrix3 PK2_GxG=O*PK2*O.Transpose();
    //cout << "   PK2 Tensor in Global Basis =    " << PK2_GxG <<endl;
    
    
    
    
    
    // sp in GxG
    Matrix3 sp_GxG=pF*PK2_GxG*pF.Transpose();
    //cout << "   Cauchy Tensor in Global Basis =    " << sp_GxG <<endl;
    
    
    // store initial stress
    Tensor *sp = mptr->GetStressTensor();
    
    // save new J
    //mptr->SetHistoryDble(J_history,J2);  // Stocking J
    
    // J as determinant of F (or sqrt root of determinant of B) normalized to residual stretch
    // Must also divide elements of B by resStretch2
    
    // Others constants
    //double J23 = pow(J, 2./3.);
    
    
    sp->xx = sp_GxG(0,0);//1./3.*E1sp*(2.*B(0,0)-B(1,1))/J23;
    sp->yy = sp_GxG(1,1);//1./3.*E1sp*(2.*B(1,1)-B(0,0))/J23;
    sp->zz = 0.;///3.*G1sp;//;*(-B(0,0)-B(1,1))/J23;
    sp->xy = sp_GxG(1,0);//G1sp*B(0,1)/J23;
    
    //cout << "   Cauchy Tensor in Global Basis in MPM=    " << sp_GxG <<endl;
    
    
    //return;
    
        
}

#pragma mark HEAnisotropic::Custom Methods

#pragma mark HEAnisotropic::Accessors

// Return the material tag
int HEAnisotropic::MaterialTag(void) const { return HEANISOTROPIC; }

// return unique, short name for this material
const char *HEAnisotropic::MaterialType(void) const { return "Hyperelastic Anisotropic"; }

// calculate wave speed in mm/sec - needs to be implemented
double HEAnisotropic::WaveSpeed(bool threeD,MPMBase *mptr) const { return 1.e-12; }

