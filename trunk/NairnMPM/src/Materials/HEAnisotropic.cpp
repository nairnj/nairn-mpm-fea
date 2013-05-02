/********************************************************************************
    HEAnisotropic.cpp
    NairnMPM
    
    Created by John Nairn, Sept 27, 2011.
    Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		HyperElastic.hpp (MaterialBase.hpp)
********************************************************************************/

#include "HEAnisotropic.hpp"
#include "MPM_Classes/MPMBase.hpp"

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
    G1 = -1.;			// required
	G2 = -1.;			// Required
    E = -1.;

}

#pragma mark HEAnisotropic::Initialization

/* This method reads properties from input file that are unique to this class.
	For details see "Creating New Material Properties" section in the
	Create_MPM_Material Wiki on nairn-fea-mem google code page
*/
// Read material properties
char *HEAnisotropic::InputMat(char *xName,int &input)
{
	// read properties for this material
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"G1")==0)
        return((char *)&G1);
    
    else if(strcmp(xName,"G2")==0)
        return((char *)&G2);
   
    else if(strcmp(xName,"E")==0)
        return((char *)&E);


    //if(strcmp(xName,"mypropertyname")==0)
    //{	input=DOUBLE_NUM;
    //    return((char *)&newproperty);
    //}
	
    return(HyperElastic::InputMat(xName,input));
}

/* This method is called after input file is read but before the new
	material is printed to the results file. If necessary, verify that
	all new properties are valid. If not return string with a description
	of the problem. If OK, must pass on to a superclass
	NOTE: np is analysis type in case that is important
*/
// verify settings and maybe some initial calculations
const char *HEAnisotropic::VerifyProperties(int np)
{
	// check properties
    
    if(G1 < 0. || G2<0.  || E<0.)
 

		return "HENisotropic material needs non-negative G1, G2, and E";


	// must call super class
	return HyperElastic::VerifyProperties(np);
}

/* Called once at beginning (by VerifyProperties() in MaterialBase). For efficiency,
	use this method to calculate new terms that are independent of the particle
	state and thus will remain constant throughout the calculation. When done
	(or before), pass on to super class (but MaterialBase and Elastic do not need it)
*/
// Constant properties used in constitutive law
void HEAnisotropic::InitialLoadMechProps(int makeSpecific,int np)
{
    // MU in specific units using initial rho
	// for MPM (units N/m^2 cm^3/g)
    E1=E;
    E2=E;
    G3=G1;
	G1sp = G1*1.0e+06/rho;			// kept at this initial value
	G2sp = G2*1.e6/rho;				// may change in some hardening laws
    G3sp = G1*1.e6/rho;				// Must introduce G0 later as input properties.
	
	// second G is not used yet
    E1sp = E*1.0e+06/rho;     // E1sp = E2sp for balanced fabric
    E2sp = E*1.0e+06/rho;
    
    // Initial Orientation of the fibers
    THETA0 = 45.;
    
	HyperElastic::InitialLoadMechProps(makeSpecific,np);
}

/* Print all mechanical properties or call parent class and print
	just the new mechanical properties. Use a format compatible with code
	that will read results file and similar to style of other materials
	(need not pass on to MaterialBase or Elastic since they print nothing)
	NOTE: This is called after VerifyProperties() and InitialLoadMechProps()
	Sometimes scaling of properties for internal units is best done here after
	they are printed.
*/
// print mechanical properties to the results
void HEAnisotropic::PrintMechanicalProperties(void)
{
    PrintProperty("E1",E,"");
    PrintProperty("E2",E,"");
    PrintProperty("G1",G1,"");
    PrintProperty("G2",G2,"");
    PrintProperty("G3",G3,"");
    cout << endl;
    
	// call superclass here if it is not Material base
	HyperElastic::PrintMechanicalProperties();
    
    
}
// Store J, which is calculated incrementally, and available for archiving
// initialize to 1
char *HEAnisotropic::InitHistoryData(void)
{
	double *p = CreateAndZeroDoubles(1);
	*p=1.;
	return (char *)p;
}

#pragma mark HEAnisotropic:Methods

/*	Apply 2D constitutive law updating all needed terms for material type. Required updates are:
		stress, strain, plastic strain (all components) (stress should be a specific stress)
		rotation strain (single angle)
		strain energy, plastic energy, and dissipated energy (dissipated needed if want to couple to conduction)
	To support thermal and solvent expansion, include their effect on strains
    If there are material-related data on the particle, update them too
	dvij are (gradient rates X time increment) to give deformation gradient change
   For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
*/
void HEAnisotropic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np)
{
    
    Matrix3 pF = mptr->GetDeformationGradientMatrix();
    
    cout << " Verif in HEAniso F =  " << pF << endl;
    cout << " Verif THETA0 =  " << THETA0 <<"    E1 =  " << E1 <<"   E2 =  " << E2 <<endl;
    cout << " Verif in HEAniso du =  " << du << endl;
    
    // store initial stress
    Tensor *sp = mptr->GetStressTensor();
    
    Tensor Btrial;
	//double detdF = IncrementDeformation(mptr,du,&Btrial,np);
    
    
	// Deformation gradients and Cauchy tensor differ in plane stress and plane strain
    // This code handles plane strain, axisymmetric, and 3D - Plane stress is blocked
	double J2 = 1.;//detdF * mptr->GetHistoryDble(J_history);
    
    // save new J
    //mptr->SetHistoryDble(J_history,J2);  // Stocking J
    
    // J as determinant of F (or sqrt root of determinant of B) normalized to residual stretch
	// Must also divide elements of B by resStretch2
	double dresStretch,resStretch = GetResidualStretch(mptr,dresStretch);
    double resStretch2 = resStretch*resStretch;
    double resStretch3 = resStretch2*resStretch;
    double J = J2/resStretch3;
	Tensor B = Btrial;
	B.xx /= resStretch2;
	B.yy /= resStretch2;
	B.zz /= resStretch2;
	B.xy /= resStretch2;
	if(np==THREED_MPM)
	{	B.xz /= resStretch2;
		B.yz /= resStretch2;
	}
	
	// Get hydrostatic stress component in subroutine
    //double dresStretch3 = dresStretch*dresStretch*dresStretch;
	//detdF /= dresStretch3;
    
    // Others constants
    double J23 = pow(J, 2./3.);
    
    
    // Checking for plastic loading
    // ============================================
    
    // Get magnitude of the deviatoric stress tensor
    // ||s|| = sqrt(s.s)
    
	// Set alpint for particle
	
    
    // these will be needed for elastic or plasti
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    
    //============================
    //  TEST
    //============================
    	
		// save on particle
        *pB = Btrial;
        
        // Get specifique stress i.e. (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
		// JAN: Just store deviatoric stress
    
        sp->xx = 1./3.*G1sp*(2.*pB->xx-pB->yy-pB->zz)/J23;
        sp->yy = 1./3.*G1sp*(2.*pB->yy-pB->xx-pB->zz)/J23;
        sp->zz = 1./3.*G1sp*(2.*pB->zz-pB->xx-pB->yy)/J23;
        sp->xy = G1sp*pB->xy/J23;
    if(np==THREED_MPM)
    {   sp->xz = G1sp*pB->xz/J23;
        sp->yz = G1sp*pB->yz/J23;
    }
    
    
        //return;
    
    
    
}
#pragma mark HEAnisotropic::Custom Methods

#pragma mark HEAnisotropic::Accessors

// Return the material tag
int HEAnisotropic::MaterialTag(void) { return HEANISOTROPIC; }

// return unique, short name for this material
const char *HEAnisotropic::MaterialType(void) { return "Hyperelastic Anisotropic"; }

// calculate wave speed in mm/sec - needs to be implemented
double HEAnisotropic::WaveSpeed(bool threeD,MPMBase *mptr) { return 1.e-12; }



//double pi = 3.14159265358979;
//double THETA0_Rad= THETA0*pi/180.;

//   INPUT SOLLOCITATION

//  Restrained Tensile
// Matrix3 F(2.,0.,0.,0.,2.,0.,0,0,1.);

//  Simple Shear
/*
 Matrix3 F_GxG(1.,2.,0.,0.,1.,0.,0,0,1.);
 
 //
 Matrix3 FT_GxG = F_GxG.Transpose();
 Matrix3 C_GxG = FT_GxG*F_GxG;
 
 // Right Stretch Tensor - Calculation with CG
 Matrix3 Eigenvals = C_GxG.Eigenvalues();
 Matrix3 U_GxG = C_GxG.RightStretch(Eigenvals);
 
 // Rotation Tensor - Calculation with FG
 Matrix3 R_GxG = F_GxG.Rotation(Eigenvals,U_GxG);
 
 //  Verification of the Polar Decomposition
 Matrix3 F_GxG_Verif = R_GxG*U_GxG;
 
 cout << "   Matrix C in Global Basis =    " << C_GxG <<endl;
 cout << "   Eigenvalues Matrix in the Principal Basis =    " << Eigenvals <<endl;
 cout << "   Right Stretch Tensor in GxG =    " << U_GxG <<endl;
 cout << "   Rotation in Basis G =    " << R_GxG <<endl;
 cout << "   Deformation Gradient in Basis GxG  =    " << F_GxG <<endl;
 cout << "   Verif of the Polar Decomposition F_Verif in GxG =    " << F_GxG_Verif <<endl;
 
 // Initial Orientation of the fibers in the Global Basis
 Matrix3 O_GxG(cos(THETA0_Rad),-sin(THETA0_Rad),0,sin(THETA0_Rad),sin(THETA0_Rad),0,0,0,1);
 
 cout << "   Verif THETA0 =    " << THETA0 <<endl;
 cout << "   Verif THETA0_Rad =    " << THETA0_Rad <<endl;
 cout << "   Orientation in Global Basis =    " << O_GxG <<endl;
 
 // Right Cauchy-Green in the Local Orientation of the fibers
 
 Matrix3 F_GxL = F_GxG*O_GxG;
 cout << "   Deformation Gradient Tensor in mixed GxL =    " << F_GxL <<endl;
 
 Matrix3 FT_LxG = F_GxL.Transpose();
 cout << "   Deformation Gradient Tensor Transpose, so in LxG =    " << FT_LxG <<endl;
 
 Matrix3 C_LxL = FT_LxG*F_GxL;
 cout << "   Right Cauchy-Green Tensor in Local Basis =    " << C_LxL <<endl;
 
 // Actual Elongations in the Meches Directions
 
 double lam1=sqrt(C_LxL(0,0));
 double lam2=sqrt(C_LxL(1,1));
 cout << "   Elongation 1 =    " << lam1 << "   Elongation 2 =    " << lam2 << endl;
 
 
 
 // Actual Rotation of the Meches
 
 double CTHETA12=C_LxL(0,1)/(lam1*lam2);
 double THETA12=acos(CTHETA12);
 double THETA12_DEG=THETA12*180./pi;
 cout << "   THETA12_DEG =    " << THETA12_DEG << endl;
 
 // PK2 in the LxL
 
 Matrix3 PK2(E1sp*(lam1*lam1-1),(G1sp*CTHETA12+G2sp*CTHETA12*CTHETA12+G3sp*CTHETA12*CTHETA12*CTHETA12)/(lam1*lam2),0,
 (G1sp*CTHETA12+G2sp*CTHETA12*CTHETA12+G3sp*CTHETA12*CTHETA12*CTHETA12)/(lam1*lam2), E2sp*(lam2*lam2-1),0,0,0,0);*/
/*
 Matrix3 PK2(E1*(lam1*lam1-1),(G1*CTHETA12+G2*CTHETA12*CTHETA12+G3*CTHETA12*CTHETA12*CTHETA12)/(lam1*lam2),0,
 (G1*CTHETA12+G2*CTHETA12*CTHETA12+G3*CTHETA12*CTHETA12*CTHETA12)/(lam1*lam2),E2*(lam2*lam2-1),0,0,0,0);
 cout << "   PK2 Tensor in Local Basis =    " << PK2 <<endl;
 */
/*
 
 if THETA12_DEG>90.
 PK2(1,2)=-PK2(1,2);
 end
 
 if THETA12_DEG==90.
 PK2(1,2)=0.;
 end
 PK2(2,1)=PK2(1,2);
 */
/*
 // PK2 in GxG
 Matrix3 PK2_GxG=O_GxG.Transpose()*PK2*O_GxG;
 cout << "   PK2 Tensor in Global Basis =    " << PK2_GxG <<endl;
 
 // sp in GxG
 Matrix3 sp_GxG=F_GxG*PK2_GxG*F_GxG.Transpose();
 cout << "   PK2 Tensor in Global Basis =    " << sp_GxG <<endl;
 // add new properties here
 */
