/********************************************************************************
    HEIsotropic.cpp
    NairnMPM
    
    Created by John Nairn, Sept 27, 2011.
    Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		HyperElastic.hpp (MaterialBase.hpp)
********************************************************************************/

#include "HEIsotropic.hpp"
#include "MPM_Classes/MPMBase.hpp"

#pragma mark HEIsotropic::Constructors and Destructors

// Constructors
HEIsotropic::HEIsotropic() {}

/* The default contructor should call a parent class constructor and
	then fill in any new initialization.
	*/
// Constructors
HEIsotropic::HEIsotropic(char *matName) : HyperElastic(matName)
{
        G1 = -1.;			// required
        Kbulk = -1.;		// required

    
    // inintialize unique properties
}

#pragma mark HEIsotropic::Initialization

// print mechanical properties output window

/* This method reads properties from input file that are unique to this class.
	For details see "Creating New Material Properties" section in the
	Create_MPM_Material Wiki on nairn-fea-mem google code page
*/
// Read material properties
char *HEIsotropic::InputMat(char *xName,int &input)
{
	// read properties for this material
    //if(strcmp(xName,"mypropertyname")==0)
    //{	input=DOUBLE_NUM;
    //    return((char *)&newproperty);
    //}
        input=DOUBLE_NUM;
        
        if(strcmp(xName,"G1")==0)
            return((char *)&G1);
        
        if(strcmp(xName,"G1")==0)
            return((char *)&G1);
    
        else if(strcmp(xName,"K")==0)
            return((char *)&Kbulk);
        
        else if(strcmp(xName,"alpha")==0)
            return((char *)&aI);
        
        else if(strcmp(xName,"beta")==0)
            return((char *)&betaI);        
	
    return(HyperElastic::InputMat(xName,input));
}

/* This method is called after input file is read but before the new
	material is printed to the results file. If necessary, verify that
	all new properties are valid. If not return string with a description
	of the problem. If OK, must pass on to a superclass
	NOTE: np is analysis type in case that is important
*/
// verify settings and maybe some initial calculations
const char *HEIsotropic::VerifyProperties(int np)
{
	// check properties

	// must call super class
    
    if(G1<0. || Kbulk < 0. )
		return "HEIsotropic Material needs non-negative G1 and K";
    
	// call super class

	return HyperElastic::VerifyProperties(np);
}

/* Called once at beginning (by VerifyProperties() in MaterialBase). For efficiency,
	use this method to calculate new terms that are independent of the particle
	state and thus will remain constant throughout the calculation. When done
	(or before), pass on to super class (but MaterialBase and Elastic do not need it)
*/
// Constant properties used in constitutive law
void HEIsotropic::InitialLoadMechProps(int makeSpecific,int np)
{

	hasMatProps=TRUE;
	
	// MU and Kbulk in Specific units using initial rho
	// for MPM (units N/m^2 cm^3/g)
	G1sp=G1*1.0e+06/rho;
	Ksp=Kbulk*1.0e+06/rho;
	
	// expansion coefficients
	CTE1 = 1.e-6*aI;
	CME1 = betaI*concSaturation;
	
	// nothing else needed from superclass
    
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
void HEIsotropic::PrintMechanicalProperties(void)
{	
	// call superclass here if it is not Material base
	HyperElastic::PrintMechanicalProperties();
	
        PrintProperty("G1",G1,"");
        PrintProperty("K",Kbulk,"");
        cout << endl;
        
        PrintProperty("a",aI,"");
        cout << endl;
    
	// add new properties here
}

// The base class history variable is cummulative equivalent plastic strain
//		(defined as dalpha = sqrt((2/3)||dep||))
// If super class needs to override this method, always save the first double
//		for the cumulative determinent increment.

char *HEIsotropic::MaterialData(void)
{
	double *p=new double;
	*p=1.;
	return (char *)p;
}


// If needed, a material can initialize particle state
// For example, ideal gas initializes to base line pressure
void HEIsotropic::SetInitialParticleState(MPMBase *mptr,int np)
{
    // get previous particle B
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    
    ZeroTensor(pB);
    pB->xx = pB->yy = pB->zz = 1.;
    
    //HyperElastic::SetInitialParticleState(mptr,np);
}

#pragma mark HEIsotropic:Methods

//===============================================================================================================
/*	Apply 2D constitutive law updating all needed terms for material type. Required updates are:
		stress, strain, plastic strain (all components) (stress should be a specific stress)
		rotation strain (single angle)
		strain energy, plastic energy, and dissipated energy (dissipated needed if want to couple to conduction)
	To support thermal and solvent expansion, include their effect on strains
    If there are material-related data on the particle, update them too
	dvij are (gradient rates X time increment) to give deformation gradient change
*/
//===============================================================================================================
void HEIsotropic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
        double delTime,int np)

{
    // get new deformation gradient and update strains and rotations
	double F[3][3];
	GetDeformationGrad(F,mptr,dvxx,dvyy,dvxy,dvyx,TRUE,FALSE);
	
    
    // get new deformation gradient
    //double F[3][3];
    //double detdF = GetDeformationGrad(F,mptr,dvxx,dvyy,dvxy,dvyx,TRUE,TRUE);

    // left Cauchy deformation tensor B = dF pB dF^T
    Tensor B = GetLeftCauchyTensor2D(mptr,dvxx,dvyy,dvxy,dvyx,TRUE);
    cout << "New" << "Bxx = " << B.xx << "  , " << "Byy = " << B.yy << "  , " << "Bxy = " << B.xy << "  , " << "Bzz = " << B.zz << endl;

    // left Cauchy deformation tensor B = F F^T
	//Tensor B;
    //B = GetLeftCauchyTensor2D(F);
    //cout << "old" << "Bxx = " << B.xx << "  , " << "Byy = " << B.yy << "  , " << "Bxy = " << B.xy << "  , " << "Bzz = " << B.zz << endl;
	
	// Deformation gradients and Cauchy tensor differ in plane stress and plane strain
	double J2;
	if(np==PLANE_STRAIN_MPM)
	{	J2 = B.xx*B.yy - B.xy*B.xy;
	}
	else
	{	
		// first we don't consider a 2D plane stress
	}
	
	// J as determinant of F (or sqrt root of determinant of B) normalized to residual stretch
//	double resStretch = GetResidualStretch(mptr);
    
    // Account for density change in specific stress
    // i.e.. Get (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
	
	// find (Cauchy stress)J/rho0 = (Kirchoff stress)/rho0


    
    // Incremental calculation of J det(F)=det(dF)det(pF)
  //  double J = detdF*mptr->GetHistoryDble();

    double J=sqrt(J2);
    double J53=pow(J,5./3.);
    //double Kse;
	//double Kterm = J*GetVolumetricTerms(J,&Kse);
    cout << "New" << "Bxx = " << B.xx << "  , " << "Byy = " << B.yy << "  , " << "Bxy = " << B.xy << "  , " << "Bzz = " << B.zz << endl;

    Tensor *sp=mptr->GetStressTensor();
    sp->xx = Ksp/2*(J-1./J)+1/3*G1sp*(2*B.xx-B.yy-B.zz)/J53;
    sp->yy = Ksp/2*(J-1./J)+1/3*G1sp*(2*B.yy-B.xx-B.zz)/J53;
    sp->zz = Ksp/2*(J-1./J)+1/3*G1sp*(2*B.zz-B.xx-B.yy)/J53;
    sp->xy = G1sp*B.xy/J53;
 
  // strain energy per unit mass (U/(rho0 V0)) and we are using
    // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
	double I1bar = (B.xx+B.yy+B.zz)/J53;
    mptr->SetStrainEnergy(0.5*(G1sp*(I1bar-3.) + Kbulk*(1/2*(J*J-1)-log(J))));
}


//===============================================================================================================
/* Apply 3D constitutive law updating all needed terms for material type. Required updates are:
		stress, strain, plastic strain (all components) (stress should be a specific stress)
		rotation strain
		strain energy, plastic energy, and dissipated energy (dissipated needed if want to couple to conduction)
	To support thermal and solvent expansion, include their effect on strains
    If there are material-related data on the particle, update them too
    dvij are (gradient rates X time increment) to give deformation gradient change
*/
//===============================================================================================================


// 3D Stress Calculation New Law
void HEIsotropic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
                              double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)

{

    
    
// left Cauchy deformation tensor B = dF pB dF^T
	Tensor B = GetLeftCauchyTensor3D(mptr,dvxx,dvyy,dvzz,dvxy,dvyx,dvxz,dvzx,dvyz,dvzy,TRUE);
    
    // get new deformation gradient
	double F[3][3];
	double detdF = GetDeformationGrad(F,mptr,dvxx,dvyy,dvzz,dvxy,dvyx,dvxz,dvzx,dvyz,dvzy,TRUE,TRUE);
    
    
    // Incremental calculation of J det(F)=det(dF)det(pF)
    double J = detdF*mptr->GetHistoryDble();
    double J53=pow(J,5./3.);

 Tensor *sp=mptr->GetStressTensor();
 sp->xx = Ksp/2*(J-1./J) + 1/3*G1sp*(2*B.xx-B.yy-B.zz)/J53;
 sp->yy = Ksp/2*(J-1./J) + 1/3*G1sp*(2*B.yy-B.xx-B.zz)/J53;
 sp->zz = Ksp/2*(J-1./J) + 1/3*G1sp*(2*B.zz-B.xx-B.yy)/J53;
 sp->xy = G1sp*B.xy/J53;
 sp->xz = G1sp*B.xz/J53;
 sp->yz = G1sp*B.yz/J53;
}
// Return normal stress term (due to bulk modulus) and twice the pressure term for strain energy.
// Each block of lines is for a different U(J).
// Any change here must also be made in 2D MPMConstLaw for the numerical solution to find B.zz in plane stress
double HEIsotropic::GetVolumetricTerms(double J,double *Kse)
{
    // This is for *Kse/2 = U(J) = (K/2)(J-1)^2
    //double Kterm = Ksp*(J-1.);
    //*Kse = Kterm*(J-1);
    //return Kterm;             // = Ksp*(J-1)
    
    // This is for for *Kse/2 = U(J) = (K/2)(ln J)^2
    // Zienkiewicz & Taylor recommend not using this one
    //double lj = log(J);
    //double Kterm =Ksp*lj;
    //*Kse = Kterm*lj;
    //return Kterm/J;           // = Ksp*(ln J)/J
    
    // This is for *Kse/2 = U(J) = (K/2)((1/2)(J^2-1) - ln J)
    // Zienkiewicz & Taylor note that stress is infinite as J->0 and J->infinity for this function, while others are not
    // Simo and Hughes also use this form (see Eq. 9.2.3)
    *Kse = Ksp*(0.5*(J*J-1.)-log(J));
    return 0.5*Ksp*(J - 1./J);      // = (Ksp/2)*(J - 1/J)
}

// Return normal stress term (due to bulk modulus) and twice the pressure term for strain energy.
// Each block of lines is for a different U(J).
// Any change here must also be made in 2D MPMConstLaw for the numerical solution to find B.zz in plane stress



#pragma mark HEIsotropic::Custom Methods

#pragma mark HEIsotropic::Accessors

// Return the material tag
int HEIsotropic::MaterialTag(void) { return HEISOTROPIC; }

// return unique, short name for this material
const char *HEIsotropic::MaterialType(void) { return "Hyperelastic isotropic"; }

// calculate wave speed in mm/sec - needs to be implemented
double HEIsotropic::WaveSpeed(bool threeD,MPMBase *mptr)

//{ return 1.e-12; }

{
    return sqrt(1.e9*(Kbulk+4.*G1/3.)/rho);
}
