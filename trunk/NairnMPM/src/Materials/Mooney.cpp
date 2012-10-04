/********************************************************************************
    Mooney.cpp
    NairnMPM
    
    Created by John Nairn on Fri Feb 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/Mooney.hpp"
#include "MPM_Classes/MPMBase.hpp"
 
#pragma mark Mooney::Constructors and Destructors

// Constructors
Mooney::Mooney()
{
}
// Constructors with arguments 
Mooney::Mooney(char *matName) : HyperElastic(matName)
{
	G1 = -1.;			// required
	G2 = 0.;			// zero is neo-Hookean
	Kbulk = -1.;		// required
}

#pragma mark Mooney::Initialization

// print mechanical properties output window
void Mooney::PrintMechanicalProperties(void)
{
	PrintProperty("G1",G1,"");
	PrintProperty("G2",G2,"");
	PrintProperty("K",Kbulk,"");
	cout << endl;
	
	PrintProperty("a",aI,"");
	cout << endl;
}
	
// Read material properties
char *Mooney::InputMat(char *xName,int &input)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"G1")==0)
        return((char *)&G1);
    
    else if(strcmp(xName,"G2")==0)
        return((char *)&G2);
    
    else if(strcmp(xName,"K")==0)
        return((char *)&Kbulk);
    
    else if(strcmp(xName,"alpha")==0)
        return((char *)&aI);
    
    else if(strcmp(xName,"beta")==0)
        return((char *)&betaI);
    
    return(HyperElastic::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *Mooney::VerifyProperties(int np)
{
    if(G1<0. || Kbulk < 0. || G2<0.)
		return "Mooney-Rivlin Hyperelastic material needs non-negative G1, G2, and K";

	// call super class
	return MaterialBase::VerifyProperties(np);
}

// Private properties used in constitutive law
void Mooney::InitialLoadMechProps(int makeSpecific,int np)
{
	hasMatProps=TRUE;
	
	// G1, G2, and Kbulk in Specific units
	// for MPM (units N/m^2 cm^3/g)
	G1sp=G1*1.0e+06/rho;
	G2sp=G2*1.0e+06/rho;
	Ksp=Kbulk*1.0e+06/rho;
	
	// expansion coefficients
	CTE1 = 1.e-6*aI;
	CME1 = betaI*concSaturation;
	
	// nothing else needed from superclass
}

#pragma mark Mooney::Methods

/* For 2D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
	Does not support thermal or moisture strains
 
    Note the incompressible rubber are not suitable for dynamic calculations because
    they have an infinite wave speed
*/
void Mooney::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
								double delTime,int np)
{
	// get new deformation gradient and update strains and rotations
	double F[3][3];
	GetDeformationGrad(F,mptr,dvxx,dvyy,dvxy,dvyx,TRUE,FALSE);
	
	// left Cauchy deformation tensor B = F F^T
	Tensor B = GetLeftCauchyTensor2D(F);
	
	// Deformation gradients and Cauchy tensor differ in plane stress and plane strain
	double J2;
	if(np==PLANE_STRAIN_MPM)
	{	J2 = B.xx*B.yy - B.xy*B.xy;
	}
	else
	{	// Find B.zz required to have zero stress in z direction
		
		// fixed arguments
		double arg = B.xx*B.yy - B.xy*B.xy;
		double arg12 = sqrt(arg);
		double arg16 = pow(arg,1./6.);
		double arg2 = B.xx+B.yy;
		
		// Newton-Rapheson starting at B.zz = 1
		// In tests finds answer in 3 or less steps
		Tensor *ep=mptr->GetStrainTensor();
		double xn16,xn12,xnp1,xn = 1.+ep->zz;
		double fx,fxp;
		int iter=1;
		while(iter<10)
		{	xn16 = pow(xn,1./6.);
			xn12 = sqrt(xn);
			fx = 3.*Ksp*xn*arg*(xn12*arg12-1.) + G1sp*(2.*xn-arg2)*xn16*arg16
				+ G2sp*(xn*arg2-2.*arg)/(xn16*arg16);
			fxp = 3.*Ksp*arg*(1.5*arg12*xn12-1.) + G1sp*xn16*arg16*(14.*xn-arg2)/(6.*xn)
				+ G2sp*(2.*arg+5.*xn*arg2)/(6.*xn16*arg16*xn);
			xnp1 = xn - fx/fxp;
			//cout << iter << ": " << xn << "," << xnp1 << "," << fabs(xn-xnp1) << endl;
			if(fabs(xn-xnp1)<1e-6)
			{	B.zz = xn;
				break;
			}
			xn = xnp1;
			iter+=1;
		}
		
		// particle strain
		F[2][2] =  sqrt(B.zz);					// 1 + dw/dz
		ep->zz = F[2][2] - 1.;
		J2 = B.zz*arg;
	}
	
	// J as determinant of F (or sqrt root of determinant of B) normalized to residual stretch
	double resStretch = GetResidualStretch(mptr);
	double J = sqrt(J2)/(resStretch*resStretch*resStretch);
    double Kse;
#ifdef CONSTANT_RHO
    // Ignore change in density in the specific stress
    // i.e., get Cauchy stress / rho0
	double J53 = pow(J, 5./3.);
	double J23 = J53/J;
	double J73 = J53*J23;
	double J43 = J23*J23;
	double Kterm = GetVolumetricTerms(J,&Kse);
#else
    // Account for density change in specific stress
    // i.e.. Get (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
	double J23 = pow(J, 2./3.);
	double J43 = J23*J23;
	double J53 = J23;                 // Using J^(2/3) = J^(5/3)/J
	double J73 = J43;                 // Using J^(4/3) = J^(7/3)/J
	double Kterm = J*GetVolumetricTerms(J,&Kse);
#endif
	
	// find (Cauchy stress)/rho0 (if CONSTANT_RHO) or (Kirchoff stress)/rho0 (if not)
	Tensor *sp=mptr->GetStressTensor();
	sp->xx = Kterm + (2*B.xx-B.yy-B.zz)*G1sp/(3.*J53)
			+ (B.xx*(B.yy+B.zz)-2*B.yy*B.zz-B.xy*B.xy)*G2sp/(3.*J73);
	sp->yy = Kterm + (2*B.yy-B.xx-B.zz)*G1sp/(3.*J53)
			+ (B.yy*(B.xx+B.zz)-2*B.xx*B.zz-B.xy*B.xy)*G2sp/(3.*J73);
	sp->xy = B.xy*G1sp/J53 + (B.zz*B.xy)*G2sp/J73;
	if(np==PLANE_STRAIN_MPM)
	{	sp->zz = Kterm + (2*B.zz-B.xx-B.yy)*G1sp/(3.*J53)
				+ (B.zz*(B.xx+B.yy)-2*B.xx*B.yy+2.*B.xy*B.xy)*G2sp/(3.*J73);
	}
	
	// strain energy (total energy divided by initial rho)
	double I1bar = (B.xx+B.yy+B.zz)/J23;
	double I2bar = 0.5*(I1bar*I1bar - (B.xx*B.xx+B.yy*B.yy+B.zz*B.zz+2.*B.xy*B.xy)/J43);
    mptr->SetStrainEnergy(0.5*(G1sp*(I1bar-3.) + G2sp*(I2bar-3.) + Kse));
}

/* For 3D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
 
    Note the incompressible rubber are not suitable for dynamic calculations because
    they have an infinite wave speed
*/
void Mooney::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
						  double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
	// get new deformation gradient
	double F[3][3];
	GetDeformationGrad(F,mptr,dvxx,dvyy,dvzz,dvxy,dvyx,dvxz,dvzx,dvyz,dvzy,TRUE,FALSE);

	// left Cauchy deformation tensor B = F F^T
	Tensor B = GetLeftCauchyTensor3D(F);
	
	// J as determinant of F (or sqrt root of determinant of B), normalized to residual stretch
	double J2 = B.xx*B.yy*B.zz + 2.*B.xy*B.xz*B.yz - B.yz*B.yz*B.xx - B.xz*B.xz*B.yy - B.xy*B.xy*B.zz;
	double resStretch = GetResidualStretch(mptr);
	double J = sqrt(J2)/(resStretch*resStretch*resStretch);
    double Kse;
#ifdef CONSTANT_RHO
    // Ignore change in density in the specific stress
    // i.e., get Cauchy stress / rho0
	double J53 = pow(J, 5./3.);
	double J23 = J53/J;
	double J43 = J23*J23;
	double J73 = J43*J;
	double Kterm = GetVolumetricTerms(J,&Kse);
#else
    // Account for density change in specific stress
    // i.e.. Get (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
	double J23 = pow(J, 2./3.);
	double J43 = J23*J23;
	double J53 = J23;                 // Using J^(2/3) = J^(5/3)/J
	double J73 = J43;                 // Using J^(4/3) = J^(7/3)/J
	double Kterm = J*GetVolumetricTerms(J,&Kse);
#endif
	
	// Find Cauchy stresses
	//double Kterm = Ksp*log(J)/J;
	Tensor *sp=mptr->GetStressTensor();
	sp->xx = Kterm + (2*B.xx-B.yy-B.zz)*G1sp/(3.*J53)
			+ (B.xx*(B.yy+B.zz)-2*B.yy*B.zz-B.xy*B.xy-B.xz*B.xz+2.*B.yz*B.yz)*G2sp/(3.*J73);
	sp->yy = Kterm + (2*B.yy-B.xx-B.zz)*G1sp/(3.*J53)
			+ (B.yy*(B.xx+B.zz)-2*B.xx*B.zz-B.xy*B.xy+2.*B.xz*B.xz-B.yz*B.yz)*G2sp/(3.*J73);
	sp->zz = Kterm + (2*B.zz-B.xx-B.yy)*G1sp/(3.*J53)
			+ (B.zz*(B.xx+B.yy)-2*B.xx*B.yy+2.*B.xy*B.xy-B.xz*B.xz-B.yz*B.yz)*G2sp/(3.*J73);
	sp->xy = B.xy*G1sp/J53 + (B.zz*B.xy-B.xz*B.yz)*G2sp/J73;
	sp->xz = B.xz*G1sp/J53 + (B.yy*B.xz-B.xy*B.yz)*G2sp/J73;
	sp->yz = B.yz*G1sp/J53 + (B.xx*B.yz-B.xy*B.xz)*G2sp/J73;
    
	// strain energy
	double I1bar = (B.xx+B.yy+B.zz)/J23;
	double I2bar = 0.5*(I1bar*I1bar - (B.xx*B.xx+B.yy*B.yy+B.zz*B.zz+2.*B.xy*B.xy+2*B.xz*B.xz+2.*B.yz*B.yz)/J43);
    mptr->SetStrainEnergy(0.5*(G1sp*(I1bar-3.) + G2sp*(I2bar-3.) + Kse));
}

// Return normal stress term (due to bulk modulus) and twice the pressure term for strain energy
// Each block of lines if for a different U(J)
double Mooney::GetVolumetricTerms(double J,double *Kse)
{
    // This is for U(J) = (K/2)(J-1)^2
    //double Kterm = Ksp*(J-1.);
    //*Kse = Kterm*(J-1);     // = Ksp*(J-1)^2
    //return Kterm;
    
    // This is for for U(J) = (K/2)(ln J)^2
    // Zienkiewicz & Taylor recommend not using this one
    //double lj = log(J);
    //double Kterm =Ksp*log(J);
    //*Kse = Kterm*lj;        // = Ksp*(ln J)^2
    //return Kterm/J;         // = Ksp*(ln J)/J
    
    // This is for U(J) = (K/2)((1/2)(J^2-1) - ln J)
    // Zienkiewicz & Taylor note that stress goes to infinite as J=0 and J->infinity for this function, while others do not
    // Simo and Hughes also use this form (see Eq. 9.2.3)
    *Kse = Ksp*(0.5*(J*J-1.)-log(J));
    return 0.5*Ksp*(J - 1./J);
}

#pragma mark Mooney::Accessors

// Return the material tag
int Mooney::MaterialTag(void) { return MOONEYRIVLIN; }

/* Calculate wave speed in mm/sec (because G in MPa and rho in g/cm^3)
	Uses sqrt((K +4G/3)/rho) which is dilational wave speed
	at low strain G = G1+G2
*/
double Mooney::WaveSpeed(bool threeD,MPMBase *mptr)
{
    return sqrt(1.e9*(Kbulk+4.*(G1+G2)/3.)/rho);
}

/* Calculate shear wave speed in mm/sec (because G1 and G2 in MPa and rho in g/cm^3)
	at low strain G = G1+G2
*/
double Mooney::ShearWaveSpeed(bool threeD,MPMBase *mptr) { return sqrt(1.e9*(G1+G2)/rho); }

// return material type
const char *Mooney::MaterialType(void) { return "Mooney-Rivlin Hyperelastic"; }


