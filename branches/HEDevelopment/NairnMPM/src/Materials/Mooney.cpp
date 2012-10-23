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

char *Mooney::MaterialData(void)
{
	double *p=new double;
	*p=1.;
	return (char *)p;
}


void Mooney::SetInitialParticleState(MPMBase *mptr,int np)
{
    // get previous particle B
    Tensor *pB = mptr->GetElasticLeftCauchyTensor();
    
    ZeroTensor(pB);
    pB->xx = pB->yy = pB->zz = 1.;
    
    //HyperElastic::SetInitialParticleState(mptr,np);
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
	
	// G1, G2, and Kbulk in Specific units using initial rho
	// for MPM (units N/m^2 cm^3/g)
	G1sp=G1*1.0e+06/rho;
	G2sp=G2*1.0e+06/rho;
	Ksp=Kbulk*1.0e+06/rho;
    cout << "Ksp haut = " << Ksp << endl;
	
	// expansion coefficients
	CTE1 = 1.e-6*aI;
	CME1 = betaI*concSaturation;
	
	// nothing else needed from superclass
    
    // If needed, a material can initialize particle state
    // For example, ideal gas initializes to base line pressure
  
}



#pragma mark Mooney::Methods

/* For 2D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy (per unit mass)
    dvij are (gradient rates X time increment) to give deformation gradient change
 
    Note that incompressible rubber is not suitable for dynamic calculations because
    they have an infinite wave speed
*/
void Mooney::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
								double delTime,int np)
{
	// get new deformation gradient and update strains and rotations
	double F[3][3];
	GetDeformationGrad(F,mptr,dvxx,dvyy,dvxy,dvyx,TRUE,FALSE);
    	    
    // get new deformation gradient
    //   double F[3][3];
    //double detdF = GetDeformationGrad(F,mptr,dvxx,dvyy,dvxy,dvyx,TRUE,TRUE);
    //cout << "detdF dans Mooney = " << detdF << endl;
    // Incremental calculation of J det(F)=det(dF)det(pF)
    //double JJ = detdF*mptr->GetHistoryDble();
    //cout << "New JJ = " << JJ << endl;

	
    Tensor B;
	// left Cauchy deformation tensor B = F F^T
    //B = GetLeftCauchyTensor2D(F);
    //cout << "Bxx = " << B.xx << "  , " << "Byy = " << B.yy << "  , " << "Bxy = " << B.xy << "  , " << "Bzz = " << B.zz << endl;
    
    // left Cauchy deformation tensor B = dF pB dF^T
    B = GetLeftCauchyTensor2D(mptr,dvxx,dvyy,dvxy,dvyx,TRUE);
    //cout << dvxx << endl;
    cout << "New" << "Bxx = " << B.xx << "  , " << "Byy = " << B.yy << "  , " << "Bxy = " << B.xy << "  , " << "Bzz = " << B.zz << endl;
	
	// Deformation gradients and Cauchy tensor differ in plane stress and plane strain
	double J2;
	if(np==PLANE_STRAIN_MPM)
	{
        J2 = B.xx*B.yy - B.xy*B.xy;
        // get new deformation gradient
        //   double F[3][3];
        //double detdF = GetDeformationGrad(F,mptr,dvxx,dvyy,dvxy,dvyx,TRUE,TRUE);
        //cout << "detdF dans Mooney = " << detdF << endl;
        // Incremental calculation of J det(F)=det(dF)det(pF)
        //double J2 = detdF*mptr->GetHistoryDble();
        //cout << "New J = " << J2 << endl;

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
        
        // solution for B.zz in xn
		while(iter<10)
		{	xn16 = pow(xn,1./6.);
			xn12 = sqrt(xn);
            
            // three options for K term in function
            // For K(J-1) term
			//fx = 3.*Ksp*xn*arg*(xn12*arg12-1.); 
			//fxp = 3.*Ksp*arg*(1.5*arg12*xn12-1.); 
            
            // For K(ln J)/J term
            //fx = 1.5*Ksp*xn12*arg12*log(xn*arg);
            //fxp = 0.75*Ksp*arg12*(2.+log(arg*xn))/xn12;
            
            // For (K/2)(J-1/J) term
            fx = 1.5*Ksp*xn12*arg12*(xn*arg-1.);            
            fxp = 0.75*Ksp*arg12*(3.*arg*xn-1.)/xn12;
           
            // N ow add the shear terms
            fx += G1sp*(2.*xn-arg2)*xn16*arg16 + G2sp*(xn*arg2-2.*arg)/(xn16*arg16);
            fxp += G1sp*xn16*arg16*(14.*xn-arg2)/(6.*xn) + G2sp*(2.*arg+5.*xn*arg2)/(6.*xn16*arg16*xn);
            
            // new prediction for solution
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
    //double J = J2;
    //cout << "old J = " << J <<","<< "old J = " << JJ << endl;
    mptr->SetHistoryDble(J);
    
    // Account for density change in specific stress
    // i.e.. Get (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
	double J23 = pow(J, 2./3.);
	double J43 = J23*J23;
	double JforG1 = J23;                   // J^(2/3) to get J*(Cauchy Stress) below
	double JforG2 = J43;                   // J^(4/3) to get J*(Cauchy Stress) below
    double Kse;
	double Kterm = J*GetVolumetricTerms(J,&Kse);
    //double Kterm=1/2*Ksp*(J-1./J);
	//cout << "Ksp = " << Ksp << endl;
	
    // find (Cauchy stress)J/rho0 = (Kirchoff stress)/rho0
	Tensor *sp=mptr->GetStressTensor();

	//cout << "old B  " <<"Bxx = " << B.xx << "Bzz = " << B.zz << ","<< "Byy = " << B.yy << endl;

    sp->xx = Kterm + (2*B.xx-B.yy-B.zz)*G1sp/(3.*JforG1)
			+ (B.xx*(B.yy+B.zz)-2*B.yy*B.zz-B.xy*B.xy)*G2sp/(3.*JforG2);
    //sp->xx = Kterm+1/3*G1sp*(2*B.xx-B.yy-B.zz)/J23;
    //sp->xx = Kterm+1/3*G1sp*(2*B.xx-B.yy-B.zz)/J23;
	sp->yy = Kterm + (2*B.yy-B.xx-B.zz)*G1sp/(3.*JforG1)
			+ (B.yy*(B.xx+B.zz)-2*B.xx*B.zz-B.xy*B.xy)*G2sp/(3.*JforG2);
	sp->xy = B.xy*G1sp/JforG1 + (B.zz*B.xy)*G2sp/JforG2;
	if(np==PLANE_STRAIN_MPM)
	{	sp->zz = Kterm + (2*B.zz-B.xx-B.yy)*G1sp/(3.*JforG1)
				+ (B.zz*(B.xx+B.yy)-2*B.xx*B.yy+2.*B.xy*B.xy)*G2sp/(3.*JforG2);
	}
    //sp->xx = 0.5*Ksp*(J*J-1)/J+1/3*G1sp*(2*B.xx-B.yy-B.zz)/J53;
    //sp->yy = 0.5*Ksp*(J-1/J)/J+1/3*G1sp*(2*B.yy-B.xx-B.zz)/J53;
    //sp->zz = 0.5*Ksp*(J-1/J)/J+1/3*G1sp*(2*B.zz-B.xx-B.yy)/J53;
    //sp->xy = G1sp*B.xy/J53;
    //cout << "Bxx = " << B.xx << ","<< "Byy = " << B.yy << endl;
	
	// strain energy per unit mass (U/(rho0 V0)) and we are using
    // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
	double I1bar = (B.xx+B.yy+B.zz)/J23;
	double I2bar = 0.5*(I1bar*I1bar - (B.xx*B.xx+B.yy*B.yy+B.zz*B.zz+2.*B.xy*B.xy)/J43);
    mptr->SetStrainEnergy(0.5*(G1sp*(I1bar-3.) + G2sp*(I2bar-3.) + Kse));
}

/* For 3D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy (per unit mass)
	dvij are (gradient rates X time increment) to give deformation gradient change
 
    Note that incompressible rubber is not suitable for dynamic calculations because
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
    
    // Account for density change in specific stress
    // i.e.. Get (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
	double J23 = pow(J, 2./3.);
	double J43 = J23*J23;
	double JforG1 = J23;                 // J^(2/3) = J^(5/3)/J to get Kirchoff stress
	double JforG2 = J43;                 // J^(4/3) = J^(7/3)/J to get Kirchoff stress
    double Kse;
	double Kterm = J*GetVolumetricTerms(J,&Kse);       // times J to get Kirchoff stress
	
	// find (Cauchy stress)J/rho0 = (Kirchoff stress)/rho0
	Tensor *sp=mptr->GetStressTensor();
	sp->xx = Kterm + (2*B.xx-B.yy-B.zz)*G1sp/(3.*JforG1)
			+ (B.xx*(B.yy+B.zz)-2*B.yy*B.zz-B.xy*B.xy-B.xz*B.xz+2.*B.yz*B.yz)*G2sp/(3.*JforG2);
	sp->yy = Kterm + (2*B.yy-B.xx-B.zz)*G1sp/(3.*JforG1)
			+ (B.yy*(B.xx+B.zz)-2*B.xx*B.zz-B.xy*B.xy+2.*B.xz*B.xz-B.yz*B.yz)*G2sp/(3.*JforG2);
	sp->zz = Kterm + (2*B.zz-B.xx-B.yy)*G1sp/(3.*JforG1)
			+ (B.zz*(B.xx+B.yy)-2*B.xx*B.yy+2.*B.xy*B.xy-B.xz*B.xz-B.yz*B.yz)*G2sp/(3.*JforG2);
	sp->xy = B.xy*G1sp/JforG1 + (B.zz*B.xy-B.xz*B.yz)*G2sp/JforG2;
	sp->xz = B.xz*G1sp/JforG1 + (B.yy*B.xz-B.xy*B.yz)*G2sp/JforG2;
	sp->yz = B.yz*G1sp/JforG1 + (B.xx*B.yz-B.xy*B.xz)*G2sp/JforG2;
    
	// strain energy per unit mass (U/(rho0 V0)) and we are using
    // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
	double I1bar = (B.xx+B.yy+B.zz)/J23;
	double I2bar = 0.5*(I1bar*I1bar - (B.xx*B.xx+B.yy*B.yy+B.zz*B.zz+2.*B.xy*B.xy+2*B.xz*B.xz+2.*B.yz*B.yz)/J43);
    mptr->SetStrainEnergy(0.5*(G1sp*(I1bar-3.) + G2sp*(I2bar-3.) + Kse));
}

// Return normal stress term (due to bulk modulus) and twice the pressure term for strain energy.
// Each block of lines is for a different U(J).
// Any change here must also be made in 2D MPMConstLaw for the numerical solution to find B.zz in plane stress
double Mooney::GetVolumetricTerms(double J,double *Kse)
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


