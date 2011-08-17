/********************************************************************************
    Mooney.cpp
    NairnMPM
    
    Created by John Nairn on Fri Feb 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/Mooney.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
 
#pragma mark Mooney::Constructors and Destructors

// Constructors
Mooney::Mooney()
{
}
// Constructors with arguments 
Mooney::Mooney(char *matName) : RubberElastic(matName)
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
    
    return(MaterialBase::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *Mooney::VerifyProperties(int np)
{
    if(G1<0. || Kbulk < 0.)
		return "Mooney-Rivlin material needs non-negative G1 and K";

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
*/
void Mooney::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
								double delTime,int np)
{
	// get new doformation gradient from current one using dF.F where dF = I + gradV * dt and F is current
	// deformation gradient (found from current strains and rotations)
	Tensor *ep=mptr->GetStrainTensor();
	TensorAntisym *wrot = mptr->GetRotationStrainTensor();
	double Fxx = 1. + ep->xx;;
	double Fxy = (ep->xy - wrot->xy)/2.;
	double Fyx = (ep->xy + wrot->xy)/2.;
	double Fyy = 1. + ep->yy;
	
	// get new 2D deformation gradient
	double F[3][3];
	F[0][0] = (1. + dvxx)*Fxx + dvxy*Fyx;		// 1 + du/dx
	F[0][1] = (1. + dvxx)*Fxy + dvxy*Fyy;		// du/dy
	//F[0][2] = 0.;								// du/dz
	F[1][0] = dvyx*Fxx + (1. + dvyy)*Fyx;		// dv/dx
	F[1][1] = dvyx*Fxy + (1. + dvyy)*Fyy;		// 1 + dv/dy
	//F[1][2] = 0.;							// dv/dz
	//F[2][0] = 0.;							// dw/dx
	//F[2][1] = 0.;							// dw/dy
	//F[2][2] = ? ;							// 1 + dw/dz (done later)
	
	// store in total strain and rotation tensors
    ep->xx = F[0][0] - 1.;
    ep->yy = F[1][1] - 1.;
    ep->xy = F[1][0] + F[0][1];
	
	// rotational strain increments
	wrot->xy = F[1][0] - F[0][1];
	
    // total residual stretch (1 + alpha dT + beta csat dConcentration)
	double resStretch = 1.0;
	double dTemp=mptr->pPreviousTemperature-thermal.reference;
	resStretch += CTE1*dTemp;
	if(DiffusionTask::active)
	{	double dConc=mptr->pPreviousConcentration-DiffusionTask::reference;
		resStretch += CME1*dConc;
	}
	
	// left Caucy deformation tensor B = F F^T
	Tensor B;
	ZeroTensor(&B);
	B.xx = F[0][0]*F[0][0] + F[0][1]*F[0][1];
	B.yy = F[1][0]*F[1][0] + F[1][1]*F[1][1];
	B.xy = F[0][0]*F[1][0] + F[0][1]*F[1][1];
	//B.zz = ? ;   done in next section
	
	// Deformation gradiates and Cauchy tensor differ in plane stress and plane strain
	double J2;
	if(np==PLANE_STRAIN_MPM)
	{	F[2][2] =  1.;					// 1 + dw/dz
		B.zz = 1.;
		J2 = B.xx*B.yy - B.xy*B.xy;
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
	
	// J as determinant of F (or sqrt root of determinant of B)
	double J = sqrt(J2)/(resStretch*resStretch*resStretch);
	double J53 = pow(J, 5./3.);
	double J73 = pow(J, 7./3.);
	
	// find Cauchy stresses
	Tensor *sp=mptr->GetStressTensor();
	sp->xx = Ksp*(J-1.) + (2*B.xx-B.yy-B.zz)*G1sp/(3.*J53)
			+ (B.xx*(B.yy+B.zz)-2*B.yy*B.zz-B.xy*B.xy)*G2sp/(3.*J73);
	sp->yy = Ksp*(J-1.) + (2*B.yy-B.xx-B.zz)*G1sp/(3.*J53)
			+ (B.yy*(B.xx+B.zz)-2*B.xx*B.zz-B.xy*B.xy)*G2sp/(3.*J73);
	sp->xy = B.xy*G1sp/J53 + (B.zz*B.xy)*G2sp/J73;
	if(np==PLANE_STRAIN_MPM)
	{	sp->zz = Ksp*(J-1.) + (2*B.zz-B.xx-B.yy)*G1sp/(3.*J53)
				+ (B.zz*(B.xx+B.yy)-2*B.xx*B.yy+2.*B.xy*B.xy)*G2sp/(3.*J73);
	}
	
	// Convert to Kirchoff Stress (time J)
	/*
	 sp->xx *= J;
	 sp->yy *= J;
	 sp->zz *= J;
	 sp->xy *= J;
	 */
	
	// Convert to nominal stress by premultiply with J F^(-1)
	// JFi is transpose of matrix of cofactors for F
	/*
	 double JFi[3][3];
	 JFi[0][0] = F[1][1]*F[2][2];
	 JFi[0][1] = -F[0][1]*F[2][2];
	 JFi[0][2] = 0.;
	 JFi[1][0] = -F[1][0]*F[2][2];
	 JFi[1][1] = F[0][0]*F[2][2];
	 JFi[1][2] = 0.;
	 JFi[2][0] = 0.;
	 JFi[2][1] = 0.;
	 JFi[2][2] = F[0][0]*F[1][1] - F[1][0]*F[0][1];
	 Tensor sp0=*sp;
	 sp->xx = JFi[0][0]*sp0.xx + JFi[0][1]*sp0.xy;
	 sp->xy = JFi[0][0]*sp0.xy + JFi[0][1]*sp0.yy;
	 sp->yy = JFi[1][0]*sp0.xy + JFi[1][1]*sp0.yy;
	 sp->zz = JFi[2][2]*sp0.zz;
	 */
	
	// strain energy
	double J23 = J53/J;
	double J43 = J23*J23;
	double I1bar = (B.xx+B.yy+B.zz)/J23;
	double I2bar = 0.5*(I1bar*I1bar - (B.xx*B.xx+B.yy*B.yy+B.zz*B.zz+2.*B.xy*B.xy)/J43);
    mptr->SetStrainEnergy(0.5*(G1sp*(I1bar-3.) + G2sp*(I2bar-3.) + Ksp*(J-1.)*(J-1.)));
}

/* For 3D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
*/
void Mooney::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
						  double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
	// get new doformation gradient from current one using dF.F where dF = I + gradV * dt and F is current
	// deformation gradient (stored int current strains and rotations)
	Tensor *ep=mptr->GetStrainTensor();
	TensorAntisym *wrot = mptr->GetRotationStrainTensor();
	double Fxx = 1. + ep->xx;;
	double Fxy = (ep->xy - wrot->xy)/2.;
	double Fxz = (ep->xz - wrot->xz)/2.;
	double Fyx = (ep->xy + wrot->xy)/2.;
	double Fyy = 1. + ep->yy;
	double Fyz = (ep->yz - wrot->yz)/2.;
	double Fzx = (ep->xz + wrot->xz)/2.;
	double Fzy = (ep->yz + wrot->yz)/2.;
	double Fzz = 1. + ep->zz;

	// get new deformation gradient
	double F[3][3];
	F[0][0] = (1. + dvxx)*Fxx + dvxy*Fyx + dvxz*Fzx;		// 1 + du/dx
	F[0][1] = (1. + dvxx)*Fxy + dvxy*Fyy + dvxz*Fzy;		// du/dy
	F[0][2] = (1. + dvxx)*Fxz + dvxy*Fyz + dvxz*Fzz;		// du/dz
	F[1][0] = dvyx*Fxx + (1. + dvyy)*Fyx + dvyz*Fzx;		// dv/dx
	F[1][1] = dvyx*Fxy + (1. + dvyy)*Fyy + dvyz*Fzy;		// 1 + dv/dy
	F[1][2] = dvyx*Fxz + (1. + dvyy)*Fyz + dvyz*Fzz;		// dv/dz
	F[2][0] = dvzx*Fxx + dvzy*Fyx + (1. + dvzz)*Fzx;		// dw/dx
	F[2][1] = dvzx*Fxy + dvzy*Fyy + (1. + dvzz)*Fzy;		// dw/dy
	F[2][2] = dvzx*Fxz + dvzy*Fyz + (1. + dvzz)*Fzz;		// 1 + dw/dz

	// store in total strain and rotation tensors
    ep->xx = F[0][0] - 1.;
    ep->yy = F[1][1] - 1.;
	ep->zz = F[2][2] - 1.;
    ep->xy = F[1][0] + F[0][1];
	ep->xz = F[2][0] + F[0][2];
	ep->yz = F[2][1] + F[1][2];
	
	// rotational strain increments
	wrot->xy = F[1][0] - F[0][1];
	wrot->xz = F[2][0] - F[0][2];
	wrot->yz = F[2][1] - F[1][2];
	
    // total residual stretch (1 + alpha dT + beta csat dConcentration)
	double resStretch = 1.0;
	double dTemp=mptr->pPreviousTemperature-thermal.reference;
	resStretch += CTE1*dTemp;
	if(DiffusionTask::active)
	{	double dConc=mptr->pPreviousConcentration-DiffusionTask::reference;
		resStretch += CME1*dConc;
	}
	
	// left Cauchy deformation tensor B = F F^T
	Tensor B;
	ZeroTensor(&B);
	int i;
	for(i=0;i<3;i++)
	{	B.xx += F[0][i]*F[0][i];
		B.yy += F[1][i]*F[1][i];
		B.zz += F[2][i]*F[2][i];
		B.xy += F[0][i]*F[1][i];
		B.xz += F[0][i]*F[2][i];
		B.yz += F[1][i]*F[2][i];
	}
	
	// J as determinant of F (or sqrt root of determinant of B)
	double J2 = B.xx*B.yy*B.zz + 2.*B.xy*B.xz*B.yz - B.yz*B.yz*B.xx - B.xz*B.xz*B.yy - B.xy*B.xy*B.zz;
	double J = sqrt(J2)/(resStretch*resStretch*resStretch);
	double J53 = pow(J, 5./3.);
	double J73 = pow(J, 7./3.);
	
	// find Cauchy stresses
	Tensor *sp=mptr->GetStressTensor();
	sp->xx = Ksp*(J-1.) + (2*B.xx-B.yy-B.zz)*G1sp/(3.*J53)
			+ (B.xx*(B.yy+B.zz)-2*B.yy*B.zz-B.xy*B.xy-B.xz*B.xz+2.*B.yz*B.yz)*G2sp/(3.*J73);
	sp->yy = Ksp*(J-1.) + (2*B.yy-B.xx-B.zz)*G1sp/(3.*J53)
			+ (B.yy*(B.xx+B.zz)-2*B.xx*B.zz-B.xy*B.xy+2.*B.xz*B.xz-B.yz*B.yz)*G2sp/(3.*J73);
	sp->zz = Ksp*(J-1.) + (2*B.zz-B.xx-B.yy)*G1sp/(3.*J53)
			+ (B.zz*(B.xx+B.yy)-2*B.xx*B.yy+2.*B.xy*B.xy-B.xz*B.xz-B.yz*B.yz)*G2sp/(3.*J73);
	sp->xy = B.xy*G1sp/J53 + (B.zz*B.xy-B.xz*B.yz)*G2sp/J73;
	sp->xz = B.xz*G1sp/J53 + (B.yy*B.xz-B.xy*B.yz)*G2sp/J73;
	sp->yz = B.yz*G1sp/J53 + (B.xx*B.yz-B.xy*B.xz)*G2sp/J73;
	
	// Convert to Kirchoff Stress (time J)
	/*
	sp->xx *= J;
	sp->yy *= J;
	sp->zz *= J;
	sp->xy *= J;
	sp->xz *= J;
	sp->yz *= J;
	*/
	
	// Convert to nominal stress by premultiply with J F^(-1)
	// JFi is transpose of matrix of cofactors for F
	/*
	double JFi[3][3];
	JFi[0][0] = F[1][1]*F[2][2] - F[2][1]*F[1][2];
	JFi[0][1] = -(F[0][1]*F[2][2] - F[2][1]*F[0][2]);
	JFi[0][2] = F[0][1]*F[1][2] - F[1][1]*F[0][2];
	JFi[1][0] = -(F[1][0]*F[2][2] - F[2][0]*F[1][2]);
	JFi[1][1] = F[0][0]*F[2][2] - F[2][0]*F[0][2];
	JFi[1][2] = -(F[0][0]*F[1][2] - F[1][0]*F[0][2]);
	JFi[2][0] = F[1][0]*F[2][1] - F[2][0]*F[1][1];
	JFi[2][1] = -(F[0][0]*F[2][1] - F[2][0]*F[0][1]);
	JFi[2][2] = F[0][0]*F[1][1] - F[1][0]*F[0][1];
	Tensor sp0=*sp;
	sp->xx = JFi[0][0]*sp0.xx + JFi[0][1]*sp0.xy + JFi[0][2]*sp0.xz;
	sp->xy = JFi[0][0]*sp0.xy + JFi[0][1]*sp0.yy + JFi[0][2]*sp0.yz;
	sp->xz = JFi[0][0]*sp0.xz + JFi[0][1]*sp0.yz + JFi[0][2]*sp0.zz;
	sp->yy = JFi[1][0]*sp0.xy + JFi[1][1]*sp0.yy + JFi[1][2]*sp0.yz;
	sp->yz = JFi[1][0]*sp0.xz + JFi[1][1]*sp0.yz + JFi[1][2]*sp0.zz;
	sp->zz = JFi[2][0]*sp0.xz + JFi[2][1]*sp0.yz + JFi[2][2]*sp0.zz;
	*/
	
	// strain energy
	double J23 = J53/J;
	double J43 = J23*J23;
	double I1bar = (B.xx+B.yy+B.zz)/J23;
	double I2bar = 0.5*(I1bar*I1bar - (B.xx*B.xx+B.yy*B.yy+B.zz*B.zz+2.*B.xy*B.xy+2*B.xz*B.xz+2.*B.yz*B.yz)/J43);
    mptr->SetStrainEnergy(0.5*(G1sp*(I1bar-3.) + G2sp*(I2bar-3.) + Ksp*(J-1.)*(J-1.)));
}

#pragma mark Mooney::Accessors

// Return the material tag
int Mooney::MaterialTag(void) { return MOONEYRIVLIN; }

/*	calculate wave speed in mm/sec (because G in MPa and rho in g/cm^3)
 Uses sqrt((K +4G/3)/rho) which is dilational wave speed
 at low strain G = G1+G2
*/
double Mooney::WaveSpeed(bool threeD)
{
    return sqrt(1.e9*(Kbulk+4.*(G1+G2)/3.)/rho);
}

/*	calculate shear wave speed in mm/sec (because G1 and G2 in MPa and rho in g/cm^3)
	at low strain G = G1+G2
*/
double Mooney::ShearWaveSpeed(bool threeD) { return sqrt(1.e9*(G1+G2)/rho); }

// return material type
const char *Mooney::MaterialType(void) { return "Mooney-Rivlin Rubber"; }


