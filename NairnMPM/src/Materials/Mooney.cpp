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
    switch(UofJOption)
    {   case J_MINUS_1_SQUARED:
            PrintProperty("U(J)",UofJOption,"[ = (K/2)(J-1)^2 ]");
            break;
            
        case LN_J_SQUARED:
            PrintProperty("U(J)",UofJOption,"[ = (K/2)(ln J)^2 ]");
            break;
            
        case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
        default:
            PrintProperty("U(J)",UofJOption,"[ = (K/2)((1/2)(J^2-1) - ln J) ]");
            break;
    }
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
    
    return(HyperElastic::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *Mooney::VerifyProperties(int np)
{
    if(G1<0. || Kbulk < 0. || G2<0.)
		return "Mooney-Rivlin Hyperelastic material needs non-negative G1, G2, and K";
    
    if(UofJOption<HALF_J_SQUARED_MINUS_1_MINUS_LN_J && UofJOption>LN_J_SQUARED)
		return "Mooney-Rivlin dilational energy (UJOption) must be 0, 1, or 2";

	// call super class
	return MaterialBase::VerifyProperties(np);
}

// Private properties used in constitutive law
void Mooney::InitialLoadMechProps(int makeSpecific,int np)
{
	// G1 and G2 in Specific units using initial rho
	// for MPM (units N/m^2 cm^3/g)
	G1sp=G1*1.0e+06/rho;
	G2sp=G2*1.0e+06/rho;
	
    // call super class for the rest
    HyperElastic::InitialLoadMechProps(makeSpecific,np);
}

// Store J, which is calculated incrementally, and available for archiving
// initialize to 1
char *Mooney::InitHistoryData(void)
{
	double *p = CreateAndZeroDoubles(1);
	*p=1.;
	return (char *)p;
}

#pragma mark Mooney::Methods


/* Take increments in strain and calculate new Particle: strains, rotation strain,
        stresses, strain energy,
    dvij are (gradient rates X time increment) to give deformation gradient change
    For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
    This material tracks pressure and stores deviatoric stress only in particle stress
        tensor
*/
void Mooney::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np)
{
	// Update strains and rotations and Left Cauchy strain
	double detDf = IncrementDeformation(mptr,du,NULL,np);
    
    // get pointer to new left Cauchy strain
    Tensor *B = mptr->GetElasticLeftCauchyTensor();
	
	// Deformation gradients and Cauchy tensor differ in plane stress
	if(np==PLANE_STRESS_MPM)
	{	// Find B->zz required to have zero stress in z direction
		
		// fixed arguments
		double arg = B->xx*B->yy - B->xy*B->xy;
		double arg12 = sqrt(arg);
		double arg16 = pow(arg,1./6.);
		double arg2 = B->xx+B->yy;
		
		// Newton-Rapheson starting at B.zz = 1
		// In tests finds answer in 3 or less steps
		Tensor *ep=mptr->GetStrainTensor();
		double xn16,xn12,xnp1,xn = 1.+ep->zz;
		double fx,fxp;
		int iter=1;
        
        // solution for B.zz in xn and J = sqrt(xn*arg)
        // Solving 0 = 3J^2 Kterm + G1(2*xn-arg2)J^(1/3) + G2*(xn*arg2 - 2*arg)/J^(1/3)
		while(iter<10)
		{	xn16 = pow(xn,1./6.);
			xn12 = sqrt(xn);
            
            switch(UofJOption)
            {   case J_MINUS_1_SQUARED:
                    // This is for Kterm = K(J-1)
                    fx = 3.*Ksp*xn*arg*(xn12*arg12-1.); 
                    fxp = 3.*Ksp*arg*(1.5*arg12*xn12-1.); 
                    break;
                    
                case LN_J_SQUARED:
                    // This is for Kterm = K(ln J)/J = (K/2)(ln J^2)/J
                    fx = 1.5*Ksp*xn12*arg12*log(xn*arg);
                    fxp = 0.75*Ksp*arg12*(2.+log(arg*xn))/xn12;
                    break;
                    
                case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
                default:
                    // This is for Kterm = (K/2)(J-1/J)
                    fx = 1.5*Ksp*xn12*arg12*(xn*arg-1.);            
                    fxp = 0.75*Ksp*arg12*(3.*arg*xn-1.)/xn12;
                    break;
            }
            
            // Now add the shear terms
            fx += G1sp*(2.*xn-arg2)*xn16*arg16 + G2sp*(xn*arg2-2.*arg)/(xn16*arg16);
            fxp += G1sp*xn16*arg16*(14.*xn-arg2)/(6.*xn) + G2sp*(2.*arg+5.*xn*arg2)/(6.*xn16*arg16*xn);
            
            // new prediction for solution
			xnp1 = xn - fx/fxp;
			//cout << iter << ": " << xn << "," << xnp1 << "," << fabs(xn-xnp1) << endl;
			if(fabs(xn-xnp1)<1e-6) break;
			xn = xnp1;
			iter+=1;
		}
        
        // Done and xn = new B->zz = Fzz^2 = dFzz*(old Bzz)*dFzz = dFzz^2*(old Bzz),
		//    and Fzz = dFzz*(old Fzz) = 1 + ep->zz
        double dFzz = sqrt(xn/B->zz);
        B->zz = xn;
		
		// particle strain ezz now known
		ep->zz = dFzz*(1.+ep->zz) - 1.;
        
        // incremental J changes
        detDf *= dFzz;
	}
    
    // Increment J and save it in history data
    double J = detDf*mptr->GetHistoryDble();
    mptr->SetHistoryDble(J);
    
    // account for residual stresses
	double resStretch = GetResidualStretch(mptr);
	J /= (resStretch*resStretch*resStretch);
    
    // update pressure
    double Kse;
	double Kterm = J*GetVolumetricTerms(J,&Kse);       // times J to get Kirchoff stress
    mptr->SetPressure(-Kterm);
    mptr->SetStrainEnergy(Kse);
	
    // Account for density change in specific stress
    // i.e.. Get (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
	double J23 = pow(J, 2./3.);
	double J43 = J23*J23;
	double JforG1 = J23;                 // J^(2/3) = J^(5/3)/J to get Kirchoff stress
	double JforG2 = J43;                 // J^(4/3) = J^(7/3)/J to get Kirchoff stress
    
	// find deviatoric (Cauchy stress)J/rho0 = deviatoric (Kirchoff stress)/rho0
	Tensor *sp=mptr->GetStressTensor();
    
	sp->xx = (2*B->xx-B->yy-B->zz)*G1sp/(3.*JforG1)
            + (B->xx*(B->yy+B->zz)-2*B->yy*B->zz-B->xy*B->xy)*G2sp/(3.*JforG2);
	sp->yy = (2*B->yy-B->xx-B->zz)*G1sp/(3.*JforG1)
            + (B->yy*(B->xx+B->zz)-2*B->xx*B->zz-B->xy*B->xy)*G2sp/(3.*JforG2);
	sp->zz = (2*B->zz-B->xx-B->yy)*G1sp/(3.*JforG1)
            + (B->zz*(B->xx+B->yy)-2*B->xx*B->yy+2.*B->xy*B->xy)*G2sp/(3.*JforG2);
	sp->xy = B->xy*G1sp/JforG1 + (B->zz*B->xy)*G2sp/JforG2;
    
    if(np==THREED_MPM)
    {   sp->xx += (2.*B->yz*B->yz-B->xz*B->xz)*G2sp/(3.*JforG2);
        sp->yy += (2.*B->xz*B->xz-B->yz*B->yz)*G2sp/(3.*JforG2);
        sp->zz -= (B->xz*B->xz+B->yz*B->yz)*G2sp/(3.*JforG2);
        sp->xy -= B->xz*B->yz*G2sp/JforG2;
	    sp->xz = B->xz*G1sp/JforG1 + (B->yy*B->xz-B->xy*B->yz)*G2sp/JforG2;
        sp->yz = B->yz*G1sp/JforG1 + (B->xx*B->yz-B->xy*B->xz)*G2sp/JforG2;
    }
    
	// strain energy per unit mass (U/(rho0 V0)) and we are using
    // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
	double I1bar = (B->xx+B->yy+B->zz)/J23;
	double I2bar = 0.5*(I1bar*I1bar - (B->xx*B->xx+B->yy*B->yy+B->zz*B->zz+2.*B->xy*B->xy+2*B->xz*B->xz+2.*B->yz*B->yz)/J43);
    mptr->AddStrainEnergy(0.5*(G1sp*(I1bar-3.) + G2sp*(I2bar-3.)));

}

#pragma mark Mooney::Accessors

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor Mooney::GetStress(Tensor *sp,double pressure)
{   Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}

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

// archive material data for this material type when requested.
double Mooney::GetHistory(int num,char *historyPtr)
{   double history=0.;
    if(num==1)
    {	double *J=(double *)historyPtr;
        history=*J;
    }
    return history;
}



