/********************************************************************************
    Mooney.cpp
    NairnMPM
    
    Created by John Nairn on Fri Feb 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/Mooney.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
 
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
    rubber = FALSE;     // not a rubber
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
        return (char *)&G1;
    
    else if(strcmp(xName,"G2")==0)
        return (char *)&G2;
    
    else if(strcmp(xName,"IdealRubber")==0)
    {   input = NOT_NUM;
        rubber = TRUE;
        return (char *)&rubber;
    }
    
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
    For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
    This material tracks pressure and stores deviatoric stress only in particle stress
        tensor
*/
void Mooney::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np)
{
	// incremental energy, store initial stress
	Tensor *sporig=mptr->GetStressTensor();
	Tensor st0 = *sporig;
	
	// Update strains and rotations and Left Cauchy strain
	double detDf = IncrementDeformation(mptr,du,NULL,np);
    
    // get pointer to new left Cauchy strain
    Tensor *B = mptr->GetElasticLeftCauchyTensor();
	
    // account for residual stresses
	// Divid J bt resStretch3 adn elements of B by resStretch2
	double dresStretch,resStretch = GetResidualStretch(mptr,dresStretch);
	double resStretch2 = resStretch*resStretch;
	double resStretch3 = resStretch2*resStretch;
	double G1eff = G1sp/resStretch2;
	double G2eff = G2sp/(resStretch2*resStretch2);

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
		double fx,fxp,J13,J0,J2;
		int iter=1;
        
        // solution for B.zz in xn and J = sqrt(xn*arg)/resStretch3 with dJ/dxn = arg/(2 sqrt(xn*arg) resStretch3)
        // Solving f=0 where f = 3J^2 Kterm + G1'(2*xn-arg2)J^(1/3) + G2'(xn*arg2 - 2*arg)/J^(1/3)
		//			where J^(1/3) = (xn*arg)^(1/6)/resStretch, G1' = G1eff, G2' = G2eff
		// df/dxn = 3 d(J^2 Kterm)/dJ dJ/dxn
		//				+ G1'(2+(2*xn-arg2)/(3J)*dJ/dxn)J^(1/3)
		//				+ G2'(arg2-(xn*arg2-2*arg)/(3J)*dJ/dxn)/J^(1/3)
		// df/dxn = d(3J^2 Kterm)/dJ dJ/dxn
		//				+ G1'((14*xn-arg2)/(6*xn))J^(1/3)
		//				+ G2'((5*xn*arg2+2*arg)/(6 xn))/J^(1/3)
		while(iter<10)
		{	xn16 = pow(xn,1./6.);
			xn12 = sqrt(xn);
			J13 = xn16*arg16/resStretch;
			J0 = xn12*arg12/resStretch3;
			J2 = J0*J0;
            
            switch(UofJOption)
            {   case J_MINUS_1_SQUARED:
                    // This is for 3J^2Kterm = 3KJ^2(J-1)
                    fx = 3.*Ksp*J2*(J0-1.);
					// d(3K(J^3-J^2))/dJ dJ/dxn = 1.5(3J^2-2J) arg12/(xn12 resStretch3)
                    fxp = 1.5*Ksp*(3.*J2-2.*J0)*arg12/(xn12*resStretch3);
                    break;
                    
                case LN_J_SQUARED:
                    // This is for 3J^2Kterm = 3J K(ln J) = 3J(K/2)(ln J^2)
                    fx = 1.5*Ksp*J0*log(J2);
					// d(1.5 J K(ln J^2))/dJ dJ/dxn = 0.75K(2 + ln J^2) arg12/(xn12 resStretch3)
                    fxp = 0.75*Ksp*arg12*(2.+log(J2))/(xn12*resStretch3);
                    break;
                    
                case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
                default:
                    // This is for 3J^2 Kterm = 3(K/2)J^2(J-1/J) = 1.5 K J(J^2-1)
                    fx = 1.5*Ksp*J0*(J2-1.);
					// d(1.5 K(J^3-J))/dJ dJ/dxn = 0.75 K (3J^2-1) arg12/(xn12 resStretch3)
                    fxp = 0.75*Ksp*arg12*(3.*J2-1.)/(xn12*resStretch3);
                    break;
            }
            
            // Now add the shear terms
            fx += G1eff*(2.*xn-arg2)*J13 + G2eff*(xn*arg2-2.*arg)/J13;
            fxp += G1eff*J13*(14.*xn-arg2)/(6.*xn) + G2eff*(2.*arg+5.*xn*arg2)/(6.*J13*xn);
            
            // new prediction for solution
			xnp1 = xn - fx/fxp;
			//cout << iter << ": " << xn << "," << xnp1 << "," << fabs(xn-xnp1) << endl;
			if(fabs(xn-xnp1)<1e-7) break;
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
	J /= resStretch3;
    detDf /= dresStretch;
	
    // update pressure
	double p0=mptr->GetPressure();
	double Kterm = J*GetVolumetricTerms(J);       // times J to get Kirchoff stress
    mptr->SetPressure(-Kterm);
	
	// incremental energy - dilational part
    double delV = 1. - 1./detDf;
    double avgP = 0.5*(p0-Kterm);
    double dilEnergy = -avgP*delV;
	
    // Account for density change in specific stress
    // i.e.. Get (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
	double J23 = pow(J, 2./3.);
	double J43 = J23*J23;
	double JforG1 = J23;				// J^(2/3) = J^(5/3)/J to get Kirchoff stress
	double JforG2 = J43;                // J^(4/3) = J^(7/3)/J to get Kirchoff stress
    
	// find deviatoric (Cauchy stress)J/rho0 = deviatoric (Kirchoff stress)/rho0
	Tensor *sp=mptr->GetStressTensor();
    
	sp->xx = (2*B->xx-B->yy-B->zz)*G1eff/(3.*JforG1)
            + (B->xx*(B->yy+B->zz)-2*B->yy*B->zz-B->xy*B->xy)*G2eff/(3.*JforG2);
	sp->yy = (2*B->yy-B->xx-B->zz)*G1eff/(3.*JforG1)
            + (B->yy*(B->xx+B->zz)-2*B->xx*B->zz-B->xy*B->xy)*G2eff/(3.*JforG2);
	sp->zz = (2*B->zz-B->xx-B->yy)*G1eff/(3.*JforG1)
            + (B->zz*(B->xx+B->yy)-2*B->xx*B->yy+2.*B->xy*B->xy)*G2eff/(3.*JforG2);
	sp->xy = B->xy*G1eff/JforG1 + (B->zz*B->xy)*G2eff/JforG2;
    
    if(np==THREED_MPM)
    {   sp->xx += (2.*B->yz*B->yz-B->xz*B->xz)*G2eff/(3.*JforG2);
        sp->yy += (2.*B->xz*B->xz-B->yz*B->yz)*G2eff/(3.*JforG2);
        sp->zz -= (B->xz*B->xz+B->yz*B->yz)*G2eff/(3.*JforG2);
        sp->xy -= B->xz*B->yz*G2eff/JforG2;
	    sp->xz = B->xz*G1eff/JforG1 + (B->yy*B->xz-B->xy*B->yz)*G2eff/JforG2;
        sp->yz = B->yz*G1eff/JforG1 + (B->xx*B->yz-B->xy*B->xz)*G2eff/JforG2;
    }
    
	// strain energy per unit mass (U/(rho0 V0)) and we are using
    // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
	// Divide B's by resStretch2 to accound for current stress free state
	//double I1bar = (B->xx+B->yy+B->zz)/(resStretch2*J23);
	//double I2bar = 0.5*(I1bar*I1bar - (B->xx*B->xx+B->yy*B->yy+B->zz*B->zz+2.*B->xy*B->xy+2*B->xz*B->xz+2.*B->yz*B->yz)
	//									/(J43*resStretch2*resStretch2));
    //mptr->AddStrainEnergy(0.5*(G1sp*(I1bar-3.) + G2sp*(I2bar-3.)));

	// incremental energy = shear energy (2D only)
    double shearEnergy = 0.5*((sp->xx+st0.xx)*du(0,0) + (sp->yy+st0.yy)*du(1,1) + (sp->zz+st0.zz)*du(2,2)+
							  (sp->xy+st0.xy)*(du(0,1)+du(1,0)))/resStretch;
    if(np==THREED_MPM)
    {   shearEnergy += 0.5*((sp->xz+st0.xz)*(du(0,2)+du(2,0)) + (sp->yz+st0.yz)*(du(1,2)+du(2,1)))/resStretch;
    }
    
    // strain energy
    double dU = dilEnergy + shearEnergy;
    mptr->AddStrainEnergy(dU);
    
    // thermodynamics depends on whether or not this is a rubber
    double dT = ConductionTask::dTemperature;
    if(rubber)
    {   // just like ideal gas
        IncrementHeatEnergy(mptr,dT,0.,dU);
    }
    else
    {   // elastic - no heating
        IncrementHeatEnergy(mptr,dT,0.,0.);
    }
    
    // the plastic energy is not otherwise used, so let's track entropy
    double dS = 0., Cv = 1000.*GetHeatCapacity(mptr);
    double Tp = mptr->pPreviousTemperature;
    if(ConductionTask::energyCoupling)
    {   double dTS = dU/Cv;
        dS = Cv*(dT/Tp - dTS/(Tp+dTS));
    }
    else
        dS = (Cv*dT - dU)/Tp;
	mptr->AddPlastEnergy(dS);
    
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



