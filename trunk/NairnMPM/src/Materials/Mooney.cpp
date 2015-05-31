/********************************************************************************
    Mooney.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Fri Feb 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/Mooney.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "System/UnitsController.hpp"

#pragma mark Mooney::Constructors and Destructors

// Constructors
Mooney::Mooney()
{
}

// Constructors with arguments 
Mooney::Mooney(char *matName) : HyperElastic(matName)
{
	G1 = -1.;			// Must enter K and G1 OR Etens and nu
	Etens = -1.;
	nu = -2.;
	G2 = 0.;			// zero is neo-Hookean
    rubber = FALSE;     // not a rubber
	J_History = 0;
}

#pragma mark Mooney::Initialization

// print mechanical properties output window
void Mooney::PrintMechanicalProperties(void) const
{
	PrintProperty("G1",G1*UnitsController::Scaling(1.e-6),"");
	PrintProperty("G2",G2*UnitsController::Scaling(1.e-6),"");
	if(nu<0.5)
		PrintProperty("K",Kbulk*UnitsController::Scaling(1.e-6),"");
	else
		cout << "K = infinite";
	cout << endl;
	
	PrintProperty("E",Etens*UnitsController::Scaling(1.e-6),"");
	PrintProperty("nu",nu,"");
	if(nu==0.5)
		cout << "incompressible";
	else
		PrintProperty("lam",nu*Etens*UnitsController::Scaling(1.e-6)/(1.+nu)/(1.-2.*nu),"");
	cout << endl;
	
	PrintProperty("a",aI,"");
	PrintProperty("gam0",gamma0,"");
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
	
	// call superclass here if it is not Material base
	HyperElastic::PrintMechanicalProperties();
}
	
// Read material properties
char *Mooney::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"G1")==0)
        return UnitsController::ScaledPtr((char *)&G1,gScaling,1.e6);
    
    else if(strcmp(xName,"G2")==0)
        return UnitsController::ScaledPtr((char *)&G2,gScaling,1.e6);
    
    else if(strcmp(xName,"E")==0)
        return UnitsController::ScaledPtr((char *)&Etens,gScaling,1.e6);
    
    else if(strcmp(xName,"nu")==0)
        return (char *)&nu;
    
    else if(strcmp(xName,"IdealRubber")==0)
    {   input = NOT_NUM;
        rubber = TRUE;
        return (char *)&rubber;
    }
    
    return(HyperElastic::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *Mooney::VerifyAndLoadProperties(int np)
{
	// must enter G1 and Kbulk OR Etens and nu
	if(G1>=0. && Kbulk>=0.)
	{	if(Etens>=0. || nu>=-1.)
			return "Mooney-Rivlin Hyperelastic material needs K and G1 OR E and nu";
		Etens = 9.*Kbulk*G1/(3.*Kbulk+G1);
		nu = (3.*Kbulk-2.*G1)/(6.*Kbulk+2.*G1);
	}
	else if(G1>=0. || Kbulk>=0. || Etens<0. || nu<-1.)
	{	return "Mooney-Rivlin Hyperelastic material needs K and G1 OR E and nu";
	}
	else
	{	// has Etens and nu
		if(nu>0.5)
			return "Mooney-Rivlin Hyperelastic Poisson's ratio must be between -1 and 1/2";
		
		// get G1+G2 equal to low strain shear modulus
		G1 = Etens/(2.*(1.+nu)) - G2;
		
		// bulk modulus, but allow membrane to be incompressible
		if(nu<0.5)
			Kbulk = Etens/(3.*(1-2.*nu));
		else if(!isMembrane())
			return "Mooney-Rivlin Hyperelastic Poisson's ratio for solid material must be less than 1/2";
		else
		{	Kbulk = 0.;	// which causes CP-Cv=0 as well
		}
	}
		
    if(G2<0.)
		return "Mooney-Rivlin Hyperelastic material needs non-negative G2";
    
    if(UofJOption<HALF_J_SQUARED_MINUS_1_MINUS_LN_J && UofJOption>LN_J_SQUARED)
		return "Mooney-Rivlin dilational energy (UJOption) must be 0, 1, or 2";

	// G1 and G2 in Specific units using initial rho (F/L^2 L^3/mass)
	G1sp=G1/rho;
	G2sp=G2/rho;
	
    // heating gamma0 (dimensionless)
    double alphaV = 3.e-6*aI;
    gamma0 = Kbulk*alphaV/(rho*heatCapacity);
	
	// call super class
	return HyperElastic::VerifyAndLoadProperties(np);
}

#pragma mark Mooney::History Data

// Store J, which is calculated incrementally, and available for archiving
// initialize to 1
char *Mooney::InitHistoryData(void)
{
	double *p = CreateAndZeroDoubles(1);
	*p=1.;
	return (char *)p;
}

// archive material data for this material type when requested.
double Mooney::GetHistory(int num,char *historyPtr) const
{   double history=0.;
    if(num==1)
    {	double *J=(double *)historyPtr;
        history=*J;
    }
    return history;
}

#pragma mark Mooney::Methods


/* Take increments in strain and calculate new Particle: strains, rotation strain,
        stresses, strain energy,
    dvij are (gradient rates X time increment) to give deformation gradient change
    For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
    This material tracks pressure and stores deviatoric stress only in particle stress
        tensor
*/
void Mooney::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// incremental energy, store initial stress
	Tensor *sporig=mptr->GetStressTensor();
	Tensor st0 = *sporig;
	
	// Update strains and rotations and Left Cauchy strain
	double detDf = IncrementDeformation(mptr,du,NULL,np);
    
    // get pointer to new left Cauchy strain
    Tensor *B = mptr->GetAltStrainTensor();
	
    // account for residual stresses
	double dresStretch,resStretch = GetResidualStretch(mptr,dresStretch,res);
	double resStretch2 = resStretch*resStretch;
	double Jres = resStretch2*resStretch;

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
		double xn16,xn12,xnp1,xn = (1.+ep->zz)*(1.+ep->zz);
		double fx,fxp,J13,J0,J2,Jeff;
		int iter=1;
        double mJ2P,mdJ2PdJ;
        
        // solution for B.zz in xn and J = sqrt(xn*arg) with dJ/dxn = arg/(2 sqrt(xn*arg))
        // Solving f=0 where f = 3J^2 Kterm + Jres*G1(2*xn-arg2)J^(1/3) + Jres*G2(xn*arg2 - 2*arg)/J^(1/3)
		//			where J^(1/3) = (xn*arg)^(1/6)
		// df/dxn = d(3J^2 Kterm)/dJ dJ/dxn
		//				+ Jres*G1((14*xn-arg2)/(6*xn))J^(1/3)
		//				+ Jres*G2((5*xn*arg2+2*arg)/(6*xn))/J^(1/3)
		while(iter<20)
		{	xn16 = pow(xn,1./6.);
			xn12 = sqrt(xn);
			J13 = xn16*arg16;
			J0 = xn12*arg12;
			J2 = J0*J0;
            Jeff = J0/Jres;
            
            // get f and df/dxn
            GetNewtonPressureTerms(Jeff, Ksp, mJ2P, mdJ2PdJ);
            fx = 3.*Jres*mJ2P + G1sp*(2.*xn-arg2)*J13 + G2sp*(xn*arg2-2.*arg)/J13;
            fxp = (1.5*J0/xn)*mdJ2PdJ + G1sp*J13*(14.*xn-arg2)/(6.*xn) + G2sp*(2.*arg+5.*xn*arg2)/(6.*J13*xn);
            
            // new prediction for solution
			xnp1 = xn - fx/fxp;
			//cout << iter << ": " << xn << "," << xnp1 << "," << fabs(xn-xnp1) << endl;
			if(fabs(xn-xnp1)<1e-10) break;
			xn = xnp1;
			iter+=1;
		}
        
        if(iter>=20) cout << "# Not enough iterations in plane stress Mooney-Rivlin material" << endl;
        
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
	double Jeff = J/Jres;
	
    // update pressure
	double p0=mptr->GetPressure();
	double Kterm = J*GetVolumetricTerms(Jeff,Ksp);       // times J to get Kirchoff stress
	
	//Kterm /= Jres;		// this uses Jeff to get Kirchoff stress
    
    // artifical viscosity
    double delV = 1. - 1./detDf;                        // total volume change
    double QAVred = 0.,AVEnergy=0.;
    if(delV<0. && artificialViscosity)
	{	QAVred = GetArtificalViscosity(delV/delTime,sqrt(Ksp*J));
        if(ConductionTask::AVHeating) AVEnergy = fabs(QAVred*delV);
    }
    double Pfinal = -Kterm + QAVred;
    
    // set the pressure
    mptr->SetPressure(Pfinal);
	
	// incremental energy - dilational part
    double avgP = 0.5*(p0+Pfinal);
    double dilEnergy = -avgP*delV;
	
	// incremental residual energy
	double delVres = 1. - 1./(dresStretch*dresStretch*dresStretch);
	double resEnergy = -avgP*delVres;
	
    // Account for density change in specific stress
    // i.e.. Get (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
	double J23 = pow(J, 2./3.);
	double J43 = J23*J23;
	double JforG1 = J23/Jres;                   // J^(5/3)/(Jres J) = J^(2/3)/Jres to get Kirchoff stress
	double JforG2 = J43/Jres;                   // J^(7/3)/(Jres J) = J^(4/3)/Jres to get Kirchoff stress
	
	//JforG1 *= Jres;			// this uses Jeff to get Kirchoff stress
	//JforG2 *= Jres;			// this uses Jeff to get Kirchoff stress
    
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
    
	// incremental work energy = shear energy
    double shearEnergy = 0.5*((sp->xx+st0.xx)*du(0,0) + (sp->yy+st0.yy)*du(1,1) + (sp->zz+st0.zz)*du(2,2)+
							  (sp->xy+st0.xy)*(du(0,1)+du(1,0)));
    if(np==THREED_MPM)
    {   shearEnergy += 0.5*((sp->xz+st0.xz)*(du(0,2)+du(2,0)) + (sp->yz+st0.yz)*(du(1,2)+du(2,1)));
    }
    
    // strain energy
    double dU = dilEnergy + shearEnergy;
    mptr->AddWorkEnergyAndResidualEnergy(dU,resEnergy);
	
	// particle isentropic temperature increment
	double Kratio;				// = rho_0 K/(rho K_0)
	switch(UofJOption)
	{   case J_MINUS_1_SQUARED:
			Kratio = Jeff;
			break;
			
		case LN_J_SQUARED:
			Kratio = (1-log(Jeff))/(Jeff*Jeff);
			break;
			
		case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
		default:
			Kratio = 0.5*(Jeff + 1./Jeff);
			break;
	}
    double dTq0 = -J*Kratio*gamma0*mptr->pPreviousTemperature*delV;
    
    // thermodynamics depends on whether or not this is a rubber
    if(rubber)
    {   // just like ideal gas
        IncrementHeatEnergy(mptr,res->dT,dTq0,dU+QAVred);
    }
    else
    {   // elastic - no heating
        IncrementHeatEnergy(mptr,res->dT,dTq0,QAVred);
    }
}

#pragma mark Mooney::Accessors

// convert J to K using isotropic method
Vector Mooney::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{	double GLS = G1+G2;
	double nuLS = (3.*Kbulk-2.*GLS)/(6.*Kbulk+2.*GLS);
	return IsotropicJToK(d,C,J0,np,nuLS,GLS);
}

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor Mooney::GetStress(Tensor *sp,double pressure) const
{   Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}

// Return the material tag
int Mooney::MaterialTag(void) const { return MOONEYRIVLIN; }

// Calculate wave speed in L/sec (because K and G in mass/(L sec^2) and rho in mass/L^3)
// Uses sqrt((K +4G/3)/rho) which is dilational wave speed at low strain G = G1+G2
double Mooney::WaveSpeed(bool threeD,MPMBase *mptr) const
{	return sqrt((Kbulk+4.*(G1+G2)/3.)/rho);
}

// Calculate shear wave speed in mm/sec (because G1 and G2 in mass/(L sec^2) and rho in mass/L^3)
double Mooney::ShearWaveSpeed(bool threeD,MPMBase *mptr) const { return sqrt((G1+G2)/rho); }

// return material type
const char *Mooney::MaterialType(void) const { return "Mooney-Rivlin Hyperelastic"; }

// if a subclass material supports artificial viscosity, override this and return TRUE
bool Mooney::SupportsArtificialViscosity(void) const { return true; }

// Get current relative volume change = J (which this material tracks)
// All subclasses must track J in J_History (or override this method)
double Mooney::GetCurrentRelativeVolume(MPMBase *mptr) const
{   return mptr->GetHistoryDble(J_History);
}


