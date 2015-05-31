/********************************************************************************
	Neohookean.cpp
	nairn-mpm-fea

	Created by John Nairn on Sun Mar 10 2014.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/Neohookean.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "System/UnitsController.hpp"

#pragma mark Neohookean::Constructors and Destructors

// Constructors
Neohookean::Neohookean()
{
}

// Constructors with arguments
Neohookean::Neohookean(char *matName) : HyperElastic(matName)
{
	G = -1.;			// Must enter K and G OR Etens and nu OR Lame and G
	Etens = -1.;
	Lame = -1.;
	nu = -2.;
	J_History = 0;
}

#pragma mark Neohookean::Initialization

// print mechanical properties output window
void Neohookean::PrintMechanicalProperties(void) const
{
	PrintProperty("G",G*UnitsController::Scaling(1.e-6),"");
	if(nu<0.5)
	{	PrintProperty("K",Kbulk*UnitsController::Scaling(1.e-6),"");
		PrintProperty("lam",Lame*UnitsController::Scaling(1.e-6),"");
	}
	else
		cout << "K = lam = infinite";
	cout << endl;
	
	PrintProperty("E",Etens*UnitsController::Scaling(1.e-6),"");
	PrintProperty("nu",nu,"");
	if(nu==0.5) cout << "incompressible";
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
char *Neohookean::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"G")==0)
		return UnitsController::ScaledPtr((char *)&G,gScaling,1.e6);
    
    else if(strcmp(xName,"Lame")==0)
		return UnitsController::ScaledPtr((char *)&Lame,gScaling,1.e6);
    
    else if(strcmp(xName,"E")==0)
		return UnitsController::ScaledPtr((char *)&Etens,gScaling,1.e6);
    
    else if(strcmp(xName,"nu")==0)
        return (char *)&nu;
    
    return(HyperElastic::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *Neohookean::VerifyAndLoadProperties(int np)
{
	// must enter G1 and Kbulk OR Etens and nu
	if(G>=0. && Kbulk>=0.)
	{	if(Etens>=0. || nu>=-1. || Lame>=0.)
			return "Neohookean Hyperelastic material needs K and G, Lame and G, OR E and nu";
		Etens = 9.*Kbulk*G/(3.*Kbulk+G);
		nu = (3.*Kbulk-2.*G)/(6.*Kbulk+2.*G);
		Lame = Kbulk - 2.*G/3.;
	}
	else if(G>=0. || Lame>=0.)
	{	if(Etens>=0. || nu>=-1. || Kbulk>=0.)
			return "Neohookean Hyperelastic material needs K and G, Lame and G, OR E and nu";
		Kbulk = Lame + 2.*G/3.;
		Etens = 9.*Kbulk*G/(3.*Kbulk+G);
		nu = (3.*Kbulk-2.*G)/(6.*Kbulk+2.*G);
	}
	else if(Etens>=0 && G<0. && Lame<0. && Kbulk<0.)
	{	// has Etens and nu
		if(nu>0.5 || nu<-1.)
			return "Neohookean Hyperelastic material needs K and G, Lame and G, OR E and nu (with nu between -1 and 1/2)";
		G = Etens/(2.*(1.+nu));
		
		// bulk modulus, but allow membrane to be incompressible
		if(nu<0.5)
		{	Kbulk = Etens/(3.*(1-2.*nu));
			Lame = nu*Etens/(1.+nu)/(1.-2.*nu);
		}
		else if(!isMembrane())
			return "Neohookean Hyperelastic Poisson's ratio for solid material must be less than 1/2";
		else
		{	Kbulk = 0.;			// whcih causes Cp-Cv=0 as well
			Lame = 0.;
		}
	}
	
    if(UofJOption<HALF_J_SQUARED_MINUS_1_MINUS_LN_J && UofJOption>LN_J_SQUARED)
		return "Neohookean dilational energy (UJOption) must be 0, 1, or 2";
	
	// Convert to specific units (F/L^2 L^3/mass)
	pr.Gsp = G/rho;
	pr.Lamesp = Lame/rho;
    pr.Ksp = pr.Lamesp + 2.*pr.Gsp/3.;
	
    // heating gamma0
    double alphaV = 3.e-6*aI;
    gamma0 = Kbulk*alphaV/(rho*heatCapacity);
	
	// call super class
	return HyperElastic::VerifyAndLoadProperties(np);
}

#pragma mark Neohookean::History Data

// Store J, which is calculated incrementally, and available for archiving
// initialize to 1
char *Neohookean::InitHistoryData(void)
{
	double *p = CreateAndZeroDoubles(1);
	*p=1.;
	return (char *)p;
}

// archive material data for this material type when requested.
double Neohookean::GetHistory(int num,char *historyPtr) const
{   double history=0.;
    if(num==1)
    {	double *J=(double *)historyPtr;
        history=*J;
    }
    return history;
}

#pragma mark Neohookean::Methods

/* Take increments in strain and calculate new Particle: strains, rotation strain,
	stresses, strain energy,
	dvij are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
	This material tracks pressure and stores deviatoric stress only in particle stress
	tensor
*/
void Neohookean::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// Update strains and rotations and Left Cauchy strain
	double detDf = IncrementDeformation(mptr,du,NULL,np);
	
	// get pointer to new left Cauchy strain
    Tensor *B = mptr->GetAltStrainTensor();
	
    // account for residual stresses
	double dresStretch,resStretch = GetResidualStretch(mptr,dresStretch,res);
	double Jres23 = resStretch*resStretch;
	double Jres = Jres23*resStretch;
	
	// Deformation gradients and Cauchy tensor differ in plane stress
	if(np==PLANE_STRESS_MPM)
	{	// Find B->zz required to have zero stress in z direction
		double arg = B->xx*B->yy - B->xy*B->xy;
		double xn;
		Tensor *ep=mptr->GetStrainTensor();
		
		switch(UofJOption)
		{   case J_MINUS_1_SQUARED:
			{	double a = pr.Lamesp*arg + pr.Gsp*pow(Jres,4./3.);
				double b = pr.Lamesp*sqrt(arg);
				xn = Jres*(b + sqrt(b*b + 4.*pr.Gsp*a))/(2.*a);
				xn *= xn;
				break;
			}
			case LN_J_SQUARED:
			{	xn = B->zz;
				double fx,fxp,xnp1,Jres23 = pow(Jres,2./3.);
				int iter=1;
				while(iter<20)
				{	// get f and df/dxn
					fx = pr.Gsp*(xn-Jres23) + 0.5*pr.Lamesp*Jres23*log(xn*arg/(Jres*Jres));
					fxp = pr.Gsp + pr.Lamesp*Jres23/(2*xn);
					
					// new prediction for solution
					xnp1 = xn - fx/fxp;
					if(fabs(xn-xnp1)<1e-10) break;
					xn = xnp1;
					iter+=1;
				}
				break;
			}
			case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
			default:
				xn = Jres*Jres*(pr.Lamesp+2.*pr.Gsp)/(pr.Lamesp*arg+2.*pr.Gsp*pow(Jres,4./3.));
				break;
		}
		
        // Done and xn = new B->zz = Fzz^2 = dFzz*(old Bzz)*dFzz = dFzz^2*(old Bzz),
		//    and Fzz = dFzz*(old Fzz) = 1 + ep->zz
        //cout << xn << "," << pow(Jres,2./3.)*pr.Lamesp*(xn*arg/(Jres*Jres)-1.) + 2.*pr.Gsp*(xn-pow(Jres,2./3.)) << endl;;
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
	
	// for incremental energy, store initial stress
	Tensor *sporig=mptr->GetStressTensor();
	Tensor st0 = *sporig;
	
	// account for residual stresses
	double Jeff = J/Jres;
	
    // update pressure
	double p0=mptr->GetPressure();
	double Pterm = J*GetVolumetricTerms(Jeff,pr.Lamesp) + Jres*pr.Gsp*((B->xx+B->yy+B->zz)/(3.*Jres23) - 1.);
	
    // artifical viscosity
    double delV = 1. - 1./detDf;                        // total volume change
    double QAVred = 0.,AVEnergy=0.;
    if(delV<0. && artificialViscosity)
	{	QAVred = GetArtificalViscosity(delV/delTime,sqrt(pr.Ksp)*J);
        if(ConductionTask::AVHeating) AVEnergy = fabs(QAVred*delV);
    }
    double Pfinal = -Pterm + QAVred;
    
    // set the pressure
    mptr->SetPressure(Pfinal);
	
	// incremental energy - dilational part
    double avgP = 0.5*(p0+Pfinal);
    double dilEnergy = -avgP*delV;
	
	// incremental residual energy
	double delVres = 1. - 1./(dresStretch*dresStretch*dresStretch);
	double resEnergy = -avgP*delVres;
	
    // Account for density change in specific stress
    // i.e.. Get (Kirchoff Stress)/rho0
	double GJeff = resStretch*pr.Gsp;		// = J*(Jres^(1/3) G/J) to get Kirchoof
    
	// find deviatoric (Cauchy stress)J/rho0 = deviatoric (Kirchoff stress)/rho0
	Tensor *sp=mptr->GetStressTensor();
    
    double I1third = (B->xx+B->yy+B->zz)/3.;
	sp->xx = GJeff*(B->xx-I1third);
	sp->yy = GJeff*(B->yy-I1third);
	sp->zz = GJeff*(B->zz-I1third);
	sp->xy = GJeff*B->xy;
    if(np==THREED_MPM)
	{	sp->xz = GJeff*B->xz;
        sp->yz = GJeff*B->yz;
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
    double Jeff1third = pow(Jeff,1./3.);
    double Gterm = pr.Gsp*(1. - Jeff1third*Jeff1third + 2./(3.*Jeff1third));
	switch(UofJOption)
	{   case J_MINUS_1_SQUARED:
			Kratio = pr.Lamesp+Jeff + Gterm;
			break;
			
		case LN_J_SQUARED:
			Kratio = pr.Lamesp*(1-log(Jeff))/(Jeff*Jeff) + Gterm;
			break;
			
		case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
		default:
			Kratio = 0.5*(Jeff + 1./Jeff);
			break;
	}
    Kratio /= pr.Ksp;
    double dTq0 = -J*Kratio*gamma0*mptr->pPreviousTemperature*delV;
    
    // thermodynamics heat and temperature
	IncrementHeatEnergy(mptr,res->dT,dTq0,QAVred);
}

#pragma mark Mooney::Accessors

// convert J to K using isotropic method
Vector Neohookean::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{	double nuLS = (3.*Kbulk-2.*G)/(6.*Kbulk+2.*G);
	return IsotropicJToK(d,C,J0,np,nuLS,G);
}

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor Neohookean::GetStress(Tensor *sp,double pressure) const
{   Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}

// Return the material tag
int Neohookean::MaterialTag(void) const { return NEOHOOKEAN; }

// Calculate wave speed in L/sec (because G in mass/(L sec^2) and rho in mass/L^3)
// Uses sqrt((K +4G/3)/rho) which is dilational wave speed at low strain
double Neohookean::WaveSpeed(bool threeD,MPMBase *mptr) const
{	return sqrt((Kbulk+4.*G/3.)/rho);
}

// Calculate shear wave speed in L/sec (because G in mass/(L sec^2) and rho in mass/L^3)
double Neohookean::ShearWaveSpeed(bool threeD,MPMBase *mptr) const { return sqrt(G/rho); }

// return material type
const char *Neohookean::MaterialType(void) const { return "Neohookean Hyperelastic"; }

// if a subclass material supports artificial viscosity, override this and return TRUE
bool Neohookean::SupportsArtificialViscosity(void) const { return true; }

// Get current relative volume change = J (which this material tracks)
// All subclasses must track J in J_History (or override this method)
double Neohookean::GetCurrentRelativeVolume(MPMBase *mptr) const
{   return mptr->GetHistoryDble(J_History);
}


