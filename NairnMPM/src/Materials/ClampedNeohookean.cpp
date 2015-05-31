/********************************************************************************
	ClampedNeohookean.hpp
	nairn-mpm-fea

	Created by John Nairn on 2/3/15.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#include "ClampedNeohookean.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"

#pragma mark ClampedNeohookean::Constructors and Destructors

// Constructors
ClampedNeohookean::ClampedNeohookean()
{
}

// Constructors with arguments
ClampedNeohookean::ClampedNeohookean(char *matName) : Neohookean(matName)
{
	critComp = 0.025;
	critTens = 0.0075;
	hardening = 10;
	elasticModel = ELASTIC_DISNEY;
}

#pragma mark ClampedNeohookean::Initialization

// Read material properties
char *ClampedNeohookean::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"CritComp")==0)
        return (char *)&critComp;
    
    else if(strcmp(xName,"CritTens")==0)
        return (char *)&critTens;
    
    else if(strcmp(xName,"xihard")==0)
        return (char *)&hardening;
    
    else if(strcmp(xName,"Elastic")==0)
	{	input = INT_NUM;
        return (char *)&elasticModel;
	}
    
    return(Neohookean::InputMaterialProperty(xName,input,gScaling));
}

// print mechanical properties output window
void ClampedNeohookean::PrintMechanicalProperties(void) const
{
	PrintProperty("Cc",critComp,"");
	PrintProperty("Tc",critTens,"");
	PrintProperty("xi",hardening,"");
	cout << endl;
	
	cout << "Elastic Model: ";
	if(elasticModel==ELASTIC_DISNEY)
		cout << "Disney" << endl;
	else
		cout << "Neo-Hookean" << endl;
	
	// call superclass here if it is not Material base
	Neohookean::PrintMechanicalProperties();
}

// verify settings and some initial calculations
const char *ClampedNeohookean::VerifyAndLoadProperties(int np)
{
	// must enter critical value
	if(critComp<0. || critTens<0.)
	{	return "Clamped Neohookean material requires both CritComp and CritTens to be non-negative";
	}
	
	// find elongation ranged (squared values)
	lamMin2 = (1.-critComp)*(1.-critComp);
	lamMax2 = (1.+critTens)*(1.+critTens);
	
	// Disney has only one pressure term
	if(elasticModel==ELASTIC_DISNEY) UofJOption = J_MINUS_1_SQUARED;
	
	// call super class
	return Neohookean::VerifyAndLoadProperties(np);
}

// plane stress not allowed in viscoelasticity
void ClampedNeohookean::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
	{	throw CommonException("Clamped Neohookean material cannot be used in plane stress MPM yet",
						  "ClampedNeohookean::ValidateForUse");
	}
	
	//call super class (why can't call super class?)
	Neohookean::ValidateForUse(np);
}

#pragma mark Neohookean::History Data

// Store J, which is calculated incrementally, and available for archiving and Jp
// initialize to 1
char *ClampedNeohookean::InitHistoryData(void)
{	double *p = CreateAndZeroDoubles(2);
	p[0] = 1.;
	p[1] = 1.;
	return (char *)p;
}

// archive material data for this material type when requested.
double ClampedNeohookean::GetHistory(int num,char *historyPtr) const
{   double history=0.;
    if(num==1 || num==2)
    {	double *J=(double *)historyPtr;
        history=J[num-1];
    }
    return history;
}

#pragma mark ClampedNeohookean::Methods

/* Take increments in strain and calculate new Particle: strains, rotation strain,
		stresses, strain energy,
	dvij are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
	This material tracks pressure and stores deviatoric stress only in particle stress tensor
 */
void ClampedNeohookean::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// Update total deformation gradient, and calculate trial B
    Tensor Btrial;
	double detDF = IncrementDeformation(mptr,du,&Btrial,np);
    
	// global J
	double J = detDF * mptr->GetHistoryDble(J_History);
	mptr->SetHistoryDble(J);
	
	// convert Btrial to matrix to get eigenvalues and check for clamping
	Matrix3 Belas(Btrial.xx,Btrial.xy,Btrial.xz,
				  Btrial.xy,Btrial.yy,Btrial.yz,
				  Btrial.xz,Btrial.yz,Btrial.zz);
	if(np!=THREED_MPM) Belas.setIs2D(true);
		
	// get Eigenvalues and Eigenvectors
	Vector lam2 = Belas.Eigenvalues();
	
	// clamp eigenvalues if needed
	bool clamped = false;
	if(lam2.x<lamMin2 || lam2.x>lamMax2 || lam2.y<lamMin2 || lam2.y>lamMax2 || lam2.z<lamMin2 || lam2.z>lamMax2)
		clamped = true;
	
	// Get Je and Jp, adjusting if clamped
	double Je,Jp;
	Matrix3 Ucol;
	if(clamped)
	{	// Find Belas = U.LAM.UT
		Ucol = Belas.Eigenvectors(lam2);
		// clamp values now
		if(lam2.x<lamMin2)
			lam2.x = lamMin2;
		else if(lam2.x>lamMax2)
			lam2.x = lamMax2;
		if(lam2.y<lamMin2)
			lam2.y = lamMin2;
		else if(lam2.y>lamMax2)
			lam2.y = lamMax2;
		if(lam2.z<lamMin2)
			lam2.z = lamMin2;
		else if(lam2.z>lamMax2)
			lam2.z = lamMax2;
		Matrix3 UcolT = Ucol.Transpose();
		Matrix3 Lam(lam2.x,0.,0.,lam2.y,lam2.z);
		Matrix3 LamUcolT = Lam*UcolT;
		Belas = Ucol*LamUcolT;
		
		// get Je and Jp
		Je = sqrt(lam2.x*lam2.y*lam2.z);
		Jp = J/Je;
		mptr->SetHistoryDble(JP_HISTORY,Jp);
	}
	else
	{	Jp = mptr->GetHistoryDble(JP_HISTORY);
		Je = J/Jp;
		if(elasticModel==ELASTIC_DISNEY)
			Ucol = Belas.Eigenvectors(lam2);
	}
	
	// store B elastic
	Tensor *sp=mptr->GetStressTensor();
    Tensor *B = mptr->GetAltStrainTensor();
	B->xx = Belas(0,0);
	B->yy = Belas(1,1);
	B->zz = Belas(2,2);
	B->xy = Belas(0,1);
	if(np==THREED_MPM)
	{	B->xz = Belas(0,2);
		B->yz = Belas(1,2);
	}
	
	// change mechanical properties by hardening
	double arg = exp(hardening*(1.-Jp));
	double altGsp = pr.Gsp*arg;
	double altLamesp = pr.Lamesp*arg;

	// account for residual stresses
	double dresStretch,resStretch = GetResidualStretch(mptr,dresStretch,res);
	double Jres23 = resStretch*resStretch;
	double Jres = Jres23*resStretch;
		
	// account for residual stresses relative to elastic J
	double Jeff = Je/Jres;
		
	// for incremental energy, store initial stress and pressure
	Tensor *sporig=mptr->GetStressTensor();
	Tensor st0 = *sporig;
	double Pfinal,p0=mptr->GetPressure();
	
	if(elasticModel==ELASTIC_DISNEY)
	{	// Use model from Disney paper
		
		// Get Cauchy stress/rho0
		double sig[3][3];
		double lam[3];
		lam[0]=sqrt(lam2.x)/resStretch;
		lam[1]=sqrt(lam2.y)/resStretch;
		lam[2]=sqrt(lam2.z)/resStretch;
		for(int i=0;i<3;i++)
		{	for(int j=i;j<3;j++)
			{	sig[i][j] = 0.;
				for(int k=0;k<3;k++)
				{	sig[i][j] += (2.*altGsp*lam[k]*(lam[k]-1)/Jeff + altLamesp*(Jeff-1))*Ucol(i,k)*Ucol(j,k);
				}
			}
		}
		
		// update pressure (*J to get Kirchoff pressure)
		Pfinal = -J*(sig[0][0]+sig[1][1]+sig[2][2])/3.;
		mptr->SetPressure(Pfinal);

		// get and set deviatoric stress
		// find eviatoric (Kirchoff stress)/rho0 = deviatoric (Cauchy stress)J/rho0
		sp->xx = J*sig[0][0]+Pfinal;
		sp->yy = J*sig[1][1]+Pfinal;
		sp->zz = J*sig[2][2]+Pfinal;
		sp->xy = J*sig[0][1];
		if(np==THREED_MPM)
		{	sp->xz = J*sig[0][2];
			sp->yz = J*sig[1][2];
		}
	}
	else
	{	// Use standard neo-Hookean law
		
		// update pressure (*J to get Kirchoff pressure)
		Pfinal = -J*(GetVolumetricTerms(Jeff,altLamesp) + (altGsp/Jeff)*((B->xx+B->yy+B->zz)/(3.*Jres23) - 1.));
		mptr->SetPressure(Pfinal);
		
		// Account for density change in specific stress
		// i.e.. Get (Kirchoff Stress)/rho0
		double GJeff = J*resStretch*altGsp/Je;		// = J*(Jres^(1/3) G/Je) to get Kirchoff
		
		// find deviatoric (Kirchoff stress)/rho0 = deviatoric (Cauchy stress)J/rho0
		double I1third = (B->xx+B->yy+B->zz)/3.;
		sp->xx = GJeff*(B->xx-I1third);
		sp->yy = GJeff*(B->yy-I1third);
		sp->zz = GJeff*(B->zz-I1third);
		sp->xy = GJeff*B->xy;
		if(np==THREED_MPM)
		{	sp->xz = GJeff*B->xz;
			sp->yz = GJeff*B->yz;
		}
	}
	
	// work and residual energies
    double delV = 1. - 1./detDF;                        // total volume change
    double avgP = 0.5*(p0+Pfinal);
    double dilEnergy = -avgP*delV;
	
	// incremental residual energy
	double delVres = 1. - 1./(dresStretch*dresStretch*dresStretch);
	double resEnergy = -avgP*delVres;

	// incremental work energy = shear energy
    double shearEnergy = 0.5*((sp->xx+st0.xx)*du(0,0) + (sp->yy+st0.yy)*du(1,1) + (sp->zz+st0.zz)*du(2,2)+
							  (sp->xy+st0.xy)*(du(0,1)+du(1,0)));
    if(np==THREED_MPM)
    {   shearEnergy += 0.5*((sp->xz+st0.xz)*(du(0,2)+du(2,0)) + (sp->yz+st0.yz)*(du(1,2)+du(2,1)));
    }
	
    // strain energy
    double dU = dilEnergy + shearEnergy;
    mptr->AddWorkEnergyAndResidualEnergy(dU,resEnergy);
	
	// thermodynamics heat and temperature
	// Should find energy dissipated by plasticity and add in third term
	IncrementHeatEnergy(mptr,res->dT,0.,0.);
}
	
#pragma mark ClampedNeohookean::Accessors

// Return the material tag
int ClampedNeohookean::MaterialTag(void) const { return CLAMPEDNEOHOOKEAN; }

// return material type
const char *ClampedNeohookean::MaterialType(void) const { return "Clamped Neohookean Hyperelastic-plastic"; }

// if a subclass material supports artificial viscosity, override this and return TRUE
bool ClampedNeohookean::SupportsArtificialViscosity(void) const { return false; }

//	calculate current wave speed in L/sec. Uses sqrt((K+4G/3)/rho) which is dilational wave speed
double ClampedNeohookean::CurrentWaveSpeed(bool threeD,MPMBase *mptr) const
{	double Jp = mptr->GetHistoryDble(JP_HISTORY);
	double arg = exp(hardening*(1.-Jp));
    return sqrt(arg*(Kbulk+4.*G/3.)/rho);
}
