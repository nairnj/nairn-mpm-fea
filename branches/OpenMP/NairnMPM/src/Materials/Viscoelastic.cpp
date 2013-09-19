/********************************************************************************
    Viscoelastic.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Feb 5 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/Viscoelastic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark Viscoelastic::Constructors and Destructors

// Constructors
Viscoelastic::Viscoelastic()
{
}

// Constructors
Viscoelastic::Viscoelastic(char *matName) : MaterialBase(matName)
{
    int i;

    ntaus=-1;
    Gk=NULL;
    tauk=NULL;
    G0=0.;
    currentGk=0;
    currentTauk=0;
    for(i=0;i<VISCO_PROPS;i++)
        read[i]=0;
}

#pragma mark Viscoelastic::Initialization

// print mechanical properties to output window
void Viscoelastic::PrintMechanicalProperties(void) const
{
    PrintProperty("K",K,"");
	PrintProperty("G0",G0,"");
	PrintProperty("ntaus",(double)ntaus,"");
    cout <<  endl;
	
	int i;
    for(i=0;i<ntaus;i++)
	{	PrintProperty("i",(double)i,"");
		PrintProperty("Gk",Gk[i],"");
		PrintProperty("tauk",1000.*tauk[i],"ms");
        cout << endl;
    }
	
	PrintProperty("a",aI,"");
    cout << endl;
}
    
// Read material properties
char *Viscoelastic::InputMat(char *xName,int &input)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"G0")==0)
        return((char *)&G0);
    
    else if(strcmp(xName,"K")==0)
    {	read[VK_PROP]=1;
        return((char *)&K);
    }
    
    else if(strcmp(xName,"alpha")==0)
    {	read[VA_PROP]=1;
        return((char *)&aI);
    }
    
    else if(strcmp(xName,"ntaus")==0)
    {	input=INT_NUM;
        return((char *)&ntaus);
    }
    
    else if(strcmp(xName,"Gk")==0)
    {	if(Gk==NULL)
        {   if(ntaus<=0)
				ThrowSAXException("Gk found before number of taus specified.");
            Gk=new double[ntaus];
        }
        currentGk++;
        if(currentGk>ntaus)
			ThrowSAXException("Too many Gk's given.");
        return((char *)&Gk[currentGk-1]);
    }
    
    else if(strcmp(xName,"tauk")==0)
    {	if(tauk==NULL)
        {   if(ntaus<=0)
				ThrowSAXException("tauk found before number of taus specified.");
            tauk=new double[ntaus];
        }
        currentTauk++;
        if(currentTauk>ntaus)
			ThrowSAXException("Too many tauk's given.");
        return((char *)&tauk[currentTauk-1]);
    }
    
    return(MaterialBase::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *Viscoelastic::VerifyAndLoadProperties(int np)
{
	// check properties
    if(currentGk<ntaus || currentTauk<ntaus)
		return "Insufficient Gk or tauk for expected number of taus.";
    
    if(ntaus<0)
		return "Number of taus was never entered.";
    
    if(!read[VK_PROP] && !read[VA_PROP])  // Oleg changed !! to &&
		return "Required bulk modulus or thermal expansion not given.";
    
    int i;
    
    // zero time shear modulus
    Ge=G0;
    dGe=0.0;
    for(i=0;i<ntaus;i++)
    {   Ge+=Gk[i];
        dGe-=Gk[i]/tauk[i];
    }
	
	// From MPa to Pa
	Ge*=1.e6;
	dGe*=1.e6;
	Ke=K*1.e6;
	
	// to absolute CTE and CME
	CTE=1.e-6*aI;
	CME=betaI*concSaturation;

	// call super class
	return MaterialBase::VerifyAndLoadProperties(np);
}

// plane stress not allowed in viscoelasticity
void Viscoelastic::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM || np==AXISYMMETRIC_MPM)
	{	throw CommonException("Viscoelastic materials require 2D plane strain or 3D MPM analysis",
							  "Viscoelastic::ValidateForUse");
	}
	
	//call super class (why can't call super class?)
	return MaterialBase::ValidateForUse(np);
}

// create and return pointer to history variables
// initialize all to zero
char *Viscoelastic::InitHistoryData(void)
{
    if(ntaus==0) return NULL;
    
    // allocate array of double pointers (3)
	int blocks = fmobj->IsThreeD() ? 6 : 3 ;
    char *p = new char[sizeof(double *)*blocks];
	
    double **h = (double **)p;
    
    // for each allocate ntaus doubles
    //	h[ij_HISTORY][0] to h[ij_HISTORY][ntaus-1] can be read
    //	in MPMConstLaw() by casting mptr->GetHistoryPtr() pointer as
    //		double **h=(double **)(mptr->GetHistoryPtr())
    h[XX_HISTORY] = new double[ntaus];
    h[YY_HISTORY] = new double[ntaus];
    h[XY_HISTORY] = new double[ntaus];
	if(blocks==6)
    {	h[ZZ_HISTORY] = new double[ntaus];
		h[XZ_HISTORY] = new double[ntaus];
		h[YZ_HISTORY] = new double[ntaus];
	}
    
    // initialize to zero
    int k;
    for(k=0;k<ntaus;k++)
    {	h[XX_HISTORY][k] = 0.;
        h[YY_HISTORY][k] = 0.;
        h[XY_HISTORY][k] = 0.;
		if(blocks==6)
		{	h[ZZ_HISTORY][k] = 0.;
			h[XZ_HISTORY][k] = 0.;
			h[YZ_HISTORY][k] = 0.;
		}
    }
    
    return p;
}

#pragma mark Viscoelastic::Methods

/* For 2D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, history variables,
		strain energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
   For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
*/
void Viscoelastic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
        double dvzz,double delTime,int np,void *properties,ResidualStrains *res) const
{
    /* ---------------------------------------------------
        Add to total strain
        G and K are in MPa and should multiplied by a factor 1.0e6
        delTime (sec)
        tau[k] (sec)
    */
	Tensor *ep = mptr->GetStrainTensor();
    ep->xx += dvxx;
    ep->yy += dvyy;
    double dgam = dvxy+dvyx;
    ep->xy += dgam;
	double dwrotxy = dvyx-dvxy;
    
    // save initial stresses
	Tensor *sp = mptr->GetStressTensor();
    Tensor st0 = *sp;

    /* ---------------------------------------------------
        find stress (Units N/m^2  cm^3/g)
		stress increament due to G(t) and corresponding strains 
    */
    
    // viscous part of stress update
    double **h =(double **)(mptr->GetHistoryPtr());
    int k;
    double dsigmaGxx = 0;
	double dsigmaGyy = 0;
	double dsigmaGxy = 0;
    for (k=0;k<ntaus;k++)
    {   dsigmaGxx += h[XX_HISTORY][k];
        dsigmaGyy += h[YY_HISTORY][k];
        dsigmaGxy += h[XY_HISTORY][k];
    }
    dsigmaGxx = 0.5*delTime*dGe*dvxx + dsigmaGxx*delTime;
    dsigmaGyy = 0.5*delTime*dGe*dvyy + dsigmaGyy*delTime;
    dsigmaGxy = 0.5*delTime*dGe*dgam + dsigmaGxy*delTime;

	// stress increment (specific stress)
    double factor23 = 2.0/3.0;
	double er = CTE*res->dT;
    if(DiffusionTask::active) er += CME*res->dC;
	double dilate = Ke*(dvxx+dvyy-3.*er);
	
    // elastic and viscous components of stree update
	double dsxxv = factor23*(2.*dsigmaGxx-dsigmaGyy)/rho;
    double dsxxe = (dilate + factor23*Ge*(2.*dvxx-dvyy))/rho;
	double dsyyv = factor23*(2.*dsigmaGyy-dsigmaGxx)/rho;
	double dsyye = (dilate + factor23*Ge*(2.*dvyy-dvxx))/rho;
	double dtxyv = dsigmaGxy/rho;
    double dtxye = Ge*dgam/rho;
	double dszzv = -factor23*(dsigmaGyy+dsigmaGxx)/rho;
	double dszze = (dilate - factor23*Ge*(dvyy+dvzz))/rho;
	Hypo2DCalculations(mptr,-dwrotxy,dsxxe+dsxxv,dsyye+dsyyv,dtxye+dtxyv);
    sp->zz += dszze+dszzv;

    // update history variables after stress has been calculated
    for(k=0;k<ntaus;k++)
    {   double tmp = exp(-delTime/tauk[k]);
		double gtk = 1.0e6*Gk[k]/tauk[k];
		h[XX_HISTORY][k] = tmp*(h[XX_HISTORY][k]-gtk*dvxx);
		h[YY_HISTORY][k] = tmp*(h[YY_HISTORY][k]-gtk*dvyy);
		h[XY_HISTORY][k] = tmp*(h[XY_HISTORY][k]-gtk*dgam);
    }
	
	// find energy from work increment
    // add total energy using total stress increment
    // Not sure correct, but does have energy + dissipated equal to total energy increment
	// energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
	dvxx -= er;
	dvyy -= er;
	dvzz -= er;
	double totalEnergy=(st0.xx+0.5*(dsxxe+dsxxv))*dvxx
							+(st0.yy+0.5*(dsyye+dsyyv))*dvyy
							+(st0.xy+0.5*(dtxye+dtxyv))*dgam
							+(st0.xy+0.5*(dszze+dszzv))*dvzz;
    mptr->AddStrainEnergy(totalEnergy);
    
	// Not sure of dissipated enery here
    // visous energy is disspated (not sure of thermo here)
    // dissipated energy using viscous stress increment only (which is negative)
	double dispEnergy = 0.5*dsxxv*dvxx + 0.5*dsyyv*dvyy
                    +0.5*dtxyv*dgam + 0.5*dszzv*dvzz;
    mptr->AddPlastEnergy(-dispEnergy);
	IncrementHeatEnergy(mptr,res->dT,0.,-dispEnergy);
}

/* For 3D MPM analysis, take increments in strain and calculate new
    Particle: nothing yet
    dvij are (gradient rates X time increment) to give deformation gradient change
*/
void Viscoelastic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
        double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np,void *properties,ResidualStrains *res) const
{
    /* ---------------------------------------------------
        Add to total strain
        G and K are in MPa and should multiplied by a factor 1.0e6
        delTime (sec)
        tau[k] (sec)
    */
	Tensor *ep=mptr->GetStrainTensor();
    ep->xx+=dvxx;
    ep->yy+=dvyy;
    ep->zz+=dvzz;
    double dgamxy=dvxy+dvyx;
    ep->xy+=dgamxy;
    double dgamxz=dvxz+dvzx;
    ep->xz+=dgamxz;
    double dgamyz=dvyz+dvzy;
    ep->yz+=dgamyz;
	
	// rotational strain increments (particle updated by Hypo3D)
	double dwrotxy=dvyx-dvxy;
	double dwrotxz=dvzx-dvxz;
	double dwrotyz=dvzy-dvyz;
	    
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;

    /* ---------------------------------------------------
        find stress (Units N/m^2  cm^3/g)
		stress increament due to G(t) and corresponding strains 
    */
    double **h =(double **)(mptr->GetHistoryPtr());
    int k;
    double dsigmaGxx=0;
	double dsigmaGyy=0;
	double dsigmaGzz=0;
	double dsigmaGxy=0;
	double dsigmaGxz=0;
	double dsigmaGyz=0;
    for (k=0;k<ntaus;k++)
    {   dsigmaGxx+=h[XX_HISTORY][k];
        dsigmaGyy+=h[YY_HISTORY][k];
        dsigmaGzz+=h[ZZ_HISTORY][k];
        dsigmaGxy+=h[XY_HISTORY][k];
        dsigmaGxz+=h[XZ_HISTORY][k];
        dsigmaGyz+=h[YZ_HISTORY][k];
    }
	double GeArg=Ge+0.5*delTime*dGe;
    dsigmaGxx = GeArg*dvxx + dsigmaGxx*delTime;
    dsigmaGyy = GeArg*dvyy + dsigmaGyy*delTime;
    dsigmaGzz = GeArg*dvzz + dsigmaGzz*delTime;
    dsigmaGxy = GeArg*dgamxy + dsigmaGxy*delTime;
    dsigmaGxz = GeArg*dgamxz + dsigmaGxz*delTime;
    dsigmaGyz = GeArg*dgamyz + dsigmaGyz*delTime;

	// stress increment (specific stress)
    double factor23=2./3.;
	double er=CTE*res->dT+CME*res->dC;
	double dilate=Ke*(dvxx+dvyy+dvzz-3.*er);
	
	double delsp[6];
	delsp[0] = (dilate+factor23*(2.*dsigmaGxx-dsigmaGyy-dsigmaGzz))/rho;
	delsp[1] = (dilate+factor23*(2.*dsigmaGyy-dsigmaGxx-dsigmaGzz))/rho;
	delsp[2] = (dilate+factor23*(2.*dsigmaGzz-dsigmaGxx-dsigmaGyy))/rho;
	delsp[3] = dsigmaGyz/rho;
	delsp[4] = dsigmaGxz/rho;
	delsp[5] = dsigmaGxy/rho;
	
	// update stress (need to make hypoelastic)
	Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,delsp);

    // update history variables after stress has been calculated
    for(k=0;k<ntaus;k++)
    {   double tmp=exp(-delTime/tauk[k]);
		double gtk=1.0e6*Gk[k]/tauk[k];
		h[XX_HISTORY][k]=tmp*(h[XX_HISTORY][k]-gtk*dvxx);
		h[YY_HISTORY][k]=tmp*(h[YY_HISTORY][k]-gtk*dvyy);
		h[ZZ_HISTORY][k]=tmp*(h[ZZ_HISTORY][k]-gtk*dvzz);
		h[XY_HISTORY][k]=tmp*(h[XY_HISTORY][k]-gtk*dgamxy);
		h[XZ_HISTORY][k]=tmp*(h[XZ_HISTORY][k]-gtk*dgamxz);
		h[YZ_HISTORY][k]=tmp*(h[YZ_HISTORY][k]-gtk*dgamyz);
    }
    
    // how much was viscous energy?
	double delspe[6];
	delspe[0] = (dilate+factor23*Ge*(2.*dvxx-dvyy-dvzz))/rho;
	delspe[1] = (dilate+factor23*Ge*(2.*dvyy-dvxx-dvzz))/rho;
	delspe[2] = (dilate+factor23*Ge*(2.*dvzz-dvxx-dvyy))/rho;
	delspe[3] = Ge*dgamyz/rho;
	delspe[4] = Ge*dgamxz/rho;
	delspe[5] = Ge*dgamxy/rho;
	
	// find energy from work increment
    // add total energy using total stress increment
    // Not sure correct, but does have energy + dissipated equal to total energy increment
	// energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
	dvxx-=er;
	dvyy-=er;
	dvzz=-er;
	double totalEnergy = (st0.xx+0.5*delsp[0])*dvxx + (st0.yy+0.5*delsp[1])*dvyy
						+ (st0.zz+0.5*delsp[2])*dvzz + (st0.yz+0.5*delsp[3])*dgamyz
						+ (st0.xz+0.5*delsp[4])*dgamxz + (st0.xy+0.5*delsp[5])*dgamxy;
    mptr->AddStrainEnergy(totalEnergy);
    
    // Not sure of thermo here
    // dissipated energy using viscous stress increment only (which is negative)
    for(k=0;k<6;k++) delsp[k] -= delspe[k];
	double dispEnergy = 0.5*delsp[0]*dvxx + 0.5*delsp[1]*dvyy + 0.5*delsp[2]*dvzz
                        + 0.5*delsp[3]*dgamyz + 0.5*delsp[4]*dgamxz + 0.5*delsp[5]*dgamxy;
    mptr->AddPlastEnergy(-dispEnergy);
 	IncrementHeatEnergy(mptr,res->dT,0.,-dispEnergy);
}

// convert J to K using isotropic method
Vector Viscoelastic::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{	double nuLS = (3.*Ke-2.*Ge)/(6.*Ke+2.*Ge);
	return IsotropicJToK(d,C,J0,np,nuLS,Ge*1.e-6);
}

#pragma mark Viscoelastic::Accessors

// return material type
const char *Viscoelastic::MaterialType(void) const { return "Viscoelastic"; }

// Return the material tag
int Viscoelastic::MaterialTag(void) const { return VISCOELASTIC; }

/* Calculate wave speed in mm/sec (because G in MPa and rho in g/cm^3)
	Uses sqrt((K +4Ge/3)/rho) which is probably the maximum wave speed possible
*/
double Viscoelastic::WaveSpeed(bool threeD,MPMBase *mptr) const { return sqrt(1.e3*(Ke + 4.*Ge/3.)/rho); }

// Should support archiving history - if it is useful
double Viscoelastic::GetHistory(int num,char *historyPtr) const { return (double)0; }

