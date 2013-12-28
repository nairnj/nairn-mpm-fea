/********************************************************************************
    BistableIsotropic.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Apr 11 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/BistableIsotropic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark BistableIsotropic::Constructors and Destructors

// Constructors
BistableIsotropic::BistableIsotropic() {}

// Constructors
BistableIsotropic::BistableIsotropic(char *matName) : IsotropicMat(matName)
{
    int i;
    
    for(i=0;i<BISTABLE_PROPS;i++)
        readbs[i]=0;
    
    dVii=0.;
    dVcrit=0.;
	transState=-1;
    rule=DILATION_RULE;
    reversible=FALSE;
	diff0=0.;
	kCond0=0;
	beta0=0.;
	betad=0.;
	a0=40.;
	ad=0.;
}

#pragma mark BistableIsotropic::Initialization

// Read material properties
char *BistableIsotropic::InputMat(char *xName,int &input)
{
    input=DOUBLE_NUM;
	
    // bulk
    if(strcmp(xName,"K0")==0)
    {	readbs[K0_PROP]=1;
        return((char *)&K0);
    }
    else if(strcmp(xName,"Kd")==0)
    {	readbs[KD_PROP]=1;
        return((char *)&Kd);
    }
	
    // shear
    else if(strcmp(xName,"G0")==0)
    {	readbs[G0_PROP]=1;
        return((char *)&G0);
    }
    else if(strcmp(xName,"Gd")==0)
    {	readbs[GD_PROP]=1;
        return((char *)&Gd);
    }
    
    // cte
    else if(strcmp(xName,"alpha0")==0)
    {	readbs[A0_PROP]=1;
        return((char *)&a0);
    }
    else if(strcmp(xName,"alphad")==0)
    {	readbs[AD_PROP]=1;
        return((char *)&ad);
    }
	
    // cme
    else if(strcmp(xName,"beta0")==0)
    {	readbs[B0_PROP]=1;
        return((char *)&beta0);
    }
    else if(strcmp(xName,"betad")==0)
    {	readbs[BD_PROP]=1;
        return((char *)&betad);
    }
	
    // diffusion
    else if(strcmp(xName,"D0")==0)
    {	readbs[DIFF0_PROP]=1;
        return((char *)&diff0);
    }
    else if(strcmp(xName,"Dd")==0)
    {	readbs[DIFFD_PROP]=1;
        return((char *)&diffd);
    }
	
    // conductivity
    else if(strcmp(xName,"kCond0")==0)
    {	readbs[KCOND0_PROP]=1;
        return((char *)&kCond0);
    }
    else if(strcmp(xName,"kCondd")==0)
    {	readbs[KCONDD_PROP]=1;
        return((char *)&kCondd);
    }
	
    // transitions and offsets
    else if(strcmp(xName,"transition")==0)
    {	readbs[TRANSITION_PROP]=1;
        input=INT_NUM;
        return((char *)&rule);
    }
    else if(strcmp(xName,"critical")==0)
    {	readbs[DVCRIT_PROP]=1;
        return((char *)&dVcrit);
    }
    else if(strcmp(xName,"DeltaVOffset")==0)
    {	readbs[DVOFF_PROP]=1;
        return((char *)&dVii);
    }
    else if(strcmp(xName,"reversible")==0)
    {	reversible=TRUE;
        input=NOT_NUM;
        return((char *)&reversible);
    }
    else if(strcmp(xName,"irreversible")==0)
    {	reversible=FALSE;
        input=NOT_NUM;
        return((char *)&reversible);
    }
	
    return MaterialBase::InputMat(xName,input);
}

// Verify properties and initial calculations
const char *BistableIsotropic::VerifyAndLoadProperties(int np)
{
    // Require initial properties, but second are optional
    //	They equal first if not provide
    if(!readbs[K0_PROP] || !readbs[G0_PROP] || !readbs[A0_PROP])
		return "Initial K0, G0, or alpha0 is missing.";
	if(!readbs[B0_PROP])
		beta0=0.;
    if(!readbs[TRANSITION_PROP] || rule<DILATION_RULE || rule>VONMISES_RULE)
		return "Phase transition rule is missing or invalid.";
    
    // if not provided, no change in property at transition
    if(!readbs[KD_PROP]) Kd=K0;
    if(!readbs[GD_PROP]) Gd=G0;
    if(!readbs[AD_PROP]) ad=a0;
    if(!readbs[BD_PROP]) betad=beta0;
    if(!readbs[DIFFD_PROP]) diffd=diff0;
    if(!readbs[KCONDD_PROP]) kCondd=kCond0;
    
    // no change in heat capacity
    
    // make conductivty specific (N mm^3/(sec-K-g))
    kCond0 *= (1000./rho);
    kCondd *= (1000./rho);
    
    // test validity of each state
    const char *err=CurrentProperties(DEFORMED_STATE,np);
    if(err!=NULL) return err;
	FillElasticProperties(&pr2,np);
	FillTransportProperties(&tr2);
	
    err=CurrentProperties(INITIAL_STATE,np);
    if(err!=NULL) return err;
	FillElasticProperties(&pr,np);
	// transport gets done in material base
	
    // convert strain rules in percent to absolute strains
    if(rule==VONMISES_RULE)
		dVcrit*=1.e6/rho;			// convert MPa to specific stress
	else
		dVcrit/=100.;
    dVii/=100.;
	
	// call super-super class (skip IsotropicMat due to conflicts and Elastic because has none)
	return MaterialBase::VerifyAndLoadProperties(np);
}

// print mechanical properties to output window
void BistableIsotropic::PrintMechanicalProperties(void) const
{
	PrintProperty("Initial:",false);
	PrintProperty("K",K0,"");
	PrintProperty("G",G0,"");
	PrintProperty("a",a0,"");
	cout << endl;
	
 	PrintProperty("Transformed:",false);
	PrintProperty("K",Kd,"");
	PrintProperty("G",Gd,"");
	PrintProperty("a",ad,"");
	cout << endl;
    
    char mline[200];
	if(rule==DILATION_RULE)
    {   sprintf(mline,"Dilation transition at dV = %g%c to V offset = %g%c",100.*dVcrit,'%',100.*dVii,'%');
	}
	else if(rule==DISTORTION_RULE)
    {   sprintf(mline,"Distortion transition at sqrt(0.5*e'ij e'ij) = %g%c",100.*dVcrit,'%');
	}
	else if(rule==VONMISES_RULE)
    {   sprintf(mline,"Distortion transition at sqrt(0.5*s'ij s'ij) = %g MPa",1.e-6*rho*dVcrit);
	}

	cout << mline << endl;
    
    if(reversible)
        cout << "Reversible" << endl;
    else
        cout << "Irreversible" << endl;
}

    
// print transport propertie to output window
void BistableIsotropic::PrintTransportProperties(void) const
{
    char mline[200];
	
	// Diffusion constants
	if(DiffusionTask::active)
	{   sprintf(mline,"D0 =%12.3g   Dd =%12.3f mm^2/sec  csat = %9.5lf",diff0,diffd,concSaturation);
		cout << mline << endl;
	    sprintf(mline,"b0 =%12.6g   bD =%12.6g 1/wt fr",beta0,betad);
		cout << mline << endl;
	}
	// Conductivity constants (Cp is also mJ/(g-K))
	if(ConductionTask::active)
	{   sprintf(mline,"k0 =%12.3g W/(m-K)  kd =%12.3g W/(m-K)  C   =%12.3g J/(kg-K)",
                            rho*kCond0/1000.,rho*kCondd/1000.,heatCapacity);
		cout << mline << endl;
	}
}

// 3D not allowed
void BistableIsotropic::ValidateForUse(int np) const
{	if(np==THREED_MPM || np==AXISYMMETRIC_MPM)
	{	throw CommonException("BistableIsotropic materials cannot do 3D or Axisymmetric MPM analysis",
							  "BistableIsotropic::ValidateForUse");
	}
	
	// call super class (why can't call super class?)
	IsotropicMat::ValidateForUse(np);
}

// calculate properties for give state
const char *BistableIsotropic::CurrentProperties(short newState,int np)
{
    double K;
	const char *err=NULL;
    
    // find properties for state
    if(newState==INITIAL_STATE)
    {	K=K0;
        G=G0;
        aI=a0;
		betaI=beta0;
		diffusionCon=diff0;
		kCond=kCond0;
    }
    else
    {	K=Kd;
        G=Gd;
        aI=ad;
		betaI=betad;
		diffusionCon=diffd;
		kCond=kCondd;
    }
    
    // analysis properties this state
    E=9.e6*K*G/(G + 3.*K);
    nu=(3.*K-2.*G)/(6.*K+2.*G);
	if(DbleEqual(E,0.0))
		return "State with zero modulus is not supported.";
    err=SetAnalysisProps(np,E,E,E,nu,nu,nu,1.e6*G,1.e6*G,1.e6*G,
			1.e-6*aI,1.e-6*aI,1.e-6*aI,betaI*concSaturation,betaI*concSaturation,betaI*concSaturation);
	
	return err;
}

// Single short to hold the current particle state
char *BistableIsotropic::InitHistoryData(void)
{
    // allocate pointer to a single short
    char *p=new char[sizeof(short)];
	
    short *h=(short *)p;
    *h=INITIAL_STATE;
	
    return p;
}

#pragma mark BistableIsotropic::Methods

// fill in transport tensors matrix if needed
void BistableIsotropic::GetTransportProps(MPMBase *mptr,int np,TransportProperties *t) const
{
	short *state=(short *)(mptr->GetHistoryPtr());     // history pointer is short * with state
	*t = *state==INITIAL_STATE ? tr : tr2;
}

/* For 2D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle,
		current state
    dvij are (gradient rates X time increment) to give deformation gradient change
   For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
*/
void BistableIsotropic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
        double dvzz,double delTime,int np,void *properties,ResidualStrains *res) const
{
    // update in current state
    short *state=(short *)(mptr->GetHistoryPtr()),transition=FALSE;
	const ElasticProperties *p = *state==INITIAL_STATE ? &pr : &pr2;
    Elastic::MPMConstLaw(mptr,dvxx,dvyy,dvxy,dvyx,dvzz,delTime,np,(void *)p,res);
	
    // Calculate critical value for transition
    double dmechV,dTrace,ds1,ds2,ds3;
	Tensor *sp=mptr->GetStressTensor();
	Tensor *ep=mptr->GetStrainTensor();
	
	switch(rule)
	{	case DILATION_RULE:
			dmechV=ep->xx+ep->yy+ep->zz;
			break;
		
		case DISTORTION_RULE:
			// find deviatoric strain inner product
			dmechV=ep->xx+ep->yy+ep->zz;
			dTrace=dmechV/3.;
			dmechV=(ep->xx-dTrace)*(ep->xx-dTrace);
			dmechV+=(ep->yy-dTrace)*(ep->yy-dTrace);
			dmechV+=(ep->zz-dTrace)*(ep->zz-dTrace);
			// 2 * (0.5 gamma) * (0.5 gamma) to true shear strain terms
			dmechV+=0.5*ep->xy*ep->xy;
			dmechV=sqrt(0.5*dmechV);
			break;
		
		case VONMISES_RULE:
			// sqrt(((#sxx-#syy)^2+(#syy-#szz)^2+(#szz-#sxx)^2+6*#sxy^2)/2)
			ds1=sp->xx - sp->yy;
			ds2=sp->xx - sp->zz;
			ds3=sp->zz - sp->yy;
			dmechV=sqrt(0.5*(ds1*ds1+ds2*ds2+ds3*ds3+6*sp->xy*sp->xy));
			break;
		
		default:
			dmechV=0.;
			break;
	}
	
	// is there a transition
    if(dmechV>=dVcrit)
    {	if(*state!=DEFORMED_STATE)
        {   // Transition to dilated state
            *state=DEFORMED_STATE;
            transition=TRUE;
        }
    }
    else if(*state==DEFORMED_STATE)
    {	if(reversible)
        {   // Transition back to undilated state (if reversible)
            *state=INITIAL_STATE;
            transition=TRUE;
        }
    }
    
    // instantaneous change in stress at constant strain (DILATION_RULE only)
    if(transition && rule==DILATION_RULE)
    {	// find changed stress by current constitutive law
		double normOffset,alphazz,betazz;
		if(*state==INITIAL_STATE)
		{	normOffset = 0.;
			alphazz = 1.e-6*a0;
			betazz = 1.e-6*beta0;
			p = &pr;
		}
		else
		{	normOffset = dVii/3.;
			alphazz = 1.e-6*ad;
			betazz = 1.e-6*betad;
			p = &pr2;
		}

        //LoadMechanicalProps(mptr,np);
		double er = p->alpha[1]*(mptr->pPreviousTemperature-thermal.reference)
						+ p->beta[1]*(mptr->pPreviousConcentration-DiffusionTask::reference);
		double erzz = alphazz*(mptr->pPreviousTemperature-thermal.reference)
						+ betazz*(mptr->pPreviousConcentration-DiffusionTask::reference);
		double exx,eyy;

		if(np==PLANE_STRAIN_MPM)
		{	// need effective offset here for (see JAN048-63)
			exx=ep->xx-normOffset*(1+nu)-er;
			eyy=ep->yy-normOffset*(1+nu)-er;
		}
		else
        {	exx=ep->xx-normOffset-er;
			eyy=ep->yy-normOffset-er;
		}
        sp->xx=p->C[1][1]*exx + p->C[1][2]*eyy;
        sp->yy=p->C[1][2]*exx + p->C[2][2]*eyy;
        sp->xy=p->C[3][3]*ep->xy;
		if(np==PLANE_STRAIN_MPM)
        {	double ezz=normOffset+erzz;			// because ezz=0 and now do not want effective properties
        	exx=ep->xx-ezz;
			eyy=ep->yy-ezz;
			sp->zz=p->C[4][1]*(exx + eyy) - p->C[4][4]*ezz;
		}
		else
			ep->zz=p->C[4][1]*(exx + eyy)+normOffset+erzz;
    }
}

#pragma mark BistableIsotropic::Accessors

// Return the material tag
int BistableIsotropic::MaterialTag(void) const { return BISTABLEISO; }

/*	calculate wave speed in mm/sec (because K,G in MPa and rho in g/cm^3)
	Uses max sqrt((K +4G/3)/rho) which is dilational wave speed
*/
double BistableIsotropic::WaveSpeed(bool threeD,MPMBase *mptr) const
{ return fmax(sqrt(1.e9*(K0+4.*G0/3.)/rho),sqrt(1.e9*(Kd+4.*Gd/3.)/rho));
}

// maximum diffusion coefficient in cm^2/sec (diff in mm^2/sec)
double BistableIsotropic::MaximumDiffusion(void) const { return max(diffd,diff0)/100.; }

// maximum diffusivity in cm^2/sec
// specific k is mJ mm^2/(sec-K-g) and Cp is mJ/(g-K) so k/C = mm^2/sec * 1e-2 = cm^2/sec
double BistableIsotropic::MaximumDiffusivity(void) const { return 0.01*max(kCondd,kCond0)/heatCapacity; }

// return material type
const char *BistableIsotropic::MaterialType(void) const { return "Bistable Isotropic"; }

// archive material data for this material type when requested.
double BistableIsotropic::GetHistory(int num,char *historyPtr) const
{
    double history=0.;
    short *state;
    
    if(num==1)
    {	state=(short *)historyPtr;
        history=(double)*state;
    }
    return history;
}
