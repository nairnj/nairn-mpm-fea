//
//  GeoMaterials.cpp
//  NairnMPM
//
//  Created by Raydel Lorenzo on 6/11/12.
//  Copyright (c) 2012 __Geotecnia-UnB__. All rights reserved.
//

#include <iostream>

/********************************************************************************
 Dependencies
 Orthotropic.hpp (TranIsotropic.hpp, MaterialBase.hpp)
 ********************************************************************************/

#include "GeoMaterials.hpp"

#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark GeoMaterials::Constructors and Destructors

// Constructors
GeoMaterials::GeoMaterials() {}

// Constructors
GeoMaterials::GeoMaterials(char *matName) : IsotropicMat(matName)
{
	readYield=FALSE;
}

#pragma mark GeoMaterials::Initialization

// Read material properties
char *GeoMaterials::InputMat(char *xName,int &input)
{
    //Ground surface level 
    if(strcmp(xName,"yref")==0)
    {   input=DOUBLE_NUM;
        readyref=true;
        return((char *)&yref);
    }
    
    //Over Consolidation Ratio
    if(strcmp(xName,"OCR")==0)
    {   input=DOUBLE_NUM;
        readOCR=true;
        return((char *)&OCR);
    }
    
    //Earth's pressure coefficient at rest
    if(strcmp(xName,"k0")==0)
    {   input=DOUBLE_NUM;
        readk0=true;
        return((char *)&k0);
    }
    
    return(IsotropicMat::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *GeoMaterials::VerifyProperties(int np)
{	
    //check properties
    if(!readyref) return " The yref is missing ";
    
    if(!readOCR) return " The OCR is missing ";
    
    if(!readk0) return " The k0 is missing ";
    
	// call super class
	return IsotropicMat::VerifyProperties(np);
}

// Private properties used in constitutive law
// For variable shear and bulk moduli, subclass can overrive
//		LoadMechanicalProps(MPMBase *mptr,int np) and set new
//		Gred and Kred
// Here gets yldred, Gred, Kred, psRed, psLr2G, and psKred
void GeoMaterials::InitialLoadMechProps(int makeSpecific,int np)
{
	// reduced prooperties
    yldred = yield*1.e6/rho;
	Gred = C66/rho;
	Kred = C33/rho - 4.*Gred/3.;	// from C33 = lambda + 2G = K + 4G/3
	
	// these are terms for plane stress calculations only
	psRed = 1./(Kred/(2.*Gred) + 2./3.);			// (1-2nu)/(1-nu) for plane stress
	psLr2G = (Kred/(2.*Gred) - 1./3.)*psRed;		// nu/(1-nu) to find ezz
	psKred = Kred*psRed;							// E/(3(1-v)) to find lambda
	
	// nothing needed from superclasses
}



// If super class needs to override this method, always save the first double
//		for the cumulative plastic strain.
char *GeoMaterials::InitHistoryData(void)
{
	double *p=new double;
	*p=0.;
	return (char *)p;
}

// If needed, a material can initialize particle state
// For example, ideal gas initializes to base line pressure
void GeoMaterials::SetInitialParticleState(MPMBase *mptr,int np)
{    
    Tensor *sp=mptr->GetStressTensor();
    
    
    sp->yy = -(yref - mptr->origpos.y) * 10.;           //multiplied by rho but divided by rho for specific stress, then without rho.
    sp->xx = (sp->yy * k0);
    sp->zz = (sp->yy * k0);
    
}



#pragma mark GeoMaterials:Methods

/* For 2D MPM analysis, take increments in strain and calculate new                         
 Particle: strains, rotation strain, plastic strain, stresses, strain energy, 
 plastic energy, dissipated energy, angle
 dvij are (gradient rates X time increment) to give deformation gradient change
 This is general analysis for isotropic plastic material with plastic volume deformation as internal variable of hardening. Subclass must define
 GetFtrial() and GetDfDsigmageo() and optionally can override more. Those methods
 require history dependent variables and rates (e.g. cum. plastic volumetric strain (dvplas).
 */
void GeoMaterials::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,double dvzz,double delTime,int np)
{
    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
	//  (note: using unreduced terms in CTE3 and CME3)
	double eres=CTE3*ConductionTask::dTemperature;
	if(DiffusionTask::active)
		eres+=CME3*DiffusionTask::dConcentration;
	
	// here dvij is total strain increment, dexxr is relative strain by subtracting off eres
    double dexxr=dvxx-eres;			// trial dexx=dvxx
    double deyyr=dvyy-eres;			// trial deyy=dvyy
	double dezzr=-eres;				// used in plane strain only, where trial dezz=0
    double dgxy=dvxy+dvyx;			// no need to substract residual strain
	
	// Trial update assuming elastic response
	
    // Elastic stress increment
	Tensor *ep=mptr->GetStrainTensor();
	Tensor *sp=mptr->GetStressTensor();
    Tensor dsigma,stkgeo,stkgeotrial,stk,st0=*sp;
    
    //Negative strees because in Geotechnics positive is compression. 
    stkgeotrial.xx = -st0.xx;
    stkgeotrial.yy = -st0.yy;
    stkgeotrial.zz = -st0.zz;
    stkgeotrial.yz = -st0.yz;
    stkgeotrial.xz = -st0.xz;
    stkgeotrial.xy = -st0.xy;
    
	if(np==PLANE_STRAIN_MPM)
		dsigma = GetElasticIncrement(mptr, np, &stkgeotrial, dexxr, deyyr, 0., dgxy, 0., 0.);
	else
        cout << " Not implemented for this material ";
    
    // Increment tensor of stress for trial
    stk.xx = st0.xx + dsigma.xx;
    stk.yy = st0.yy + dsigma.yy;
	stk.zz = st0.zz + dsigma.zz;
    stk.yz = st0.yz + dsigma.yz;
    stk.xz = st0.xz + dsigma.xz;
    stk.xy = st0.xy + dsigma.xy;
    
    // Negative Stress Tensor to use in the GeoMaterials constituve law, becouse compresion is positive for calculation. But in the final he gave positive tensor
    stkgeo.xx = -stk.xx;
    stkgeo.yy = -stk.yy;
    stkgeo.zz = -stk.zz;
    stkgeo.yz = -stk.yz;
    stkgeo.xz = -stk.xz;
    stkgeo.xy = -stk.xy;
   
    UpdateDefVolPlas(mptr);                                         //update volumetric plastic strain
    ftrial = GetFtrial(mptr,np,&stkgeo);                           //Call subclass and Calcule f-> in funtion of material

	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
        
		// Add input strain increment to elastic strain on particle
		ep->xx += dvxx;
		ep->yy += dvyy;
		ep->xy += dgxy;
		double dwrotxy=dvyx-dvxy;
        
		// increment stress (Units N/m^2  cm^3/g)
		Hypo2DCalculations(mptr,-dwrotxy,dsigma.xx,dsigma.yy,dsigma.xy);
		
		// out of plane
			sp->zz=stk.zz;

		
		// strain energy (by midpoint rule)
		mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr + (st0.yy+sp->yy)*deyyr
                                   + (st0.xy+sp->xy)*dgxy));
		if(np==PLANE_STRAIN_MPM)
			mptr->AddStrainEnergy(0.5*(st0.zz+sp->zz)*dezzr);
		
		// give material chance to update history variables that change in elastic updates
		ElasticUpdateFinished(mptr,np,delTime);
		
		return;
    }
	
	// Return Derivative of F with sigma for actual state GEO
    Tensor dfdsigma = GetDfDsigmaGeo(mptr,np,&stkgeo);
    
    
    // Return Plastic Multiplicator with sigma for actual state GEO
    double multplas = GetMultPlast(mptr,np, &stkgeo, dexxr,deyyr,0.,dgxy,0.,0.);
    	    
    // Plastic strain increments on particle
    double dexxp = multplas * dfdsigma.xx;
    double deyyp = multplas * dfdsigma.yy;
    double dezzp = multplas * dfdsigma.zz;
    double dgxyp = multplas * dfdsigma.xy;     
	
    //delta volumetric plastic strain GEO
    ddvplas = dexxp + deyyp + dezzp;
    
    Tensor *eplast=mptr->GetPlasticStrainTensor();
    eplast->xx += dexxp;
    eplast->yy += deyyp;
    eplast->xy += dgxyp;
    eplast->zz += dezzp;
    
    // Elastic strain increments on particle
    ep->xx += (dvxx-dexxp);
    ep->yy += (dvyy-deyyp);
    dgxy -= dgxyp;
    ep->xy += dgxy;
    ep->zz -= dezzp;
	
	// rotational strain
	double dwrotxy=dvyx-dvxy;
	
    // Elastic strain increment minus the residual terms by now subtracting plastic parts
    dexxr -= dexxp;
    deyyr -= deyyp;
	//dgxy -= dgxyp;			// done above
	//dezzr -= dezzp;			// plane strain only done above
    
	// calculate stress correction due to plastic deformation and increment particle stresses
    Tensor dsigmap = GetElasticIncrement(mptr, np, &stkgeo, dexxp, deyyp, dezzp, dgxyp, 0., 0.);
	
    dsigma.xx -= dsigmap.xx;
    dsigma.yy -= dsigmap.yy;
    dsigma.zz -= dsigmap.zz;
    dsigma.xy -= dsigmap.xy;

    sp->zz += dsigma.zz;
    
	Hypo2DCalculations(mptr,-dwrotxy,dsigma.xx,dsigma.yy,dsigma.xy);
	
    // Elastic energy density
    mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr
                               + (st0.yy+sp->yy)*deyyr
                               + (st0.xy+sp->xy)*dgxy));
    
    // Plastic energy increment
	double dispEnergy=0.5*((st0.xx+sp->xx)*dexxp
                           + (st0.yy+sp->yy)*deyyp
                           + (st0.xy+sp->xy)*dgxyp);
	
	if(np==PLANE_STRAIN_MPM)
    {	mptr->AddStrainEnergy(0.5*(st0.zz+sp->zz)*dezzr);
		dispEnergy += 0.5*(st0.zz+sp->zz)*dezzp;
	}
    
	// add plastic energy to the particle
	mptr->AddDispEnergy(dispEnergy);
    mptr->AddPlastEnergy(dispEnergy);

    // update internal variables GEO
    StoragePlasticInternal(mptr);
}

/* For 3D MPM analysis, take increments in strain and calculate new
 Particle: strains, rotation strain, stresses, strain energy, angle
 dvij are (gradient rates X time increment) to give deformation gradient change
 */
void GeoMaterials::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
                                double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
	//  (note: using unreduced terms in CTE3 and CME3)
	double eres=CTE3*ConductionTask::dTemperature;
	if(DiffusionTask::active)
		eres+=CME3*DiffusionTask::dConcentration;
	
	// here dvij is total strain increment, dexxr is relative strain by subtracting off eres
    double dexxr=dvxx;
    double deyyr=dvyy;
	double dezzr=dvzz;				// use in plane strain only
    double dgxy=dvxy+dvyx;
    double dgxz=dvxz+dvzx;
	double dgyz=dvyz+dvzy;
	
	// rotational strain increments (particle updated by Hypo3D)
	double dwrotxy=dvyx-dvxy;
	double dwrotxz=dvzx-dvxz;
	double dwrotyz=dvzy-dvyz;
	
    // Elastic stress increment
	Tensor *ep=mptr->GetStrainTensor();
	Tensor *sp=mptr->GetStressTensor();
    Tensor stk,st0=*sp;
    double dsig[6];
    Tensor stkgeotrial, stkgeo;          // this tensor is just used in the Methods of Geomaterial because compresion is positive
    // Negative Stress Tensor to use in the GeoMaterials constituve law, because compresion is positive for calculation. But in the final he gave positive tensor for traction   
    stkgeotrial.xx = -st0.xx;
    stkgeotrial.yy = -st0.yy;
    stkgeotrial.zz = -st0.zz;
    stkgeotrial.yz = -st0.yz;
    stkgeotrial.xz = -st0.xz;
    stkgeotrial.xy = -st0.xy;
    
    Tensor dsigma = GetElasticIncrement(mptr, np, &stkgeotrial, dexxr,deyyr,dezzr,dgxy,dgyz,dgxz);
    
    // this array is just for method Hypoeslastic
    dsig[XX] = dsigma.xx;
    dsig[YY] = dsigma.yy;
    dsig[ZZ] = dsigma.zz;
    dsig[YZ] = dsigma.yz;
    dsig[XZ] = dsigma.xz;
    dsig[XY] = dsigma.xy;
    
    // Increment tensor of stress for trial
    stk.xx = st0.xx + dsigma.xx;
    stk.yy = st0.yy + dsigma.yy;
	stk.zz = st0.zz + dsigma.zz;
    stk.yz = st0.yz + dsigma.yz;
    stk.xz = st0.xz + dsigma.xz;
    stk.xy = st0.xy + dsigma.xy;
        
    // Negative Stress Tensor to use in the GeoMaterials constituve law, becouse compresion is positive for calculation. But in the final he gave positive tensor    
    stkgeo.xx = -stk.xx;
    stkgeo.yy = -stk.yy;
    stkgeo.zz = -stk.zz;
    stkgeo.yz = -stk.yz;
    stkgeo.xz = -stk.xz;
    stkgeo.xy = -stk.xy;
    
    //Update Plastic volumetric strain
    UpdateDefVolPlas(mptr);
    
    //get yield surface GEO
    ftrial =  GetFtrial(mptr,np,&stkgeo);
    
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
        
		// Add to total strain
		ep->xx+=dvxx;
		ep->yy+=dvyy;
		ep->zz+=dvzz;
		ep->xy+=dgxy;
		ep->xz+=dgxz;
		ep->yz+=dgyz;
		
		// update stress (need to make hypoelastic)
		Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,dsig);
		
		// strain energy
		mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr
                                   + (st0.yy+sp->yy)*deyyr
                                   + (st0.zz+sp->zz)*dezzr
                                   + (st0.yz+sp->yz)*dgyz
                                   + (st0.xz+sp->xz)*dgxz
                                   + (st0.xy+sp->xy)*dgxy));
        
		// give material chance to update history variables that change in elastic updates
		ElasticUpdateFinished(mptr,np,delTime);
		return;
	}
	
    // Find direction of plastic strain and lambda for this plastic state
    // Return Derivative of F with sigma for actual state GEO
    Tensor dfdsigma = GetDfDsigmaGeo(mptr,np,&stkgeo);
        
    // Return Plastic Multiplicator with sigma for actual state GEO
    double multplas = GetMultPlast(mptr, np, &stkgeo, dexxr,deyyr,dezzr,dgxy,dgyz,dgxz);
        
	// Now have multplas, finish update on this particle    
    // Plastic strain increments on particle GEO
    double dexxp = multplas * dfdsigma.xx;
    double deyyp = multplas * dfdsigma.yy;
	double dezzp = multplas * dfdsigma.zz;
    double dgxyp = multplas * dfdsigma.xy;
    double dgxzp = multplas * dfdsigma.xz;
    double dgyzp = multplas * dfdsigma.yz;
    
    //delta volumetric plastic strain GEO
    ddvplas = dexxp + deyyp + dezzp;

    
	Tensor *eplast=mptr->GetPlasticStrainTensor();
	eplast->xx += dexxp;
    eplast->yy += deyyp;
	eplast->zz += dezzp;
    eplast->xy += dgxyp;
	eplast->xz += dgxzp;
    eplast->yz += dgyzp;

    // Elastic strain increments on particle
    ep->xx += (dvxx-dexxp);
    ep->yy += (dvyy-deyyp);
    ep->zz += (dvzz-dezzp);
    dgxy -= dgxyp;
    ep->xy += dgxy;
    dgxz -= dgxzp;
    ep->xz += dgxz;
    dgyz -= dgyzp;
	ep->yz += dgyz;
	
    // Elastic strain increment minus the residual terms by now subtracting plastic parts
    dexxr -= dexxp;
    deyyr -= deyyp;
	dezzr -= dezzp;				// plain strain only
	//dgxy, dgxz, dgyz done above
    
	// increment particle stresses
    // decrease trial stress because plastic strain
    Tensor dsigmap = GetElasticIncrement(mptr, np, &stkgeo, dexxp, deyyp, dezzp, dgxyp, dgyzp, dgxzp);
    
	dsigma.xx -= dsigmap.xx;
	dsigma.yy -= dsigmap.yy;
	dsigma.zz -= dsigmap.zz;
	dsigma.yz -= dsigmap.yz;
	dsigma.xz -= dsigmap.xz;
	dsigma.xy -= dsigmap.xy;
    
    //this array is just for method hypoeslastic
    dsig[XX] = dsigma.xx;
    dsig[YY] = dsigma.yy;
    dsig[ZZ] = dsigma.zz;
    dsig[YZ] = dsigma.yz;
    dsig[XZ] = dsigma.xz;
    dsig[XY] = dsigma.xy;
    
    
 /*   Tensor sigmac = CorrectDrift(mptr, np, &dsigma, &stkgeotrial);            // Correcting for drift
    
    dsig[XX] = sigmac.xx - stkgeotrial.xx;
    dsig[YY] = sigmac.yy - stkgeotrial.yy;
    dsig[ZZ] = sigmac.zz - stkgeotrial.zz;
    dsig[XY] = sigmac.xy - stkgeotrial.xy;
    dsig[YZ] = sigmac.yz - stkgeotrial.yz;
    dsig[XZ] = sigmac.xz - stkgeotrial.xz;
    
    cout <<"dsig " << dsig[XX];*/
    
	Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,dsig);
	
    // Elastic energy density
	mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr
                               + (st0.yy+sp->yy)*deyyr
                               + (st0.zz+sp->zz)*dezzr
                               + (st0.yz+sp->yz)*dgyz
                               + (st0.xz+sp->xz)*dgxz
                               + (st0.xy+sp->xy)*dgxy));
    
    // Plastic energy increment
	double dispEnergy=0.5*(0.5*((st0.xx+sp->xx)*dexxp
                                + (st0.yy+sp->yy)*deyyp
                                + (st0.zz+sp->zz)*dezzp
                                + (st0.yz+sp->yz)*dgyzp
                                + (st0.xz+sp->xz)*dgxzp
                                + (st0.xy+sp->xy)*dgxyp));
    
	// add plastic energy to the particle
	mptr->AddDispEnergy(dispEnergy);
    mptr->AddPlastEnergy(dispEnergy);

                     
    // update internal variables GEO
    StoragePlasticInternal(mptr);
    
}

#pragma mark GeoMaterials::Custom Methods

/* Get internal variables*/


//Actual plastic volumetric strain  GEO
void GeoMaterials::UpdateDefVolPlas(MPMBase * mptr)
{
    dvplas=mptr->GetHistoryDble();
    ddvplas=0.;
    
}


//Storage delta plastic volumetric strain  GEO
void GeoMaterials::StoragePlasticInternal(MPMBase * mptr)
{

    mptr->SetHistoryDble(dvplas+ddvplas);                    //storage in 0 position the volumetric plastic strain
    
}


// material can override if history variable changes during elastic update (e.g., history dependent on plastic strain rate now zero)
void GeoMaterials::ElasticUpdateFinished(MPMBase *mptr,int np,double delTime)
{
}

#pragma mark GeoMaterials::Accessors



