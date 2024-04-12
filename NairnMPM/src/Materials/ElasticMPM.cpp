/********************************************************************************
    ElasticMPM.cpp - more Elastic for MPM code
    nairn-mpm-fea
    
    Created by John Nairn on Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/Elastic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Elements/ElementBase.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/MPMWarnings.hpp"

// for numerical softening methods
#include "SofteningLaw.hpp"

#pragma mark Elastic::Initialization

// print any properties common to all MPM material types
void Elastic::PrintCommonProperties(void) const
{
	if(useLargeRotation)
		cout << "Large rotation method for hypoelastic materials" << endl;
	
	MaterialBase::PrintCommonProperties();
}
	
#pragma mark Elastic::Methods

/* Take increments in strain and calculate new Particle: strains, rotation strain,
		stresses, strain energy,
	dvij are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
	This method only called by TransIsotropic and subclasses
 */
void Elastic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 dv,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{	if(useLargeRotation)
		LRConstitutiveLaw(mptr,dv,delTime,np,properties,res);
	else
	{	// increment deformation gradient
		HypoIncrementDeformation(mptr,dv);

		if(np==THREED_MPM)
		{   SRConstitutiveLaw3D(mptr,dv,delTime,np,properties,res);
		}
		else
		{   SRConstitutiveLaw2D(mptr,dv,delTime,np,properties,res);
		}
	}
}

// return pointer to elastic properties (subclass might store them different in properties
ElasticProperties *Elastic::GetElasticPropertiesPointer(void *properties) const { return (ElasticProperties *)properties; }

#pragma mark Elastic::Methods (Large Rotation)

// Entry point for large rotation
void Elastic::LRConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// get rotations and initial rotation R0
	Matrix3 dR,Rnm1,Rtot;
	Matrix3 R0 = mptr->GetInitialRotation();
	
	// get incremental strain and rotation
	Matrix3 de = LRGetStrainIncrement(INITIALMATERIAL,mptr,du,&dR,&R0,&Rnm1,&Rtot);
	
	// Get rotation to material axes for state n-1 in n (or tot)
	Matrix3 Rtotnm1 = Rnm1*R0;
	Rtot *= R0;
	if(np==THREED_MPM) mptr->SetRtot(Rtot);
	
	// cast pointer to material-specific data
	ElasticProperties *p = GetElasticPropertiesPointer(properties);

    // residual strains (thermal and moisture) in material axes
	double exxr,eyyr,ezzr;
	if(np==THREED_MPM)
	{	exxr = p->alpha[0]*res->dT;
		eyyr = p->alpha[1]*res->dT;
		ezzr = p->alpha[2]*res->dT;
		if(fmobj->HasFluidTransport())
		{	exxr += p->beta[0]*res->dC;
			eyyr += p->beta[1]*res->dC;
			ezzr += p->beta[2]*res->dC;
		}
	}
	else
	{	exxr = p->alpha[1]*res->dT;
		eyyr = p->alpha[2]*res->dT;
		ezzr = p->alpha[4]*res->dT;
		if(fmobj->HasFluidTransport())
		{	exxr += p->beta[1]*res->dC;
			eyyr += p->beta[2]*res->dC;
			ezzr += p->beta[4]*res->dC;
		}
	}
	Matrix3 er = Matrix3(exxr,0.,0.,eyyr,ezzr);
	
	// finish up
	LRElasticConstitutiveLaw(mptr,de,er,Rtot,dR,Rtotnm1,np,properties,res);
}

// Once stress deformation has been decomposed, finish calculations in the material axis system
// When stress are found, they are rotated back to the global axes (using Rtot and dR)
// Similarly, strain increments are rotated back to find work energy (done in global system)
void Elastic::LRElasticConstitutiveLaw(MPMBase *mptr,Matrix3 &de,Matrix3 &er,Matrix3 &Rtot,Matrix3 &dR,
									   Matrix3 &Rnm1tot,int np,void *properties,ResidualStrains *res) const
{
	// effective strains
	double dvxxeff = de(0,0)-er(0,0);
	double dvyyeff = de(1,1)-er(1,1);
	double dvzzeff = de(2,2)-er(2,2);
	double dgamxy = 2.*de(0,1);
	
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
 
	// stress increments
	// cast pointer to material-specific data
	ElasticProperties *p = GetElasticPropertiesPointer(properties);
	Tensor dsig;
	if(np==THREED_MPM)
	{	double dgamyz = 2.*de(1,2);
		double dgamxz = 2.*de(0,2);
		
		// get dsigma
		dsig.xx = p->C[0][0] * dvxxeff + p->C[0][1] * dvyyeff + p->C[0][2] * dvzzeff;
		dsig.yy = p->C[1][0] * dvxxeff + p->C[1][1] * dvyyeff + p->C[1][2] * dvzzeff;
		dsig.zz = p->C[2][0] * dvxxeff + p->C[2][1] * dvyyeff + p->C[2][2] * dvzzeff;
		dsig.yz = p->C[3][3]*dgamyz;
		dsig.xz = p->C[4][4]*dgamxz;
		dsig.xy = p->C[5][5]*dgamxy;

		// update sigma = dR signm1 dRT + Rtot dsigma RtotT
		dsig = Rtot.RVoightRT(&dsig, true, false);
		*sp = dR.RVoightRT(sp, true, false);
		AddTensor(sp, &dsig);
		
		// stresses are in global coordinates so need to rotate strain and residual
		// strain to get work energy increment per unit mass (dU/(rho0 V0))
		Matrix3 derot = de.RMRT(Rtot);
		Matrix3 errot = er.RMRT(Rtot);
		mptr->AddWorkEnergyAndResidualEnergy(sp->xx*derot(0,0) + sp->yy*derot(1,1) + sp->zz*derot(2,2)
											  + 2.*(sp->yz*derot(1,2) + sp->xz*derot(0,2) + sp->xy*derot(0,1)),
											 sp->xx*errot(0,0) + sp->yy*errot(1,1) + sp->zz*errot(2,2)
											  + 2.*(sp->yz*errot(1,2) + sp->xz*errot(0,2) + sp->xy*errot(0,1)) );
	}
	else
	{	// find stress increment
		// this does xx, yy, and xy only. zz done later if needed
		if(np==AXISYMMETRIC_MPM)
		{	dsig.xx = p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff + p->C[4][1]*dvzzeff ;
			dsig.yy = p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff + p->C[4][2]*dvzzeff;
		}
		else
		{	dsig.xx = p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff;
			dsig.yy = p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff;
		}
		dsig.xy = p->C[3][3]*dgamxy;
		dsig.zz = 0.;
		
		// update sigma = dR signm1 dRT + Rtot dsigma RtotT
		dsig = Rtot.RVoightRT(&dsig, true, true);
		*sp = dR.RVoightRT(sp, true, true);
		AddTensor(sp, &dsig);
		
		// stresses are in global coordinate so need to rotate strain and residual
		// strain to get work energy increment per unit mass (dU/(rho0 V0))
		double ezzr = er(2,2);
		Matrix3 derot = de.RMRT(Rtot);
		Matrix3 errot = er.RMRT(Rtot);
		
		// work and residual strain energy increments and sigma or F in z direction
		double workEnergy = sp->xx*derot(0,0) + sp->yy*derot(1,1) + 2.0*sp->xy*derot(0,1);
		double resEnergy = sp->xx*errot(0,0) + sp->yy*errot(1,1) + 2.*sp->xy*errot(0,1);
		if(np==PLANE_STRAIN_MPM)
		{	// need to add back terms to get from reduced cte to actual cte
			sp->zz += p->C[4][1]*(dvxxeff+p->alpha[5]*ezzr) + p->C[4][2]*(dvyyeff+p->alpha[6]*ezzr) - p->C[4][4]*ezzr;
			
			// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			resEnergy += sp->zz*ezzr;
		}
		else if(np==PLANE_STRESS_MPM)
		{	// zz deformation
			de.set(2,2,p->C[4][1]*dvxxeff + p->C[4][2]*dvyyeff + ezzr);
			mptr->IncrementDeformationGradientZZ(de(2,2));
		}
		else
		{	// axisymmetric hoop stress
			sp->zz += p->C[4][1]*dvxxeff + p->C[4][2]*dvyyeff + p->C[4][4]*dvzzeff;
			
			// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			workEnergy += sp->zz*de(2,2);
			resEnergy += sp->zz*ezzr;
		}
		mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
	}
    
    // no heat energy (should get dTq0)
    //IncrementHeatEnergy(mptr,dTq0,0.);
	
#ifdef POROELASTICITY
	UndrainedPressIncrement(mptr,de(0,0),de(1,1),de(2,2));
#endif
}

#pragma mark Elastic::Methods (Small Rotation)

/* For 2D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
	Assumes linear elastic, uses hypoelastic correction
   For Axisymmetric MPM, x->R, y->Z, x->theta, and dvzz is change in hoop strain
	(i.e., du/r on particle and dvzz will be zero if not axisymmetric)
*/
void Elastic::SRConstitutiveLaw2D(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// Find strain increments
	double dvxx = du(0,0);
	double dvyy = du(1,1);
	double dgam = du(1,0)+du(0,1);
	double dwrotxy = du(1,0)-du(0,1);
	
	// cast pointer to material-specific data
	ElasticProperties *p = (ElasticProperties *)properties;
	
    // residual strains (thermal and moisture)
	double erxx = p->alpha[1]*res->dT;
	double eryy = p->alpha[2]*res->dT;
	double erxy = p->alpha[3]*res->dT;
	double erzz = CTE3*res->dT;
	if(fmobj->HasFluidTransport())
	{	erxx += p->beta[1]*res->dC;
		eryy += p->beta[2]*res->dC;
		erxy += p->beta[3]*res->dC;
		erzz += CME3*res->dC;
	}
	
    // moisture and thermal strain and temperature change
	//   (when diffusion, conduction, OR thermal ramp active)
    double dvxxeff = dvxx - erxx;
    double dvyyeff = dvyy - eryy;
	double dgameff = dgam - erxy;
    
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;
	
	// find stress (Units N/m^2  mm^3/g)
	// this does xx, yy, ans xy only. zz do later if needed
	double c1,c2,c3;
	if(np==AXISYMMETRIC_MPM)
	{	// hoop stress affect on RR, ZZ, and RZ stresses
		double dvzzeff = du(2,2) - erzz;
		c1=p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff + p->C[4][1]*dvzzeff + p->C[1][3]*dgameff;
		c2=p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff + p->C[4][2]*dvzzeff + p->C[2][3]*dgameff;
		c3=p->C[1][3]*dvxxeff + p->C[2][3]*dvyyeff + p->C[4][3]*dvzzeff + p->C[3][3]*dgameff;
	}
	else
    {	c1=p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff + p->C[1][3]*dgameff;
		c2=p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff + p->C[2][3]*dgameff;
		c3=p->C[1][3]*dvxxeff + p->C[2][3]*dvyyeff + p->C[3][3]*dgameff;
	}
	Hypo2DCalculations(mptr,dwrotxy,c1,c2,c3);
    
	// work and resdiaul strain energy increments
	double workEnergy = 0.5*((st0.xx+sp->xx)*dvxx + (st0.yy+sp->yy)*dvyy + (st0.xy+sp->xy)*dgam);
	double resEnergy = 0.5*((st0.xx+sp->xx)*erxx + (st0.yy+sp->yy)*eryy + (st0.xy+sp->xy)*erxy);
	if(np==PLANE_STRAIN_MPM)
	{	// need to add back terms to get from reduced cte to actual cte
		sp->zz += p->C[4][1]*(dvxxeff+p->alpha[5]*erzz)+p->C[4][2]*(dvyyeff+p->alpha[6]*erzz)
					+p->C[4][3]*(dgameff+p->alpha[7]*erzz)-p->C[4][4]*erzz;
		
		// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
		resEnergy += 0.5*(st0.zz+sp->zz)*erzz;
	}
	else if(np==PLANE_STRESS_MPM)
	{	// zz deformation
		mptr->IncrementDeformationGradientZZ(p->C[4][1]*dvxxeff+p->C[4][2]*dvyyeff+p->C[4][3]*dgameff+erzz);
	}
	else
	{	// axisymmetric hoop stress
		sp->zz += p->C[4][1]*dvxxeff + p->C[4][2]*dvyyeff + p->C[4][4]*(du(2,2) - erzz) + p->C[4][3]*dgameff;
		
		// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
		workEnergy += 0.5*(st0.zz+sp->zz)*du(2,2);
		resEnergy += 0.5*(st0.zz+sp->zz)*erzz;
	}
	mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
	
	// no more heat energy (should get dTq0)
	//IncrementHeatEnergy(mptr,dTq0,0.);
	
	// poroelasticity - none because never in SR mode
}

/* For 3D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    duij are (gradient rates X time increment) to give deformation gradient change
    Assumes linear elastic, uses hypoelastic correction
*/
void Elastic::SRConstitutiveLaw3D(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// strain increments
	double dvxx = du(0,0);
	double dvyy = du(1,1);
	double dvzz = du(2,2);
	
	// engineering shear strain icrements
    double dgamxy = du(0,1)+du(1,0);
    double dgamxz = du(0,2)+du(2,0);
    double dgamyz = du(1,2)+du(2,1);
	
	// rotational strain increments
	double dwrotxy = du(1,0)-du(0,1);
	double dwrotxz = du(2,0)-du(0,2);
	double dwrotyz = du(2,1)-du(1,2);

	// cast pointer to material-specific data
	ElasticProperties *p = (ElasticProperties *)properties;
	
    // residual strains (thermal and moisture) (isotropic only)
	double exxr = p->alpha[0]*res->dT;
	double eyyr = p->alpha[1]*res->dT;
	double ezzr = p->alpha[2]*res->dT;
	double eyzr = p->alpha[3]*res->dT;
	double exzr = p->alpha[4]*res->dT;
	double exyr = p->alpha[5]*res->dT;
	if(fmobj->HasFluidTransport())
	{	exxr += p->beta[0]*res->dC;
		eyyr += p->beta[1]*res->dC;
		ezzr += p->beta[2]*res->dC;
		eyzr += p->beta[3]*res->dC;
		exzr += p->beta[4]*res->dC;
		exyr += p->beta[5]*res->dC;
	}
	
	// effective strains
	double dvxxeff = dvxx-exxr;
	double dvyyeff = dvyy-eyyr;
	double dvzzeff = dvzz-ezzr;
	double dgamyzeff = dgamyz-eyzr;
	double dgamxzeff = dgamxz-exzr;
	double dgamxyeff = dgamxy-exyr;
	
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;
	
	// stress increments
	double delsp[6];
	int i;
	for(i=0;i<6;i++)
    {   delsp[i] = p->C[i][0]*dvxxeff + p->C[i][1]*dvyyeff + p->C[i][2]*dvzzeff
                    + p->C[i][3]*dgamyzeff + p->C[i][4]*dgamxzeff + p->C[i][5]*dgamxyeff;
    }
	
	// update stress (need to make hypoelastic)
	Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,delsp);
	
	// work energy increment per unit mass (dU/(rho0 V0))
    mptr->AddWorkEnergyAndResidualEnergy(
						0.5*((st0.xx+sp->xx)*dvxx + (st0.yy+sp->yy)*dvyy
                            + (st0.zz+sp->zz)*dvzz  + (st0.yz+sp->yz)*dgamyz
                            + (st0.xz+sp->xz)*dgamxz + (st0.xy+sp->xy)*dgamxy),
						0.5*((st0.xx+sp->xx)*exxr + (st0.yy+sp->yy)*eyyr
							 + (st0.zz+sp->zz)*ezzr  + (st0.yz+sp->yz)*eyzr
							 + (st0.xz+sp->xz)*exzr + (st0.xy+sp->xy)*exyr));
    
	// no more heat energy (should get dTq0)
	//IncrementHeatEnergy(mptr,dTq0,0.);

	// poroelasticity - none because never in SR mode
}

#pragma mark Elastic:Methods

// From thermodyanamics Cp-Cv = C [a].[a] T/rho (and rotation independent)
// If orthotropic = (C11*a11+C12*a22+C13*a33, C21*a11+C22*a22+C23*a33, C31*a11+C32*a22+C33*a33).(a11,a22,a33)T/rho
//                = (C11*a11^2+(C12+C21)*a22*a11+(C13_C31)*a33*a11+C22*a22^2+(C23+C32)*a33*a22+C33*a33^2)T/rho
// Cadota in nJ/(g-K^2) so output in nJ/(g-K)
double Elastic::GetCpMinusCv(MPMBase *mptr) const
{   return mptr!=NULL ? Cadota*mptr->pPreviousTemperature : Cadota*thermal.reference;
}

// Increment small strain material deformation gradient using F(n) = (I + grad u)F(n-1)
void Elastic::HypoIncrementDeformation(MPMBase *mptr,Matrix3 du) const
{
	// get incremental deformation gradient
	const Matrix3 dF = du.Exponential(1);
	
	// current deformation gradient
	Matrix3 pF = mptr->GetDeformationGradientMatrix();
	
	// new deformation matrix
	const Matrix3 F = dF*pF;
    mptr->SetDeformationGradientMatrix(F);
}

#pragma mark Elastic:Softening Methods

// For softening materials find area that the plane with given normal vector makes with
// box define by the particle volume. Return that area divided by particle volume
// Note that if particle have some rotation (either initial or during simulation), that
//   norm is vector in the unrotated particle coordinates
double AreaOverVolume3D(Vector *,double,double,double); // function prototype
double Elastic::GetAcOverVp(int np,MPMBase *mptr,Vector *norm) const
{
	double AcOverVp;

	// get particle size (relative to grid size)
	Vector lp;
	mptr->GetDimensionlessSize(lp);
	
	// get cell size (which may differ in x, y, and z)
	Vector grid = theElements[mptr->ElemID()]->GetDeltaBox();
	
	// undeformed particle size in x and y directions
	double dx = lp.x*grid.x;
	double dy = lp.y*grid.y;

	if(np==THREED_MPM)
	{	
		double dz = lp.z*grid.z;
		AcOverVp = AreaOverVolume3D(norm,dx,dy,dz);
	}
	else
	{	// in 2D, sin(thetap) is sin of angle of particle diagonal with x axis
		//      and n = (sin(thetac),-cos(thetac)) where theta c is angle of crack plane with x axis
		// Can use absolute values
		//      if thetac>thetap, intersects top and bottom Ac/Vp = (dy/sin(thetac)) / (dx*dy)
		//         otherwise, intersects left and right Ac/Vp = (dx/cos(thetac) / (dx*dy)
		double sinthetap = dy/sqrt(dx*dx+dy*dy);
		if(fabs(norm->x)>sinthetap)
			AcOverVp = 1./(dx*fabs(norm->x));
		else
			AcOverVp = 1./(dy*fabs(norm->y));
	}

	return AcOverVp;
}

// Soften in one direction (normal, one shear, or coupled shear)
//		den is strain increment in that direction
//		history - softening variables for this material with elements
//					deltaParam - maximum cracking strain
//					dParam - damage parameter
//		sigma - strength this property
//		scale - FM scaling = rG Ac/(Vp rS sigmai) including relative strength and toughness
//		softLaw - softening law this direction
//		Es - stiffness (e.g., C11,C55, or C66)
//		T0 - current traction in that direction (sign adjusted if was shear)
//		sigmaAlpha - for other variables on strength
// 		scaleAlpha - for other variables on strength
// Return false if no damage
//		call must find elastic cracking strain increment and watch for contact
// Return true if damage and also do the following
//		Set deCrack to the cracking strain increment
//		Set dispEnergy to increment in energy dissipation
//		Set criticalStrain to true if decohesion occured
bool Elastic::SoftenAxis(MPMBase *mptr,double den,double *history,int deltaParam,int dParam,
                         double sigma,double scale,SofteningLaw *softLaw,
						 double Es,double T0,double sigmaAlpha,double scaleAlpha,
						 double &deCrack,double &dispEnergy,bool &criticalStrain) const
{
	// no damage if negative here (shear accounts for this by caller) and no change in alpha variables
	if(den<0. && scaleAlpha<0.) return false;
	
	// get damage parameter, trial stress increment and stress
	double d = history[dParam];
	double dTTrial = Es*(1-d)*den;
	double TTrial = T0 + dTTrial;
	
	// Elastic or damaging?
	double delta = history[deltaParam];
	double fval = softLaw->GetFFxn(delta,scale);
	double FI = sigma*fval;
	
	// look for alpha variables
 	if(sigmaAlpha>0.)
    {   // switch to delta+ddeltaElastic, alpha+dAlpha on damage surface
        double ddeltaElastic = softLaw->GetDDeltaElastic(delta,sigmaAlpha,scaleAlpha,d,Es*(1-d));
		
		// update strength and delta
		sigma = sigmaAlpha;
		scale = scaleAlpha;
		delta += ddeltaElastic;
		
		// evaluate softening law at new delta and alpha+dalpha
		fval = softLaw->GetFFxn(delta,scale);
		FI = sigma*fval;
	}
	
	// If no damage, update delta (if needed) and return false
	if(TTrial <= FI)
	{	history[deltaParam] = delta;			// in case changed
		return false;
	}
	
	// Find portion  of the step that causes damage
	// Note that for alpha terms, FI is F(delta+ddeltaElastic,alpha+dAlpha)
	// The elastic strain part is de(1) = dee-de(2) (with or without alpha dependence)
	double den2 = (TTrial - FI)/(Es*(1.-d));
	
	// increment in delta strain (note: returns < 0 if decohesion)
	// note that when sigmaAlpha>0, sigma, scale, and delta are
	//			changed to sigmaAlpha, scaleAlpha, and delta+ddeltaElastic
	double en0 = sigma/Es;				 // initiation strain
	double ddelta = softLaw->GetDDelta(den2,en0,delta,d,scale);
	
	// check if failed
	double deltaPrevious = delta;			// = delta + ddeltaelastic if occurred
	if(ddelta<0.)
	{	// Decohesion: set delta=deltaMax, d=1, and set crackingStrain flag true
		// call will get cracking strain in post-failure update code
		delta = softLaw->GetDeltaMax(scale);
		criticalStrain = true;
		history[dParam] = 1.;
	}
	else
	{	// damage propagation
		delta += ddelta;
		fval = softLaw->GetFFxn(delta,scale);
		history[dParam] = delta/(delta+en0*fval);
		
		// total cracking strain increment is some of elastic part and damage part
		// deCrack = d*(de(1)+de(elastic)) + ddelta
		deCrack = d*(den-den2) + ddelta;
		
	}
	
	// dissipated energy increment per unit mass using dOmega = (1/2) phi (ddelta-ddeltaElastic)
	// note that deltaPrevious = deltai + ddeltaelastic, so deltaf-deltaPrevious = ddeltaf-ddelti-ddeltaElastic
	double dissipated = 0.5*sigma*softLaw->GetPhiFxn(deltaPrevious,scale)*(delta-deltaPrevious);
	
	// This is prior discrete calculation. It is identical to above for linear softening
	//double dissipated = sigma*(softLaw->GetGToDelta(delta,scale) - softLaw->GetGToDelta(deltaPrevious,scale));
	dispEnergy += dissipated;
	
	// new delta paramater
	history[deltaParam] = delta;
    
	return true;
}

// Calculate incremental cracking strains and incremental stress in post failure state
// str is current stress rotated into crack axis system, etr is current cracking strain in global coordinate
// Rtot rotates to crack axis system, C11 is normal direction stiffness
// den, dgamxy, and dgamxz are total increments
void Elastic::PostFailureUpdate(double &decxx,double &dgcxy,double &dgcxz,Tensor *dsig,Tensor *str,Tensor *etr,
								Matrix3 Rtot,double C11,double den,double Gxy,double dgamxy,double Gxz,double dgamxz,
								bool is2D,double coef) const
{
	bool inContact = false;
	
	// normal direction is opening
	if(den>0.)
	{	if(etr->xx > 0.)
		{	// opened crack stays open
			dsig->xx = -str->xx;			// set stress to zero
			decxx = den - dsig->xx/C11;		// add to cracking strain
		}
		else
		{	if(decxx>=0.)
			{	// closed crack opens
				dsig->xx = -str->xx;			// set stress to zero
				decxx = den - dsig->xx/C11;		// add to cracking strain
			}
			else
			{	// closed crack stays closed
				decxx = 0.;				// no cracking strain
				dsig->xx = C11*den;		// elastic stress change
				inContact = true;
			}
		}
	}
	
	// normal direction is closing
	else
	{	if(etr->xx > 0.)
		{	if(etr->xx+den > 0.)
			{	// opened crack stays open
				dsig->xx = -str->xx;			// set stress to zero
				decxx = den - dsig->xx/C11;		// add to cracking strain
			}
			else
			{	// opened crack closes
				decxx = -etr->xx;				// cracking strain to zero
				dsig->xx = C11*(den - decxx);	// elastic compression stress
				inContact = true;
			}
		}
		else
		{	// closed crack remains closed
			decxx = 0.;							// stay at zero cracking strain
			dsig->xx = C11*den;					// add to elastic compression
			inContact = true;
		}
	}

	if(inContact && coef>0.)
	{	// normal compression
		double N = -(str->xx + dsig->xx);
		
		// to test (phi=1 is stick and phi=-tauxy/(Gxy*dgamxy) is frictionless)
		double phi=1.,stick=coef*N,txyTest=str->xy+Gxy*dgamxy;
		
		if(is2D)
		{	// slip or stick
			if(txyTest > stick)
				phi = (stick - str->xy)/(Gxy*dgamxy);
			else if(txyTest<-stick)
				phi = (-stick - str->xy)/(Gxy*dgamxy);
			
			// stress update
			dsig->xy = phi*Gxy*dgamxy;
			dsig->xz = -str->xz;
			
			// cracking strain
			dgcxy = (1.-phi)*dgamxy;
			dgcxz = dgamxz;
					
			// TO DO: frictional heating = coef*N*|dgcxy|*vol(p)
		}
		else
		{	double txzTest=str->xz+Gxz*dgamxz,stick2=stick*stick;
			double mag2 = txyTest*txyTest + txzTest*txzTest;
			if(mag2>stick2)
			{	double c = str->xy*str->xy + str->xz*str->xz - stick2;
				double gdxy = Gxy*dgamxy,gdxz = Gxz*dgamxz;
				double a = gdxy*gdxy + gdxz*gdxz;
				double b = 2.*(str->xy*gdxy + str->xz*gdxz);
				double r1, r2;
				if(RealQuadraticRoots(a,b,c,r1,r2))
				{	if(c<=0.)
					{	// c<0 means started inside friction circle and went out, will get 1 positive root at crossing
						phi = fmax(r1,r2);
					}
					else if(r1>=0.)
					{	// c>0 and two positive roots, take smallest or closest to circle
						phi = fmin(r1,r2);
					}
					else
					{	// c>0 and two negative roots, take largest or closest to circle
						phi = fmax(r1,r2);
					}
					
					// stress update
					dsig->xy = phi*Gxy*dgamxy;
					dsig->xz = phi*Gxz*dgamxz;
					
					// cracking strain
					dgcxy = (1.-phi)*dgamxy;
					dgcxz = (1.-phi)*dgamxz;
				}
				else
				{	// two imaginary roots means c>0 (started outside the circle)
					// and increment line does not intersect the circle
					// convert to radial return
					phi = sqrt(stick2/mag2);
					
					// stress update
					dsig->xy = phi*txyTest - str->xy;
					dsig->xz = phi*txzTest - str->xz;
					
					// shear cracking strains
					dgcxy = dgamxy - dsig->xy/Gxy;
					dgcxz = dgamxz - dsig->xz/Gxz;
				}

				// TO DO: frictional heating = coef*N*sqrt(dgcxy^2+dgcxz^2)*vol(p)
			}
			else
			{	// stick or set Phi to 1
				dsig->xy = Gxy*dgamxy;
				dsig->xz = Gxz*dgamxz;
				
				// cracking strain
				dgcxy = 0.;
				dgcxz = 0.;
			}
		}
	}
	else
	{	// zero shear stress
		dsig->xy = -str->xy;
		dsig->xz = -str->xz;
	
		// shear cracking strains
		dgcxy = dgamxy - dsig->xy/Gxy;
		dgcxz = dgamxz - dsig->xz/Gxz;
	}
}

// Update cracking strain by rotation (decxx, dgcxy, dgcxz) to crack axis
// system and them adding to alt strain tensor
// Materials that store cracking strain elsewhere must override
void Elastic::UpdateCrackingStrain(int np,Tensor *ecrack,double decxx,double dgcxy,double dgcxz,Matrix3 Rtot,double *soft) const
{
	if(np==THREED_MPM)
	{	// rotate from crack axis to global axes
		Tensor dec = MakeTensor(decxx,0.,0.,0.,dgcxz,dgcxy);
		dec=Rtot.RVoightRT(&dec,false,false);
		AddTensor(ecrack, &dec);
	}
	else
	{	// rotate from crack axis to global axes
		Tensor dec = MakeTensor2D(decxx, 0., 0., dgcxy);
		dec = Rtot.RVoightRT(&dec,false,true);
		ecrack->xx += dec.xx;
		ecrack->yy += dec.yy;
		ecrack->xy += dec.xy;
	}
}

// Rotate stress increment in crack axis system to global axes
// and then add to sp (str has sp rotated already by dR)
void Elastic::UpdateCrackingStress(int np,Tensor *sp,Tensor *dsig,Tensor *str,Matrix3 Rtot) const
{
	if(np==THREED_MPM)
	{	// rotate stress increment from crack axis to global axes
		*dsig = Rtot.RVoightRT(dsig,true,false);
		
		// add up
		sp->xx = str->xx + dsig->xx;
		sp->yy = str->yy + dsig->yy;
		sp->zz = str->zz + dsig->zz;
		sp->xy = str->xy + dsig->xy;
		sp->xz = str->xz + dsig->xz;
		sp->yz = str->yz + dsig->yz;
	}
	else
	{	// rotate stress increment from crack axis to global axes
		*dsig = Rtot.RVoightRT(dsig,true,true);
		
		// add up
		sp->xx = str->xx + dsig->xx;
		sp->yy = str->yy + dsig->yy;
		sp->zz = str->zz + dsig->zz;
		sp->xy = str->xy + dsig->xy;
	}
}
