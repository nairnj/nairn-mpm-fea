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

// for numerical softening methods
#include "SofteningLaw.hpp"
#define NEWTON_MAX 20
static bool exceedNewtonMax = false;

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
    
    // no more heat energy (should get dTq0)
    //IncrementHeatEnergy(mptr,dTq0,0.);
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
//		dStr is strain increment in that direction
//		history - softening variables for this material with elements
//					maxStr - maximum cracking strain
//					dParam - damage parameter
//		sigma - strength this property
//		scale - FM scaling
//		softLaw - softening law this direction
//		stiff - stiffness (e.g., C11,C55, or C66)
//		T0 - current traction in that direction (adjust for coupled shear, and sign for shear)
//		str0 - initiation strain this direction (adjust for coupled shear, and sign for shear)
// Last four to set cracking strain(s), add to dissipated energy, and flag failure
//		For coupled shear, deAlt!=NULL
// return false if no damage (dStr<0 or Updated Traction < Current Strength) and caller
//		must handle elastic update with contact
bool Elastic::SoftenAxis(double dStr,double *history,int maxStr,int dParam,double sigma,double scale,SofteningLaw *softLaw,double stiff,double T0,double str0,double &deCrack,double *deAlt,double &dispEnergy,bool &criticalStrain) const
{
	// no damage if negative here (shear accounts for this by caller)
	if(dStr<0.) return false;
	
	// get parameter and elastic increment
	double delta = history[maxStr];
	double fval = softLaw->GetFFxn(delta,scale);
	double dam = history[dParam];			// could also calculate: dam = delta/(delta+str0*fval);

	// Assume damage occurting (dStr>0)
	
	// trial stress increment
	double dTTrial = stiff*(1-dam)*dStr;
	double currentStrength = sigma*fval;
	
	// Elastic or damaging?
	double TTrial = T0 + dTTrial;
	if(TTrial <= currentStrength) return false;
	
	// damage propagation strain
	double dam2 = (TTrial - currentStrength)/(stiff*(1.-dam));
	
	// increment in max crack opening strain
	double ddelta = softLaw->GetDDelta(dam2,str0,delta,scale);
	
	// check if failed
	if(ddelta<0.)
	{	// Make this step a post failure step by setting ddelta = dStr
		// cracking strain gets set to same in post failure update
		//ddelta = deltaMax - delta; - not needed but gives (delta += ddelta) = deltaMax
		delta = softLaw->GetDeltaMax(scale);
		criticalStrain = true;
		history[dParam] = 1.;
	}
	else
	{	// damage propagation
		delta += ddelta;
		fval = softLaw->GetFFxn(delta,scale);
		history[dParam] = delta/(delta+str0*fval);
		
		if(deAlt==NULL)
			deCrack = dam*(dStr-dam2) + ddelta;
		else
		{	// get each component from ds*dgxy^(1) + d(ds*gxy^(2)) - see notes for details
			// load these terms temporarily into cracking strains, caller must decode as shown
			deCrack = delta/(delta+str0*fval);
			*deAlt = (deCrack-dam)/(stiff*(1.-dam));
			//dgcxy = deCrack*dgamxy + str(0,1)*(*deAlt);
			//dgcxz = deCrack*dgamxz + str(0,2)*(*deAlt);
		}
	}
	
	// dissipated energy increment per unit mass
	// here delta is new value and history[maxStr] is previous one
	dispEnergy += sigma*(softLaw->GetGToDelta(delta,scale) - softLaw->GetGToDelta(history[maxStr],scale));
	
	// new max strain
	history[maxStr] = delta;

	return true;
}

// Soften two shear components coupled by elliptical failure surface
bool Elastic::SoftenTwoAxes(double *history, double Txy0, double Gxy, double dgxy,
	double sigxy, double scxy, SofteningLaw *softxy, int dxyParam, int deltaxyParam,
	double Txz0, double Gxz, double dgxz, double sigxz, double scxz, SofteningLaw *softxz,
	int dxzParam, int deltaxzParam,double &dgcxy,double &dgcxz,double &dispEnergy,bool &criticalStrain) const
{
	// Work in positive stress quadrant (caller must fix sign on output dgcxy and dgcxz)
	if(Txy0<0.)
	{	Txy0 = -Txy0;
		dgxy = -dgxy;
	}
	if(Txz0<0.)
	{	Txz0 = -Txz0;
		dgxz = -dgxz;
	}
	
	// no damage
	if(dgxy<=0. && dgxz<=0.) return false;
	
	// get parameters
	double dxy = history[dxyParam];
	double delxy = history[deltaxyParam];
	double dTrixy = Gxy*(1 - dxy)*dgxy;
	double dxz = history[dxzParam];
	double delxz = history[deltaxzParam];
	double dTrixz = Gxz*(1 - dxz)*dgxz;

	// see if elastic
	double fxy = softxy->GetFFxn(delxy, scxy);
	double fxz = softxz->GetFFxn(delxz, scxz);
	double s1 = Txy0 / (sigxy*fxy);
	double s2 = Txz0 / (sigxz*fxz);
	double sd1 = dTrixy / (sigxy*fxy);
	double sd2 = dTrixz / (sigxz*fxz);
	if ((s1+sd1)*(s1+sd1) + (s2+sd2)*(s2+sd2) < 1.) return false;
	
	// common terms
	double gxy0 = sigxy/Gxy, gxz0 = sigxz/Gxz;
	
	// final answers
	double x1=0.,x2=0.;
	double dgxy1,dgxz1,dgxy2,dgxz2;
	
	// unaxial damage (and can only be one)
	if(dgxy<=0. || dgxz<=0.)
	{	// find which direction
		if(dgxz<=0.)
		{	// Damage only to xy direction
			double T2 = sqrt(1-(s2+sd2)*(s2+sd2));
			
			// partition elastic and damage
			double phi = (sigxy*fxy*T2-Txy0)/dTrixy;
			if(phi>1.) cout << "# phi>1 when dgxz<0" << endl;
			
			//  "1" is elastic strain increment and "2" is damage strain increment
			dgxy1 = phi*dgxy;
			dgxy2 = dgxy - dgxy1;
			dgxz1 = dgxz;		// all elastic
			
			// increment in max crack opening strain and check if failed
			x1 = softxy->GetDDelta(dgxy2,T2*gxy0,delxy,scxy);
			if(x1<0.)
			{	// Make large enough to fail, failure picked up below
				x1 = softxy->GetDeltaMax(scxy)+.1;
				fxy = 0.;
			}
			else
				fxy = softxy->GetFFxn(delxy+x1, scxy);
		}
		else
		{	// Damage only to xz direction
			double T2 = sqrt(1.-(s1+sd1)*(s1+sd1));
			
			// partition elastic and damage
			double phi = (sigxz*fxz*T2-Txz0)/dTrixz;
			if(phi>1.) cout << "# phi>1 when dgxy<0" << endl;
			
			//  "1" is elastic strain increment and "2" is damage strain increment
			dgxz1 = phi*dgxz;
			dgxz2 = dgxz - dgxz1;
			dgxy1 = dgxy;		// all elastic
			
			// increment in max crack opening strain
			x2 = softxy->GetDDelta(dgxz2,T2*gxz0,delxz,scxz);
			if(x2<0.)
			{	// Make large enough to fail, failure picked up below
				x2 = softxz->GetDeltaMax(scxz)+.1;
				fxz = 0.;
			}
			else
				fxz = softxz->GetFFxn(delxz+x2, scxy);
		}
	}
	
	// Biaxial Damage (using Newton's Method)
	else
	{	// solve for elastic fraction phi from a phi^2 + b phi + c = 0
		// c=0 if currently on failure surfaces (or phi=0), othersize c<0
		// Note that c>0 is likely small and is roundoff error from previous
		//   step moving to failure surface. Treat it as phi=0
		double phi = 0.;
		double c = s1*s1 + s2*s2 - 1.;
		if(c<0.)
		{	double a = sd1*sd1 + sd2*sd2;
			double b = 2.*(s1*sd1 + s2*sd2);
			
			// for b>0, two roots -2c/barg and -barg/2a
			// here c<0 for first is positive, a>0 so second is negative; we take the first
			double barg = b+sqrt(b*b-4.*a*c);
			phi = -2.*c/barg;
			if(phi>1.) cout << "# Found invalid phi = " << phi << endl;
		}
		
		//  "1" is elastic strain increment and "2" is damage strain increment
		dgxy1 = phi*dgxy;
		dgxz1 = phi*dgxz;
		dgxy2 = dgxy - dgxy1;
		dgxz2 = dgxz - dgxz1;
		
		// constant terms
		double T1 = Txy0/Gxy + (1 - dxy)*dgxy1 + dgxy2;
		double T2 = Txz0/Gxz + (1 - dxz)*dgxz1 + dgxz2;
		double gxy02 = gxy0*gxy0, gxz02 = gxz0*gxz0 ;
		
		// scaled stresses
		double nsigxy2 = (sigxy/Gxy)*(sigxy/Gxz);
		double nsigxz2 = (sigxz/Gxy)*(sigxz/Gxz);

//#define TRACK_NEWTON
		// cannot exceed maximum cracking opening strain (or decohesion if it does)
		double x1max = softxy->GetDeltaMax(scxy)-history[deltaxyParam];
		double x2max = softxz->GetDeltaMax(scxz)-history[deltaxzParam];
		double x1min = 0.;
		double x2min = 0.;
		double x1Converge = 1.e-7*(x1max+history[deltaxyParam]-x1min);
		double x2Converge = 1.e-7*(x2max+history[deltaxzParam]-x2min);
		int hitMax = 0;
		
		// initial guess (one step faster to start here instead of x1=x2=0).
		// linear steps assuming other axis does not change
		double fpxy = softxy->GetFpFxn(delxy, scxy);
		double fpxz = softxz->GetFpFxn(delxz, scxz);
		//double trialTau = sqrt(1-(s2+sd2)*(s2+sd2));
		//x1 = fmin(x1max,fmax(0.,dgxy2/(1.+trialTau*gxy0*fpxy)));
		//trialTau = sqrt(1-(s1+sd1)*(s1+sd1));
		//x2 = fmin(x2max,fmax(0.,dgxz2/(1.+trialTau*gxz0*fpxz)));
		x1 = fmin(0.5*x1max,dgxy2);
		x2 = fmin(0.5*x2max,dgxz2);
#ifdef TRACK_NEWTON
		double minCS = 0.9*softxy->GetDeltaMax(scxy);
		if(history[deltaxyParam]>minCS || history[deltaxzParam]>minCS)
			cout << "Start biaxial Newton ("<<dgxy2<<","<<x1<<","<<x1max << ")("<<dgxz2<<","<<x2<<","<<x2max<<")" << endl;
#endif
		fxy = softxy->GetFFxn(delxy+x1, scxy);
		fxz = softxz->GetFFxn(delxz+x2, scxz);
		if(!softxy->IsLinear()) fpxy = softxy->GetFpFxn(delxy+x1, scxy);
		if(!softxz->IsLinear()) fpxz = softxz->GetFpFxn(delxz+x2, scxz);
		
		// Newton's Method loop
		int nstep=1;
		while(true)
		{	// current fxy, fxz, fpxy, and fpxz at delxy+x1 and delxz+x2 are evaluated
			// get function values
			double fxy2 = fxy*fxy;
			double fxz2 = fxz*fxz;
			double dT1 = (T1 - x1)/gxy0;
			double dT2 = (T2 - x2)/gxz0;
			double f1 = fxz2*dT1*dT1 + fxy2*dT2*dT2 - fxy2*fxz2;
			double f2 = fxz2*(T1 - x1)*(dxz*dgxz2 - x2)/nsigxy2
							- fxy2*(T2 - x2)*(dxy*dgxy2 - x1)/nsigxz2;
#ifdef TRACK_NEWTON
			if(history[deltaxyParam]>minCS || history[deltaxzParam]>minCS)
				cout << "..." << nstep << ". (x1,x2)=("<<x1<<","<<x2<<") (f1,f2)=("<<f1<<","<<f2<<")" << endl;
#endif

			// get Jacobian
			double J11 = (-2.*fxz2/gxy02)*(T1 - x1) + 2.*fxy*fpxy*(dT2*dT2-fxz2);
			double J12 = (-2.*fxy2/gxz02)*(T2 - x2) + 2.*fxz*fpxz*(dT1*dT1-fxy2);
			double J21 = -fxz2*(dxz*dgxz2-x2)/nsigxy2
						+ fxy*(T2-x2)*(fxy - 2.*(dxy*dgxy2-x1)*fpxy)/nsigxz2;
			double J22 = fxy2*(dxy*dgxy2-x1)/nsigxz2
						- fxz*(T1-x1)*(fxz - 2.*(dxz*dgxz2-x2)*fpxz)/nsigxy2;
			
			// invert
			double detJ = J11*J22 - J12*J21;
			
			// increment x(k+1) = x(k) - J^{-1} (f1(k),f2(k))
			
			double detJ_inv = detJ / (0.00000001 + detJ*detJ); // regularized inverse
			double dx1 = (J22*f1 - J12*f2) *detJ_inv;
			double dx2 = (-J21*f1 + J11*f2) * detJ_inv;
			x1 -= dx1;
			x2 -= dx2;
			// limit to range (not sure best approach here)
			if(x1<x1min)
				x1 = x1min;
			else if(x1>x1max)
			{	x1 = x1max;
				hitMax++;
			}
			if(x2<x2min)
				x2 = x2min;
			if(x2>x2max)
			{	x2 = x2max;
				hitMax++;
			}
			
			// new softening values
			fxy = softxy->GetFFxn(delxy+x1, scxy);
			fxz = softxz->GetFFxn(delxz+x2, scxz);
			
			// is it converged
			//if((fabs(f1) < 1.e-14 && fabs(f2) < 1.e-14) || (nstep>NEWTON_MAX))
			if((fabs(dx1) < x1Converge && fabs(dx2) < x2Converge) || hitMax>1 || (nstep>NEWTON_MAX))
			{	if(nstep>NEWTON_MAX && !exceedNewtonMax)
				{
#pragma omp critical (output)
					{	cout << "# Warning: Newton's Method Exceeded " << NEWTON_MAX << " steps for biaxial shear softening material" << endl;
						cout << "#    Final (f1,f2) = (" << f1 << "," << f2 << ") (fxy,fxz) = (" << fxy << "," << fxz << ")" << endl;
						exceedNewtonMax = true;
					}
				}
				break;
			}
			
			// new derivatives
			if(!softxy->IsLinear()) fpxy = softxy->GetFpFxn(delxy+x1, scxy);
			if(!softxz->IsLinear()) fpxz = softxz->GetFpFxn(delxz+x2, scxz);
			nstep++;
		}
#ifdef TRACK_NEWTON
		if(history[deltaxyParam]>minCS || history[deltaxzParam]>minCS)
			cout << "...End Newton ("<<x1<<","<<x2<<")" << endl;
#endif
	}

	// check if failed in either direction
	if(softxy->HasFailed(delxy,x1,scxy))
	{	softxz->HasFailed(delxz,x2,scxz);
		history[deltaxyParam] += x1;
		history[deltaxzParam] += x2;
		criticalStrain = true;
	}
	else if(softxz->HasFailed(delxz,x2,scxz))
	{	history[deltaxyParam] += x1;
		history[deltaxzParam] += x2;
		criticalStrain = true;
	}
	else
	{	// damage propagation
		history[deltaxyParam] += x1;
		history[deltaxzParam] += x2;
		history[dxyParam] = history[deltaxyParam]/(history[deltaxyParam]+gxy0*fxy);
		history[dxzParam] = history[deltaxzParam]/(history[deltaxzParam]+gxz0*fxz);
		
		// absolute value of cracking strain
		dgcxy = dxy*dgxy1 + x1;
		dgcxz = dxz*dgxz1 + x2;
	}
	
	// dissipated energy increment per unit mass
	// Here history[] is new value and delij is previous value
	dispEnergy += sigxy*(softxy->GetGToDelta(history[deltaxyParam],scxy)
						 - softxy->GetGToDelta(delxy,scxy));
	dispEnergy += sigxz*(softxz->GetGToDelta(history[deltaxzParam],scxz)
						 - softxz->GetGToDelta(delxz,scxz));
	
	// compare to increment
	return true;
}

// Calculate incremental cracking strains and incremental stress in post failure state
// str is current stress rotated into crack axis system, etr is current cracking strain in global coordinate
// Rtot rotates to crack axis system, C11 is normal direction stiffness
// den, dgamxy, and dgamxz are total increments
void Elastic::PostFailureUpdate(double &decxx,double &dgcxy,double &dgcxz,Tensor *dsig,Tensor *str,Tensor *etr,
								Matrix3 Rtot,double C11,double den,double dgamxy,double dgamxz,bool is2D) const
{
	// normal direction is opening
	if(den>0.)
	{	if(etr->xx > 0.)
		{	// opened crack stays open
			decxx = den;				// add to cracking strain
			dsig->xx = -str->xx;		// set stress to zero
		}
		else
		{	decxx = den+str->xx/C11;	// cracking strain subtracting strain to get zero stress
			if(decxx>=0.)
			{	// closed crack opens
				dsig->xx = -str->xx;	// set stress to zero
			}
			else
			{	// closed crack stays closed
				decxx = 0.;				// no cracking strain
				dsig->xx = C11*den;		// elastic stress change
			}
		}
	}
	
	// normal direction is closing
	else
	{	if(etr->xx > 0.)
		{	if(etr->xx+den>0.)
			{	// opened crack stays open
				decxx = den;			// add to cracking strain
				dsig->xx = -str->xx;	// set stress to zero
			}
			else
			{	// opened crack closes
				decxx = -etr->xx;				// cracking strain to zero
				dsig->xx = C11*(den - decxx);	// elastic compression stress
			}
		}
		else
		{	// closed crack remains closed
			decxx = 0.;							// stay at zero cracking strain
			dsig->xx = C11*den;					// add to elastic compression
		}
	}
	
	// zero shear stress
	dsig->xy = -str->xy;
	dsig->xz = -str->xz;
	
	// shear cracking strains
	dgcxy = dgamxy;
	dgcxz = dgamxz;
}

// Update cracking strain and total stress
// Here str is current particle stress rotated by dR to updated configuration
// Return cracking strain in updated configuration in ecrack
void Elastic::UpdateCrackingStrainStress(int np,Tensor *ecrack,Tensor *sp,double decxx,double dgcxy,double dgcxz,
										 Tensor *dsig,Tensor *str,Matrix3 Rtot) const
{
	if(np==THREED_MPM)
	{	// rotate from crack axis to global axes
		Tensor dec = MakeTensor(decxx,0.,0.,0.,dgcxz,dgcxy);
		dec=Rtot.RVoightRT(&dec,false,false);
		AddTensor(ecrack, &dec);
		
		// rotate stress increment from crack axis to global axes
		*dsig = Rtot.RVoightRT(dsig,true,false);
		
		// add up
		sp->xx = str->xx + dsig->xx ;
		sp->yy = str->yy + dsig->yy;
		sp->zz = str->zz + dsig->zz;
		sp->xy = str->xy + dsig->xy;
		sp->xz = str->xz + dsig->xz;
		sp->yz = str->yz + dsig->yz;
	}
	else
	{	// rotate from crack axis to global axes
		Tensor dec = MakeTensor2D(decxx, 0., 0., dgcxy);
		dec = Rtot.RVoightRT(&dec,false,true);
		ecrack->xx += dec.xx;
		ecrack->yy += dec.yy;
		ecrack->xy += dec.xy;
		
		// rotate stress increment from crack axis to global axes
		*dsig = Rtot.RVoightRT(dsig,true,true);
		
		// add up
		sp->xx = str->xx + dsig->xx;
		sp->yy = str->yy + dsig->yy;
		sp->zz = str->zz + dsig->zz;
		sp->xy = str->xy + dsig->xy;
	}
}
