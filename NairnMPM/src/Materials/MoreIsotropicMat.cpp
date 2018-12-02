/********************************************************************************
	MoreIsotropicMat.cpp
	nairn-mpm-fea
 
	Created by John Nairn on 3/30/2015.
	opyright (c) 2015 John A. Nairn, All rights reserved.
 
	MPM-specific code for istropic, low-strain materials
********************************************************************************/

#include "stdafx.h"
#include "Materials/IsotropicMat.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"

#pragma mark IsotropicMat::Methods

// Isotropic material can use read-only initial properties
void *IsotropicMat::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{	return (void *)&pr;
}

// convert J to K using isotropic method
Vector IsotropicMat::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{	return IsotropicJToK(d,C,J0,np,nu,G);
}

/* Take increments in strain and calculate new Particle: strains, rotation strain,
	stresses, strain energy,
	du are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
 */
void IsotropicMat::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{	if(useLargeRotation)
		LRConstitutiveLaw(mptr,du,delTime,np,properties,res);
	else
	{	// increment deformation gradient
		HypoIncrementDeformation(mptr,du);

		if(np==THREED_MPM)
			SRConstitutiveLaw3D(mptr,du,delTime,np,properties,res);
		else
			SRConstitutiveLaw2D(mptr,du,delTime,np,properties,res);
	}
}
		
#pragma mark IsotropicMat::Methods (Large Rotation)

/* Take increments in strain and calculate new Particle: strains, rotation strain,
 stresses, strain energy,
 dvij are (gradient rates X time increment) to give deformation gradient change
 For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
 */
void IsotropicMat::LRConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
    // get incremental strain and rotation
    Matrix3 dR;
    Matrix3 de = LRGetStrainIncrement(CURRENT_CONFIGURATION,mptr,du,&dR,NULL,NULL,NULL);
	
	// cast pointer to material-specific data
	ElasticProperties *p = (ElasticProperties *)properties;
	
	// save initial stresses
	Tensor *sp = mptr->GetStressTensor();
	
	// residual strains (thermal and moisture)
	double eres = CTE1*res->dT;
	double ezzres = CTE3*res->dT;
	if(DiffusionTask::HasFluidTransport())
	{	eres += CME1*res->dC;
		ezzres += CME3*res->dC;
	}
	if(np==PLANE_STRESS_MPM || np==PLANE_STRAIN_MPM)
	{	// for generalized plane stress or strain, for isotropic = -nu*dszz/E or -nu*dezz
		eres += p->C[5][1]*res->doopse;
	}
	
	// effective strains
	double dvxxeff = de(0,0)-eres;
	double dvyyeff = de(1,1)-eres;
	double dgamxy = de(0,1)+de(1,0);
	double dVoverV = de(0,0)+de(1,1);
	
	// stress increments
	Tensor delsp;
	if(np==THREED_MPM)
	{	// more effective strains
		double dvzzeff = de(2,2)-eres;
		double dgamyz = de(1,2)+de(2,1);
		double dgamxz = de(0,2)+de(2,0);
		dVoverV += de(2,2);
		
		delsp.xx = p->C[0][0]*dvxxeff + p->C[0][1]*dvyyeff + p->C[0][2]*dvzzeff;
		delsp.yy = p->C[1][0]*dvxxeff + p->C[1][1]*dvyyeff + p->C[1][2]*dvzzeff;
		delsp.zz = p->C[2][0]*dvxxeff + p->C[2][1]*dvyyeff + p->C[2][2]*dvzzeff;
		delsp.yz = p->C[3][3]*dgamyz;
		delsp.xz = p->C[4][4]*dgamxz;
		delsp.xy = p->C[5][5]*dgamxy;
		
		// increment stress: rotate previous stress and add new one
		*sp = dR.RVoightRT(sp,true,false);
		AddTensor(sp,&delsp);
		
		// work energy increment per unit mass (dU/(rho0 V0))
		mptr->AddWorkEnergyAndResidualEnergy(sp->xx*de(0,0) + sp->yy*de(1,1) + sp->zz*de(2,2)
						+ sp->yz*(de(1,2)+de(2,1)) + sp->xz*(de(0,2)+de(2,0)) + sp->xy*(de(0,1)+de(1,0)),
							(sp->xx + sp->yy + sp->zz)*eres);
	}
	
	else
	{	if(np==AXISYMMETRIC_MPM)
		{	// hoop stress affect on RR, ZZ, and RZ stresses
			double dvzzeff = de(2,2) - eres;
			delsp.xx = p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff + p->C[4][1]*dvzzeff;
			delsp.yy = p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff + p->C[4][2]*dvzzeff;
			delsp.xy = p->C[3][3]*dgamxy;
		}
		else
		{	delsp.xx = p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff;
			delsp.yy = p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff;
			delsp.xy = p->C[3][3]*dgamxy;
		}
		
		// increment stress: rotate previous stress and add in-plane components
		*sp = dR.RVoightRT(sp,true,true);
		sp->xx += delsp.xx;
		sp->yy += delsp.yy;
		sp->xy += delsp.xy;
		
		// work and resdidual strain energy increments
		double workEnergy = sp->xx*de(0,0) + sp->yy*de(1,1) + sp->xy*(de(0,1)+de(1,0));
		double resEnergy = (sp->xx + sp->yy)*ezzres;
		if(np==PLANE_STRAIN_MPM)
		{	// for generalized plane strain, get zz stress (need to use actual cte)
			sp->zz += p->C[4][1]*(de(0,0)-ezzres) + p->C[4][2]*(de(1,1)-ezzres) + p->C[4][4]*(res->doopse-ezzres);
			
			// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			workEnergy += sp->zz*res->doopse;
			resEnergy += sp->zz*ezzres;
			dVoverV += res->doopse;
		}
		else if(np==PLANE_STRESS_MPM)
		{	// for generalized plane stress, get zz deformation (need to use actual cte)
			double dezz = res->doopse/p->C[4][4] + p->C[4][1]*(de(0,0)-ezzres)
								+ p->C[4][2]*(de(1,1)-ezzres) + ezzres;
			mptr->IncrementDeformationGradientZZ(dezz);
			
			// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			workEnergy += sp->zz*dezz;
			resEnergy += sp->zz*ezzres;
			dVoverV += dezz;
		}
		else
		{	// axisymmetric hoop stress
			sp->zz += p->C[4][1]*dvxxeff + p->C[4][2]*dvyyeff + p->C[4][4]*(de(2,2) - eres);
			
			// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			workEnergy += sp->zz*de(2,2);
			resEnergy += sp->zz*eres;
			dVoverV += de(2,2);
		}
		mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
	}
	
	// Isoentropic temperature rise = -(K 3 alpha T)/(rho Cv) (dV/V) = - gamma0 T (dV/V)
	double dTq0 = -gamma0*mptr->pPreviousTemperature*dVoverV;
    
    // track heat energy
    IncrementHeatEnergy(mptr,dTq0,0.);
}

#pragma mark IsotropicMat::Methods (Small Rotation)

/* For 2D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
	Assumes linear elastic and istropic, uses hypoelastic correction
	For Axisymmetric MPM, x->R, y->Z, x->theta, and dvzz is change in hoop strain
		(i.e., du/r on particle and dvzz will be zero if not axisymmetric)
*/
void IsotropicMat::SRConstitutiveLaw2D(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// Strain increments to local variables
    double dvxx = du(0,0);
	double dvyy = du(1,1);
	double dvzz = du(2,2);
    double dgam = du(0,1)+du(1,0);
	double dwrotxy = du(1,0)-du(0,1);
	
	// cast pointer to material-specific data
	ElasticProperties *p = (ElasticProperties *)properties;
	
    // residual strains (thermal and moisture)
	double eres = CTE1*res->dT;
	double ezzres = CTE3*res->dT;
	if(DiffusionTask::HasFluidTransport())
	{	eres += CME1*res->dC;
		ezzres += CME3*res->dC;
	}
	if(np==PLANE_STRESS_MPM || np==PLANE_STRAIN_MPM)
	{	// for generalized plane stress or strain, for isotropic = -nu*dszz/E or -nu*dezz
		eres += p->C[5][1]*res->doopse;
	}
	
    // moisture and thermal strain and temperature change
	//   (when diffusion, conduction, OR thermal ramp active)
    double dvxxeff = dvxx - eres;
    double dvyyeff = dvyy - eres;
	double dVoverV = dvxx + dvyy;
    
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;
	
	// find stress (Units N/m^2  mm^3/g)
	// this does xx, yy, ans xy only. zz do later if needed
	double c1,c2,c3;
	if(np==AXISYMMETRIC_MPM)
	{	// hoop stress affect on RR, ZZ, and RZ stresses
		double dvzzeff = dvzz - eres;
		c1 = p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff + p->C[4][1]*dvzzeff;
		c2 = p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff + p->C[4][2]*dvzzeff;
		c3 = p->C[3][3]*dgam;
	}
	else
    {	c1 = p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff;
		c2 = p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff;
		c3 = p->C[3][3]*dgam;
	}
	Hypo2DCalculations(mptr,dwrotxy,c1,c2,c3);
    
	// work and residual strain energy increments
	double workEnergy = 0.5*((st0.xx+sp->xx)*dvxx + (st0.yy+sp->yy)*dvyy + (st0.xy+sp->xy)*dgam);
	double resEnergy = 0.5*(st0.xx+sp->xx+st0.yy+sp->yy)*ezzres;
	if(np==PLANE_STRAIN_MPM)
	{	// for generalized plane strain, get zz stress (need to use actual cte)
		sp->zz += p->C[4][1]*(dvxx-ezzres) + p->C[4][2]*(dvyy-ezzres) + p->C[4][4]*(res->doopse-ezzres);
		
		// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
		workEnergy += 0.5*(st0.zz+sp->zz)*res->doopse;
		resEnergy += 0.5*(st0.zz+sp->zz)*ezzres;
		dVoverV += res->doopse;
	}
	else if(np==PLANE_STRESS_MPM)
	{	// for generalized plane stress, get zz deformation (need to use actual cte)
		double dezz = res->doopse/p->C[4][4] + p->C[4][1]*(dvxx-ezzres)
							+ p->C[4][2]*(dvyy-ezzres) + ezzres;
		mptr->IncrementDeformationGradientZZ(dezz);
		
		// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
		workEnergy += 0.5*(st0.zz+sp->zz)*dezz;
		resEnergy += 0.5*(st0.zz+sp->zz)*ezzres;
		dVoverV += dezz;
	}
	else
	{	// axisymmetric hoop stress
		sp->zz += p->C[4][1]*dvxxeff + p->C[4][2]*dvyyeff + p->C[4][4]*(dvzz - eres);
		
		// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
		workEnergy += 0.5*(st0.zz+sp->zz)*dvzz;
		resEnergy += 0.5*(st0.zz+sp->zz)*eres;
		dVoverV += dvzz;
	}
	mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
	
	// Isoentropic temperature rise = -(K 3 alpha T)/(rho Cv) (dV/V) = - gamma0 T (dV/V)
	double dTq0 = -gamma0*mptr->pPreviousTemperature*dVoverV;
    
    // track heat energy
    IncrementHeatEnergy(mptr,dTq0,0.);
}

/* For 3D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
	Assumes linear elastic and isotropic, uses hypoelastic correction
*/
void IsotropicMat::SRConstitutiveLaw3D(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
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
	double eres = CTE3*res->dT;
	if(DiffusionTask::HasFluidTransport())
		eres += CME3*res->dC;
	
	// effective strains
	double dvxxeff = dvxx-eres;
	double dvyyeff = dvyy-eres;
	double dvzzeff = dvzz-eres;
	double dVoverV = dvxx + dvyy + dvzz;
	
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;
	
	// stress increments
	double delsp[6];
	delsp[0] = p->C[0][0]*dvxxeff + p->C[0][1]*dvyyeff + p->C[0][2]*dvzzeff;
	delsp[1] = p->C[1][0]*dvxxeff + p->C[1][1]*dvyyeff + p->C[1][2]*dvzzeff;
	delsp[2] = p->C[2][0]*dvxxeff + p->C[2][1]*dvyyeff + p->C[2][2]*dvzzeff;
	delsp[3] = p->C[3][3]*dgamyz;
	delsp[4] = p->C[4][4]*dgamxz;
	delsp[5] = p->C[5][5]*dgamxy;
	
	// update stress (need to make hypoelastic)
	Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,delsp);
	
	// work energy increment per unit mass (dU/(rho0 V0))
    mptr->AddWorkEnergyAndResidualEnergy(0.5*((st0.xx+sp->xx)*dvxx + (st0.yy+sp->yy)*dvyy
											  + (st0.zz+sp->zz)*dvzz  + (st0.yz+sp->yz)*dgamyz
											  + (st0.xz+sp->xz)*dgamxz + (st0.xy+sp->xy)*dgamxy),
										 0.5*(st0.xx+sp->xx + st0.yy+sp->yy + st0.zz+sp->zz)*eres);
    
	// Isoentropic temperature rise = -(K 3 alpha T)/(rho Cv) (dV/V) = - gamma0 T (dV/V)
	double dTq0 = -gamma0*mptr->pPreviousTemperature*dVoverV;
    
    // track heat energy
    IncrementHeatEnergy(mptr,dTq0,0.);
}

#pragma mark IsotropicMat::Accessors

// Get magnitude of the deviatoric stress tensor when input is a deviatoric stress
// ||s|| = sqrt(s.s) = sqrt(2J2) where J2 = (1/2)s.s
// In 2D 2J2 = sx^2 + sy^2 + sz^2 + 2*txy^2
// In 3D 2J2 = sx^2 + sy^2 + sz^2 + 2*txy^2 + 2*txz^2 + 2*tyz^2
double IsotropicMat::GetMagnitudeSFromDev(Tensor *st,int np) const
{
	double s,t;
	
	switch(np)
	{   case THREED_MPM:
			s = st->xx*st->xx + st->yy*st->yy + st->zz*st->zz;
			t = st->xy*st->xy + st->xz*st->xz + st->yz*st->yz;
			break;
			
		default:
			s = st->xx*st->xx + st->yy*st->yy + st->zz*st->zz;
			t = st->xy*st->xy;
			break;
	}
	return sqrt(s+t+t);
}

// Get magnitude of the deviatoric stress tensor when input is a total stress
double IsotropicMat::GetMagnitudeSFromTotal(Tensor *st,int np) const
{
	double totalP = -(st->xx + st->yy + st->zz)/3.;
	Tensor sdev = np==THREED_MPM ?
		MakeTensor(st->xx+totalP,st->yy+totalP,st->zz+totalP,st->yz,st->xz,st->xy) :
		MakeTensor2D(st->xx+totalP,st->yy+totalP,st->zz+totalP,st->xy);
	return GetMagnitudeSFromDev(&sdev,np);
}

// Calculate wave speed in L/sec
// Because G in F/L^2 = mass/(L sec^2) and rho in mass/L^3, G/rho in L^2/sec^2
// Uses sqrt((K +4G/3)/rho) which is dilational wave speed
// Identity also: K + 4G/3 = Lambda + 2G = 2G(1-nu)/(1-2 nu)
double IsotropicMat::WaveSpeed(bool threeD,MPMBase *mptr) const
{	return sqrt(2.*G*(1.-nu)/(rho*(1.-2.*nu)));
}

// Calculate shear wave speed in L/sec
double IsotropicMat::ShearWaveSpeed(bool threeD,MPMBase *mptr,int offset) const { return sqrt(G/rho); }



