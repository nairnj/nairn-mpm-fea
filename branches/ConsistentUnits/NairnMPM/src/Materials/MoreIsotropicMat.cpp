/********************************************************************************
	MoreIsotropicMat.cpp
	nairn-mpm-fea
 
	Created by John Nairn on 3/30/2015.
	opyright (c) 2015 John A. Nairn, All rights reserved.
 
	MPM-specific code for istropic, low-strain materials
********************************************************************************/

#include "Materials/IsotropicMat.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"

#pragma mark IsotropicMat::Methods

// Isotropic material can use read-only initial properties
void *IsotropicMat::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
{	return (void *)&pr;
}

// convert J to K using isotropic method
Vector IsotropicMat::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{	return IsotropicJToK(d,C,J0,np,nu,G);
}

/* For 2D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
	Assumes linear elastic and istropic, uses hypoelastic correction
	For Axisymmetric MPM, x->R, y->Z, x->theta, and dvzz is change in hoop strain
		(i.e., du/r on particle and dvzz will be zero if not axisymmetric)
*/
void IsotropicMat::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,double dvzz,
                          double delTime,int np,void *properties,ResidualStrains *res) const
{
	// cast pointer to material-specific data
	ElasticProperties *p = (ElasticProperties *)properties;
	
	// Add to total strain
	Tensor *ep=mptr->GetStrainTensor();
    ep->xx+=dvxx;
    ep->yy+=dvyy;
    double dgam=dvxy+dvyx;
    ep->xy+=dgam;
	double dwrotxy=dvyx-dvxy;
	
    // residual strains (isotropic thermal and moisture)
	double eres = CTE1*res->dT;
	if(DiffusionTask::active)
		eres += CME1*res->dC;
	
    // moisture and thermal strain and temperature change
	//   (when diffusion, conduction, OR thermal ramp active)
    double dvxxeff = dvxx - eres;
    double dvyyeff = dvyy - eres;
    
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;
	
	// find stress (Units N/m^2  mm^3/g)
	// this does xx, yy, ans xy only. zz do later if needed
	double c1,c2,c3;
	if(np==AXISYMMETRIC_MPM)
	{	// axisymmetric strain
		ep->zz += dvzz;
		
		// hoop stress affect on RR, ZZ, and RZ stresses
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
	Hypo2DCalculations(mptr,-dwrotxy,c1,c2,c3);
    
	// work and resdiaul strain energy increments
	double workEnergy = 0.5*((st0.xx+sp->xx)*dvxx + (st0.yy+sp->yy)*dvyy + (st0.xy+sp->xy)*dgam);
	double resEnergy = 0.5*(st0.xx+sp->xx + st0.yy+sp->yy)*eres;
	if(np==PLANE_STRAIN_MPM)
	{	// need to add back terms to get from reduced cte to actual cte
		sp->zz += p->C[4][1]*(dvxx+p->alpha[5]*eres)+p->C[4][2]*(dvyy+p->alpha[6]*eres)-p->C[4][4]*eres;
		
		// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
		resEnergy += 0.5*(st0.zz+sp->zz)*eres;
	}
	else if(np==PLANE_STRESS_MPM)
	{	ep->zz += p->C[4][1]*dvxxeff+p->C[4][2]*dvyyeff;
	}
	else
	{	// axisymmetric hoop stress
		sp->zz += p->C[4][1]*dvxxeff + p->C[4][2]*dvyyeff + p->C[4][4]*(dvzz - eres);
		
		// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
		workEnergy += 0.5*(st0.zz+sp->zz)*dvzz;
		resEnergy += 0.5*(st0.zz+sp->zz)*eres;
	}
	mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
	
    // track heat energy
    IncrementHeatEnergy(mptr,res->dT,0.,0.);
}

/* For 3D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
	Assumes linear elastic and isotropic, uses hypoelastic correction
*/
void IsotropicMat::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
						  double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// cast pointer to material-specific data
	ElasticProperties *p = (ElasticProperties *)properties;
	
	// Add to total strain
	Tensor *ep = mptr->GetStrainTensor();
    ep->xx += dvxx;
    ep->yy += dvyy;
    ep->zz += dvzz;
    double dgamxy = dvxy+dvyx;
    ep->xy += dgamxy;
    double dgamxz = dvxz+dvzx;
    ep->xz += dgamxz;
    double dgamyz = dvyz+dvzy;
    ep->yz += dgamyz;
	
	// rotational strain increments (particle updated by Hypo3D)
	double dwrotxy = dvyx-dvxy;
	double dwrotxz = dvzx-dvxz;
	double dwrotyz = dvzy-dvyz;
	
    // residual strains (thermal and moisture) (isotropic only)
	double eres = CTE1*res->dT;
	if(DiffusionTask::active)
		eres += CME1*res->dC;
	
	// effective strains
	double dvxxeff = dvxx-eres;
	double dvyyeff = dvyy-eres;
	double dvzzeff = dvzz-eres;
	
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
    
    // track heat energy
    IncrementHeatEnergy(mptr,res->dT,0.,0.);
	
}

#ifdef USE_PSEUDOHYPERELASTIC
/* Take increments in strain and calculate new Particle: strains, rotation strain,
		stresses, strain energy,
	dvij are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
*/
void IsotropicMat::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// current previous deformation gradient and stretch
	const Matrix3 pFnm1 = mptr->GetDeformationGradientMatrix();
	Matrix3 Vnm1 = pFnm1.LeftDecompose(NULL,NULL);
	
    // get incremental deformation gradient and decompose it
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
    Matrix3 dR;
    Matrix3 dV = dF.LeftDecompose(&dR,NULL);
	
	// get strain increments de = (dV-I) dR Vnma dRT
	dV(0,0) -= 1.;
	dV(1,1) -= 1.;
	dV(2,2) -= 1.;
	Matrix3 de = dV*Vnm1.RMRT(dR);
	
	// Update total deformation gradient
	Matrix3 pF = dF*pFnm1;
	mptr->SetDeformationGradientMatrix(pF);
	
	// cast pointer to material-specific data
	ElasticProperties *p = (ElasticProperties *)properties;
	
	// save initial stresses
	Tensor *sp = mptr->GetStressTensor();
	Tensor st0 = *sp;
	
	// stress increments
	double delsp[6];
	if(np==THREED_MPM)
	{	// residual strains (thermal and moisture)
		double eres = CTE1*res->dT;
		if(DiffusionTask::active)
			eres += CME1*res->dC;
		
		// effective strains
		double dvxxeff = de(0,0)-eres;
		double dvyyeff = de(1,1)-eres;
		double dvzzeff = de(2,2)-eres;
		double dgamyz = de(1,2)+de(2,1);
		double dgamxz = de(0,2)+de(2,0);
		double dgamxy = de(0,1)+de(1,0);
		
		for(int i=0;i<3;i++)
		{   delsp[i] = p->C[i][0]*dvxxeff + p->C[i][1]*dvyyeff + p->C[i][2]*dvzzeff;
		}
		delsp[3] = p->C[3][3]*dgamyz;
		delsp[4] = p->C[4][4]*dgamxz;
		delsp[5] = p->C[5][5]*dgamxy;
		
		// incremental rotate of prior strain
		Matrix3 stn(sp->xx,sp->xy,sp->xz,sp->xy,sp->yy,sp->yz,sp->xz,sp->yz,sp->zz);
		Matrix3 str = stn.RMRT(dR);
		
		sp->xx = str(0,0)+delsp[0];
		sp->yy = str(1,1)+delsp[1];
		sp->zz = str(2,2)+delsp[2];
		sp->yz = str(1,2)+delsp[3];
		sp->xz = str(0,2)+delsp[4];
		sp->xy = str(0,1)+delsp[5];
		
		// work energy increment per unit mass (dU/(rho0 V0))
		mptr->AddWorkEnergyAndResidualEnergy( 0.5*((st0.xx+sp->xx)*du(0,0) + (st0.yy+sp->yy)*du(1,1)
												   + (st0.zz+sp->zz)*du(2,2)  + (st0.yz+sp->yz)*(du(1,2)+du(2,1))
												   + (st0.xz+sp->xz)*(du(0,2)+du(2,0)) + (st0.xy+sp->xy)*(du(0,1)+du(1,0))),
											 0.5*(st0.xx+sp->xx + st0.yy+sp->yy + st0.zz+sp->zz)*eres);
	}
	
	else
	{	// residual strains (thermal and moisture)
		double eres = CTE1*res->dT;
		if(DiffusionTask::active)
			eres += CME1*res->dC;
		
		// effective strains
		double dvxxeff = de(0,0)-eres;
		double dvyyeff = de(1,1)-eres;
		double dgamxy = de(0,1)+de(1,0);
		
		if(np==AXISYMMETRIC_MPM)
		{	// hoop stress affect on RR, ZZ, and RZ stresses
			double dvzzeff = de(2,2) - eres;
			delsp[0] = p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff + p->C[4][1]*dvzzeff;
			delsp[1] = p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff + p->C[4][2]*dvzzeff;
			delsp[5] = p->C[3][3]*dgamxy;
		}
		else
		{	delsp[0] = p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff;
			delsp[1] = p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff;
			delsp[5] = p->C[3][3]*dgamxy;
		}
		
		// rotate previous stress
		Matrix3 stn(sp->xx,sp->xy,sp->xy,sp->yy,sp->zz);
		Matrix3 str = stn.RMRT(dR);
		
		// update in plane stress
		sp->xx = str(0,0)+delsp[0];
		sp->yy = str(1,1)+delsp[1];
		sp->xy = str(0,1)+delsp[5];
		
		// work and resdidual strain energy increments
		double workEnergy = 0.5*((st0.xx+sp->xx)*de(0,0) + (st0.yy+sp->yy)*de(1,1) + (st0.xy+sp->xy)*dgamxy);
		double resEnergy = 0.5*(st0.xx+sp->xx + st0.yy+sp->yy)*eres;
		if(np==PLANE_STRAIN_MPM)
		{	// need to add back terms to get from reduced cte to actual cte
			sp->zz += p->C[4][1]*(de(0,0)+p->alpha[5]*eres)+p->C[4][2]*(de(1,1)+p->alpha[6]*eres)-p->C[4][4]*eres;
			
			// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			resEnergy += 0.5*(st0.zz+sp->zz)*eres;
		}
		else if(np==PLANE_STRESS_MPM)
		{	Tensor *ep=mptr->GetStrainTensor();
			ep->zz += p->C[4][1]*dvxxeff+p->C[4][2]*dvyyeff+eres;
		}
		else
		{	// axisymmetric hoop stress
			sp->zz += p->C[4][1]*dvxxeff + p->C[4][2]*dvyyeff + p->C[4][4]*(de(2,2) - eres);
			
			// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			workEnergy += 0.5*(st0.zz+sp->zz)*de(2,2);
			resEnergy += 0.5*(st0.zz+sp->zz)*eres;
		}
		mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
	}
    
    // track heat energy
    IncrementHeatEnergy(mptr,res->dT,0.,0.);
}
#endif

#pragma mark IsotropicMat::Accessors

// Calculate wave speed in L/sec
// Because G in F/L^2 = mass/(L sec^2) and rho in mass/L^3, G/rho in L^2/sec^2
// Uses sqrt((K +4G/3)/rho) which is dilational wave speed
// Identity also: K + 4G/3 = Lambda + 2G = 2G(1-nu)/(1-2 nu)
double IsotropicMat::WaveSpeed(bool threeD,MPMBase *mptr) const
{	return sqrt(2.*G*(1.-nu)/(rho*(1.-2.*nu)));
}

// Calculate shear wave speed in L/sec
double IsotropicMat::ShearWaveSpeed(bool threeD,MPMBase *mptr) const { return sqrt(G/rho); }



