/********************************************************************************
    ElasticMPM.cpp - more Elastic for MPM code
    NairnMPM
    
    Created by John Nairn on Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/Elastic.hpp"
#include "MPM_Classes/MPMBase.hpp"
//#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Global_Quantities/ThermalRamp.hpp"

#pragma mark Elastic::Methods

/* For 2D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
	Assumes linear elastic, uses hypoelastic correction
   For Axisymmetric MPM, x->R, y->Z, x->theta, and dvzz is change in hoop strain
	(i.e., du/r on particle and dvzz will be zero if not axisymmetric)
*/
void Elastic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,double dvzz,
                          double delTime,int np,void *properties,ResidualStrains *res)
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
	
    // residual strains (thermal and moisture)
	double erxx=p->alpha[1]*res->dT;
	double eryy=p->alpha[2]*res->dT;
	double erxy=p->alpha[3]*res->dT;
	double erzz=CTE3*res->dT;
	if(DiffusionTask::active)
	{	erxx+=p->beta[1]*res->dC;
		eryy+=p->beta[2]*res->dC;
		erxy+=p->beta[3]*res->dC;
		erzz+=CME3*res->dC;
	}
	
    // moisture and thermal strain and temperature change
	//   (when diffusion, conduction, OR thermal ramp active)
    dvxx -= erxx;
    dvyy -= eryy;
	dgam -= erxy;
    
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;
	
	// find stress (Units N/m^2  cm^3/g)
	// this does xx, yy, ans xy only. zz do later if needed
	double c1,c2,c3;
	if(np==AXISYMMETRIC_MPM)
	{	// axisymmetric strain
		ep->zz += dvzz;
		
		// hoop stress affect on RR, ZZ, and RZ stresses
		dvzz -= erzz;
		c1=p->C[1][1]*dvxx + p->C[1][2]*dvyy + p->C[4][1]*dvzz + p->C[1][3]*dgam;
		c2=p->C[1][2]*dvxx + p->C[2][2]*dvyy + p->C[4][2]*dvzz + p->C[2][3]*dgam;
		c3=p->C[1][3]*dvxx + p->C[2][3]*dvyy + p->C[4][3]*dvzz + p->C[3][3]*dgam;
	}
	else
    {	c1=p->C[1][1]*dvxx + p->C[1][2]*dvyy + p->C[1][3]*dgam;
		c2=p->C[1][2]*dvxx + p->C[2][2]*dvyy + p->C[2][3]*dgam;
		c3=p->C[1][3]*dvxx + p->C[2][3]*dvyy + p->C[3][3]*dgam;
	}
	Hypo2DCalculations(mptr,-dwrotxy,c1,c2,c3);
    
	// out of plane stress or strain and energy incrment
	if(np==PLANE_STRAIN_MPM)
	{	// need to add back terms to get from reduced cte to actual cte
		sp->zz += p->C[4][1]*(dvxx+p->alpha[5]*erzz)+p->C[4][2]*(dvyy+p->alpha[6]*erzz)
					+p->C[4][3]*(dgam+p->alpha[7]*erzz)-p->C[4][4]*erzz;
		
		// strain energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (uJ/g)
		mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dvxx + (st0.yy+sp->yy)*dvyy
								   + (st0.xy+sp->xy)*dgam)-0.5*(st0.zz+sp->zz)*erzz);
	}
	else if(np==PLANE_STRESS_MPM)
	{	ep->zz += p->C[4][1]*dvxx+p->C[4][2]*dvyy+p->C[4][3]*dgam+erzz;
		
		// strain energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (uJ/g)
		mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dvxx + (st0.yy+sp->yy)*dvyy
								   + (st0.xy+sp->xy)*dgam));
	}
	else
	{	// axisymmetric hoop stress
		sp->zz += p->C[4][1]*dvxx + p->C[4][2]*dvyy + p->C[4][4]*dvzz + p->C[4][3]*dgam;
		
		// strain energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (uJ/g)
		mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dvxx + (st0.yy+sp->yy)*dvyy
								   + (st0.xy+sp->xy)*dgam) + (st0.zz+sp->zz)*dvzz);
	}
	
    // track heat energy
    IncrementHeatEnergy(mptr,res->dT,0.,0.);
}

/* For 3D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
    Assumes linear elastic, uses hypoelastic correction
*/
void Elastic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
        double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np,void *properties,ResidualStrains *res)
{
	// cast pointer to material-specific data
	ElasticProperties *p = (ElasticProperties *)properties;
	
	// Add to total strain
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
	
    // residual strains (thermal and moisture) (isotropic only)
	dvxx-=p->alpha[0]*res->dT;
	dvyy-=p->alpha[1]*res->dT;
	dvzz-=p->alpha[2]*res->dT;
	dgamyz-=p->alpha[3]*res->dT;
	dgamxz-=p->alpha[4]*res->dT;
	dgamxy-=p->alpha[5]*res->dT;
	if(DiffusionTask::active)
	{	dvxx-=p->beta[0]*res->dC;
		dvyy-=p->beta[1]*res->dC;
		dvzz-=p->beta[2]*res->dC;
		dgamyz-=p->beta[3]*res->dC;
		dgamxz-=p->beta[4]*res->dC;
		dgamxy-=p->beta[5]*res->dC;
	}
	
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;
	
	// stress increments
	double delsp[6];
	int i;
	for(i=0;i<6;i++)
		delsp[i]=p->C[i][0]*dvxx + p->C[i][1]*dvyy + p->C[i][2]*dvzz + p->C[i][3]*dgamyz + p->C[i][4]*dgamxz + p->C[i][5]*dgamxy;
	
	// update stress (need to make hypoelastic)
	Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,delsp);
	
	// strain energy increment per unit mass (dU/(rho0 V0))
    mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dvxx
                            + (st0.yy+sp->yy)*dvyy
                            + (st0.zz+sp->zz)*dvzz
                            + (st0.yz+sp->yz)*dgamyz
                            + (st0.xz+sp->xz)*dgamxz
                            + (st0.xy+sp->xy)*dgamxy));
    
    // track heat energy
    IncrementHeatEnergy(mptr,res->dT,0.,0.);

}
