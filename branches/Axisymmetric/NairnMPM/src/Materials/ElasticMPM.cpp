/******************************************************************************** 
    ElasticMPM.cpp - more Elastic for MPM code
    NairnMPM
    
    Created by John Nairn on Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/Elastic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Global_Quantities/ThermalRamp.hpp"

#pragma mark Elastic::Methods

/* For 2D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
    Assumes linear elastic, uses hypoelastic correction
*/
void Elastic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
        double delTime,int np)
{
	// Add to total strain
	Tensor *ep=mptr->GetStrainTensor();
    ep->xx+=dvxx;
    ep->yy+=dvyy;
    double dgam=dvxy+dvyx;
    ep->xy+=dgam;
	double dwrotxy=dvyx-dvxy;
	
    // residual strains (thermal and moisture)
	double erxx=me0[1]*ConductionTask::dTemperature;
	double eryy=me0[2]*ConductionTask::dTemperature;
	double erxy=me0[3]*ConductionTask::dTemperature;
	double erzz=CTE3*ConductionTask::dTemperature;
	if(DiffusionTask::active)
	{	erxx+=mc0[1]*DiffusionTask::dConcentration;
		eryy+=mc0[2]*DiffusionTask::dConcentration;
		erxy+=mc0[3]*DiffusionTask::dConcentration;
		erzz+=CME3*DiffusionTask::dConcentration;
	}
	
    // thermal strain and temperature change (if conduction OR thermal ramp active)
    dvxx-=erxx;
    dvyy-=eryy;
	dgam-=erxy;
    
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;
	
    /* ---------------------------------------------------
        find stress (Units N/m^2  cm^3/g)
    */
    double c1=mdm[1][1]*dvxx + mdm[1][2]*dvyy + mdm[1][3]*dgam;
    double c2=mdm[1][2]*dvxx + mdm[2][2]*dvyy + mdm[2][3]*dgam;
    double c3=mdm[1][3]*dvxx + mdm[2][3]*dvyy + mdm[3][3]*dgam;
	Hypo2DCalculations(mptr,-dwrotxy,c1,c2,c3);
    
	// out of plane stress or strain
	if(np==PLANE_STRAIN_MPM)
	{	// need to add back terms to get from reduced cte to actual cte
		sp->zz+=mdm[4][1]*(dvxx+me0[5]*erzz)+mdm[4][2]*(dvyy+me0[6]*erzz)
								+mdm[4][3]*(dgam+me0[7]*erzz)-mdm[4][4]*erzz;
		
	}
	else
		ep->zz+=mdm[4][1]*dvxx+mdm[4][2]*dvyy+mdm[4][3]*dgam+erzz;
	
    // strain energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule)
    mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dvxx
                            + (st0.yy+sp->yy)*dvyy
                            + (st0.xy+sp->xy)*dgam));
	if(np==PLANE_STRAIN_MPM)
	{	mptr->AddStrainEnergy(-0.5*(st0.zz+sp->zz)*erzz);
	}
}

/* For Axismmetric MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
    Assumes linear elastic, uses hypoelastic correction
    Here x -> r, y -> z, and z -> theta directions
*/
void Elastic::MPMConstLaw(MPMBase *mptr,double dvrr,double dvzz,double dvrz,double dvzr,double dvtt,
                          double delTime,int np)
{
	// Add to total strain
	Tensor *ep=mptr->GetStrainTensor();
    ep->xx += dvrr;
    ep->yy += dvzz;
    double dgam = dvrz+dvzr;
    ep->xy += dgam;
	double dwrotrz = dvzr-dvrz;
	
    // residual strains (thermal and moisture)
	double errr = me0[1]*ConductionTask::dTemperature;
	double erzz = me0[2]*ConductionTask::dTemperature;
	double errz = me0[3]*ConductionTask::dTemperature;
	double ertt = CTE3*ConductionTask::dTemperature;
	if(DiffusionTask::active)
	{	errr += mc0[1]*DiffusionTask::dConcentration;
		erzz += mc0[2]*DiffusionTask::dConcentration;
		errz += mc0[3]*DiffusionTask::dConcentration;
		ertt += CME3*DiffusionTask::dConcentration;
	}
	
    // thermal strain and temperature change (if conduction OR thermal ramp active)
    dvrr -= errr;
    dvzz -= erzz;
	dgam -= errz;
    dvtt -= ertt;
    
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;
	
    /* ---------------------------------------------------
        find stress (Units N/m^2  cm^3/g)
    */
    double c1 = mdm[1][1]*dvrr + mdm[1][2]*dvzz + mdm[1][3]*dgam + mdm[4][1]*dvtt;
    double c2 = mdm[1][2]*dvrr + mdm[2][2]*dvzz + mdm[2][3]*dgam + mdm[4][2]*dvtt;
    double c3 = mdm[1][3]*dvrr + mdm[2][3]*dvzz + mdm[3][3]*dgam + mdm[4][3]*dvtt;
	Hypo2DCalculations(mptr,-dwrotrz,c1,c2,c3);
    
	// hoop stress
    sp->zz += mdm[4][1]*dvrr + mdm[4][2]*dvzz + mdm[4][3]*dgam + mdm[4][4]*dvtt;
	
    // strain energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule)
    mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dvrr + (st0.yy+sp->yy)*dvzz
                               + (st0.xy+sp->xy)*dgam) + (st0.zz+sp->zz)*dvtt);
}

/* For 3D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
    Assumes linear elastic, uses hypoelastic correction
*/
void Elastic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
        double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
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
	dvxx-=me0[0]*ConductionTask::dTemperature;
	dvyy-=me0[1]*ConductionTask::dTemperature;
	dvzz-=me0[2]*ConductionTask::dTemperature;
	dgamyz-=me0[3]*ConductionTask::dTemperature;
	dgamxz-=me0[4]*ConductionTask::dTemperature;
	dgamxy-=me0[5]*ConductionTask::dTemperature;
	if(DiffusionTask::active)
	{	dvxx-=mc0[0]*DiffusionTask::dConcentration;
		dvyy-=mc0[1]*DiffusionTask::dConcentration;
		dvzz-=mc0[2]*DiffusionTask::dConcentration;
		dgamyz-=mc0[3]*DiffusionTask::dConcentration;
		dgamxz-=mc0[4]*DiffusionTask::dConcentration;
		dgamxy-=mc0[5]*DiffusionTask::dConcentration;
	}
	
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;
	
	// stress increments
	double delsp[6];
	int i;
	for(i=0;i<6;i++)
		delsp[i]=mdm[i][0]*dvxx + mdm[i][1]*dvyy + mdm[i][2]*dvzz + mdm[i][3]*dgamyz + mdm[i][4]*dgamxz + mdm[i][5]*dgamxy;
	
	// update stress (need to make hypoelastic)
	Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,delsp);
	
	// strain energy increment per unit mass (dU/(rho0 V0))
    mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dvxx
                            + (st0.yy+sp->yy)*dvyy
                            + (st0.zz+sp->zz)*dvzz
                            + (st0.yz+sp->yz)*dgamyz
                            + (st0.xz+sp->xz)*dgamxz
                            + (st0.xy+sp->xy)*dgamxy));
}
