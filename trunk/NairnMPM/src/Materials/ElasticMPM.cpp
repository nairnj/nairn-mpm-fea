/********************************************************************************
    ElasticMPM.cpp - more Elastic for MPM code
    nairn-mpm-fea
    
    Created by John Nairn on Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/Elastic.hpp"
#include "MPM_Classes/MPMBase.hpp"
//#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Global_Quantities/ThermalRamp.hpp"

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
void Elastic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 dv,double delTime,int np,void *properties,ResidualStrains *res) const
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
	// get previous deformation gradient
	Matrix3 pFnm1 = mptr->GetDeformationGradientMatrix();
	
    // get incremental deformation gradient and decompose it
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
    Matrix3 dR;
    Matrix3 dU = dF.RightDecompose(&dR,NULL);

	// get pinitial rotation R0
	Matrix3 R0 = mptr->GetInitialRotation();
	
#ifdef TRACK_RTOT
	// read previous rotation
	Matrix3 *Rtotnm1 = mptr->GetRtotPtr();
	
	// get incremental strain
	dU(0,0) -= 1.;
	dU(1,1) -= 1.;
	dU(2,2) -= 1.;
	Matrix3 RtotTUmI = (*Rtotnm1).Transpose()*dU;
	Matrix3 Fnm1tot = pFnm1*R0;
	Matrix3 de = RtotTUmI*Fnm1tot;
	
	// get total rotation and update
	Matrix3 Rtot = dR*(*Rtotnm1);
	mptr->SetRtot(Rtot);
#else
	// get previous rotation and stretch
    Matrix3 Rnm1;
	Matrix3 Unm1 = pFnm1.RightDecompose(&Rnm1,NULL);
	
	// get strain increments de = R0T.[(Rnm1T.dU.Rnm1 - I).Unm1].R0
	Matrix3 dUrot = dU.RTMR(Rnm1);
	dUrot(0,0) -= 1.;
	dUrot(1,1) -= 1.;
	dUrot(2,2) -= 1.;
	dUrot *= Unm1;
	
	// apply initial rotation to get strain increment in the material coordinates
	Matrix3 de = dUrot.RTMR(R0);
	Matrix3 Rtotnm1M3 = Rnm1*R0;
	Matrix3 *Rtotnm1 = &Rtotnm1M3;
	
	// get total rotation
	Matrix3 Rtot = dR*Rtotnm1M3;
#endif
	
	// Update total deformation gradient
	Matrix3 pF = dF*pFnm1;
	mptr->SetDeformationGradientMatrix(pF);
	
	// cast pointer to material-specific data
	ElasticProperties *p = GetElasticPropertiesPointer(properties);

    // residual strains (thermal and moisture) in material axes
	double exxr,eyyr,ezzr;
	if(np==THREED_MPM)
	{	exxr = p->alpha[0]*res->dT;
		eyyr = p->alpha[1]*res->dT;
		ezzr = p->alpha[2]*res->dT;
		if(DiffusionTask::active)
		{	exxr += p->beta[0]*res->dC;
			eyyr += p->beta[1]*res->dC;
			ezzr += p->beta[2]*res->dC;
		}
	}
	else
	{	exxr = p->alpha[1]*res->dT;
		eyyr = p->alpha[2]*res->dT;
		ezzr = p->alpha[4]*res->dT;
		if(DiffusionTask::active)
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
// Similar, strain increments are rotated back to find work energy (done in global system)
void Elastic::LRElasticConstitutiveLaw(MPMBase *mptr,Matrix3 de,Matrix3 er,Matrix3 Rtot,Matrix3 dR,
									   Matrix3 *Rnmatot,int np,void *properties,ResidualStrains *res) const
{
	// effective strains
	double dvxxeff = de(0,0)-er(0,0);
	double dvyyeff = de(1,1)-er(1,1);
	double dvzzeff = de(2,2)-er(2,2);
	double dgamxy = 2.*de(0,1);
	
    // save initial stresses
	Tensor *sp=mptr->GetStressTensor();
    Tensor st0=*sp;
	
	// stress increments
	// cast pointer to material-specific data
	ElasticProperties *p = GetElasticPropertiesPointer(properties);
	if(np==THREED_MPM)
	{	double dgamyz = 2.*de(1,2);
		double dgamxz = 2.*de(0,2);
		
		// update sigma = dR signm1 dRT + Rtot dsigma RtotT
		double dsigyz = p->C[3][3]*dgamyz;
		double dsigxz = p->C[4][4]*dgamxz;
		double dsigxy = p->C[5][5]*dgamxy;
		Matrix3 dsig(p->C[0][0]*dvxxeff + p->C[0][1]*dvyyeff + p->C[0][2]*dvzzeff, dsigxy, dsigxz,
					 dsigxy, p->C[1][0]*dvxxeff + p->C[1][1]*dvyyeff + p->C[1][2]*dvzzeff, dsigyz,
					 dsigxz, dsigyz, p->C[2][0]*dvxxeff + p->C[2][1]*dvyyeff + p->C[2][2]*dvzzeff);
		Matrix3 dsigrot = dsig.RMRT(Rtot);
		Matrix3 stn(sp->xx,sp->xy,sp->xz,sp->xy,sp->yy,sp->yz,sp->xz,sp->yz,sp->zz);
		Matrix3 str = stn.RMRT(dR);
		sp->xx = str(0,0) + dsigrot(0,0);
		sp->yy = str(1,1) + dsigrot(1,1);
		sp->zz = str(2,2) + dsigrot(2,2);
		sp->xy = str(0,1) + dsigrot(0,1);
		sp->xz = str(0,2) + dsigrot(0,2);
		sp->yz = str(1,2) + dsigrot(1,2);
		
		// stresses are in global coordinates so need to rotate strain and residual
		// strain to get work energy increment per unit mass (dU/(rho0 V0))
		Matrix3 derot = de.RMRT(Rtot);
		Matrix3 errot = er.RMRT(Rtot);
		mptr->AddWorkEnergyAndResidualEnergy(0.5*((st0.xx+sp->xx)*derot(0,0) + (st0.yy+sp->yy)*derot(1,1)
												  + (st0.zz+sp->zz)*derot(2,2)) + (st0.yz+sp->yz)*derot(1,2)
												  + (st0.xz+sp->xz)*derot(0,2) + (st0.xy+sp->xy)*derot(0,1),
											 0.5*((st0.xx+sp->xx)*errot(0,0) + (st0.yy+sp->yy)*errot(1,1)
												  + (st0.zz+sp->zz)*errot(2,2)) + (st0.yz+sp->yz)*errot(1,2)
												  + (st0.xz+sp->xz)*errot(0,2) + (st0.xy+sp->xy)*errot(0,1));
	}
	else
	{	// find stress increment
		// this does xx, yy, and xy only. zz done later if needed
		double dsxx,dsyy;
		if(np==AXISYMMETRIC_MPM)
		{	dsxx = p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff + p->C[4][1]*dvzzeff ;
			dsyy = p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff + p->C[4][2]*dvzzeff;
		}
		else
		{	dsxx = p->C[1][1]*dvxxeff + p->C[1][2]*dvyyeff;
			dsyy = p->C[1][2]*dvxxeff + p->C[2][2]*dvyyeff;
		}
		double dsxy = p->C[3][3]*dgamxy;
		
		// update sigma = dR signm1 dRT + Rtot dsigma RtotT
		Matrix3 dsig(dsxx,dsxy,dsxy,dsyy,0.);
		Matrix3 dsigrot = dsig.RMRT(Rtot);
		Matrix3 stn(sp->xx,sp->xy,sp->xy,sp->yy,sp->zz);
		Matrix3 str = stn.RMRT(dR);
		sp->xx = str(0,0) + dsigrot(0,0);
		sp->yy = str(1,1) + dsigrot(1,1);
		sp->xy = str(0,1) + dsigrot(0,1);
		
		// stresses are in global coordinate so need to rotate strain and residual
		// strain to get work energy increment per unit mass (dU/(rho0 V0))
		double ezzr = er(2,2);
		Matrix3 derot = de.RMRT(Rtot);
		Matrix3 errot = er.RMRT(Rtot);
		
		// work and residual strain energy increments and sigma or F in z direction
		double workEnergy = 0.5*((st0.xx+sp->xx)*derot(0,0) + (st0.yy+sp->yy)*derot(1,1)) + (st0.xy+sp->xy)*derot(0,1);
		double resEnergy = 0.5*((st0.xx+sp->xx)*errot(0,0) + (st0.yy+sp->yy)*errot(1,1)) + (st0.xy+sp->xy)*errot(0,1);
		if(np==PLANE_STRAIN_MPM)
		{	// need to add back terms to get from reduced cte to actual cte
			sp->zz += p->C[4][1]*(de(0,0)+p->alpha[5]*ezzr) + p->C[4][2]*(de(1,1)+p->alpha[6]*ezzr) - p->C[4][4]*ezzr;
			
			// extra residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			resEnergy += 0.5*(st0.zz+sp->zz)*ezzr;
		}
		else if(np==PLANE_STRESS_MPM)
		{	// zz deformation
			mptr->IncrementDeformationGradientZZ(p->C[4][1]*dvxxeff + p->C[4][2]*dvyyeff + ezzr);
		}
		else
		{	// axisymmetric hoop stress
			sp->zz += p->C[4][1]*dvxxeff + p->C[4][2]*dvyyeff + p->C[4][4]*dvzzeff;
			
			// extra work and residual energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
			workEnergy += 0.5*(st0.zz+sp->zz)*de(2,2);
			resEnergy += 0.5*(st0.zz+sp->zz)*ezzr;
		}
		mptr->AddWorkEnergyAndResidualEnergy(workEnergy, resEnergy);
		
	}
    
    // track heat energy
    IncrementHeatEnergy(mptr,res->dT,0.,0.);
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
	if(DiffusionTask::active)
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
	Hypo2DCalculations(mptr,dwrotxy,dvxx+dvyy,c1,c2,c3);
    
	// work and resdiaul strain energy increments
	double workEnergy = 0.5*((st0.xx+sp->xx)*dvxx + (st0.yy+sp->yy)*dvyy + (st0.xy+sp->xy)*dgam);
	double resEnergy = 0.5*((st0.xx+sp->xx)*erxx + (st0.yy+sp->yy)*eryy + (st0.xy+sp->xy)*erxy);
	if(np==PLANE_STRAIN_MPM)
	{	// need to add back terms to get from reduced cte to actual cte
		sp->zz += p->C[4][1]*(dvxx+p->alpha[5]*erzz)+p->C[4][2]*(dvyy+p->alpha[6]*erzz)
					+p->C[4][3]*(dgam+p->alpha[7]*erzz)-p->C[4][4]*erzz;
		
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
	
    // track heat energy
    IncrementHeatEnergy(mptr,res->dT,0.,0.);
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
	if(DiffusionTask::active)
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
    
    // track heat energy
    IncrementHeatEnergy(mptr,res->dT,0.,0.);

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

