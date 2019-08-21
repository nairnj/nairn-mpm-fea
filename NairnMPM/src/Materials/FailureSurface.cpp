/********************************************************************************
	FailureSurface.cpp
	nairn-mpm-fea

	Created by John Nairn, June 26, 2015.
	Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		FailureSurface.hpp
********************************************************************************/

#include "stdafx.h"
#include "Materials/FailureSurface.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"

#pragma mark FailureSurface::Constructors and Destructors

// Constructors
FailureSurface::FailureSurface() {}

// Main Constructor
FailureSurface::FailureSurface(MaterialBase *pair)
{
	criticalNormal = 1.e20;
	criticalShear = 1.e20;
	parent = pair;
}

#pragma mark FailureSurface::Initialize

// Read hardening law properties
char *FailureSurface::InputInitationProperty(char *xName,int &input,double &gScaling)
{
    // principle stress failure loads
    if(strcmp(xName,"sigmac")==0)
    {   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalNormal,gScaling,1.e6);
    }
    
    else if(strcmp(xName,"tauc")==0)
    {   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalShear,gScaling,1.e6);
    }
	
    // is not a hardening law property
    return NULL;
}

// verify settings and some initial calculations
const char *FailureSurface::VerifyAndLoadProperties(int np)
{
	// reduced failure stresses
	double rho = parent->GetRho(NULL);
    sigmacRed = criticalNormal/rho;
    taucRed = criticalShear/rho;
	
	return NULL;
}

// print just initiation properties to output window
void FailureSurface::PrintInitiationProperties(void) const
{
    cout << GetInitiationLawName() << endl;
    MaterialBase::PrintProperty("sigc",criticalNormal*UnitsController::Scaling(1.e-6),"");
    MaterialBase::PrintProperty("tauc",criticalShear*UnitsController::Scaling(1.e-6),"");
    cout << endl;
}

#pragma mark FailureSurface::Methods

// check stress state in initial axis system to see if failed
// If failed, fill the normal vector, and return failure mode
// If not failed, return 0
int FailureSurface::ShouldInitiateFailure(Tensor *str,Vector *normal,int np,double relStrength) const
{
	// initialize
	int failureMode = NO_FAILURE;
	
	// get principal stress
	if(np==THREED_MPM)
	{	// get principle stress
		Matrix3 strmat = TensorToMatrix(str,true);
		Vector lam = strmat.Eigenvalues();
		double sigma1 = fmax(lam.x,lam.y);
		sigma1 = fmax(sigma1,lam.z);
		double maxShear = fabs(0.5*(lam.x-lam.y));
		maxShear = fmax(maxShear,fabs(0.5*(lam.x-lam.z)));
		maxShear = fmax(maxShear,fabs(0.5*(lam.y-lam.z)));
		
		// check maximum principle stress and maximum shear stress
		if(sigma1>=sigmacRed*relStrength || maxShear>=taucRed*relStrength)
		{	// get eigenvectors
			Matrix3 R = strmat.Eigenvectors(lam);
			
			// get largest one first
			if(lam.y>lam.z)
			{	if(lam.y>lam.x)
				{	R.SwapColumns(0,1);
					double temp = lam.x;
					lam.x=lam.y;
					lam.y=temp;
				}
			}
			else if(lam.z>lam.x)
			{	R.SwapColumns(0,2);
				double temp = lam.x;
				lam.x=lam.z;
				lam.z=temp;
			}
			
			// flip z if determinant is < 0
			if(R.determinant()<0.)
			{	R(0,2) = -R(0,2);
				R(1,2) = -R(1,2);
				R(2,2) = -R(2,2);
			}
			
			// if shear failure, rotate +/-45 about axis with the middle value
			// must rotate about the principle axis using axis and angle rotation matrix
			// For details see https://en.wikipedia.org/wiki/Rotation_matrix
			if(sigma1 < sigmacRed*relStrength)
			{	failureMode = SHEAR_FAILURE;
				
				double c45 = 0.5*sqrt(2.);				// cos(45deg)
				double s45 = c45;						// sin(45deg)
				double ux,uy,uz;
				if(lam.y>lam.z)
				{	// rotate ccw by +45 about principle y axis
					ux = R(0,1);
					uy = R(1,1);
					uz = R(2,1);
				}
				else
				{	// rotate ccw by -45 about prinviple z axis
					ux = R(0,2);
					uy = R(1,2);
					uz = R(2,2);
					s45 = -s45;							// sin(-45deg)
				}
				Matrix3 R45 = Matrix3(c45+ux*ux*(1.-c45),ux*uy*(1.-c45)-uz*s45,ux*uz*(1.-c45)+uy*s45,
									  uy*ux*(1.-c45)+uz*s45,c45+uy*uy*(1.-c45),uy*uz*(1.-c45)-ux*s45,
									  uz*ux*(1.-c45)-uy*s45,uz*uy*(1.-c45)+ux*s45,c45+uz*uz*(1.-c45));
				R = R45*R;
			}
			else
				failureMode = TENSILE_FAILURE;
			
			// decode to Euler angles
			// see https://en.wikipedia.org/wiki/Euler_angles - but it has the equation wrong
			// There are two solutions (where second angle is one way or the other) for ZYZ Euler angles
			//   and either can be used for softening material to reconstruct the ZYZ R matrix
			// normal to crack will be first column on reconstructed R matrix
			if(R(2,2)<-1.)
				normal->y = acos(-1.);
			else if(R(2,2)>1.)
				normal->y = acos(1.);
			else
				normal->y = acos(R(2,2));				// beta  = arccos(Z3)
			if(R(1,2)==0. && R(0,2)==0.)
				normal->x = 0.;
			else
				normal->x = atan2(R(1,2),R(0,2));		// alpha = atan2(Z2,Z1)
			if(R(2,1)==0. && R(2,0)==0.)
				normal->z = 0.;
			else
				normal->z = atan2(R(2,1),-R(2,0));		// gamma = atan2(Y3,-X3)
#ifdef CHECK_NAN
			if(normal->x!=normal->x || normal->y!=normal->y || normal->z!=normal->z)
			{
#pragma omp critical (output)
				{	cout << "\n# FailureSurface::ShouldInitiateFailure: bad 3D normals" << endl;
					cout << "# Normals = (" << normal->x << "," << normal->y << "," << normal->z << ")" << endl;
					cout << "# From R = " << R << ", Mode = " << failureMode << endl;
				}
			}
#endif
			return failureMode;
		}
	}
	else
	{	// get principle stress
		Vector lam = TensorEigenvalues(str,true,true);
		double sigma1 = fmax(lam.x,lam.y);
		double maxShear = fabs(0.5*(lam.x-lam.y));
		
		// check maximum principle stress and maximum shear stress
		if(sigma1>=sigmacRed*relStrength || maxShear>=taucRed*relStrength)
		{	// find angle to maximum principle stress
			double ssum = 0.5*(str->xx+str->yy);
			double sdif = 0.5*(str->xx-str->yy);
			double twoTheta;
			if(DbleEqual(str->xx,str->yy))
				twoTheta = PI_CONSTANT/2.;
			else
				twoTheta = atan(str->xy/sdif);
			
			// the correct angle is theta or theta+90
			double cs = cos(twoTheta);
			double sn = sin(twoTheta);
			double sigmaxp = ssum + sdif*cs + str->xy*sn;
			double sigmayp = ssum - sdif*cs - str->xy*sn;
			double theta = fabs(sigmaxp-sigma1)<fabs(sigmayp-sigma1) ?
			0.5*twoTheta : 0.5*(twoTheta+PI_CONSTANT);
			
			// if shear failure add -45 ccw rotation about principle z axis
			if(sigma1 < sigmacRed*relStrength)
			{	failureMode = SHEAR_FAILURE;
				theta -= PI_CONSTANT/4.;
			}
			else
				failureMode = TENSILE_FAILURE;
			
			// found normal (theta is ccw rotation angle from x axis to crack normal)
			normal->x = cos(theta);
			normal->y = sin(theta);
			normal->z = 0.;
			return failureMode;
		}
	}
	
	// has not failed yet
	return failureMode;
}

#pragma mark FailureSurface::Accessors

// reduced mode I stength
double FailureSurface::sigmaI(void) const { return sigmacRed; }

// reduced mode II stength
double FailureSurface::sigmaII(void) const { return taucRed; }

// initiation law name
const char *FailureSurface::GetInitiationLawName(void) const { return "Isotropic principal stresses"; }

