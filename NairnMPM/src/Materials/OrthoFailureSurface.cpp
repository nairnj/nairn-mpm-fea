/********************************************************************************
	OrthoFailureSurface.cpp
	nairn-mpm-fea

	Created by John Nairn, 2/17/20.
	Copyright (c) 2020 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/OrthoFailureSurface.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark OrthoFailureSurface::Constructors and Destructors

// Constructors
OrthoFailureSurface::OrthoFailureSurface() {}

OrthoFailureSurface::OrthoFailureSurface(MaterialBase *pair) : TIFailureSurface(pair)
{
	// criticalNormal, criticalShear, criticalAxialNormal, criticalAxialShear,
 	// and criticalTransverseShear in superclass are for 5 strength properties
	
	// These new values are needed for remaining properties
	criticalZZNormal = 1.e40;
	criticalYZ_ZShear = 1.e40;
	criticalXZ_XShear = 1.e40;
	criticalXZ_ZShear = 1.e40;
}

#pragma mark OrthoFailureSurface::Initialize

// Read hardening law properties
char *OrthoFailureSurface::InputInitationProperty(char *xName,int &input,double &gScaling)
{
	// axial xx strength
	if(strcmp(xName,"sigmaXXc")==0)
	{   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalAxialNormal,gScaling,1.e6);
	}
	
	// axial yy strength
	else if(strcmp(xName,"sigmaYYc")==0)
	{   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalNormal,gScaling,1.e6);
	}
	
	// axial yy strength
	else if(strcmp(xName,"sigmaZZc")==0)
	{   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalZZNormal,gScaling,1.e6);
	}
	
	// axial shear strength
	else if(strcmp(xName,"tauXY-Xc")==0)
	{   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalAxialShear,gScaling,1.e6);
	}
	
	// transverse shear strength
	else if(strcmp(xName,"tauXY-Yc")==0)
	{   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalTransverseShear,gScaling,1.e6);
	}
	
	// axial shear strength
	else if(strcmp(xName,"tauXZ-Xc")==0)
	{   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalXZ_XShear,gScaling,1.e6);
	}
	
	// transverse shear strength
	else if(strcmp(xName,"tauXZ-Zc")==0)
	{   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalXZ_ZShear,gScaling,1.e6);
	}

	// axial shear strength
	else if(strcmp(xName,"tauYZ-Yc")==0)
	{   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalShear,gScaling,1.e6);
	}
	
	// transverse shear strength
	else if(strcmp(xName,"tauYZ-Zc")==0)
	{   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalYZ_ZShear,gScaling,1.e6);
	}
	
	// is not a hardening law property
	return NULL;
}

// Convert properties to reduces stresses
const char *OrthoFailureSurface::VerifyAndLoadProperties(int np)
{
	// reduced failure stresses
	double rho = parent->GetRho(NULL);
	
	// normal direction xx, yy, zz
	sigmacARed = criticalAxialNormal/rho;
	sigmacRed = criticalNormal/rho;
	sigmaZZRed = criticalZZNormal/rho;
	
	// tau_{xy} values
	taucARed = criticalAxialShear/rho;
	taucTRed = criticalTransverseShear/rho;
	
	// tau_{yz} values
	taucRed = criticalShear/rho;
	taucYZ_ZRed = criticalYZ_ZShear/rho;

	// tau_{xz} values
	taucXZ_XRed = criticalXZ_ZShear/rho;
	taucXZ_ZRed = criticalXZ_ZShear/rho;
	
	return NULL;
}

// Swap y an z axes
void OrthoFailureSurface::RemapProperties(int swap)
{
    if(swap==1)
    {   MaterialBase::pswap(criticalAxialNormal,criticalZZNormal);     // XX with ZZ
        MaterialBase::pswap(criticalXZ_XShear,criticalXZ_ZShear);      // XZX with ZXZ/XZZ
        MaterialBase::pswap(criticalTransverseShear,criticalShear);    // XYY with ZYY/YZY
        MaterialBase::pswap(criticalAxialShear,criticalYZ_ZShear);     // XYX with ZYZ/YZZ
    }
    else if(swap>1)
    {   MaterialBase::pswap(criticalNormal,criticalZZNormal);           // YY with ZZ
        MaterialBase::pswap(criticalShear,criticalYZ_ZShear);           // YZY with ZYZ/YZZ
        MaterialBase::pswap(criticalTransverseShear,criticalXZ_ZShear); // XYY with XZZ
        MaterialBase::pswap(criticalAxialShear,criticalXZ_XShear);      // XYX with XZX
    }
}

// print just initiation properties to output window
void OrthoFailureSurface::PrintInitiationProperties(void) const
{
	cout << GetInitiationLawName() << endl;
	MaterialBase::PrintProperty("sigXXc",criticalAxialNormal*UnitsController::Scaling(1.e-6),"");
	MaterialBase::PrintProperty("tauYZ-Yc",criticalShear*UnitsController::Scaling(1.e-6),"");
	MaterialBase::PrintProperty("tauYZ-Zc",criticalYZ_ZShear*UnitsController::Scaling(1.e-6),"");
	cout << endl;
	MaterialBase::PrintProperty("sigYYc",criticalNormal*UnitsController::Scaling(1.e-6),"");
	MaterialBase::PrintProperty("tauXZ-Xc",criticalXZ_XShear*UnitsController::Scaling(1.e-6),"");
	MaterialBase::PrintProperty("tauXZ-Zc",criticalXZ_ZShear*UnitsController::Scaling(1.e-6),"");
	cout << endl;
	MaterialBase::PrintProperty("sigZZc",criticalZZNormal*UnitsController::Scaling(1.e-6),"");
	MaterialBase::PrintProperty("tauXY-Xc",criticalAxialShear*UnitsController::Scaling(1.e-6),"");
	MaterialBase::PrintProperty("tauXY-Yc",criticalTransverseShear*UnitsController::Scaling(1.e-6),"");
	cout << endl;
}

#pragma mark OrthoFailureSurface::Methods

// Note that this failure surface only used for transversely isotropic materials.
// On call, the stress is stress in the material axis system, but axial direction may be along Y or Z axes
//		for 2D it will have only x-y plane stresses (and not zz stress)
// For 3D, axial direction will always be along z direction
// If failed, find the normal vector in the same material axis system, but for
//     3D, find Euler angles in normal, change relStrength to
//          (relStrength)/(failure_stress)
//     for scaling when using smoothed stresses, and return failure mode
// If not failed return 0 (NO_FAILURE)
int OrthoFailureSurface::ShouldInitiateFailure(Tensor *str,Vector *normal,int np,
                                               double &relStrength,GenADaMVariables *alpha) const
{
	int failureMode = NO_FAILURE;
	ZeroVector(normal);

	if(np==THREED_MPM)
	{	// 3D: get ZYZ rotation in normal
		
		// sigma_{xx}
		if(str->xx > sigmaXX()*relStrength)
		{	// axial failure by tension in x direction
#ifdef SMOOTHED_STRESS
			relStrength /= str->xx;
#endif
			return EA_FAILURE;
		}

		// sigma_{yy}
		if(str->yy > sigmaXX()*relStrength)
		{	// axial failure by tension in y direction
			normal->x = PI_CONSTANT/2.;   // swaps x and y
#ifdef SMOOTHED_STRESS
			relStrength /= str->xx;
#endif
			return TENSILE_FAILURE;
		}
		
		// sigma_{zz}
		if(str->zz > sigmaZZ()*relStrength)
		{	// axial failure by tension in x direction
			normal->y = PI_CONSTANT/2.;	  // swaps x and z
#ifdef SMOOTHED_STRESS
			relStrength /= str->xx;
#endif
			return EZZ_FAILURE;
		}
		
		// tau_{xy}
		if(fabs(str->xy)>tauXYX()*relStrength || fabs(str->xy)>tauXYY()*relStrength)
		{	if(tauXYX() < tauXYY())
			{	// fails by XYX with normal in y direction, otherwise normal in x direction
				normal->x = PI_CONSTANT/2.;   // swaps x and y
			}
#ifdef SMOOTHED_STRESS
			relStrength /= fabs(str->xy);
#endif
			return GA_FAILURE;
		}
		
		// tau_{xz}
		if(fabs(str->xz)>tauXZX()*relStrength || fabs(str->xz)>tauXZZ()*relStrength)
		{	// axial failure by shear
			if(tauXZX() < tauXZZ())
			{	// fails by XYX with normal in z direction, otherwise normal in x direction
				normal->y = PI_CONSTANT/2.;	  // swaps z and x
			}
#ifdef SMOOTHED_STRESS
			relStrength /= fabs(str->xz);
#endif
			return GXZ_FAILURE;
		}
		
		// tau_{yz}
		if(fabs(str->yz)>tauYZY()*relStrength || fabs(str->yz)>tauYZZ()*relStrength)
		{	// axial failure by shear
			if(tauYZY() <= tauYZZ())
			{	// y crack is lower, normal in z direction
				normal->y = PI_CONSTANT/2.;	  // swaps z and x
			}
			else
			{	// z crack is lower, normal in y direction
				normal->x = PI_CONSTANT/2.;   // swaps x and y
			}
#ifdef SMOOTHED_STRESS
			relStrength /= fabs(str->yz);
#endif
			return SHEAR_FAILURE;
		}
	}
	
	else
	{	// 2D: get cos(theta),sin(theta) in normal
		normal->z = 0.;
		
		// sigmna_{xx}
		if(str->xx > sigmaXX()*relStrength)
		{	// axial failure by tension in x direction
			normal->x = 1.;
	#ifdef SMOOTHED_STRESS
			relStrength /= str->xx;
	#endif
			return EA_FAILURE;
		}
	
		// sigmna_{yy}
		if(str->yy > sigmaYY()*relStrength)
		{	// axial failure by tension in y direction
			normal->y = 1.;			// swaps y and x
	#ifdef SMOOTHED_STRESS
			relStrength /= str->xx;
	#endif
			return TENSILE_FAILURE;
		}
		
		// tau_{xy}
		if(fabs(str->xy)>tauXYX()*relStrength || fabs(str->xy)>tauXYY()*relStrength)
		{	// axial failure by shear
			if(tauXYX() < tauXYY())
			{	// x crack is lower, normal in y direction
				normal->y = 1.;
			}
			else
			{	// y crack is lower, normal in x direction
				normal->x = 1.;
			}
#ifdef SMOOTHED_STRESS
			relStrength /= fabs(str->xy);
#endif
			return GA_FAILURE;
		}
	}
	
	// if not failed, return false
	return failureMode;
}

#pragma mark OrthoFailureSurface::Accessors

// reduced yy tensile strength
double OrthoFailureSurface::sigmaXX(void) const { return sigmacARed; }
double OrthoFailureSurface::sigmaYY(void) const { return sigmacRed; }
double OrthoFailureSurface::sigmaZZ(void) const { return sigmaZZRed; }

// reduced shear strength - get minimum or maximum of the two
double OrthoFailureSurface::tauXYX(void) const { return taucARed; }
double OrthoFailureSurface::tauXYY(void) const { return taucTRed; }
double OrthoFailureSurface::tauXZX(void) const { return taucXZ_XRed; }
double OrthoFailureSurface::tauXZZ(void) const { return taucXZ_ZRed; }
double OrthoFailureSurface::tauYZY(void) const { return taucRed; }
double OrthoFailureSurface::tauYZZ(void) const { return taucYZ_ZRed; }

// initiation law name
const char *OrthoFailureSurface::GetInitiationLawName(void) const { return "Orthotropic failure"; }

