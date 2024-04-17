/********************************************************************************
	TIFailureSurface.cpp
	nairn-mpm-fea

	Created by John Nairn, 4 Aug 2016.
	Copyright (c) 2016 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/TIFailureSurface.hpp"
#include "Materials/MaterialBase.hpp"
#include "Materials/TransIsotropic.hpp"
#include "System/UnitsController.hpp"

#pragma mark TIFailureSurface::Constructors and Destructors

// Constructors
TIFailureSurface::TIFailureSurface() {}

TIFailureSurface::TIFailureSurface(MaterialBase *pair) : FailureSurface(pair)
{
	// criticalNormal and criticalShear in superclass are used to model
	// damage in the isotropic plane and crack normal will be perpendicular
	// to axial direciton of the material (parent is set in superclass too)
	// shear is like rolling shear in wood
	
	// These new values are for fiber breakage with crack normal in
	// axial direction, and two types of shear (like axial and transverse shear in wood)
	criticalAxialNormal = 1.e40;
	criticalAxialShear = 1.e40;
	criticalTransverseShear = 1.e40;
}

#pragma mark TIFailureSurface::Initialize

// Read hardening law properties
char *TIFailureSurface::InputInitationProperty(char *xName,int &input,double &gScaling)
{
    // axial strength
    if(strcmp(xName,"sigmacA")==0)
    {   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalAxialNormal,gScaling,1.e6);
    }
    
	// axial shear strength
    else if(strcmp(xName,"taucA")==0)
    {   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalAxialShear,gScaling,1.e6);
    }
	
	// transverse shear strength
	else if(strcmp(xName,"taucT")==0)
	{   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&criticalTransverseShear,gScaling,1.e6);
	}
	
	// is not a hardening law property
	return FailureSurface::InputInitationProperty(xName,input,gScaling);
}

// Convert properties to reduces stresses
const char *TIFailureSurface::VerifyAndLoadProperties(int np)
{
	// reduced failure stresses
	double rho = parent->GetRho(NULL);
	
	// failure in isotropic plane uses super class variables
    sigmacRed = criticalNormal/rho;
    taucRed = criticalShear/rho;
	
	// these are for failure in axial direciton
	sigmacARed = criticalAxialNormal/rho;
	taucARed = criticalAxialShear/rho;
	taucTRed = criticalTransverseShear/rho;
	
	// this initiation law does not support pressure dependences
	tauhRed = -1.;
	
	return NULL;
}

// Remap properties of anisotropic failure surfaces
void TIFailureSurface::RemapProperties(int swap) {}

// print just initiation properties to output window
void TIFailureSurface::PrintInitiationProperties(void) const
{
    cout << GetInitiationLawName() << endl;
	MaterialBase::PrintProperty("sigcA",criticalAxialNormal*UnitsController::Scaling(1.e-6),"");
    MaterialBase::PrintProperty("sigcT",criticalNormal*UnitsController::Scaling(1.e-6),"");
	cout << endl;
	MaterialBase::PrintProperty("taucA",criticalAxialShear*UnitsController::Scaling(1.e-6),"");
    MaterialBase::PrintProperty("taucT",criticalTransverseShear*UnitsController::Scaling(1.e-6),"");
	MaterialBase::PrintProperty("taucRS",criticalShear*UnitsController::Scaling(1.e-6),"");
    cout << endl;
}

#pragma mark TIFailureSurface::Methods

// Note that this failure surface only used for transversely isotropic materials.
// On call, the stress is stress in the material axis system, but axial direction may be along Y or Z axes
//		for 2D it will have only x-y plane stresses (and not zz stress)
// For 3D, axial direction will always be along z direction
// If failed, find the normal vector in the same material axis system, but for
//     3D, find Euler angles in normal, change relStrength to
//          (relStrength)/(failure_stress)
//     for scaling when using smoothed stresses, and return failure mode
// If not failed return 0 (NO_FAILURE)
int TIFailureSurface::ShouldInitiateFailure(Tensor *str,Vector *normal,int np,double &relStrength,GenADaMVariables *alpha) const
{
	int failureMode = NO_FAILURE;
	
	if(np==THREED_MPM)
	{	// in 3D, material axial direction always in the z direction
		// Two possible normals: n.y = pi/2 (softens axial modulus and GA using min(taucT,taucA)
		//    n.x = theta (softens transverse modulus, GA using min(taucT,taucA), and GT
		
		// check axial direction tensile file (e.g., fiber breakage)
		if(str->zz>sigmacARed*relStrength)
		{	// ZYZ rotation = (0,90,0) (crack x will be in original axial (or z) direction)
			// Use Dx
			normal->x = 0.;
			normal->y = PI_CONSTANT/2.;   // swaps z and x
			normal->z = 0.;
#ifdef SMOOTHED_STRESS
			relStrength /= str->zz;
#endif
			return EA_FAILURE;
		}
		
		// check axial shear (shear failure across fiber with same normal as tensile failure)
		double maxAxialShear = sqrt(str->xz*str->xz+str->yz*str->yz);
		if(maxAxialShear > taucARed*relStrength || maxAxialShear > taucTRed*relStrength)
		{	// ZYZ rotation = (0,90,0) (crack x will be axial direction)
			if(taucTRed < taucARed)
			{	// fiber direction is weaker (crack x will in original axial (or z) direction)
				// e.g., wood if transverse shear was lower than axial shear (it is not)
				// Use Dx
				normal->x = 0.;
				normal->y = PI_CONSTANT/2.;		// swaps z and x
				normal->z = 0.;
			}
			else
			{	// transverse direction is weaker
				// Material n = (txz,tyz,0)/maxAxialShear
				// normal->x is ccw rotation angle from x axis to crack normal
				// Use Dz
				normal->x = asin(str->xz/maxAxialShear);
				normal->y = 0.;
				normal->z = 0.;
			}
#ifdef SMOOTHED_STRESS
			relStrength /= maxAxialShear;
#endif
			return GA_FAILURE;
		}
		
		// get transverse plane max principle and max shear
		double ssum = 0.5*(str->xx+str->yy);
		double sdif = 0.5*(str->xx-str->yy);
		double maxShear = sqrt(sdif*sdif + str->xy*str->xy);
		double sigma1 = ssum + maxShear;
		
		if(sigma1>=sigmacRed*relStrength || maxShear>=taucRed*relStrength)
		{	// find angle to maximum principle stress
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
#ifdef SMOOTHED_STRESS
				relStrength /= maxShear;
#endif
				theta -= PI_CONSTANT/4.;
			}
			else
			{	failureMode = TENSILE_FAILURE;
#ifdef SMOOTHED_STRESS
				relStrength /= sigma1;
#endif
			}

			// ZYZ rotation (axial direction will reamain in z direction)
			// found normal (theta is ccw rotation angle from x axis to crack normal)
			// Use Dz
			normal->x = theta;
			normal->y = 0.;
			normal->z = 0.;
			return failureMode;
		}
	}
	else
	{	// get axial direction
		int axialAxis = ((TransIsotropic *)parent)->AxialDirection();
		
		if(axialAxis == AXIAL_Y)
		{	// Three types of failure when axial direction is along y axis
			if(str->yy > sigmacARed*relStrength)
			{	// axial failure by tension. Use Dx
				normal->x = 0.;
				normal->y = 1.;
				normal->z = 0.;
#ifdef SMOOTHED_STRESS
				relStrength /= str->yy;
#endif
				return EA_FAILURE;
			}
			else if(fabs(str->xy)>taucARed*relStrength || fabs(str->xy)>taucTRed*relStrength)
			{	// axial failure by shear
				if(taucTRed < taucARed)
				{	// transverse shear is lower
					// Use Dx
					normal->x = 0.;
					normal->y = 1.;
					normal->z = 0.;
				}
				else
				{	// axial shear is lower
					// Use Dy
					normal->x = 1.;
					normal->y = 0.;
					normal->z = 0.;
				}
#ifdef SMOOTHED_STRESS
				relStrength /= fabs(str->xy);
#endif
				return GA_FAILURE;
			}
			else if(str->xx > sigmacRed*relStrength)
			{	// transverse tensile failure
				// can't fail in transverse shear because tauxz=0
				// TIQUERY - can't fail by str(2,2) failure because could not model it
				// Use Dy
				normal->x = 1.;
				normal->y = 0.;
				normal->z = 0.;
#ifdef SMOOTHED_STRESS
				relStrength /= str->xx;
#endif
				return TENSILE_FAILURE;
			}
		}
		else
		{	// isotropic failure when axial direction along z axis
			// parent class handles it and returns normal in material axis system if failed
			// Use Dz if get TENSILE_FAILURE or SHEAR_FAILURE
			return FailureSurface::ShouldInitiateFailure(str,normal,np,relStrength,NULL);
		}
	}
	
	// if not failed, return false
	return failureMode;
}

#pragma mark TIFailureSurface::Accessors

// reduced axial tensile strength
double TIFailureSurface::sigmaA(void) const { return sigmacARed; }

// reduced axial shear stength
double TIFailureSurface::tauA(void) const { return taucARed; }

// reduced transverse shear stength
double TIFailureSurface::tauT(void) const { return taucTRed; }

// normally only need minimum of these two
double TIFailureSurface::tauMin(void) const { return fmin(taucARed,taucTRed); }

// initiation law name
const char *TIFailureSurface::GetInitiationLawName(void) const { return "Transversely isotropic failure"; }

