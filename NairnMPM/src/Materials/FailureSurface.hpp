/********************************************************************************
	FailureSurface.hpp
	nairn-mpm-fea

	Created by John Nairn, June 26, 2015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
 
	Base class for all failure surfaces used in softening materials

	Dependencies
		none
********************************************************************************/

#ifndef _FailureSurface_
#define _FailureSurface_

#define MAXPRINCIPALSTRESS 1

// TENSILE_FAILURE is max principle stress for isotropic or max normal
//    stress in transverse direction for TransIso where
// SHEAR_FAILURE is max shear stress or isotropic or max shear in transverse
//	  direcetion for TransIso where
// EA_FAILURE is tensile failure in axial direction for TransIso
// GA_FAILURE is axial shear failure for TransIso
enum { NO_FAILURE=0,TENSILE_FAILURE,SHEAR_FAILURE,EA_FAILURE,GA_FAILURE };

class MaterialBase;

class FailureSurface
{
	public:
		// constructors and destructors
		FailureSurface();
		FailureSurface(MaterialBase *);
 	
		// initialize
		virtual char *InputInitationProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintInitiationProperties(void) const;
	
		// methods
		virtual int ShouldInitiateFailure(Tensor *,Vector *,int,double) const;
	
		// accessors
		virtual double sigmaI(void) const;
		virtual double sigmaII(void) const;
		virtual const char *GetInitiationLawName(void) const;

	protected:
		double criticalNormal,sigmacRed;
		double criticalShear,taucRed;
		MaterialBase *parent;

};

#endif
