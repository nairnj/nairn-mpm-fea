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

// TENSILE_FAILURE is max principle stress
// SHEAR_FAILURE is max shear stress
// (note these must be 51 to 149
#define NO_FAILURE 0
#define TENSILE_FAILURE 90
#define SHEAR_FAILURE 110

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
