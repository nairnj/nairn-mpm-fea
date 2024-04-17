/********************************************************************************
	TIFailureSurface.hpp
	nairn-mpm-fea

	Created by John Nairn, 4 AUG 2016.
	Copyright (c) 2016 John A. Nairn, All rights reserved.

	Dependencies
		FailureSurface.hpp
********************************************************************************/

#ifndef _TIFailureSurface_
#define _TIFailureSurface_

#define TIFAILURESURFACE 2

#include "Materials/FailureSurface.hpp"

// new failure modes (note TransIsoSoftening might change them)
#define EA_FAILURE 95
#define GA_FAILURE 115

class MaterialBase;

class TIFailureSurface : public FailureSurface
{
	public:
		// constructors and destructors
		TIFailureSurface();
		TIFailureSurface(MaterialBase *);
 	
		// initialize
		virtual char *InputInitationProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
        virtual void RemapProperties(int);
		virtual void PrintInitiationProperties(void) const;
	
		// methods
		virtual int ShouldInitiateFailure(Tensor *,Vector *,int,double &,GenADaMVariables *) const;
	
		// accessors
		virtual double sigmaA(void) const;
		virtual double tauA(void) const;
		virtual double tauT(void) const;
		virtual double tauMin(void) const;
	
		virtual const char *GetInitiationLawName(void) const;

	protected:
		double criticalAxialNormal,sigmacARed;
		double criticalAxialShear,taucARed;
		double criticalTransverseShear,taucTRed;
};

#endif
