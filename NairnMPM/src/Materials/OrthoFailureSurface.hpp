/********************************************************************************
	OrthoFailureSurface.hpp
	nairn-mpm-fea

	Created by John Nairn, 2/17/20.
	Copyright (c) 2020 John A. Nairn, All rights reserved.

	Dependencies
	TIFailureSurface,FailureSurface.hpp
********************************************************************************/

#ifndef _OrthoFailureSurface_
#define _OrthoFailureSurface_

#define ORTHOFAILURESURFACE 3

#include "Materials/TIFailureSurface.hpp"

class MaterialBase;

// new failure modes (note TransIsoSoftening might change them)
#define EZZ_FAILURE 85
#define GXZ_FAILURE 120

class OrthoFailureSurface : public TIFailureSurface
{
	public:
		// constructors and destructors
		OrthoFailureSurface();
		OrthoFailureSurface(MaterialBase *);
	
		// initialize
		virtual char *InputInitationProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
        virtual void RemapProperties(int);
		virtual void PrintInitiationProperties(void) const;
	
		// methods
		virtual int ShouldInitiateFailure(Tensor *,Vector *,int,double &,GenADaMVariables *) const;
	
		// accessors
		virtual double sigmaXX(void) const;
		virtual double sigmaYY(void) const;
		virtual double sigmaZZ(void) const;
		virtual double tauXYX(void) const;
		virtual double tauXYY(void) const;
		virtual double tauXZX(void) const;
		virtual double tauXZZ(void) const;
		virtual double tauYZY(void) const;
		virtual double tauYZZ(void) const;
		virtual const char *GetInitiationLawName(void) const;
	
	protected:
		// criticalXXNormal in criticalAxialNormal and sigmacARed
		// criticalYYNormal in criticalNormal and sigmacRed
		double criticalZZNormal,sigmaZZRed;
		// criticalXY_XNormal in criticalAxialShear and taucARed
		// criticalXY_YNormal in criticalTransverseShear and taucTRed
		// criticalYZ_YNormal in criticalShear and taucRed
		double criticalYZ_ZShear,taucYZ_ZRed;
		double criticalXZ_XShear,taucXZ_XRed;
		double criticalXZ_ZShear,taucXZ_ZRed;
};

#endif
