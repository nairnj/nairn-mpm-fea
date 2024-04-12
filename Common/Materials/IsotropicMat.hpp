/********************************************************************************
    IsotropicMat.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 31 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.

	Dependencies
		Elastic.hpp (MaterialBase.hpp)
********************************************************************************/

#ifndef ISOTROPIC

#define ISOTROPIC 1

#include "Materials/Elastic.hpp"

enum {E_PROP=0,G_PROP,NU_PROP,ISO_PROPS};

class IsotropicMat : public Elastic
{
    public:
        
        // constructors and destructors
        IsotropicMat(char *,int);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
        
#ifdef MPM_CODE
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
        virtual Vector ConvertJToK(Vector,Vector,Vector,int);
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
		virtual void LRConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
		virtual void SRConstitutiveLaw2D(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
		virtual void SRConstitutiveLaw3D(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
#endif
		
		// accessors
		virtual const char *MaterialType(void) const;
#ifdef MPM_CODE
		virtual double GetMagnitudeSFromDev(Tensor *,int) const;
		virtual double GetMagnitudeSFromTotal(Tensor *,int) const;
		virtual double WaveSpeed(bool,MPMBase *) const;
        virtual double ShearWaveSpeed(bool,MPMBase *,int) const;
#endif

    protected:
		double E,nu,G,aI,gamma0;
        char read[ISO_PROPS];
};

#endif

