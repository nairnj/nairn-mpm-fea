/********************************************************************************
    IsotropicMat.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 31 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.

	Dependencies
		Elastic.hpp
********************************************************************************/

#ifndef ISOTROPIC

#define ISOTROPIC 1

#include "Materials/Elastic.hpp"

enum {E_PROP=0,G_PROP,NU_PROP,ISO_PROPS};

class IsotropicMat : public Elastic
{
    public:
        double E,nu,G,aI;
        
        // constructors and destructors
        IsotropicMat();
        IsotropicMat(char *);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
        
		virtual void FillElasticProperties(ElasticProperties *,int);
#ifdef MPM_CODE
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *) const;
        virtual Vector ConvertJToK(Vector,Vector,Vector,int);
#endif
		
		// accessors
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
#ifdef MPM_CODE
		virtual double WaveSpeed(bool,MPMBase *) const;
        virtual double ShearWaveSpeed(bool,MPMBase *) const;
#endif

    protected:
        char read[ISO_PROPS];
};

#endif

