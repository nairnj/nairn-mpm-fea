/********************************************************************************
    IsotropicMat.hpp
    NairnMPM
    
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
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyProperties(int);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
        
#ifdef MPM_CODE
        // convert J-integrals into stress intensity factors (YJG)
        virtual Vector ConvertJToK(Vector,Vector,Vector,int);
#endif
		
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
#ifdef MPM_CODE
		virtual void InitialLoadTransProps(void);
		virtual double WaveSpeed(bool,MPMBase *);
        virtual double ShearWaveSpeed(bool,MPMBase *);
#endif

    protected:
        char read[ISO_PROPS];
};

#endif

