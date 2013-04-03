/********************************************************************************
    Orthotropic.hpp
    NairnMPM
    
    Created by John Nairn on Tues Jan 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.

	Dependencies
		TranIsotropic.hpp (MaterialBase.hpp)
********************************************************************************/

#ifndef ORTHO

#define ORTHO 4

#include "Materials/TransIsotropic.hpp"

// for property flags, see TransIsotropic.hpp

class Orthotropic : public TransIsotropic
{
    public:
        double Ex,Ey,Ez,Gxy,Gyz,Gxz,ax,ay,az,betax,betay,betaz;
        double nuxy,nuyx,nuxz,nuzx,nuzy,nuyz;
        
        // constructors and destructors
        Orthotropic();
        Orthotropic(char *);
        
        // initilize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintMechanicalProperties(void) const;
#ifdef MPM_CODE
		virtual void PrintTransportProperties(void) const;
#endif
		
		// accessors
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
#ifdef MPM_CODE
        virtual double WaveSpeed(bool,MPMBase *) const;
		virtual double GetDiffZ(void);
		virtual double GetKcondZ(void);
#endif

	protected:
#ifdef MPM_CODE
		double Dz,kcondz;
#endif
	
};

#endif
