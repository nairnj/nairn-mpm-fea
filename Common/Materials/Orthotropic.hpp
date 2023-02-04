/********************************************************************************
    Orthotropic.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Jan 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.

	Dependencies
		TranIsotropic.hpp (Elastic.hpp MaterialBase.hpp)
********************************************************************************/

#ifndef ORTHO

#define ORTHO 4

#include "Materials/TransIsotropic.hpp"

// for property flags, see TransIsotropic.hpp

class Orthotropic : public TransIsotropic
{
    public:
        
        // constructors and destructors
        Orthotropic(char *,int);
        
        // initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintMechanicalProperties(void) const;
#ifdef MPM_CODE
		virtual void PrintTransportProperties(void) const;
#endif
		
		// accessors
		virtual const char *MaterialType(void) const;
#ifdef MPM_CODE
        virtual double WaveSpeed(bool,MPMBase *) const;
		virtual double GetDiffZ(void) const;
		virtual double GetKcondZ(void) const;
#endif

	protected:
		double Ex,Ey,Ez,Gxy,Gyz,Gxz,ax,ay,az,betax,betay,betaz;
		double nuxy,nuyx,nuxz,nuzx,nuzy,nuyz;
	
#ifdef MPM_CODE
		double Dz,kCondz;
#ifdef POROELASTICITY
		double Darcyz,alphazPE;
#endif
#endif
	
};

#endif
