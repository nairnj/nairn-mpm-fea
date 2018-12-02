/********************************************************************************
    NewMaterial.hpp
    nairn-mpm-fea
    
    Created by John Nairn, July 13, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef NEWMATERIAL

#define NEWMATERIAL 100

#include "Materials/MaterialBase.hpp"

class NewMaterial : public MaterialBase
{
    public:
		// unique properties
		double newproperty;
        
        // constructors and destructors
		NewMaterial(char *,int);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
	
		// contitutive law methods
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
				
		// accessors
		virtual const char *MaterialType(void) const;
		virtual double WaveSpeed(bool,MPMBase *) const;
		
    protected:

};

#endif

