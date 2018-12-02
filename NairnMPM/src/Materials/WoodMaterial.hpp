/********************************************************************************
	WoodMaterial.hpp
	nairn-mpm-fea
 
	Created by John Nairn on 6/4/10.
	Copyright 2010 Oregon State University. All rights reserved.
 
	Dependencies
		HillPlastic.hpp (AnisoPlasticity.hpp Orthotropic.hpp TranIsotropic.hpp
			Elastic.hpp MaterialBase.hpp)
********************************************************************************/

#ifndef WOODMATERIAL

#define WOODMATERIAL 19

#include "Materials/HillPlastic.hpp"

class WoodMaterial : public HillPlastic
{
	public:
		// unique properties
		
		// constructors and destructors
		WoodMaterial(char *,int);
		
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
		
		// methods
		void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
	
		// accessors
		virtual const char *MaterialType(void) const;
		
	protected:
	double tempC1,tempC2;
	
};

#endif


