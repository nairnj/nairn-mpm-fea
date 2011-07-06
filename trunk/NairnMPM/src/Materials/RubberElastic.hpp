/********************************************************************************
    RubberElastic.hpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
 
    Isotropic, hyperelastic materials with subclasses to implement various
	theories

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef _RUBBERELASTIC_

#define _RUBBERELASTIC_

#include "Materials/MaterialBase.hpp"

class RubberElastic : public MaterialBase
{
    public:
		double aI;				// thermal expansion isotropic
		// double beta;			// moisture expansion isotopic (in base material)
        
        // constructors and destructors
        RubberElastic();
        RubberElastic(char *);
        
};

#endif

