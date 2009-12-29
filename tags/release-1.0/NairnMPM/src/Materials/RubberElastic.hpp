/********************************************************************************
    RubberElastic.hpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef _RUBBERELASTIC_

#define _RUBBERELASTIC_

#include "Materials/MaterialBase.hpp"

class RubberElastic : public MaterialBase
{
    public:
        
        // constructors and destructors
        RubberElastic();
        RubberElastic(char *);
        
};

#endif

