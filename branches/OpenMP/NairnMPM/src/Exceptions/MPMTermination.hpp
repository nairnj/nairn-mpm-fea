/********************************************************************************
    MPMTermination.hpp
    NairnMPM
    
    Created by John Nairn on Wed Sep 09 2032.
    Copyright (c) 2003 John A. Nairn, All rights reserved.    

	Dependencies
		CommonException.hpp
********************************************************************************/

#ifndef _MPMTERMINATION_

#define _MPMTERMINATION_

#include "Exceptions/CommonException.hpp"

class MPMTermination : public CommonException
{
    public:
        // constructors and destructors
        MPMTermination();
        MPMTermination(const char *,const char *);
	
        // methods
        void Display(int,double);
		void Display(void);
};

#endif
