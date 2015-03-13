/********************************************************************************
    DugdaleMat.hpp
    NairnMPM
    
    Created by John Nairn on Wed Sep 04 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.

	Dependencies
		VonMisesHardening.hpp
********************************************************************************/

#ifndef _DUGDALEMAT_

#define _DUGDALEMAT_

#define DUGDALE 6

#include "Materials/VonMisesHardening.hpp"

class DugdaleMat : public VonMisesHardening
{
    public:
        // constructors and destructors
        DugdaleMat();
        DugdaleMat(char *);
        
        // methods
		virtual double GetF(MPMBase *,Tensor *,int);
		virtual void GetDfDsigma(MPMBase *,Tensor *,int);
		virtual char *MaterialType(void);
};

#endif
