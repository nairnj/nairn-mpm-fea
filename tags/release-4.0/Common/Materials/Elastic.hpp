/********************************************************************************
    Elastic.hpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef _ELASTIC_

#define _ELASTIC_

#include "Materials/MaterialBase.hpp"

class Elastic : public MaterialBase
{
    public:
        
        // constructors and destructors
        Elastic();
        Elastic(char *);
		
		// initialize
        
		// methods
#ifdef MPM_CODE
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int);
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
#else
        virtual double GetStressStrainZZ(double,double,double,double,double,int);
#endif

	protected:
		double prop1,prop2;
#ifdef MPM_CODE
        double mdm[6][6];
        double me0[8],mc0[6];
		int hasTransProps;			// flag set TRUE whenever transport props have been evaluated
#else
        double prop3;
#endif

		// methods
        virtual const char *SetAnalysisProps(int,double,double,double,double,double,
                        double,double,double,double,double,double,double,double,double,double);
};

#endif
