/********************************************************************************
    Elastic.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef _ELASTIC_

#define _ELASTIC_

#include "Materials/MaterialBase.hpp"

#ifdef MPM_CODE
// The full stiffness matrix in C
// alpha and beta are thermal and moisture expansion coefficients
// although some elements are used for other things
typedef struct {
	double C[6][6];
	double alpha[8];
	double beta[6];
} ElasticProperties;
#endif

class Elastic : public MaterialBase
{
    public:
        
        // constructors and destructors
        Elastic();
        Elastic(char *);
		
		// initialize
        
		// methods
#ifdef MPM_CODE
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int,void *,ResidualStrains *) const;
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int,void *,ResidualStrains *) const;
#else
        virtual double GetStressStrainZZ(double,double,double,double,double,int);
#endif

	protected:
		double prop1,prop2;
#ifdef MPM_CODE
		ElasticProperties pr;
		double Cadota;
#else
        double prop3;
#endif

		// methods
        virtual const char *SetAnalysisProps(int,double,double,double,double,double,
                        double,double,double,double,double,double,double,double,double,double);
};

#endif
