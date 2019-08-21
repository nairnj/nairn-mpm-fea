/********************************************************************************
	MatPointAS.hpp
	nairn-mpm-fea

	Created by John Nairn on 10/23/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.

	Dependencies
		MatPoint2D.hpp, MPMBase.hpp
 ********************************************************************************/

#ifndef _MATPOINTAS_

#define _MATPOINTAS_

#include "MPM_Classes/MatPoint2D.hpp"

class MaterialBase;

class MatPointAS : public MatPoint2D
{
	public:
	
		// constructors and destructors
		MatPointAS();
		MatPointAS(int,int,double,double);
	
		// methods
        virtual void UpdateStrain(double,int,int,void *,int);
		virtual void PerformConstitutiveLaw(Matrix3,double,int,void *,ResidualStrains *,Tensor *);
		virtual void GetFintPlusFext(Vector *,double,double,double,double);
        virtual void SetOrigin(Vector *);
		virtual double GetVolume(int);
        virtual double GetUnscaledVolume(void);
        virtual void GetCPDINodesAndWeights(int);
};

#endif
