/********************************************************************************
    HyperElastic.hpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
 
    Isotropic, hyperelastic materials with subclasses to implement various
	theories

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef _HYPERELASTIC_

#define _HYPERELASTIC_

/* When this is defined, subclasses should find (Cauchy Stress)/rho0, i.e.,
    ignore changes in density due to large deformation. The base material
    class will always return 1 for current relative volume.
   When it is not defined, subclasses should find (Kirchoff Stress)/rho0, i.e.,
    find J*(Cauchy Stress)/rho0 = (Cauchy Stress)/rho. The relative volume
    method should return J = det F = lam1*lam2*lam3
   I think the appropriate choice is to not define it and get Kirchoff stress
	but, this assumes GIMP shape function is divided by actual Vp and not
	initial Vp0. Hence when using CPDI, the "not defined" choice seems right.
    When using uGIMP, it might be better to define it, but then CPDI would
    be wrong and uGIMP won't work to high deformation anyway. A possible
    improvement for uGIMP would to divide shape function by J, hence
	dividing by J Vp0 = Vp.
*/
//#define CONSTANT_RHO

#include "Materials/MaterialBase.hpp"

class HyperElastic : public MaterialBase
{
    public:
		double aI;				// thermal expansion isotropic
		// double beta;			// moisture expansion isotopic (in base material)
        
        // constructors and destructors
        HyperElastic();
        HyperElastic(char *);
        
		// Methods (make virtual if any subclass needs them)
		double GetDeformationGrad(double F[][3],MPMBase *,double,double,double,double,bool,bool);
		double GetDeformationGrad(double F[][3],MPMBase *,double,double,double,double,double,double,double,double,double,bool,bool);
		Tensor GetLeftCauchyTensor2D(double F[][3]);
		Tensor GetLeftCauchyTensor3D(double F[][3]);
		double GetResidualStretch(MPMBase *);
        void ConvertToNominalStress2D(MPMBase *,double F[][3]);
    
        // Accessors
#ifndef CONSTANT_RHO
        virtual double GetCurrentRelativeVolume(MPMBase *);
#endif
};

#endif

