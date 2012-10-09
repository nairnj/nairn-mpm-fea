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
    method will return J = det F = lam1*lam2*lam3
   I think when done correctly, they give the same result when converted to
    Cauchy Stress. Default is to define the constant because it is slightly
    more efficient
   Logic: GIMP development implies Cauchy stress. In finite (or uniform) GIMP
    used here, stress is multiplied by Vp0 to get shape function integral
    normalized by Vp0. In MPM this is converted to mp/rho0 to get mass-weighted
    integrals. Thus we want mass-weighted average of (Cauchy Stress)/rho0.
    Note that this might change if uniform GIMP was changed to another method.
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
        Tensor GetLeftCauchyTensor2D(MPMBase *mptr,double,double,double,double,bool);
		Tensor GetLeftCauchyTensor3D(double F[][3]);
        Tensor GetLeftCauchyTensor3D(MPMBase *mptr,double,double,double,double,double,double,double,double,double,bool);
		double GetResidualStretch(MPMBase *);
        void ConvertToNominalStress2D(MPMBase *,double F[][3]);
    
        // Accessors
#ifndef CONSTANT_RHO
        virtual double GetCurrentRelativeVolume(MPMBase *);
#endif
};

#endif

