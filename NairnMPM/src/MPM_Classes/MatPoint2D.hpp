/********************************************************************************
    MatPoint2D.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		MPMBase.hpp
********************************************************************************/

#ifndef _MATPOINT2D_

#define _MATPOINT2D_

#include "MPM_Classes/MPMBase.hpp"

class MaterialBase;

class MatPoint2D : public MPMBase
{
    public:
        // constructors and destructors
        MatPoint2D();
        MatPoint2D(int,int,double,double);
        
        // methods
		virtual void SetOrigin(Vector *);
        virtual void SetPosition(Vector *);
        virtual void SetVelocity(Vector *);
        virtual void SetStress(Tensor *);
        virtual double thickness(void);
		virtual void UpdateStrain(double,int,int,void *,int);
		virtual void PerformConstitutiveLaw(Matrix3,double,int,void *,ResidualStrains *);
		virtual void GetFintPlusFext(Vector *,double,double,double,double);
		virtual void MovePosition(double,Vector *,Vector *,double);
		virtual void MoveVelocity(double,Vector *);
		virtual void MovePosition(double);
		virtual void SetVelocitySpeed(double);
		virtual void AddTemperatureGradient(int,Vector *);
		virtual double FCond(int,double,double,double,TransportProperties *);
		virtual void AddConcentrationGradient(Vector *);
		virtual double FDiff(double,double,double,TransportProperties *);
		virtual double KineticEnergy(void);
		virtual Matrix3 GetDeformationGradientMatrix(void) const;
		virtual void SetDeformationGradientMatrix(Matrix3);
		virtual Matrix3 GetDisplacementGradientMatrix(void) const;
        virtual Matrix3 GetElasticLeftCauchyMatrix(void);
        virtual void GetDeformationGradient(double F[][3]) const;
        virtual double GetRelativeVolume(void);
		virtual double GetVolume(int);
        virtual void GetSemiSideVectors(Vector *,Vector *,Vector *) const;
        virtual void GetUndeformedSemiSides(double *,double *,double *) const;
		virtual void GetCPDINodesAndWeights(int);
		virtual double GetTractionInfo(int,int,int *,Vector *,Vector *,int *);
		virtual Matrix3 GetInitialRotation(void);
		virtual Matrix3 GetElasticBiotStrain(void);
    
    protected:
        double thick;
};

#endif
