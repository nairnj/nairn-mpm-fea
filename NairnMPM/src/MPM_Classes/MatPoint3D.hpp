/********************************************************************************
    MatPoint3D.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 7/21/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
	
	Dependencies
		MPMBase.hpp
********************************************************************************/

#ifndef _MATPOINT3D_

#define _MATPOINT3D_

#include "MPM_Classes/MPMBase.hpp"

class MaterialBase;

class MatPoint3D : public MPMBase
{
    public:
        // constructors and destructors
        MatPoint3D();
        MatPoint3D(int,int,double);
        
        // methods
		virtual void SetOrigin(Vector *);
        virtual void SetPosition(Vector *);
        virtual void SetVelocity(Vector *);
        virtual double thickness(void);
		virtual void UpdateStrain(double,int,int,void *,int,bool);
        virtual Matrix3 ExtraVelocityGradient(void);
		virtual void PerformConstitutiveLaw(Matrix3,double,int,void *,ResidualStrains *,Tensor *);
		virtual void GetFintPlusFext(Vector *,double,double,double,double);
		virtual void MoveParticle(GridToParticleExtrap *);
		virtual void MovePosition(double);
		virtual void SetVelocitySpeed(double);
		virtual void AddTemperatureGradient(int,Vector *);
		virtual void AddConcentrationGradient(int,Vector *);
		virtual double FCond(int,double,double,double,TransportProperties *);
		virtual double FDiff(double,double,double,TransportProperties *,int);
		virtual double KineticEnergy(void);
		virtual Matrix3 GetDeformationGradientMatrix(void) const;
		virtual void SetDeformationGradientMatrix(Matrix3);
		virtual Matrix3 GetDisplacementGradientMatrix(void) const;
		virtual Matrix3 GetDisplacementGradientForJ(const MaterialBase *);
       	virtual Matrix3 GetElasticLeftCauchyMatrix(void);
        virtual void GetDeformationGradient(double F[][3]) const;
        virtual double GetRelativeVolume(void);
		virtual double GetVolume(int);
        virtual void GetSemiSideVectors(Vector *,Vector *,Vector *) const;
		virtual void ScaleSemiSideVectorsForCPDI(Vector *,Vector *,Vector *) const;
		virtual double GetDeformedRadius(Vector *) const;
        virtual void GetUndeformedSemiSides(double *,double *,double *) const;
		virtual Vector GetParticleSize(void) const;
		virtual double GetMinParticleLength(void) const;
		virtual void GetCPDINodesAndWeights(int);
		virtual Vector GetSurfaceInfo(int,int,int *,Vector *,Vector *,int *,double *);
		virtual Matrix3 GetInitialRotation(void);
		virtual Matrix3 GetElasticBiotStrain(void);
};

#endif
