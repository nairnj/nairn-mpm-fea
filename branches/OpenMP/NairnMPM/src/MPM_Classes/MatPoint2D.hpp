/********************************************************************************
    MatPoint2D.hpp
    NairnMPM
    
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
        virtual double thickness(void);
		virtual void UpdateStrain(double,int,int,void *,int);
		virtual void GetFintPlusFext(Vector *,double,double,double,double);
		virtual void MovePosition(double,Vector *);
		virtual void MoveVelocity(double,double,Vector *);
		virtual void SetVelocitySpeed(double);
		virtual void AddTemperatureGradient(Vector *);
		virtual double FCond(double,double,double,TransportProperties *);
		virtual void AddConcentrationGradient(Vector *);
		virtual double FDiff(double,double,double,TransportProperties *);
		virtual double KineticEnergy(void);
		virtual Matrix3 GetDeformationGradientMatrix(void);
        virtual Matrix3 GetElasticLeftCauchyMatrix(void);
        virtual void GetDeformationGradient(double F[][3]);
        virtual double GetRelativeVolume(void);
		virtual double GetVolume(bool);
		virtual void GetCPDINodesAndWeights(int);
		virtual double GetTractionInfo(int,int,int *,Vector *,Vector *,int *);
    
    protected:
        double thick;
};

#endif
