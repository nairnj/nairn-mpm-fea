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

class MatPoint2D : public MPMBase
{
    public:
		// used in set up
        double thick;
        
        // constructors and destructors
        MatPoint2D();
        MatPoint2D(int,int,double,double);
        
        // methods
		virtual void SetOrigin(Vector *);
        virtual void SetPosition(Vector *);
        virtual void SetVelocity(Vector *);
        virtual double thickness(void);
		virtual void SetDilatedVolume(void);
		virtual void UpdateStrain(double,int,int);
		virtual Vector Fint(double,double,double);
		virtual Vector Fext(double fni);
		virtual void MovePosition(double,Vector *);
		virtual void MoveVelocity(double,double,Vector *);
		virtual void SetVelocitySpeed(double);
		virtual void AddTemperatureGradient(void);
		virtual void AddTemperatureGradient(Vector *);
		virtual double FCond(double,double,double);
		virtual void AddConcentrationGradient(void);
		virtual void AddConcentrationGradient(Vector *);
		virtual double FDiff(double,double,double);
};

#endif
