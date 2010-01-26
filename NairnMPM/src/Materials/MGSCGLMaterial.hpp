/********************************************************************************
    MGSCGLMaterial.hpp
    NairnMPM
    
    Created by John Nairn, 11/12/2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Isoplasticity.hpp MaterialBase.hpp
********************************************************************************/

#ifndef MGSCGLMATERIAL

#define MGSCGLMATERIAL 17

#include "Materials/IsoPlasticity.hpp"

class MGSCGLMaterial : public IsoPlasticity
{
    public:
		// unique properties
		double gamma0,GPp,GTp,beta,nhard,yieldMax;
		double C0,S1,S2;
        
        // constructors and destructors
		MGSCGLMaterial();
		MGSCGLMaterial(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyProperties(int);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
		virtual void PrintYieldProperties(void);
		virtual void PrintTransportProperties(void);
		
		// methods
		virtual double GetPressureChange(MPMBase *,double &,int);
		virtual double GetYield(MPMBase *,int,double);
 		virtual double GetKPrime(MPMBase *,int,double);
		virtual double GetK2Prime(MPMBase *,double,double);
				
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
        virtual double WaveSpeed(bool);
		
    protected:
		double k1,k2,k3,gRhoCv,GPpred,G0red,yldMaxred;

};

#endif

