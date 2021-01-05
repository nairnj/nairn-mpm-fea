/********************************************************************************
	TrilinearTraction.hpp
	nairn-mpm-fea

	Created by John Nairn on 9/17/10.
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		CohesiveZone.hpp, TractionLaw.hpp, MaterialBase.hpp
 ********************************************************************************/

#ifndef TRILINEARTRACTIONMATERIAL

#define TRILINEARTRACTIONMATERIAL 20

#include "Materials/ExponentialTraction.hpp"

class TrilinearTraction : public ExponentialTraction
{
    public:
		
		// constructors and destructors
        TrilinearTraction(char *,int);
		
		// methods
		virtual char *InputTractionLawProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;

        // Core functions
        virtual double Strength(int,double);
        virtual double WorkEnergy(int,double);
        virtual double DissipatedEnergy(int,double);
        virtual double GetDFromDelta(int,double);
        virtual double GetDeltaFromD(int,double);

		// Trilinear methods (don't override)
        void PrintTriLinearModel(const char *,double,double,double,double,double,double,double) const;
        double GetTLDeltaPrimeFromD(double,double,double,double,double,double,double,bool);
        const char *SetTLTractionLaw(double &,double &,double &,double &,double &,double &,double &,
                                     bool &,double &,double &);
    
        // accessors (override)
		virtual const char *MaterialType(void) const;
		
	protected:
		double sI2,uI2;
		double sII2,uII2;
		bool break1is2I,break1is2II;
		double JI_1c,JI_2c,JII_1c,JII_2c;
        double DIbreak,DIIbreak;
        double phiII2,phiI2;
        double RI2,RII2;
};

#endif



