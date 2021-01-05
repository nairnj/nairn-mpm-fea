/********************************************************************************
    CohesiveZone.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 3/21/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		TractionLaw.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef TRIANGULARTRACTIONMATERIAL

#define TRIANGULARTRACTIONMATERIAL 12

#include "Materials/TractionLaw.hpp"

enum { CZ_DELN=0,CZ_DELT };

class CohesiveZone : public TractionLaw
{
    public:
	
		// constructors and destructors
        CohesiveZone(char *,int);
		
		// methods
		virtual char *InputTractionLawProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
	
		// history data
		virtual char *InitHistoryData(char *);
	
		// const methods
		virtual void PrintMechanicalProperties(void) const;
 
		// the traction law
		virtual void CrackTractionLaw(CrackSegment *,double,double,Vector *,Vector *,double);
		virtual double CrackWorkEnergy(CrackSegment *,double,double);
		virtual void CrackDissipatedEnergy(CrackSegment *,double &,double &);
    
        // Core functions
        virtual double Strength(int,double);
        virtual double WorkEnergy(int,double);
        virtual double DissipatedEnergy(int,double);
        virtual double StrengthPrime(int,double);
        virtual double GetDFromDelta(int,double);
        virtual double GetDeltaFromD(int,double);

        // Sawtooth methods (don't override)
        const char *SetSawToothTractionLaw(double &,double &,double &,double &,double &,double &,double &);
        void PrintSawToothModel(const char *,double,double,double,double,double,double) const;
        double GetSawToothDeltaPrimeFromD(double,double,double);
    
        // for Exponential law (don't override)
        const char *SetExponentialTractionLaw(double &,double &,double &,double &,double &,double &);

		// accessors
		virtual const char *MaterialType(void) const;
		
	protected:
		double kI1,kII1;
		double umidI,umidII;
        double phiI1,phiII1;
        double RI1,RII1;
};

#endif


