/********************************************************************************
    TractionLaw.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Feb 22 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef _TRACTIONMATERIAL_

#define _TRACTIONMATERIAL_

#include "Materials/MaterialBase.hpp"

class CrackSegment;

class TractionLaw : public MaterialBase
{
    public:
	
		// constructors and destructors
        TractionLaw(char *,int);
		
		// methods - subclass should override InputContactProperty only
		virtual char *InputMaterialProperty(char *,int &,double &);
		virtual char *InputTractionLawProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
        virtual char *PackLabel(char *,const char *,const char *,const char *) const;
		virtual void PrintTransportProperties(void) const;
		virtual void ReportDebond(double,CrackSegment *,double,double);
		virtual void CalculateTimeFunction();
		
        // prototypes for traction law methods (subclasses must override)
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
        virtual double DissipationRate(int,double);
        virtual double RatioFunction(int,double);

        // cubic law methods (don't override)
        const char *SetCubicTractionLaw(double &,double &,double &,double &);
        void PrintCubicModel(const char *,double,double,double,double) const;

		// accessors
		virtual const char *MaterialType(void) const;
		virtual double WaveSpeed(bool,MPMBase *) const;
		virtual int MaterialStyle(void) const;
        virtual int NumberOfHistoryDoubles(void) const;

	protected:
		double stress1,stress2;
        int numTractionHistory;
	
 };

#endif


