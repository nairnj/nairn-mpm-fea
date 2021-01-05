/********************************************************************************
    Exponential.hpp
    nairn-mpm-fea

    Created by John Nairn on 12/22/2020.
    Copyright (c) 2020 John A. Nairn, All rights reserved.

    Dependencies
        CohesiveZone.hpp, TractionLaw.hpp
********************************************************************************/

#ifndef EXPONENTIALTRACTIONMATERIAL

#define EXPONENTIALTRACTIONMATERIAL 34

#include "Materials/TractionLaw.hpp"
#include "Materials/CohesiveZone.hpp"

class ExponentialTraction : public CohesiveZone
{
    public:
    
        // constructors and destructors
        ExponentialTraction(char *,int);

        // initialization
        virtual char *InputTractionLawProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintMechanicalProperties(void) const;

        // Core functions
        virtual double Strength(int,double);
        virtual double WorkEnergy(int,double);
        virtual double DissipatedEnergy(int,double);
        virtual double GetDFromDelta(int,double);
        virtual double GetDeltaFromD(int,double);
        virtual double DissipationRate(int,double);
        virtual double RatioFunction(int,double);

        // accessors
        virtual const char *MaterialType(void) const;
    
    protected:
        double alphaI,alphaII;
        double expalphaI,romexpalphaI;
        double expalphaII,romexpalphaII;
};

#endif
