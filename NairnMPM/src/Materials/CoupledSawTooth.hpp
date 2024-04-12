/********************************************************************************
    CoupledSawTooth.hpp
    nairn-mpm-fea

    Created by John Nairn on 8/30/12.
    Copyright (c) 2012 John A. Nairn, All rights reserved.

    Dependencies
        CohesiveZone.hpp, TractionLaw.hpp
********************************************************************************/

#ifndef COUPLEDSAWTOOTHMATERIAL

#define COUPLEDSAWTOOTHMATERIAL 23

#include "Materials/TractionLaw.hpp"
#include "Materials/CohesiveZone.hpp"

enum { HOG_LAMBDA=0,HOG_DW,HOG_UN,HOG_UT };

class CoupledSawTooth : public CohesiveZone
{
    public:
	
        // constructors and destructors
        CoupledSawTooth(char *,int);

		// initialization
        virtual char *InputTractionLawProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintMechanicalProperties(void) const;

		// history data
		virtual char *InitHistoryData(char *,MPMBase *);
	
        // the traction law
        virtual void CrackTractionLaw(CrackSegment *,double,double,Vector *,Vector *,double);
		virtual double CrackWorkEnergy(CrackSegment *,double,double);
		virtual void CrackDissipatedEnergy(CrackSegment *,double &,double &);
    
        // accessors
        virtual const char *MaterialType(void) const;
	
	protected:
		double unpbar,utpbar,unpbar2,utpbar2;
        double alpha,expalpha,romexpalpha;
        int model;
};

#endif
