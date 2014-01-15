/********************************************************************************
    CoupledSawTooth.hpp
    nairn-mpm-fea

    Created by John Nairn on 8/30/12.
    Copyright (c) 2102 John A. Nairn, All rights reserved.

    Dependencies
        CohesiveZone.hpp, TractionLaw.hpp
********************************************************************************/

#ifndef COUPLEDSAWTOOTHMATERIAL

#define COUPLEDSAWTOOTHMATERIAL 23

#include "Materials/TractionLaw.hpp"
#include "Materials/CohesiveZone.hpp"

class CoupledSawTooth : public CohesiveZone
{
    public:
	
        // constructors and destructors
        CoupledSawTooth(char *);

        // methods
        virtual const char *VerifyAndLoadProperties(int);
		virtual char *InitHistoryData(void);
	
		// const methods
        virtual void PrintMechanicalProperties(void) const;
    
        // the traction law
        virtual void CrackTractionLaw(CrackSegment *,double,double,double,double,double);
        virtual double CrackTractionEnergy(CrackSegment *,double,double,bool);
    
        // accessors
        virtual const char *MaterialType(void) const;
        virtual int MaterialTag() const;
};

#endif
