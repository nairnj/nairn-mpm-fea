/********************************************************************************
    CohesiveZone.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 3/21/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		TractionLaw.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef COHESIVEZONEMATERIAL

#define COHESIVEZONEMATERIAL 12

#include "Materials/TractionLaw.hpp"

class CohesiveZone : public TractionLaw
{
    public:
	
		// constructors and destructors
        CohesiveZone(char *);
		
		// methods
		virtual char *InputMat(char *,int &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual char *InitHistoryData(void);
	
		// const methods
		virtual void PrintMechanicalProperties(void) const;
	
		// the traction law
		virtual void CrackTractionLaw(CrackSegment *,double,double,double,double,double);
		virtual double CrackTractionEnergy(CrackSegment *,double,double,bool);
		
		// accessors
		virtual const char *SetTractionLaw(double &,double &,double &,double &,double &);
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
		
	protected:
		double kI1,kII1,sIc,sIIc;
		double umidI,umidII;
 };

#endif


