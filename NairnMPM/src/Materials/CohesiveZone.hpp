/********************************************************************************
    CohesiveZone.hpp
    NairnMPM
    
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
		virtual const char *VerifyProperties(int);
		virtual void PrintMechanicalProperties(void);
		virtual char *MaterialData();
	
		// the traction law
		virtual void CrackTractionLaw(CrackSegment *,double,double,double,double,double);
		virtual double CrackTractionEnergy(CrackSegment *,double,double,bool);
		
		// accessors
		virtual const char *SetTractionLaw(double &,double &,double &,double &,double &);
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		
	protected:
		double kI1,kII1;
		double umidI,umidII;
 };

#endif


