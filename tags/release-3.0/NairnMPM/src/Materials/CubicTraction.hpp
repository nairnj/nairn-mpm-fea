/********************************************************************************
    CubicTraction.hpp
    NairnMPM
    
    Created by John Nairn on 4/6/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		TractionLaw.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef CUBICTRACTIONMATERIAL

#define CUBICTRACTIONMATERIAL 14

#include "Materials/TractionLaw.hpp"

class CubicTraction : public TractionLaw
{
    public:
	
		// constructors and destructors
        CubicTraction(char *);
		
		// methods
		virtual const char *VerifyProperties(int);
		virtual void PrintMechanicalProperties(void);
		char *MaterialData(void);
	
		// the traction law
		virtual void CrackTractionLaw(CrackSegment *,double,double,double,double,double);
		virtual double CrackTractionEnergy(CrackSegment *,double,double,bool);
		
		// accessors
		virtual const char *MaterialType(void);
		virtual const char *SetTractionLaw(double &,double &,double &,double &);
		virtual int MaterialTag();
		
	protected:
		double kI1,kII1;
};

#endif



