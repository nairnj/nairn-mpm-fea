/********************************************************************************
	TrilinearTraction.hpp
	NairnMPM

	Created by John Nairn on 9/17/10.
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		CohesiveZone.hpp, TractionLaw.hpp, MaterialBase.hpp
 ********************************************************************************/

#ifndef TRILINEARTRACTIONMATERIAL

#define TRILINEARTRACTIONMATERIAL 20

#include "Materials/CohesiveZone.hpp"

class TrilinearTraction : public CohesiveZone
{
    public:
		
		// constructors and destructors
        TrilinearTraction(char *);
		
		// methods
		virtual char *InputMat(char *,int &);
		virtual const char *VerifyProperties(int);
		virtual void PrintMechanicalProperties(void);
		
		// the traction law
		virtual void CrackTractionLaw(CrackSegment *,double,double,double,double,double);
		virtual double CrackTractionEnergy(CrackSegment *,double,double,bool);
		
		// accessors
		virtual const char *SetTLTractionLaw(double &,double &,double &,double &,double &,double &,double &);
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		
	protected:
		double sI2,uI2;
		double sII2,uII2;
		double sIc2,sIIc2;
		bool break1is2I,break1is2II;
	
};

#endif



