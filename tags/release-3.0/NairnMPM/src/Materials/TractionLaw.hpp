/********************************************************************************
    TractionLaw.hpp
    NairnMPM
    
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
        TractionLaw(char *);
		
		// methods
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyProperties(int);
		virtual void PrintTransportProperties(void);
		virtual void ReportDebond(double,CrackSegment *,double);
		
		// no need for constitutive laws
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,int);
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
  		
        // prototypes for traction law methods (subclasses must override)
		virtual void CrackTractionLaw(CrackSegment *,double,double,double,double,double);
		virtual double CrackTractionEnergy(CrackSegment *,double,double,bool);
		
		// accessors
		virtual const char *MaterialType(void);
		virtual double WaveSpeed(bool);
		virtual bool isTractionLaw(void);
	
	protected:
		double stress1,stress2;
	
 };

#endif


