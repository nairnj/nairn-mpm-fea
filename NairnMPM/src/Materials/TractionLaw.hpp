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
        TractionLaw(char *);
		
		// methods
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintTransportProperties(void) const;
		virtual void ReportDebond(double,CrackSegment *,double,double);
		virtual void CalculateTimeFunction();
		
        // prototypes for traction law methods (subclasses must override)
		virtual void CrackTractionLaw(CrackSegment *,double,double,double,double,double);
		virtual double CrackTractionEnergy(CrackSegment *,double,double,bool);
		
		// accessors
		virtual const char *MaterialType(void) const;
		virtual double WaveSpeed(bool,MPMBase *) const;
		virtual int MaterialStyle(void) const;
	
	protected:
		double stress1,stress2;
	
 };

#endif


