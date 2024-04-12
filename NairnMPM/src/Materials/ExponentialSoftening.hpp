/********************************************************************************
	ExponentialSoftening.hpp
	nairn-mpm-fea
	 
	Created by John Nairn, Dec 25, 2016.
	Copyright (c) 2016 John A. Nairn, All rights reserved.
	 
	Dependencies
		SofteningLaw.hpp
********************************************************************************/

#ifndef _ExponentialSoftening_
#define _ExponentialSoftening_

#include "Materials/SofteningLaw.hpp"

class ExponentialSoftening : public SofteningLaw
{
	public:
		// required methods
		virtual const char *GetSofteningLawName(void) const;
		virtual double GetFFxn(double,double) const;
		virtual double GetDeltaMax(double) const;
		virtual double GetFpFxn(double,double) const;
		virtual double GetGToDelta(double,double) const;
		virtual double GetGoverGc(double,double) const;
		virtual double GetEtaStability(void) const;
		virtual double GetPhiFxn(double,double) const;
		virtual double GetRdFxn(double,double,double) const;
	
		// optional methods
		virtual double GetRelFFxn(double,double,double) const;
		virtual double GetRelFpFxn(double,double,double) const;
		virtual double GetDeltaFromDamage(double,double,double,double);
		virtual double GetDDelta(double,double,double,double,double) const;

};

#endif
