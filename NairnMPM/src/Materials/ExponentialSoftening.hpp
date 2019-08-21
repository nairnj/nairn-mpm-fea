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
		// methods
		virtual double GetFFxn(double,double) const;
		virtual double GetFpFxn(double,double) const;
		virtual double GetRelFFxn(double,double,double) const;
		virtual double GetRelFpFxn(double,double,double) const;
		virtual double GetGToDelta(double,double) const;
		virtual double GetGoverGc(double,double) const;
		virtual double GetMaxSlope(double) const;
	
		// accessors
		virtual const char *GetSofteningLawName(void) const;
		virtual double GetDeltaMax(double) const;
};

#endif
