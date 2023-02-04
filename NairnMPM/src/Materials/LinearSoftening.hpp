/********************************************************************************
	LinearSoftening.hpp
	nairn-mpm-fea

	Created by John Nairn, Jan 21 2017.
	Copyright (c) 2017 John A. Nairn, All rights reserved.

	Dependencies
		SofteningLaw.hpp
********************************************************************************/

#ifndef _LinearSoftening_
#define _LinearSoftening_

#include "Materials/SofteningLaw.hpp"

class LinearSoftening : public SofteningLaw
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
		virtual double GetDeltaFromDamage(double,double,double,double);
		virtual double GetDDelta(double,double,double,double,double) const;
        virtual double GetDDeltaElastic(double,double,double,double,double) const;
	
		// optoinal accessor
		virtual bool IsLinear(void) const;
};

#endif
