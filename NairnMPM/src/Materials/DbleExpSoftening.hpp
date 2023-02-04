/********************************************************************************
	DbleExpSoftening.hpp
	nairn-mpm-fea
	 
	Created by John Nairn, 6 July 2021.
	Copyright (c) 2021 John A. Nairn, All rights reserved.
	 
	Dependencies
		SofteningLaw.hpp
********************************************************************************/

#ifndef _DbleExpSoftening_
#define _DbleExpSoftening_

#include "Materials/SofteningLaw.hpp"

class DbleExpSoftening : public SofteningLaw
{
	public:
		// constructors and destructors
		DbleExpSoftening();

		// initialize
		virtual char *InputSofteningProperty(char *,int &,double &);
		virtual void PrintSofteningProperties(double);
	
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
	
	protected:
		double alpha,beta;
		double kterm,ab,pkRatio,fmaxRatio;
		double oneOver1minusBeta;
};

#endif
