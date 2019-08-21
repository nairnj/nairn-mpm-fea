/********************************************************************************
	SofteningLaw.hpp
	nairn-mpm-fea

	Created by John Nairn, June 26, 2015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Base class for all softening laws

	Dependencies
	none
********************************************************************************/

#ifndef _SofteningLaw_
#define _SofteningLaw_

class SofteningLaw
{
	public:
		// constructors and destructors
		SofteningLaw();
        virtual ~SofteningLaw();
	
		// initialize
		virtual char *InputSofteningProperty(char *,int &,double &);
		virtual void PrintSofteningProperties(double);
	
		// methods must be overridden
		virtual double GetFFxn(double,double) const = 0;
		virtual double GetFpFxn(double,double) const = 0;
		virtual double GetGToDelta(double,double) const = 0;
		virtual double GetGoverGc(double,double) const = 0;
		virtual double GetMaxSlope(double) const = 0;
	
		// optional override
		virtual double GetRelFFxn(double,double,double) const;
		virtual double GetRelFpFxn(double,double,double) const;
		virtual double GetDDelta(double,double,double,double) const;
		virtual bool HasFailed(double,double &,double) const;
		virtual double GetDeltaFromDamage(double,double,double);

		// accessors must be overridden
		virtual const char *GetSofteningLawName(void) const =0;
		virtual double GetDeltaMax(double) const = 0;
	
		// accessors optional override
		virtual double GetGc(void) const;
		virtual bool IsLinear(void) const;
	
	protected:
		double Gc;
		double minFdelta;
};

#endif
