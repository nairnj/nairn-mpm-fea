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
	
		// name for the law
		virtual const char *GetSofteningLawName(void) const =0;

		// get f(delta,s) (dimensionless softening law)
		virtual double GetFFxn(double,double) const = 0;
	
		// maximum value of delta
		virtual double GetDeltaMax(double) const = 0;

		// get df(delta,s)/ddelta
		virtual double GetFpFxn(double,double) const = 0;
	
		// Get Omega(delta,s)/sigma0 in notes = int_0^delta f(u) du - 0.5 delta f(delta)
		virtual double GetGToDelta(double,double) const;
	
		// Get Omega(delta,s)/Omega(deltamax,s) (i.e., fracture of energy released)
		virtual double GetGoverGc(double,double) const;
	
		// the eta stability factor
		// Get max(-f'[delta,s])*sGc but no law depends on s
		virtual double GetEtaStability(void) const = 0;

		// Get Phi functino
		virtual double GetPhiFxn(double,double) const = 0;

		// Get R functino
		virtual double GetRdFxn(double,double,double) const = 0;
	
		// General methods to override if if a better approach
		virtual double GetRelFFxn(double,double,double) const;
		virtual double GetRelFpFxn(double,double,double) const;
		virtual double GetDeltaFromDamage(double,double,double,double);
		virtual double GetDDelta(double,double,double,double,double) const;
        virtual double GetDDeltaElastic(double,double,double,double,double) const;


		// general accessors with optional override
		virtual double GetGc(void) const;
		virtual bool IsLinear(void) const;
		virtual bool HasFailed(double,double &,double) const;

	protected:
		double Gc;
		double minFdelta,unscaledDeltaMax;
};

#endif
