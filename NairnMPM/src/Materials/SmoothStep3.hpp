/********************************************************************************
	SmoothStep3.hpp
	nairn-mpm-fea
 
	Created by John Nairn, 9/4/17.
	Copyright (c) 2017 John A. Nairn, All rights reserved.
 
	Dependencies
        SofteningLaw.hpp
 ********************************************************************************/

#ifndef _SmoothStep3Softening_
#define _SmoothStep3Softening_

#include "Materials/SofteningLaw.hpp"

// Defining this term uses analytical solution for ddelta, and it seems to work.
// If comment out, uses Newton's method
// The efficiency difference is small, but analyical solution looks smoother
#define SS3_ANALYTICAL
#define CUBEROOT3 1.442249570307408
#define CUBEROOT32 2.080083823051904

class SmoothStep3 : public SofteningLaw
{
    public:
		// constructors and destructors
		SmoothStep3();
	
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

		// optional methods
#ifdef SS3_ANALYTICAL
        virtual double GetDDelta(double,double,double,double,double) const;
#endif
 	
	protected:
		double k;
		double k2,kterm,k6;
		bool stepOnly;
};

#endif
