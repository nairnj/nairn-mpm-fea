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

class SmoothStep3 : public SofteningLaw
{
    public:
        // methods
        virtual double GetFFxn(double,double) const;
        virtual double GetFpFxn(double,double) const;
        virtual double GetGToDelta(double,double) const;
        virtual double GetGoverGc(double,double) const;
        virtual double GetMaxSlope(double) const;
#ifdef SS3_ANALYTICAL
        virtual double GetDDelta(double,double,double,double) const;
#endif
		virtual double GetRdFxn(double,double,double) const;
		virtual double GetPhiFxn(double,double) const;
	
        // accessors
        virtual const char *GetSofteningLawName(void) const;
        virtual double GetDeltaMax(double) const;
		virtual double GetEtaStability(void) const;
};

#endif
