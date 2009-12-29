/********************************************************************************
    BodyForce.hpp
    NairnMPM
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.  
	
	Header for gravity and damping
		Global bodyFrc

	Dependencies
		none
********************************************************************************/

#ifndef _BODYFORCE_

#define _BODYFORCE_

class ROperation;

class BodyForce
{
    public:
		double gforcex,gforcey;		// gravity forces
		double damping;				// external damping
		double dampingCoefficient;	// 1/Q in Nose-Hoover feedback damping
		short useFeedback;
		
        // constructors and destructors
        BodyForce();
		virtual ~BodyForce();
    
        // methods
		void Activate(void);
		short GetGravity(double *,double *);
		double GetDamping(void);
		void Output(void);
		double GetAlpha(void);
		void TrackAlpha(void);
		void TrackAlpha(Vector *,double);
		void UpdateAlpha(double,double);
		void SetTargetFunction(char *);
	
	private:
		short gravity;
		double alpha;
		double kineticEnergy,totalMass;
		ROperation *function;
		static double varTime;
};

extern BodyForce bodyFrc;

#endif
