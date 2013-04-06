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
class MPMBase;
class NodalPoint;

class BodyForce
{
    public:
		double damping;				// external damping
	bool useDamping;
		double dampingCoefficient;	// 1/Q in Nose-Hoover feedback damping
		short useFeedback;
		bool useGridFeedback;
        Vector gforce;              // gravity forces
		
        // constructors and destructors
        BodyForce();
		virtual ~BodyForce();
    
        // methods
		void Activate(void);
		void AddGravity(double,double,Vector *);
		double GetDamping(double);
		void Output(void);
		double GetAlpha(void);
		void UpdateAlpha(double,double);
		void SetTargetFunction(char *);
        void SetMaxAlpha(double);
		void SetGridDampingFunction(char *);
	
	private:
        ROperation *function;
        ROperation *gridfunction;
		bool gravity;               // true if gravity turned on
		double alpha,maxAlpha;
		static double varTime;
};

extern BodyForce bodyFrc;

#endif
