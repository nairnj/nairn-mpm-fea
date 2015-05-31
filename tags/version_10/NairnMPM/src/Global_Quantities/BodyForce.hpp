/********************************************************************************
    BodyForce.hpp
    nairn-mpm-fea
    
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
        bool useDamping;            // true when used
        double dampingCoefficient;	// 1/Q in Nose-Hoover feedback damping
        bool useFeedback;           // true when used
    
        double pdamping;                // same for particle damping
        bool usePDamping;
        double pdampingCoefficient;
        bool usePFeedback;
    
        bool useGridFeedback;       // use grid kinetic energy to evolve feedback terms
        Vector gforce;              // gravity forces
		
        // constructors and destructors
        BodyForce();
		virtual ~BodyForce();
    
        // methods
		void Activate(void);
		void AddGravity(Vector *,double,double);
		double GetDamping(double);
        double GetParticleDamping(double);
        double GetNonPICDamping(double);
        double GetNonPICParticleDamping(double);
		void Output(void);
		void UpdateAlpha(double,double);
		void SetTargetFunction(char *,bool);
        void SetMaxAlpha(double,bool);
		void SetGridDampingFunction(char *,bool);
        void SetFractionPIC(double);
	
	private:
        double alpha,maxAlpha;
        ROperation *function;
        ROperation *gridfunction;
    
        double palpha,maxPAlpha;
        ROperation *pgridfunction;
        ROperation *pfunction;
    
        double fractionPIC;
        bool usePICDamping;
    
		bool gravity;               // true if gravity turned on
		static double varTime;
};

extern BodyForce bodyFrc;

#endif
