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

class Expression;
class MPMBase;
class NodalPoint;

#define VSTAR_NOT_USED 0
#define VSTAR_NO_CONTACT 1
#define VSTAR_WITH_CONTACT 2
#define FMPM_WITH_CONTACT 2
#include "System/MPMPrefix.hpp"

enum { GRIDBC_COMBINED=0,GRIDBC_VELOCITY_ONLY,GRIDBC_LUMPED_ONLY };

class BodyForce
{
    public:
		double damping;				// external damping
		bool useDamping;            // true when grid damping (alpha, feedback)
		double dampingCoefficient;	// 1/Q in Nose-Hoover feedback damping
		bool useFeedback;           // true when grid feedback damping is on (useDamping true too)
        double pdamping;            // same for particle damping
        bool usePDamping;           // true when constant particle or particle feedback damping is on
        double pdampingCoefficient; // 1/Q in Nose-Hoover feedback damping
        bool usePFeedback;          // true when particle feedback damping is on (usePDamping true too)
    
		bool useGridFeedback;       // use grid kinetic energy to evolve feedback terms
        Vector gforce;              // gravity forces
		bool gravity;               // true if gravity turned on
		bool hasGridBodyForce;		// true if has any grid body force functions
		
        // constructors and destructors
        BodyForce();
		virtual ~BodyForce();
    
        // methods
		void Activate(void);
		bool HasGridDampingForces(void);
		void GetGridBodyForce(Vector *,Vector *,double);
		double GetParticleDamping(double);
		double GetGridDamping(double);
		void Output(void);
		void UpdateAlpha(double,double);
		void SetTargetFunction(char *,bool);
        void SetMaxAlpha(double,bool);
		void SetGridDampingFunction(char *,bool);
		void SetGridBodyForceFunction(char *,int);
	
		void SetXPICOrder(int);
		int GetXPICOrder(void);
		int UsingVstar(void);
		void SetUsingVstar(int);
		int XPICVectors(void);
		void SetXPICVectors(int);
		void SetUsingFMPM(bool);
		bool UsingFMPM(void);
		int GridBCOption(void);
		void SetGridBCOption(int);

	private:
		double alpha,maxAlpha;
		Expression *function;
		Expression *gridfunction;
		double palpha,maxPAlpha;
		Expression *pgridfunction;
		Expression *pfunction;
		int XPICOrder;				// XPIC oder (0 FLIP, 1 or higher for FMPM(k))
		int isUsingVstar;			// if current time step uses vstar (0=0, 1=no contact, 2=when contact)
		int xpicVectors;			// 3 or 1 as needed
		bool usingFMPM;				// True if FMPM (otherwise XPIC or FLIP)
		int gridBCOption;			// velocity BCs on grid for FMPM and XPIC k>1

		Expression *gridBodyForceFunction[3];
};

extern BodyForce bodyFrc;

#endif
