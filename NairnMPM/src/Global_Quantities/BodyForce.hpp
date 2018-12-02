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

class BodyForce
{
    public:
		double damping;				// external damping
		bool useDamping;            // true when constant grid or grid feedback damping is on
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
		void GetGridBodyForce(Vector *,Vector *,double);
		double GetDamping(double);
		double GetPICDamping(void);
        double GetParticleDamping(double);
        double GetNonPICDamping(double);
        double GetNonPICParticleDamping(double);
		void Output(void);
		void UpdateAlpha(double,double);
		void SetTargetFunction(char *,bool);
        void SetMaxAlpha(double,bool);
		void SetGridDampingFunction(char *,bool);
		void SetGridBodyForceFunction(char *,int);
		void SetFractionPIC(double);
		double GetFractionPIC(void);
		bool IsUsingPICDamping(void);

	private:
		double alpha,maxAlpha;
		Expression *function;
		Expression *gridfunction;
	
		double palpha,maxPAlpha;
		Expression *pgridfunction;
		Expression *pfunction;

		double fractionPIC;         // (1-beta) in my notes or alpha(PIC) = fractionPIC/dt
		bool usePICDamping;         // true when PIC damping activated (damping and pdamping need not be true)
		int XPICOrder;				// XPIC oder (1=normal PIC or 2 or higher for extended PIC, always 1 in NairnMPM)

		Expression *gridBodyForceFunction[3];
};

extern BodyForce bodyFrc;

#endif
