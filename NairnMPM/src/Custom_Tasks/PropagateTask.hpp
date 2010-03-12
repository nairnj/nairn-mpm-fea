/********************************************************************************
    PropagateTask.hpp
    NairnMPM
    
    Created by John Nairn on Mon Aug 18 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.

	Dependencies
		CustomTask.hpp
********************************************************************************/

#ifndef _PROPAGATETASK_

#define _PROPAGATETASK_

#include "Custom_Tasks/CustomTask.hpp"

class PropagateTask : public CustomTask
{
    public:
        double nextPropTime;
        int doPropCalcs,theResult;
        double totalPlastic,totalPotential;
        
		// static variables
		static double cellsPerPropagationStep;
	
        // constructors and destructors
        PropagateTask();
                
        // standard methods
		virtual const char *TaskName(void);
        virtual CustomTask *Initialize(void);
        virtual CustomTask *PrepareForStep(bool &);
        virtual CustomTask *BeginExtrapolations(void);
        virtual CustomTask *ParticleExtrapolation(MPMBase *);
        virtual CustomTask *StepCalculation(void);
        
        // special methods
        void ArrestGrowth(int);
        bool Arrested(void);
        
    private:
        bool arrested;
		
};

extern PropagateTask *propagateTask;

#endif
