/********************************************************************************
	CarnotCycle.hpp
	nairn-mpm-fea

	Created by John Nairn on 9/24/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.

	Dependencies
		CustomTask.hpp
********************************************************************************/

#ifndef _CARNOTCYCLE_

#define _CARNOTCYCLE_

#include "Custom_Tasks/CustomTask.hpp"

class CarnotCycle : public CustomTask
{
	public:
    
		// constructors and destructors
		CarnotCycle();
    
		// standard methods
		virtual const char *TaskName(void);
		virtual char *InputParam(char *,int &,double &);
		virtual CustomTask *Initialize(void);
	
		virtual CustomTask *StepCalculation(void);
    
	private:
		double V1rel,V2rel,V3rel;
		double T2;
		int carnotStep;
    
};

#endif

