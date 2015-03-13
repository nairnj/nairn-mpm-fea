/********************************************************************************
    ReverseLoad.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Aug 20 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.

	Dependencies
		CustomTask.hpp
********************************************************************************/

#ifndef _REVERSELOADTASK_

#define _REVERSELOADTASK_

#include "Custom_Tasks/CustomTask.hpp"

enum { REVERSE=0,HOLD,NOCHANGE,ABORT };
enum { CHECKING_PHASE=0,HOLDING_PHASE,REVERSED_PHASE };

class ReverseLoad : public CustomTask
{
    public:
        int crackNum,style;
        double finalLength,holdTime;
        int reversed;
		int quantity,subcode,whichMat;
		char quant[64];
        
        // constructors and destructors
        ReverseLoad();
        
        // standard methods
		virtual const char *TaskName(void);
        virtual char *InputParam(char *,int &);
        virtual CustomTask *Initialize(void);
	
        virtual CustomTask *FinishForStep(void);

    private:
        double finalTime,endHoldTime;
};

#endif
