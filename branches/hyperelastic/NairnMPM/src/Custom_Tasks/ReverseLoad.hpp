/********************************************************************************
    ReverseLoad.hpp
    NairnMPM
    
    Created by John Nairn on Wed Aug 20 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.

	Dependencies
		CustomTask.hpp
********************************************************************************/

#ifndef _REVERSELOADTASK_

#define _REVERSELOADTASK_

#include "Custom_Tasks/CustomTask.hpp"

enum { REVERSE=0,HOLD,NOCHANGE,ABORT };

class ReverseLoad : public CustomTask
{
    public:
        int crackNum,style;
        double finalLength,finalLoad;
        bool reversed;
        
        // constructors and destructors
        ReverseLoad();
        
        // standard methods
		virtual const char *TaskName(void);
        virtual char *InputParam(char *,int &);
        virtual CustomTask *Initialize(void);
        virtual CustomTask *FinishForStep(void);

    private:
};

#endif
