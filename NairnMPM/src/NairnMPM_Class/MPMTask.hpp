/********************************************************************************
	 MPMTask.hpp
	 nairn-mpm-fea
	 
	 Created by John Nairn on July 22, 2010
	 Copyright (c) 2010 John A. Nairn, All rights reserved.
	 
	 Dependencies
		CommonTask
********************************************************************************/

#ifndef _MPMTASK_

#define _MPMTASK_

class NodalPoint;

#include "System/CommonTask.hpp"

class MPMTask : public CommonTask
{
	public:
	
		// constructor
		MPMTask(const char *);
	
		// methods
#ifdef LOG_PROGRESS
		void WriteLogFile(void);
		void WriteLogFile(const char *);
#endif
	
		void WriteProfileResults(int,double,double);
		void TrackTimes(double,double);
	
        // class methods
        static int GetPatchNumber(void);
        static NodalPoint *GetNodePointer(int,int);
		static int GetNumberOfThreads(void);
    
	protected:
	
};

#endif
