/********************************************************************************
	HistoryArchive.hpp
	nairn-mpm-fea

	Created by John Nairn on 10/26/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.

	Dependencies
		CustomTask.hpp
********************************************************************************/

#ifndef _HISTORYARCHIVETASK_

#define _HISTORYARCHIVETASK_

#include "Custom_Tasks/CustomTask.hpp"

class HistoryArchive : public CustomTask
{
public:
	
	// constructors and destructors
	HistoryArchive();
	
	// standard methods
	virtual const char *TaskName(void);
	virtual char *InputParam(char *,int &,double &);
	virtual CustomTask *Initialize(void);
	
	virtual CustomTask *PrepareForStep(bool &);
	virtual CustomTask *StepCalculation(void);
	
private:
	vector< int > quantity;
	double customArchiveTime,nextCustomArchiveTime;
	bool doHistoryExport;			// flag to export the file this step
	
};

#endif

