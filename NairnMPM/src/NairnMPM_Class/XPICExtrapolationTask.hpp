/********************************************************************************
	XPICExtrapolationTask.hpp
	nairn-mpm-fea

	Created by John Nairn on June 30, 2016
	Copyright (c) 2016 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _XPICEXTRAPOLATIONTASK_

#define _XPICEXTRAPOLATIONTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class NodalPoint;
class GridPatch;
class MPMBase;

class XPICExtrapolationTask : public MPMTask
{
	public:
	
		// constructor
		XPICExtrapolationTask(const char *);
	
		// required methods
		virtual bool Execute(int);

		// custom methods
		virtual int GetXPICOrder(void);
		virtual void InitializeXPICData(NodalPoint *,double,int);
		virtual void InitializeXPICData(GridPatch *,int);
		virtual void XPICDoubleLoop(MPMBase *,int,int *,double *,int,double,double);
		virtual void ReduceXPICData(GridPatch *,int);
		virtual void UpdateXStar(NodalPoint *,double,int,int,double);
	
	protected:
	
};

extern XPICExtrapolationTask *XPICMechanicsTask;

#endif
