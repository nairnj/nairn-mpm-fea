/********************************************************************************
	XPICExtrapolationTaskTO.hpp
	nairn-mpm-fea
 
	Created by John Nairn on Feb 1, 2018
	Copyright (c) 2018 John A. Nairn, All rights reserved.
 
	Dependencies
 		MPMTask, CommonTask
 ********************************************************************************/

#ifndef _XPICEXTRAPOLATIONTASKTO_

#define _XPICEXTRAPOLATIONTASKTO_

#include "NairnMPM_Class/XPICExtrapolationTask.hpp"

class XPICExtrapolationTaskTO : public XPICExtrapolationTask
{
	public:
	
		// constructor
		XPICExtrapolationTaskTO(const char *);

		// custom methods
		virtual int GetXPICOrder(void);
		virtual void InitializeXPICData(NodalPoint *,double,int);
		virtual void InitializeXPICData(GridPatch *,int);
		virtual bool XPICDoubleLoopNeedsGradients(void);
		virtual void XPICDoubleLoop(MPMBase *,int,int *,double *,int,double,double,double *,double *,double *);
		virtual void ReduceXPICData(GridPatch *,int);
		virtual void UpdateXStar(NodalPoint *,double,int,int,double);
		virtual void CopyXStar(NodalPoint *);
	
	protected:
	
};

extern XPICExtrapolationTaskTO *XPICTransportTask;

#endif
