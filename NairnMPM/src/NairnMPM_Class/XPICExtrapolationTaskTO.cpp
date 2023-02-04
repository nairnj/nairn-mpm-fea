/********************************************************************************
	XPICExtrapolationTaskTO.cpp
	nairn-mpm-fea
 
	Created by John Nairn on Feb 1, 2018
	Copyright (c) 2018 John A. Nairn, All rights reserved.
 
	The task only activated with XPIC is on.
 
	The tasks are:
	--------------
	* Find Tstar on all nodes
 	* Use to change particle transport values
 ********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/XPICExtrapolationTaskTO.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Patches/GridPatch.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"

// class globals
XPICExtrapolationTaskTO *XPICTransportTask=NULL;

#pragma mark CONSTRUCTORS

XPICExtrapolationTaskTO::XPICExtrapolationTaskTO(const char *name) : XPICExtrapolationTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get order to transport calculations
int XPICExtrapolationTaskTO::GetXPICOrder(void)
{	// Return 0 for FLIP or FMPM(1), otherwise the order
	return TransportTask::XPICOrder > 1 ? TransportTask::XPICOrder : 0;
}

// Initialize real node gTstar, gTnext, and gTprev
void XPICExtrapolationTaskTO::InitializeXPICData(NodalPoint *ndptr,double timestep,int xpicOption)
{	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
	{	if(nextTransport->IsUsingTransportXPIC())
			nextTransport = nextTransport->InitializeForXPIC(ndptr,timestep,xpicOption);
		else
			nextTransport = nextTransport->GetNextTransportTask();
	}
}

// gTnext needs to be zero on ghost nodes too
void XPICExtrapolationTaskTO::InitializeXPICData(GridPatch *patchPtr,int xpicOption)
{	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
	{	if(nextTransport->IsUsingTransportXPIC())
			patchPtr->InitializeForXPICTransport(nextTransport,xpicOption);
		nextTransport = nextTransport->GetNextTransportTask();
	}
}

// need gradients in XPIC double loop
bool XPICExtrapolationTaskTO::XPICDoubleLoopNeedsGradients(void) { return false; }

// Double XPIC loop for transport tasks
void XPICExtrapolationTaskTO::XPICDoubleLoop(MPMBase *mpmptr,int matfld,int *nds,double *fn,int pn,double scale,
											 double scaleContact,double *xDeriv,double *yDeriv,double *zDeriv)
{
	// number of nodes
	int numnds = nds[0];
	
	// loop over transport tasks
	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
	{	if(nextTransport->IsUsingTransportXPIC())
		{	// add to gTnext for each node
			for(int i=1;i<=numnds;i++)
			{	// get ci from the real node
				TransportField *gTransi = nextTransport->GetTransportFieldPtr(nd[nds[i]]);
				double ci = gTransi->gVCT;
				
				// get particle Vp*CTp*(shape function)/ci
				double VpCTpWeight = nextTransport->GetVpCTp(mpmptr)*fn[i]/ci;
				
				// get node pointer to increment (which may now be ghost node)
				NodalPoint *ndptri = GetNodePointer(pn,nds[i]);
				gTransi = nextTransport->GetTransportFieldPtr(ndptri);
				
				// loop over nodes again
				for(int j=1;j<=numnds;j++)
				{	// read gTprev from real node j and add to node i (which may be ghost)
					TransportField *gTransj = nextTransport->GetTransportFieldPtr(nd[nds[j]]);
					gTransi->gTnext += scale*VpCTpWeight*fn[j]*gTransj->gTprev;
				}
			}
		}
		
		// next task
		nextTransport = nextTransport->GetNextTransportTask();
	}
}

// copy gTnext to real node and zero it
void XPICExtrapolationTaskTO::ReduceXPICData(GridPatch *patchPtr,int k)
{	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
	{	if(nextTransport->IsUsingTransportXPIC())
			patchPtr->XPICReductionTransport(k,nextTransport);
		nextTransport = nextTransport->GetNextTransportTask();
	}
}

// Update gTstar and reset gTprev and GTnext
void XPICExtrapolationTaskTO::UpdateXStar(NodalPoint *ndptr,double timestep,int m,int k,double vsign)
{	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
	{	if(nextTransport->IsUsingTransportXPIC())
		{	// add to gTstar
			TransportField *gTrans = nextTransport->GetTransportFieldPtr(ndptr);
			gTrans->gTstar += vsign*gTrans->gTnext;
			
			// copy next to prev and zero
			gTrans->gTprev = gTrans->gTnext;
			gTrans->gTnext = 0.;
		}
		
		// next transport
		nextTransport = nextTransport->GetNextTransportTask();
	}
}

// Copy gTstar to gTvalue
void XPICExtrapolationTaskTO::CopyXStar(NodalPoint *ndptr)
{	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
	{	if(nextTransport->IsUsingTransportXPIC() && TransportTask::XPICOrder>1)
		{	// swap gTstar to gTvalue (swap option may not be needed in the end)
			TransportField *gTrans = nextTransport->GetTransportFieldPtr(ndptr);
			double tmp = gTrans->gTValue;
			
			// Now gTValue = theta(k) = ck^{-1} tau (from paper)
			gTrans->gTValue = gTrans->gTstar;
			
			// cTstar stores lumped result theta^(L+) on grid if ever needed
			gTrans->gTstar = tmp;
		}
		// next transport
		nextTransport = nextTransport->GetNextTransportTask();
	}
}
