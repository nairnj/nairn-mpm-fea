/********************************************************************************
	UpdateMomentaTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	This task updates momenta on the nodes using
 
		pk(i+1) = pk(i) + ftot(i)*dt
 
	for each velocity field on each node.
 
	Once get new momenta, check for material contact and then crack
	contact. If either changes, change force too to keep consistent with
	momentum change. The material contact checks all nodes. The crack contact
	looks only at nodes known to have cracks
 
	Finally, update analog of momenta (genearalized forces) in all
	transport tasks
********************************************************************************/

#include "NairnMPM_Class/UpdateMomentaTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Cracks/CrackNode.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark CONSTRUCTORS

UpdateMomentaTask::UpdateMomentaTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Update grid momenta and transport rates
void UpdateMomentaTask::Execute(void)
{
	CommonException *umErr = NULL;
	
#pragma omp parallel for
	for(int i=1;i<=nnodes;i++)
	{	NodalPoint *ndptr = nd[i];
		ndptr->UpdateMomentaOnNode(timestep);
		
		// get grid transport rates (update transport properties when particle state updated)
		// do first so both material and crack contact will have actual rates
		TransportTask *nextTransport=transportTasks;
		while(nextTransport!=NULL)
			nextTransport=nextTransport->TransportRates(ndptr,timestep);

		// material contact
		if(fmobj->multiMaterialMode)
		{	try
			{	ndptr->MaterialContactOnNode(timestep,UPDATE_MOMENTUM_CALL,NULL,NULL);
			}
			catch(CommonException err)
			{	if(umErr==NULL)
				{
#pragma omp critical (error)
					umErr = new CommonException(err);
				}
			}
		}
		
	}
	
	// throw error now
	if(umErr!=NULL) throw *umErr;
		
	// adjust momenta and forces for crack contact on known nodes
	CrackNode::CrackContactTask4(timestep);
}
	
