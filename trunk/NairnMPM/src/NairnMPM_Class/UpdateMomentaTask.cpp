/********************************************************************************
	UpdateMomentaTask.cpp
	NairnMPM

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
 
	Input Variables
		mvf[]->pk, ftot, mass, disp, volume, volumeGrad
		cvf[]->norm
		nd[]->fcond, gRhoVp, fdiff, gVolume
 
	Output Variables
		mvf[]->pk
		nd[]->fcond, fdiff
********************************************************************************/

#include "NairnMPM_Class/UpdateMomentaTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Cracks/CrackNode.hpp"

#pragma mark CONSTRUCTORS

UpdateMomentaTask::UpdateMomentaTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Update grid momenta and transport rates
void UpdateMomentaTask::Execute(void)
{
#ifdef _PROFILE_TASKS_
	double beginTime=fmobj->CPUTime();
#endif

	// Update momenta on all nodes
	int i;
	for(i=1;i<=nnodes;i++)
		nd[i]->UpdateMomentaOnNode(timestep);
	
	// adjust momenta and forces for multimaterial contact
	NodalPoint::MaterialContact(fmobj->multiMaterialMode,TRUE,timestep);
	
	// adjust momenta and forces for crack contact on known nodes
	CrackNode::CrackContactTask4(timestep);
	
	// get grid transport rates (update transport properties when particle state updated)
	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport=nextTransport->TransportRates(timestep);
	
#ifdef _PROFILE_TASKS_
	totalTaskTime+=fmobj->CPUTime()-beginTime;
#endif
}
	
