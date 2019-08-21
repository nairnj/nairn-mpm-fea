/********************************************************************************
	MoveCracksTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	The tasks are:
	--------------
	* Decode time functions for crack tractions
	* Move crack surfaces
	* Move crack planes (using surface or CM as requested)
	* Update crack tractinos
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/MoveCracksTask.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Exceptions/CommonException.hpp"
#include "Materials/TractionLaw.hpp"

#pragma mark CONSTRUCTORS

MoveCracksTask::MoveCracksTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Run this task to move cracks and update tactions
// throws CommonException() if crack surface or position has left the grid
void MoveCracksTask::Execute(int taskOption)
{
	CommonException *mcErr = NULL;
    
    int numCracksPerProc = (int)((double)numberOfCracks/(double)(fmobj->GetNumberOfProcessors()));
    
    // uncomment this line to prevent parallelization of this task
    //numCracksPerProc = 0;

	// update crack tractions
	if(fmobj->hasTractionCracks)
	{	// evaluate time dependent functions outside the parallel loop
		// and only once per traction law
		for(int i=0;i<nmat;i++)
		{	if(theMaterials[i]->MaterialStyle()==TRACTION_MAT)
			{	((TractionLaw *)theMaterials[i])->CalculateTimeFunction();
			}
		}
	}
	
	// move crack surface in their local velocity fields
#pragma omp parallel if(numCracksPerProc>1)
	{
#pragma omp for
		for(int cn=0;cn<numberOfCracks;cn++)
		{   try
			{	if(!crackList[cn]->MoveCrack(ABOVE_CRACK))
				{	char errMsg[100];
					sprintf(errMsg,"Crack No. %d surface (above) has left the grid.",cn+1);
					throw CommonException(errMsg,"MoveCracksTask::Execute");
				}
				if(!crackList[cn]->MoveCrack(BELOW_CRACK))
				{	char errMsg[100];
					sprintf(errMsg,"Crack No. %d surface (below) has left the grid.",cn+1);
					throw CommonException(errMsg,"MoveCracksTask::Execute");
				}
			}
			catch(CommonException& err)
			{	if(mcErr==NULL)
				{
#pragma omp critical (error)
					mcErr = new CommonException(err);
				}
			}
			catch(std::bad_alloc&)
			{	if(mcErr==NULL)
				{
#pragma omp critical (error)
					mcErr = new CommonException("Memory error","MoveCracksTask::Execute");
				}
			}
			catch(...)
			{	if(mcErr==NULL)
				{
#pragma omp critical (error)
					mcErr = new CommonException("Memory error","MoveCracksTask::Execute");
				}
			}
		}
	}

	// throw any errors
	if(mcErr!=NULL) throw *mcErr;
	
	// get center of mass if moving crack planes that way
	if(!contact.GetMoveOnlySurfaces())
	{
#pragma omp parallel for
		for(int i=1;i<=nnodes;i++)
		{	nd[i]->GetCenterOfMassVelocity(NULL,false);
		}
	}

	// Move crack plane by one of two methods. When moving only surface, the plane will move to average
	// of the two surfaces found above, otherwise moves in the just-calculated center or mass field
	// After moving, crack plane is checked for crossing, if that option is activated
	// Also calculate tractions
#pragma omp parallel if(numCracksPerProc>1)
	{
#pragma omp for
		for(int cn=0;cn<numberOfCracks;cn++)
		{   try
			{   if(!crackList[cn]->MoveCrack())
				{	char errMsg[100];
					sprintf(errMsg,"Crack No. %d position or surface has left the grid.",cn+1);
					throw CommonException(errMsg,"MoveCracksTask::Execute");
				}
				
				// tractions
				crackList[cn]->UpdateCrackTractions();
			}
			catch(CommonException& err)
			{	if(mcErr==NULL)
				{
#pragma omp critical (error)
					mcErr = new CommonException(err);
				}
			}
			catch(std::bad_alloc&)
			{	if(mcErr==NULL)
				{
#pragma omp critical (error)
					mcErr = new CommonException("Memory error","MoveCracksTask::Execute");
				}
			}
			catch(...)
			{	if(mcErr==NULL)
				{
#pragma omp critical (error)
					mcErr = new CommonException("Memory error","MoveCracksTask::Execute");
				}
			}
		}
	}

	// throw any errors
	if(mcErr!=NULL) throw *mcErr;
	
}
	
