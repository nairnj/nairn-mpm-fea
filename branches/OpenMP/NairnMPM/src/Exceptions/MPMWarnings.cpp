/********************************************************************************
    MPMWarnings.hpp
    NairnMPM
    
    Created by John Nairn on Tues Feb 08 2005.
    Copyright (c) 2005 John A. Nairn, All rights reserved.    
********************************************************************************/

#include "Exceptions/MPMWarnings.hpp"
#include "System/ArchiveData.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"

// Single global warnings object
MPMWarnings warnings;

/*******************************************************************
	MPMWarnings: Constructors and Destructors
*******************************************************************/

// Constructors
MPMWarnings::MPMWarnings()
{
}

// add a warning
// not thread safe due to push_back()
int MPMWarnings::CreateWarning(const char *message,int abortStep,int maxIDs)
{
	WarningsData *newWarning = new WarningsData;
	newWarning->numSteps=0;
	newWarning->msg=new char[strlen(message)+1];
	newWarning->maxIssues=abortStep;
	strcpy(newWarning->msg,message);
	if(maxIDs==0)
		newWarning->issuedIDs = NULL;
	else
	{	newWarning->issuedIDs = (int *)malloc((maxIDs+1)*sizeof(int));
		newWarning->issuedIDs[0] = 0;
	}
	warningSet.push_back(newWarning);
	return warningSet.size()-1;
}

/*******************************************************************
	MPMWarnings: Methods
*******************************************************************/

// Activate all warnings for the next step
void MPMWarnings::BeginStep(void)
{	unsigned i;
	for(i=0;i<warningSet.size();i++)
	{	warningSet[i]->thisStep=TRUE;
	}
}

// issue a warning of a certain type
// but, only count once per step, and display message only on the first step this warning occurs
// return SILENT_WARNING, GAVE_WARNING, or REACHED_MAX_WARNINGS
// not thread safe due to push_back()
int MPMWarnings::Issue(int warnKind,int theID)
{	int warnResult=SILENT_WARNING;
	
	// exit if invalid warning ID
	if(warnKind<0 || warnKind>=(int)warningSet.size()) return warnResult;
	
	// check the warning
	WarningsData *warn=warningSet[warnKind];
	
	// will output warning if first time or if first time for this ID for warnding that use IDs
	bool newID=FALSE;
	if(theID>0 && warn->issuedIDs!=NULL)
	{	int i=0;
		newID=TRUE;
		while(warn->issuedIDs[i]!=0)
		{	if(theID==warn->issuedIDs[i])
			{	newID=FALSE;
				break;
			}
			i++;
		}
		if(newID)
		{	// warning - when create warning, better save room for all possible IDs
			warn->issuedIDs[i]=theID;
			warn->issuedIDs[i+1]=0;
		}
	}
	
	// increment number of steps (first time on each step)
	if(warn->thisStep) warn->numSteps++;
	
	// if first time in analysis, print message and force archiving
	if((warn->numSteps==1 && warn->thisStep) || newID)
	{   archiver->ForceArchiving();
		warn->firstStep=fmobj->mstep;
		cout << "# " << warn->msg;
		if(theID>0) cout << " (ID=" << theID << ")";
		cout << endl;
		warnResult=GAVE_WARNING;
	}
	
	// deactivate for this step
	warn->thisStep=FALSE;
	
	// abort if reach maximum issues
	if(warn->numSteps>=warn->maxIssues && warn->maxIssues>0)
		warnResult=REACHED_MAX_WARNINGS;
	return warnResult;
}

// Report all warning results
void MPMWarnings::Report(void)
{
	unsigned i;
	bool hasTitle=FALSE;
	
	// check for any warning
	for(i=0;i<warningSet.size();i++)
	{	WarningsData *warn=warningSet[i];
	
		// continue if never occurred
		if(warningSet[i]->numSteps==0) continue;
	
		// print warning section title once
		if(!hasTitle)
		{	PrintSection("ANALYSIS WARNINGS");
			hasTitle=TRUE;
		}
		
		// this warning
		cout << "The warning '" << warn->msg << "' occured on " << warn->numSteps << " time steps." << endl;
        cout << "   First such warning warning was on step " << warn->firstStep << endl;
        cout << endl;
	}
}
