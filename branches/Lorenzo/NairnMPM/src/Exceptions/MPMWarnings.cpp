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
int MPMWarnings::CreateWarning(const char *message,int abortStep)
{
	WarningsData *newWarning = new WarningsData;
	newWarning->numSteps=0;
	newWarning->msg=new char[strlen(message)+1];
	newWarning->maxIssues=abortStep;
	strcpy(newWarning->msg,message);
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
int MPMWarnings::Issue(int warnKind,int theID)
{	int warnResult=SILENT_WARNING;
	
	// exit if invalid warning ID
	if(warnKind<0 || warnKind>=(int)warningSet.size()) return warnResult;
	
	// check the warning
	WarningsData *warn=warningSet[warnKind];
	
	// will output warning if first time or if first time for this ID
	bool newID=FALSE;
	if(theID>0)
	{	int i;
		newID=TRUE;
		for(i=0;i<(int)warn->issuedIDs.size();i++)
		{	if(theID==warn->issuedIDs[i])
			{	newID=FALSE;
				break;
			}
		}
		if(newID) warn->issuedIDs.push_back(theID);
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
