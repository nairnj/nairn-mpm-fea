/********************************************************************************
    MPMWarnings.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 08 2005.
    Copyright (c) 2005 John A. Nairn, All rights reserved.    
********************************************************************************/

#include "Exceptions/MPMWarnings.hpp"
#include "System/ArchiveData.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MPMTask.hpp"

// Single global warnings object
MPMWarnings warnings;

#pragma mark MPMWarnings: Constructors and Destructors

// Constructors
MPMWarnings::MPMWarnings()
{	lastWarning = NULL;
	numWarnings = 0;
	warningSet = NULL;
}

// add a warning message
// abortStep = maximum allowed before calculations stopped, use -1 to never stop
// maxIDs = if warning occurs several places with different IDs, this is maximum ID that is used
// throws std::bad_alloc
int MPMWarnings::CreateWarning(const char *message,int abortStep,int maxIDs)
{
	// do not duplicate previous warning
	WarningsData *nextWarning = lastWarning;
	int i = numWarnings-1;
	while(nextWarning!=NULL)
	{	if(strcmp(message,nextWarning->msg)==0)
		{	// warning already exists, so reuse it
			return i;
		}
		nextWarning = (WarningsData *)(nextWarning->prevWarning);
		i--;
	}
	
	// it is a new warning
	WarningsData *newWarning = new WarningsData;
	newWarning->numSteps=0;
	newWarning->numWarnings=0;
	newWarning->maxIssues=abortStep;
	newWarning->msg=new char[strlen(message)+1];
	strcpy(newWarning->msg,message);
	newWarning->maxIDs=maxIDs;
	if(maxIDs<=0)
		newWarning->issuedIDs = NULL;
	else
	{	newWarning->issuedIDs = new int[maxIDs+1];
		newWarning->issuedIDs[0] = 0;
	}
	newWarning->prevWarning = (void *)lastWarning;
	
	lastWarning = newWarning;
	numWarnings++;
	return numWarnings-1;
}

// assemble are created warnings into a list
void MPMWarnings::CreateWarningList(void)
{
	// no warnings
	if(numWarnings==0) return;
	
	// memory
	warningSet = new WarningsData *[numWarnings];
	
	WarningsData *nextWarning = lastWarning;
	int i = numWarnings-1;
	while(nextWarning!=NULL)
	{	warningSet[i] = nextWarning;
		nextWarning = (WarningsData *)(warningSet[i]->prevWarning);
		i--;
	}
}

#pragma mark MPMWarnings: Methods

// Activate all warnings for the next step
void MPMWarnings::BeginStep(void)
{	unsigned i;
	for(i=0;i<numWarnings;i++)
	{	warningSet[i]->thisStep = true;
	}
}

// issue a warning of a certain type
// but, only count once per step, and display message only on the first step this warning occurs
// return SILENT_WARNING, GAVE_WARNING, or REACHED_MAX_WARNINGS
// not thread safe due to push_back()
int MPMWarnings::Issue(int warnKind,int theID) { return Issue(warnKind,theID,NULL); }
int MPMWarnings::Issue(int warnKind,int theID,char *comment)
{	int warnResult=SILENT_WARNING;
	
	// exit if invalid warning ID
	if(warnKind<0 || warnKind>=numWarnings) return warnResult;
	
	// check the warning
	WarningsData *warn=warningSet[warnKind];
	
	// will output warning if first time or if first time for this ID for warning that use IDs
	bool newID = false, tooManyIDs = false;
	if(theID>0 && warn->issuedIDs!=NULL)
	{	int i=0;
		newID = true;
		while(warn->issuedIDs[i]!=0)
		{	if(theID==warn->issuedIDs[i])
			{	newID = false;
				break;
			}
			i++;
		}
		if(newID)
		{	// warning - when create warning, better save room for all possible IDs
			if(i<warn->maxIDs)
			{	warn->issuedIDs[i] = theID;
				warn->issuedIDs[i+1] = 0;
			}
			else
				tooManyIDs = true;
		}
	}
	
	// increment number of steps (first time on each step)
	if(warn->thisStep) warn->numSteps++;
    
    // increment total number of warnings
    warn->numWarnings++;
	
	// if first time in analysis, print message and force archiving
	if((warn->numSteps==1 && warn->thisStep) || newID)
	{	warn->firstStep=fmobj->mstep;
		warnResult=GAVE_WARNING;
		if(MPMTask::GetNumberOfThreads()>0)
		{
#pragma omp critical (output)
			{	archiver->ForceArchiving();
				cout << "# " << warn->msg;
				if(theID>0) cout << " (ID=" << theID << ")";
				if(tooManyIDs) cout << " (exceeded maximum allowed IDs)";
				cout << endl;
				if(comment!=NULL) cout << "#  " << comment << endl;
			}
		}
		else
		{	archiver->ForceArchiving();
			cout << "# " << warn->msg;
			if(theID>0) cout << " (ID=" << theID << ")";
			if(tooManyIDs) cout << " (exceeded maximum allowed IDs)";
			cout << endl;
			if(comment!=NULL) cout << "#  " << comment << endl;
		}
	}
	
	// deactivate for this step
	warn->thisStep = false;
	
	// abort if reach maximum issues or encoutered too many IDs
	if((warn->numWarnings>=warn->maxIssues && warn->maxIssues>0) || tooManyIDs)
		warnResult=REACHED_MAX_WARNINGS;
    
    // return result
	return warnResult;
}

// Report all warning results
void MPMWarnings::Report(void)
{
	unsigned i;
	bool hasTitle=FALSE;
	
	// check for any warning
	for(i=0;i<numWarnings;i++)
	{	WarningsData *warn=warningSet[i];
	
		// continue if never occurred
		if(warn->numSteps==0) continue;
	
		// print warning section title once
		if(!hasTitle)
		{	PrintSection("ANALYSIS WARNINGS");
			hasTitle=TRUE;
		}
		
		// this warning
		cout << "The warning '" << warn->msg << "' occurred on " << warn->numSteps << " time steps." << endl;
        cout << "   First such warning warning was on step " << warn->firstStep << endl;
        cout << "   Total number of warnings was " << warn->numWarnings << endl;
        cout << endl;
	}
}
