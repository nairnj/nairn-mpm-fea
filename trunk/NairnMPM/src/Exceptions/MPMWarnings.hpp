/********************************************************************************
    MPMWarnings.hpp
    NairnMPM
    
    Created by John Nairn on Tues Feb 08 2005.
    Copyright (c) 2005 John A. Nairn, All rights reserved.    

	Dependencies
		none
********************************************************************************/

#ifndef _MPMWARNINGS_

#define _MPMWARNINGS_

enum { SILENT_WARNING=0, GAVE_WARNING, REACHED_MAX_WARNINGS };

// WarningsData structure
typedef struct {
	// variables (changed in MPM time step)
	short thisStep;				// was there a warning this step yet?
	long numSteps;				// how many steps had this warning
	vector< int > issuedIDs;	// list of IDS for this warning
	long firstStep;				// first step with this warning
	
	// constants (not changed in MPM time step)
	long maxIssues;				// abort if reach this number of steps with this warning
	char *msg;					// warning message
} WarningsData;

class MPMWarnings
{
    public:
        // constructors and destructors
        MPMWarnings();
		int CreateWarning(const char *,long);
		
		// methods
		void BeginStep(void);
		int Issue(int,int);
		void Report(void);
	
	private:
		// constants (not changed in MPM time step)
		vector< WarningsData * > warningSet;
};

extern MPMWarnings warnings;
    
#endif
