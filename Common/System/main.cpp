/*********************************************************************
    main.cpp
    Nairn Research Group MPM and FEA Code
    
    Created by John Nairn on Mon Nov 19 2001.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
*********************************************************************/

#include "stdafx.h"
#ifdef MPM_CODE
	#include "NairnMPM_Class/NairnMPM.hpp"
#else
	#include "NairnFEA_Class/NairnFEA.hpp"
#endif
#include "Exceptions/CommonException.hpp"

/*********************************************************************
    Main entry point for FEA and MPM code
		1. Create main analysis object
		2. Read command line parameters
		3. Read input file (path in last parameter)
		4. Do the analysis
*********************************************************************/

int main(int argc,const char *argv[])
{
    int retval,parmInd;
	unsigned int optInd;
	bool abort=FALSE;
	bool useWorkingDir=FALSE;
    int numProcs=1;
	
	// ---------------------------------------------
    // 1. Create main analysis object
#ifdef MPM_CODE
	fmobj=new NairnMPM();
#else
	fmobj=new NairnFEA();
#endif
	
	// Initialize things for random numbers
	InitRandom(0);
    
	// ---------------------------------------------
    // 2. Check command line and extract arguments.
    if(argc<2)
    {	fmobj->Usage();
        return NoInputErr;
    }
	
    // Check for options
    for(parmInd=1;parmInd<argc && argv[parmInd][0]=='-';parmInd++)
	{	// each option in the argument
        unsigned arglen = (int)strlen(argv[parmInd]);
		for(optInd=1;optInd<arglen;optInd++)
		{	// Help request
			if(argv[parmInd][optInd]=='H')
			{	fmobj->Usage();
				return(noErr);
			}
			else if(argv[parmInd][optInd]=='h')
			{	fmobj->Usage();
				return(noErr);
			}

#ifdef MPM_CODE
			else if(argv[parmInd][optInd]=='r')
				fmobj->SetReverseBytes(TRUE);
#endif
				
			else if(argv[parmInd][optInd]=='v')
				fmobj->SetValidate((bool)TRUE);
			
			else if(argv[parmInd][optInd]=='a')
				abort=TRUE;
				
			else if(argv[parmInd][optInd]=='w')
				useWorkingDir=TRUE;
			
			else if(argv[parmInd][optInd]=='n' && optInd==arglen-2 && argv[parmInd][optInd+1]=='p' && parmInd<argc-2)
            {   // get number in next argument
                parmInd++;
                sscanf(argv[parmInd],"%d",&numProcs);
                break;
            }
 			
			else
			{   cerr << "\nUnknown " << fmobj->CodeName() << " option '" << argv[parmInd][optInd]
					 << "' was used.\n";
				return NoInputErr;
			}
		}
    }
    
    //  Last parameter must be the file name
    if(parmInd+1!=argc)
    {	fmobj->Usage();
        return NoInputErr;
    }
    
    // set number of processors (default 1)
#ifdef _OPENMP
	// pick number of processors, but no more than number available
    // if set to 0, will set to maximum number
	int maxThreads = omp_get_max_threads();
	int numProcsAvail = omp_get_num_procs();
    int maxProcs = maxThreads>numProcsAvail ? maxThreads : numProcsAvail;
    if(numProcs>0)
    {   if(numProcs > maxProcs) numProcs = maxProcs;
        omp_set_num_threads(numProcs);
    }
    else
        numProcs = maxProcs;
#else
	// zero means omp is disabled
    numProcs = 0;
#endif
    fmobj->SetNumberOfProcessors(numProcs);
	
	//-------------------------------------------------------------
    // 3. Read the input file, exceptions handled in ReadFile()
	retval=fmobj->ReadFile(argv[parmInd],useWorkingDir);
    if(retval!=noErr) return retval;
	
	// ------------------------------------------------------------
    // 4. Main analysis block
	//		Catches CommonException(), char * exceptions, std::bad__alloc,
	//		exception, anything else
    try
	{	// Start analysis
		fmobj->StartAnalysis(abort);
    }
    
    catch(CommonException& err)
	{   cout << "Warning: pipe timing may prevent error details from appearing below" << endl;
#ifdef MPM_CODE
    	err.Display(fmobj->mstep,mtime);
#else
    	err.Display();
#endif
        return err.ErrorCode();
    }
	
    catch(CommonException* err)
    {   cout << "Warning: pipe timing may prevent error details from appearing below" << endl;
#ifdef MPM_CODE
        err->Display(fmobj->mstep,mtime);
#else
        err->Display();
#endif
        return err->ErrorCode();
    }
	catch(const char *errMsg)
    {   // send to output results and error pipe
        cout << "\n" << errMsg << endl;
		cerr << "\n" << errMsg << endl;
		return AnalysisErr;
	}
    
	catch(std::bad_alloc& ba)
	{	// memory error most likely in new command
        cout << "Memory error: " << ba.what() << endl;
        cerr << "Memory error: " << ba.what() << endl;
		return AnalysisErr;
	}
	
    catch(exception& e)
    {   // send to output results and error pipe
        cout << "Standard exception: " << e.what() << endl;
        cerr << "Standard exception: " << e.what() << endl;
		return AnalysisErr;
    }
	
	catch(...)
    {   // send to output results and error pipe
        cout << "\nMain analysis block in main exited with unknown exception" << endl;
		cerr << "\nMain analysis block in main exited with unknown exception" << endl;
		return AnalysisErr;
	}
    
    return noErr;
}
