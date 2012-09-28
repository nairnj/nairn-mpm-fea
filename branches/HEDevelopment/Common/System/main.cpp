/*********************************************************************
    main.cpp
    Nairn Research Group MPM and FEA Code
    
    Created by John Nairn on Mon Nov 19 2001.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
*********************************************************************/

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
	
	// ---------------------------------------------
    // 1. Create main analysis object
#ifdef MPM_CODE
	fmobj=new NairnMPM();
#else
	fmobj=new NairnFEA();
#endif
    
	// ---------------------------------------------
    // 2. Check command line and extract arguments.
    if(argc<2)
    {	fmobj->Usage();
        return NoInputErr;
    }
	
    // Check for options
    for(parmInd=1;parmInd<argc && argv[parmInd][0]=='-';parmInd++)
	{	// each option in the argument
		for(optInd=1;optInd<strlen(argv[parmInd]);optInd++)
		{	// Help request
			if(argv[parmInd][optInd]=='H')
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
    
	//-------------------------------------------------------------
    // 3. Read the input file, exceptions handled in ReadFile()
	retval=fmobj->ReadFile(argv[parmInd]);
    if(retval!=noErr) return retval;
	
	// ------------------------------------------------------------
    // 4. Main analysis block
	//		Catches CommonException() cand char * exceptions
    try
	{	// Start analysis
		fmobj->StartAnalysis(abort);
    }
    
    catch(CommonException err)
	{   
#ifdef MPM_CODE
    	err.Display(fmobj->mstep,mtime);
#else
    	err.Display();
#endif
        return err.ErrorCode();
    }
	
	catch(const char *errMsg)
	{	cerr << "\n" << errMsg << endl;
		return AnalysisErr;
	}
	
	catch(...)
	{	cerr << "\nStartAnalysis() block in main exited with unknown exception" << endl;
		return AnalysisErr;
	}
    
    return noErr;
}

