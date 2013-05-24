/********************************************************************************
    CommonException.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 07 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.    

	Dependencies
		none
********************************************************************************/

#ifndef _COMMONEXCEPTION_

#define _COMMONEXCEPTION_

// return error codes
enum { noErr=0,NoInputErr,XercesInitErr,ReadFileErr,MPMErr,MPMTerm,AnalysisErr };

class CommonException
{
    public:
        
        // constructors and destructors
        CommonException();
        CommonException(const char *,const char *);
        
        // methods
        void Display(int,double);
        void Display(void);
        int ErrorCode(void);
		void SetErrorCode(int);
		char *Message(void);
		char *InCode(void);
	
    protected:
        char *msg,*code;
        int errID;
};

#endif
