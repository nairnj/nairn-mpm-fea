/********************************************************************************
    CommonException.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 07 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
 
    On error, code should throw CommonException (not new CommonException)
    Catch blocks can catch(CommonException& err)
    Only main and NairnMPM had catch(CommonException* err) because many
        MPMTasks might through pointer instead of the the object itself
********************************************************************************/

#include "stdafx.h"
#include "Exceptions/CommonException.hpp"

#pragma mark CommonException: Constructors and Destructors

// Constructors
CommonException::CommonException()
{
}

// Copy Constructor
// new CommonException(CommonException err) is implied by compiler

// Constructors - inCode without the  "()" because it is added in Display
CommonException::CommonException(const char *errMsg,const char *inCode)
{
    msg=new char[strlen(errMsg)+1];
    strcpy(msg,errMsg);
    
    code=new char[strlen(inCode)+1];
    strcpy(code,inCode);
    
    errID=MPMErr;
}

// Constructors - inCode without the  "()" because it is added in Display and give error ID
CommonException::CommonException(const char *errMsg,const char *inCode,int errCode)
{
	msg=new char[strlen(errMsg)+1];
	strcpy(msg,errMsg);
	
	code=new char[strlen(inCode)+1];
	strcpy(code,inCode);
	
	errID=errCode;
}

#pragma mark CommonException: Methods

// display message and step info
void CommonException::Display(int mstep,double mtime)
{   // in both results window and in error pipe
	cout << msg << "\nIn Subroutine: " << code << "()" << endl;
    if(mstep>0)
        cout << "Step Number: " << mstep << " Current Time: " << mtime << endl;
	
	if(errID!=noErr)
	{	cerr << "\n" << msg << "\nIn Subroutine: " << code << "()" << endl;
		if(mstep>0)
			cerr << "Step Number: " << mstep << " Current Time: " << mtime << endl;
	}
}

// display message only
void CommonException::Display(void)
{	// in both results window and in error pipe
    cout << msg << "\nIn Subroutine: " << code << "()" << endl;
	if(errID!=noErr)
    	cerr << "\n" << msg << "\nIn Subroutine: " << code << "()" << endl;
}

// get the error code
int CommonException::ErrorCode(void) { return errID; }
void CommonException::SetErrorCode(int err) { errID = err; }

// get just the message
char *CommonException::Message(void) { return msg; }

// get just the message
char *CommonException::InCode(void) { return code; }

