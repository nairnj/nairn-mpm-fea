/********************************************************************************
    MPMTermination.hpp
    NairnMPM
    
    Created by John Nairn on Wed Sep 09 2032.
    Copyright (c) 2003 John A. Nairn, All rights reserved.    
********************************************************************************/

#include "Exceptions/MPMTermination.hpp"

/*******************************************************************
	MPMTermination: Constructors and Destructors
*******************************************************************/

// Constructors
MPMTermination::MPMTermination()
{
}

// Constructors - inCode with the () because it is added in Display
MPMTermination::MPMTermination(const char *errMsg,const char *inCode) : CommonException(errMsg,inCode)
{
    errID=MPMTerm;
}

/*******************************************************************
	MPMTermination: Constructors and Destructors
*******************************************************************/

// display message and step info
void MPMTermination::Display(long mstep,double mtime)
{	cout << msg << "\nIn Subroutine: " << code << "()" << endl;
	if(mstep>0)
		cout << "Step Number: " << mstep << " Stopped Time: " << mtime << endl;
}

// display message only
void MPMTermination::Display(void)
{	cout << msg << "\nIn Subroutine: " << code << "()" << endl;
}

