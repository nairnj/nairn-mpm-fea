/********************************************************************************
    CrackController.hpp
    NairnFEA
    
    Created by John Nairn on 6/27/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		ParseController.hpp
********************************************************************************/

#ifndef _CRACKCONTROLLER_

#define _CRACKCONTROLLER_

#include "Read_XML/ParseController.hpp"

class CrackHeader;
class CrackSegment;

class CrackController : public ParseController
{
    public:
	
		// methods
		void AddCrack(CrackHeader *);
		int AddSegment(CrackSegment *);
        void FinishCrack(void);
};

extern CrackController *crackCtrl;

#endif
