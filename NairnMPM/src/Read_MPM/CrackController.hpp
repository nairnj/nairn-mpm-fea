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
class NodalPoint;
class ElementBase;

class CrackController : public ParseController
{
    public:
		// contructor
		CrackController();
	
		// methods
		void AddCrack(CrackHeader *,bool);
		int AddSegment(CrackSegment *,bool);
		bool FinishCrack(void);
        int SetCracksArray(void);
    
    protected:

};

extern CrackController *crackCtrl;

#endif
