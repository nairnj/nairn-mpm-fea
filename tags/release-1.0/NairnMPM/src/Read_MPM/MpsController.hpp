/********************************************************************************
    MpsController.hpp
    NairnFEA
    
    Created by John Nairn on 6/27/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		ParseController.hpp
********************************************************************************/

#ifndef _MPSCONTROLLER_

#define _MPSCONTROLLER_

#include "Read_XML/ParseController.hpp"

class MPMBase;

class MpsController : public ParseController
{
    public:
	
		// methods
		void AddMaterialPoint(MPMBase *,double,double);
		int SetPtOrVel(char *,Vector *);
		void SetPtMass(double);
		int SetMPArray(void);
};

extern MpsController *mpCtrl;

#endif
