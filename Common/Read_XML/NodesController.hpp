/********************************************************************************
    NodesController.hpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		ParseController.hpp
********************************************************************************/

#ifndef _NODESCONTROLLER_

#define _NODESCONTROLLER_

#include "Read_XML/ParseController.hpp"

class NodesController : public ParseController
{
    public:
	
		// methods
#ifdef MPM_CODE
		void AddNode(double,double,double);
#else
		void AddNode(double,double,double,double);
#endif
		int SetNodeArray(double *,double *,double *,double *,double *,double *);
		int SetNodeArray(int *);
		int NextNodeNumber(void);
#ifdef FEA_CODE
		void MidPoint(int *,int,Vector *);
#endif
};

extern NodesController *theNodes;

#endif
