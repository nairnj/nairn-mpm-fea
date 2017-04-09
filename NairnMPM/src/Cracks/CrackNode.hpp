/********************************************************************************
    CrackNode.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 2/24/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.

	Dependencies
		NodalPoint.hpp
********************************************************************************/

#ifndef _CRACKNODE_

#define _CRACKNODE_

#include "Nodes/NodalPoint.hpp"

class CrackNode
{
    public:
		// variables (changed in MPM time step)
		static CrackNode *currentCNode;

        // constructors and destructors
        CrackNode(NodalPoint *,CrackNode *);
	
		// common methods
		void SetPrevNode(CrackNode *);
		CrackNode *GetPrevNode(void);
	
		// interface and contact methods
		virtual CrackNode *NodalCrackContact(void);
		virtual CrackNode *NodalCrackContactAndForces(double);

		// class methods
		static void RemoveCrackNodes(void);
		static void CrackContactTask4(double);
		static void ContactOnKnownNodes(void);
		
	protected:
		// variables (changed in MPM time step)
		NodalPoint *theNode;
        CrackNode *prevNode;
};

#endif


