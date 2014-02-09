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
		
		// methods
		CrackNode *NodalContact(void);
		CrackNode *NodalContactAndForces(double);
		void SetPrevBC(CrackNode *);
		CrackNode *GetPrevBC(void);
		CrackNode *InterfaceForce(void);
        
		// class methods
		static void RemoveCrackNodes(void);
		static void CrackContactTask4(double);
		static void ContactOnKnownNodes(void);
		static void InterfaceOnKnownNodes(void);
		
	private:
		// variables (changed in MPM time step)
        Vector pk[MAX_FIELDS_FOR_CRACKS];
		NodalPoint *theNode;
        CrackNode *prevBC;
};

#endif


