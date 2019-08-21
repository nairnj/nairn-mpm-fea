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

class NodalPoint;

class CrackNode
{
    public:
		// list of crack nodes
		static vector< CrackNode * > crackContactNodes;

        // constructors and destructors
		CrackNode(NodalPoint *,int,CrackNode *);
	
		// common methods
		void SetPrevNode(CrackNode *);
		CrackNode *GetPrevNode(void);
	
		// interface and contact methods
		virtual void NodalCrackContact(double,int);

		// class methods
		static bool ContactOnKnownNodes(double,int);
		static void ReleaseContactNodes(void);
	
	protected:
		// variables (changed in MPM time step)
		NodalPoint *theNode;
        CrackNode *prevNode;
		int hasFlags;
};

#endif


