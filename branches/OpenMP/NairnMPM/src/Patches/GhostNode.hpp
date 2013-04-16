/********************************************************************************
	GridNode.hpp
	NairnMPM

	Created by John Nairn on 4/11/13.
	Copyright (c) 2013 John A. Nairn, All rights reserved.

	Dependencies
		none
********************************************************************************/

#ifndef _GHOSTNODE_

#define _GHOSTNODE_

class NodalPoint;

class GhostNode
{
	public:
    
		// constructors and destructors
		GhostNode(int,int,bool,bool,bool);
	
		// methods
		void InitializeForTimeStep();
		void InitializationReduction(void);
		void MassAndMomentumReduction(void);
		void GridForcesReduction(void);
	
		// accessors
		NodalPoint *GetNodePointer(void);
	
	private:
		NodalPoint *ghost;
		NodalPoint *real;
};

#endif
