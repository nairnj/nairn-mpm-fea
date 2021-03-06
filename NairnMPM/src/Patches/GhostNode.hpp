/********************************************************************************
	GridNode.hpp
	nairn-mpm-fea

	Created by John Nairn on 4/11/13.
	Copyright (c) 2013 John A. Nairn, All rights reserved.

	Dependencies
		none
********************************************************************************/

#ifndef _GHOSTNODE_

#define _GHOSTNODE_

class NodalPoint;
class TransportTask;

class GhostNode
{
	public:
    
		// constructors and destructors
		GhostNode(int,int,bool,bool);
		GhostNode(int,int,int,bool,bool,bool);
	
		// methods
		void InitializeForTimeStep();
		void InitializationReduction(void);
		void MassAndMomentumReduction(void);
        void MassAndMomentumReductionLast(void);
        void RezeroNodeTask6(double);
		void GridForcesReduction(void);
        void ZeroDisp(void);
        void JKTaskReduction(void);
        void DeleteDisp(void);
		void XPICSupport(int,int,NodalPoint *,double,int,int,double);
        void InitializeForXPICTransport(TransportTask *,int);
        void XPICReductionTransport(int,TransportTask *);

		// accessors
		NodalPoint *GetNodePointer(void);
		NodalPoint *GetGhostNodePointer(void);
		NodalPoint *GetRealNodePointer(void);
	
	private:
		NodalPoint *ghost;
		NodalPoint *real;
};

#endif
