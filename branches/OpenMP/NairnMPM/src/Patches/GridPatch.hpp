/********************************************************************************
    GridPatch.hpp
    nairn-mpm-fea

    Created by John Nairn on 4/10/13.
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        none
********************************************************************************/

#ifndef _GRIDPATCH_

#define _GRIDPATCH_

class MPMBase;
class GhostNode;
class NodalPoint;

enum { FIRST_NONRIGID=0,FIRST_RIGID_CONTACT,FIRST_RIGID_BC };

class GridPatch
{
    public:
        static int ghostRows;
    
        // constructors and destructors
        GridPatch(int,int,int,int,int,int);
		bool CreateGhostNodes(void);
    
		// methods
		void InitializeForTimeStep(void);
		void InitializationReduction(void);
		void MassAndMomentumReduction(void);
        void MassAndMomentumReductionLast(void);
        void RezeroNodeTask6(double);
		void GridForcesReduction(void);
        void ZeroDisp(void);
        void JKTaskReduction(void);
        void DeleteDisp(void);
		void AddParticle(MPMBase *);
		void RemoveParticleAfter(MPMBase *,MPMBase *);
	
		// accessors
		MPMBase *GetFirstBlockPointer(int);
		NodalPoint *GetNodePointer(int);
        NodalPoint *GetNodePointer(int,bool);
	
    private:
        int x0,x1,y0,y1,z0,z1;					// element ranges (0-based row, col, rank)
		int xn,yn,zn;							// node count in each direction
		int numGhosts;							// total number of GhostNodes in 0 based list
		GhostNode **ghosts;						// 0 based points to the ghost nodes
        MPMBase *firstNR;                       // first non-rigid particle
        MPMBase *firstRC;                       // first rigid contact particle
        MPMBase *firstRBC;                      // first rigid BC particle
		int interiorRow;
		int fullRow;
		int basePartial;
		int baseTop;
		int interiorRank;
		int fullRank;
        int baseInterior;
        int baseApex;
};

extern GridPatch **patches;

#endif
