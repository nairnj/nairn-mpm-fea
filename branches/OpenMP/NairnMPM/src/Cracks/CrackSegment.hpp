/********************************************************************************
    CrackSegment.hpp
    NairnMPM
    
    Created by John Nairn on Wed Apr3 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		CrackHeader.hpp
********************************************************************************/

#ifndef _CRACKSEGMENT_

#define _CRACKSEGMENT_

#include "Cracks/CrackHeader.hpp"

class CrackLeaf;

enum { STATIONARY=0, PROPAGATING, ARRESTING, ARRESTED, SLOWLYPROPAGATING,
        BEGINPROPAGATING};

// Crack Segment Class
class CrackSegment
{
    public:
        double x,y,origx,origy;
        double surfx[2],surfy[2];
        int planeInElem,surfInElem[2];
        CrackSegment *nextSeg,*prevSeg;
        Vector Jint,sif,tract;
        int tipMatnum;					// crack tip material (not a traction law)
        double potential[3],plastic[3],clength[3];
        int steadyState;
        double speed,theGrowth;
		int crackIncrements;
		double release,absorb,propagationJ;
		short heating;
		double heatRate,heatEndTime;

#ifdef HIERARCHICAL_CRACKS
        double cnear[4],cfar[4];
        CrackLeaf *parent;
#endif
    
        // constructors
        CrackSegment();
        CrackSegment(double,double,int,int);
        
        // methods
        void FillArchive(char *,int);
        int FindElement(void);
        int FindElement(short);
        void MovePosition(double,double);
        void MovePosition(void);
        bool MoveSurfacePosition(short,double,double,bool,double);
		int CheckSurfaces(void);
		bool MoveToPlane(int,double,double,bool,double);
		void CollapseSurfaces(void);
		void StartCrackTipHeating(double,double);
		double HeatRate(void);
		double ForwardArea(double,double,int);
		Vector SlightlyMoved(int);
		int MatID(void);
		void SetMatID(int);
		void AddTractionForceSeg(CrackHeader *);
		double AddTractionForceSegSide(CrackHeader *,int,double);
		void FindCrackTipMaterial(void);
		void UpdateTractions(CrackHeader *);
		Vector GetTangential(double *);
		double TractionEnergy(Vector *,int,bool);
		double SegmentTractionEnergy(bool);
		Vector FTract(double);
		void SetHistoryData(char *p);
		char *GetHistoryData(void);

#ifdef HIERARCHICAL_CRACKS
        void CreateSegmentExtents(bool);
#endif
				
	private:
		Vector cFtract;				// traction law force
		double dxPlane,dyPlane,aboveMass;
		bool hadAboveNodes;
		bool planeMove;
		int matnum;						// 1-based material ID for traction law
        char *historyData;				// history dependent traction law data
    
};

#endif

