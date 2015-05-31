/********************************************************************************
    CrackSegment.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Apr3 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		CrackHeader.hpp
********************************************************************************/

#ifndef _CRACKSEGMENT_

#define _CRACKSEGMENT_

#include "Cracks/CrackHeader.hpp"

class CrackLeaf;

enum { STATIONARY=0, PROPAGATING, ARRESTING, ARRESTED };

// Crack Segment Class
class CrackSegment
{
    public:
        double x,y,origx,origy;
        int planeInElem,surfInElem[2];
        double surfx[2],surfy[2];
        CrackSegment *nextSeg,*prevSeg;
        Vector Jint,sif,tract;
        int tipMatnum;					// crack tip material (not a traction law)
        int steadyState;
        double speed,theGrowth;
		double propagationJ;
		short heating;
		double heatRate,heatEndTime;
        double cnear[4],cfar[4];
        CrackLeaf *parent;
    
        // constructors
        CrackSegment();
        CrackSegment(double,double,int,int);
        
        // methods
        void FillArchive(char *,int);
        int FindElement(void);
        int FindElement(short);
        void MovePosition(double,double);
        void MovePosition(void);
		void MovePositionToMidpoint(void);
        bool MoveSurfacePosition(short,double,double,bool);
		int CheckSurfaces(void);
		bool MoveToPlane(int,double,double,bool,double);
		void CollapseSurfaces(void);
		void StartCrackTipHeating(double,double);
		double HeatRate(void);
		double ForwardArea(double,double,int);
		Vector SlightlyMovedIfNotMovedYet(int);
		int MatID(void);
		void SetMatID(int);
		void AddTractionForceSeg(CrackHeader *);
		double AddTractionForceSegSide(CrackHeader *,int,double);
		void FindCrackTipMaterial(int);
		void UpdateTractions(CrackHeader *);
		Vector GetTangential(double *);
		double TractionEnergy(Vector *,int,bool,CrackSegment **);
		double SegmentTractionEnergy(bool);
		Vector FTract(double);
		void SetHistoryData(char *p);
		char *GetHistoryData(void);
        void CreateSegmentExtents(bool);
				
	private:
		Vector cFtract;				// traction law force
		double dxPlane,dyPlane;
		bool hadAboveNodes;
		bool planeMove;
		int matnum;						// 1-based material ID for traction law
        char *historyData;				// history dependent traction law data
    
};

#endif

