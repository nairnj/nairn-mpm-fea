/********************************************************************************
    CrackSegment.hpp
    NairnMPM
    
    Created by John Nairn on Wed Apr3 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _CRACKSEGMENT_

#define _CRACKSEGMENT_

class CrackHeader;

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
        
        // constructors
        CrackSegment();
        CrackSegment(double,double,int,int);
        
        // methods
        void FillArchive(char *,long);
        int FindElement(void);
        int FindElement(short);
        void MovePosition(double,double);
        void MovePosition(void);
        void MovePosition(short,double,double,bool,double);
		int CheckSurfaces(void);
		bool MoveToPlane(int,double,double,bool,double);
		void StartCrackTipHeating(double,double);
		double HeatRate(void);
		double ForwardArea(double,double,int);
		Vector SlightlyMoved(int);
		int MatID(void);
		void SetMatID(int);
		void AddTractionFext(CrackHeader *);
		void AddTractionFext(CrackHeader *,int,double);
		void UpdateTractions(CrackHeader *);
		Vector GetTangential(double *);
		double TractionEnergy(Vector *,int,bool);
		double SegmentTractionEnergy(bool);
		Vector FTract(double);
		void SetHistoryData(char *p);
		char *GetHistoryData(void);
				
	private:
		Vector cFtract;				// traction law force
		double dxPlane,dyPlane,aboveMass;
		bool hadAboveNodes;
		bool planeMove;
		int matnum;						// 1-based material ID for traction law
        char *historyData;				// history dependent traction law data
};

#endif

