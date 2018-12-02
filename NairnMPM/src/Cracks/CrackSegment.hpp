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
		Vector cp,surf[2],orig;			// crack plane, surface, and original position
		Vector cpVel,surfVel[2];		// crack plane and surface velocities
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
        CrackSegment(Vector *,int,int);
	
		// initialize
		void FindInitialElement(void);
	
        // methods
        void FillArchive(char *,int);
        int FindElement(void);
        int FindElement(short);
        virtual void MovePosition(void);
		virtual void MovePositionToMidpoint(void);
		void MovePosition(Vector *,Vector *,double,double);
		bool MoveSurfacePosition(short,Vector *,Vector *,double,bool);
		int CheckSurfaces(void);
		bool MoveToPlane(int,double,double,bool,double);
		bool CollapseSurfaces(void);
		void StartCrackTipHeating(double,double);
		double HeatRate(void);
		virtual Vector SlightlyMovedIfNotMovedYet(int);
		int MatID(void);
		void SetMatID(int);
		void AddTractionForceSeg(CrackHeader *);
		void AddTractionForceSegSide(CrackHeader *,int,double);
		void FindCrackTipMaterial(int);
		virtual void UpdateTractions(CrackHeader *);
		virtual double GetNormalAndTangent(CrackHeader *,Vector *,Vector *,double &,double &) const;
		double TractionEnergy(Vector *,int,bool,CrackSegment **);
		double SegmentTractionEnergy(bool);
		virtual Vector FTract(double);
		void SetHistoryData(char *p);
		char *GetHistoryData(void);
        void CreateSegmentExtents(bool);
		int planeElemID(void) const;
		int surfaceElemID(int side) const;
	
	protected:
		Vector cFtract;				// traction law force
		Vector dPlane;			// crack plane movement
		bool hadAboveNodes;
		bool planeMove;
		int matnum;						// 1-based material ID for traction law
        char *historyData;				// history dependent traction law data
		int planeInElem,surfInElem[2];	// 1-based element numbers
	
};

#endif

