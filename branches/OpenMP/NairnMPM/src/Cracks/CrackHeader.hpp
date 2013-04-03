/********************************************************************************
    CrackHeader.hpp
    NairnMPM
    
    Created by John Nairn on Wed Apr3 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _CRACKHEADER_

#define _CRACKHEADER_

/* Define method to find crack tip direction and COD.
	If _LINEAR_INTERPOLATION_ is defined then finds crack direction from direction
		of last crack segment and finds COD at adjacent segment. To get normalized
		cod compares it to crack direction
	If _LINEAR_INTERPOLATION_ is not defined will use cubic splines and there are
		two options. Defining _CUBIC_INTERPOLATION_ is relaxed cubic spline through
		the last 4 crack points. COD can be taken anywhere from tip to 3 point.
		Normalized COD is from tangent direction to spline at that point.
		If _CUBIC_INTERPOLATION_ is not defined then draws 3 Bezier curves
		with crack points as control points. COD is similar by from different
		curves.
*/
//#define _CUBIC_INTERPOLATION_
//#define _LINEAR_INTERPOLATION_

class CrackSegment;
class ContourPoint;
class CrackLeaf;

#define START_OF_CRACK 0
#define END_OF_CRACK 1
#define EXTERIOR_CRACK -2

/* The HIERARCHICAL_CRACKS constants determines the method used to screen
    cracks for intersections in the CrackCross() method and eliminate those
    that do not intersect crack extents, which are tracked.
        If defined, cracks are hierarchical structures with each branch having
            its own extents
        If not defined, track one global extent. See EXTENT_NORMALS to pick
            the shape of the crack bounding region
*/
#define HIERARCHICAL_CRACKS

// Determines full crack bounding region. Must be 2, 4, or 6, for box, octagon, dodecahedron
// Setting is irrelevant when HIERARCHICAL_CRACKS is defined
#define EXTENT_NORMALS 4

class CrackHeader : public LinkedObject
{
    public:
        CrackSegment *firstSeg,*lastSeg;
		static double codLocation;				// 0 to 2
		static int codInterval;					// 0 or 1
		static double bezArg[4];
		static double bezDer[4];
		static int warnThreeFields;
		static int warnNodeOnCrack;
		static int warnThreeCracks;
		
        // constructors and destructors
        CrackHeader();
        ~CrackHeader();
		void PreliminaryCrackCalcs(void);
        
        // methods
        short add(CrackSegment *);
        short add(CrackSegment *, int);
        void Archive(ofstream &);
        short CrackCross(double,double,double,double,Vector *);
        short MoveCrack(void);
        short MoveCrack(short);
		void UpdateCrackTractions(void);
        int Count(void);
        double Length(void);
        int NumberOfSegments(void);
		void CrackTipHeating(void);
		void SetFixedCrack(int);
		void SetFriction(double);
		void SetDn(double);
		void SetDnc(double);
		void SetDt(double);
		void SetContact(double,double,double,double);
		void Output(void);
		void Describe(void);

#ifdef HIERARCHICAL_CRACKS
        void CreateHierarchy(void);
        void MoveHierarchy(void);
        void ExtendHierarchy(CrackSegment *);
        short CrackCrossLeaf(CrackLeaf *,double,double,double,double,Vector *,short);
        short CrackCrossOneSegment(CrackSegment *,double,double,double,double,Vector *,short);
        short FlatCrackCross(double,double,double,double,Vector *);
        void CFFlatCrossing(double,double,double,double,Vector *,short *,int,int);
#else
        void CreateExtents(double,double);
        void CheckExtents(double,double);
#endif
        
		void SetNumber(int);
		int GetNumber(void);
		void SetThickness(double);
		double GetThickness(void);
		double *GetThicknessPtr(void);
    
        // calculate J-integral (YJG)
        int inMat;                        // material ID of crack
        bool SegmentsCross(ContourPoint *,Vector &,Vector &,Vector *);
        void JIntegral(void);      	  // J-Integral calculation
		CrackSegment *GetCrackTip(int);
		CrackSegment *GetAdjToCrackTip(int);
		int GetWhichTip(CrackSegment *);
		int CriterionNeeds(void);
		void GetCOD(CrackSegment *,Vector &,bool);
		void CrackTipAndDirection(int,CrackSegment **,Vector &);
		void TractionFext(void);
		bool GetHasTractionLaws(void);
		void GetInitialDirection(CrackSegment *,Vector &);
		void InterpolatePosition(int,CrackSegment **,Vector &,bool);
		bool GetAllowAlternate(int);
		void SetAllowAlternate(int,bool);
       
        // propagation
        CrackSegment *Propagate(Vector &,int,int);
		
		// class methods
		static void ContactConditions(int);
		static void SetCodLocation(double);
        static double Triangle(double,double,double,double,double,double);

    private:
        int numberSegments;
		int fixedCrack,number;
		double crackFriction,crackDn,crackDnc,crackDt;
		bool customContact,hasTractionLaws;
		Vector initialDirection[2];
		bool allowAlternate[2];
		double thickness;						// 2D tractions and crack-tip heating

#ifdef HIERARCHICAL_CRACKS
        CrackLeaf *rootLeaf;
#else
        double cnear[EXTENT_NORMALS],cfar[EXTENT_NORMALS];
#endif
        
};

extern CrackHeader *firstCrack;
extern int JGridSize,JContourType,JTerms;

#endif
