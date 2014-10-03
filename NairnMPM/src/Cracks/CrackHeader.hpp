/********************************************************************************
    CrackHeader.hpp
    nairn-mpm-fea
    
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

/* New method for j integral to handle multiple cracks
   Search for GRID_JTERMS to see changes to allow grid options
*/
#define MCJ_INTEGRAL
#define MCJ_HIERCONTOURCROSS

// Debugging
//#define CONTOUR_PARTS
//#define PRINT_CROSS_STATUS

class CrackSegment;
class ContourPoint;
class CrackLeaf;
class NodalPoint;

#define START_OF_CRACK 0
#define END_OF_CRACK 1
#define EXTERIOR_CRACK -2

class CrackHeader : public LinkedObject
{
    public:
        CrackSegment *firstSeg,*lastSeg;
		static double codLocation;				// 0 to 2
		static int codInterval;					// 0 or 1
		static double bezArg[4];
		static double bezDer[4];
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
        short MoveCrack(void);
        short MoveCrack(short);
		void UpdateCrackTractions(void);
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
		bool NodeNearTip(NodalPoint *,double);

        bool CreateHierarchy(void);
        void MoveHierarchy(void);
        void ExtendHierarchy(CrackSegment *);
        short CrackCross(double,double,double,double,Vector *) const;
        short CrackCrossLeaf(CrackLeaf *,double,double,double,double,Vector *,short) const;
        short CrackCrossOneSegment(CrackSegment *,double,double,double,double,Vector *,short) const;
        CrackSegment *ContourCrossCrack(ContourPoint *,Vector *) const;
#ifdef MCJ_HIERCONTOURCROSS
        CrackSegment *ContourCrossLeaf(CrackLeaf *,double,double,double,double,Vector *,int) const;
        bool SegmentsCross(CrackSegment *,double,double,double,double,Vector *,int) const;
#else
        bool SegmentsCross(ContourPoint *,Vector &,Vector &,Vector *) const;
#endif
    
        // These two are not used, but left here for testing of hierarchical cracks
        short FlatCrackCrossTest(double,double,double,double,Vector *);
        void CFFlatCrossing(double,double,double,double,Vector *,short *,int,int);
    
		void SetNumber(int);
		int GetNumber(void);
		void SetThickness(double);
		double GetThickness(void);
		double *GetThicknessPtr(void);
    
        // calculate J-integral (YJG)
        void JIntegral(void);      	  // J-Integral calculation
        void PrintContour(ContourPoint *,ContourPoint *,Vector &);
		CrackSegment *GetCrackTip(int);
		CrackSegment *GetAdjToCrackTip(int);
		int GetWhichTip(CrackSegment *);
		int CriterionNeeds(bool &);
		void GetCOD(CrackSegment *,Vector &,bool);
		void CrackTipAndDirection(int,CrackSegment **,Vector &);
		void AddTractionForce(void);
		bool GetHasTractionLaws(void);
		void GetInitialDirection(CrackSegment *,Vector &);
		void InterpolatePosition(int,CrackSegment **,Vector &,bool);
		bool GetAllowAlternate(int);
		void SetAllowAlternate(int,bool);
       
        // propagation
        CrackSegment *Propagate(Vector &,int,int);
		
		// class methods
		static void SetCodLocation(double);
        static double Triangle(double,double,double,double,double,double);
        static bool LineIsInExtents(double,double,double,double,double *,double *);

    private:
        int numberSegments;
		int fixedCrack,number;
		double crackFriction,crackDn,crackDnc,crackDt;
		bool customContact,hasTractionLaws;
		Vector initialDirection[2];
		bool allowAlternate[2];
		double thickness;						// 2D tractions and crack-tip heating
        CrackLeaf *rootLeaf;
        
};

extern CrackHeader *firstCrack;
// GRID_JTERMS
extern int JGridSize,JContourType,JTerms,JGridEnergy;

extern CrackHeader **crackList;
extern int numberOfCracks;


#endif
