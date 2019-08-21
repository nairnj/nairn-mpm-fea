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
		with crack points as control points. COD is similar but from different
		curves.
*/
//#define _CUBIC_INTERPOLATION_
//#define _LINEAR_INTERPOLATION_

// Debugging
//#define CONTOUR_PARTS
//#define PRINT_CROSS_STATUS
//#define JTERM_SUMMARY
#define JDEBUG_STEP 100

// to do axisymmetric J by Broberg method (any other number used Bergkvist and Huong method)
#define AXISYM_BROBERG_J 1

class CrackSegment;
class ContourPoint;
class CrackLeaf;
class NodalPoint;
class ParseController;

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
		virtual ~CrackHeader();
	
		// crack setup
		bool addSegment(CrackSegment *,bool);
		bool addSegmentTip(CrackSegment *, int);
		virtual void PreliminaryCrackCalcs(double,bool);
	
		// crack hierarchy
		virtual bool CreateHierarchy(void);
		virtual void MoveHierarchy(void);
		void ExtendHierarchy(CrackSegment *);
	
        // methods
        void Archive(ofstream &);
        short MoveCrack(void);
        short MoveCrack(short);
		void UpdateCrackTractions(void);
		void CrackTipHeating(void);
		void AddTractionForce(void);

		// crack crossing
		virtual short CrackCross(Vector *,Vector *,Vector *,int) const;
        short CrackCrossLeaf(CrackLeaf *,Vector *,Vector *,Vector *,short) const;
        short CrackCrossOneSegment(CrackSegment *,Vector *,Vector *,Vector *,short) const;
		short CrackCrossLeafOnce(CrackLeaf *,Vector *,Vector *,CrackSegment **) const;

		// Methods added for 3D
		virtual void Update_NodeList_Edges_Normal();
	
        // calculate J-integral (YJG)
        virtual void JIntegral(void);      	  // J-Integral calculation
        void PrintContour(ContourPoint *,ContourPoint *,Vector &);
		void GetCOD(CrackSegment *,Vector &,bool);
		void CrackTipAndDirection(int,CrackSegment **,Vector &);
		void InterpolatePosition(int,CrackSegment **,Vector &,bool);
	
		// J countour crossing in 2D
		CrackSegment *ContourCrossCrack(ContourPoint *,Vector *) const;
		CrackSegment *ContourCrossLeaf(CrackLeaf *,double,double,double,double,Vector *,int) const;
		bool SegmentsCross(CrackSegment *,double,double,double,double,Vector *,int) const;
	
		// propagation
        CrackSegment *Propagate(Vector &,int,int);
		
		// crossing while propagating in 2D
		double AdjustGrowForCrossing(Vector *,CrackSegment *,double,Vector *);
		short CrackCrossOnce(double,double,double,double,CrackSegment **) const;
	
		// Accessors
		double Length(void) const;
		int NumberOfSegments(void) const;
		int NumberOfFacets(void) const;
		void SetNumber(int);
		int GetNumber(void) const;
		void SetThickness(double);
		double GetThickness(void) const;
		double *GetThicknessPtr(void);
		virtual bool IsThreeD(void) const;
		bool GetAllowAlternate(int) const;
		void SetAllowAlternate(int,bool);
		void SetFixedCrack(int);
		void SetContactLawID(int);
		void Output(void);
		int CriterionNeeds(void) const;
		virtual void Describe(void) const;
		bool GetHasTractionLaws(void) const;
		void GetInitialDirection(CrackSegment *,Vector &) const;
		int GetWhichTip(CrackSegment *) const;
		CrackSegment *GetCrackTip(int) const;
		CrackSegment *GetAdjToCrackTip(int) const;
	
		// class methods
		static void SetCodLocation(double);
		static double Triangle(double,double,double,double,double,double);
        static double Triangle(Vector *,Vector *,Vector *);
		static bool LineIsInExtents(double,double,double,double,double *,double *);
		static bool LineIsInExtents(Vector *,Vector *,double *,double *);

    protected:
        int numberSegments;
		int numberFacets;						// only used for 3D cracks, otherwise 0
		int fixedCrack,number,customContactLawID;
		bool hasTractionLaws;
		Vector initialDirection[2];
		bool allowAlternate[2];
		double thickness;						// 2D tractions and crack-tip heating
        CrackLeaf *rootLeaf;
		ParseController *crossedCracks;
    
};

extern CrackHeader *firstCrack;
// GRID_JTERMS
extern int JGridSize,JContourType,JTerms,JGridEnergy;

extern CrackHeader **crackList;
extern int numberOfCracks;


#endif
