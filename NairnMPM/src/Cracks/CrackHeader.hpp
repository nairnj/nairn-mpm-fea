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
		static int warnThreeFields;
		static int warnNodeOnCrack;
		static int warnThreeCracks;
		double thickness;						// 2D tractions and crack-tip heating
		
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
		void CrackTipHeating(void);
		void SetFixedCrack(int);
		void SetFriction(double);
		void SetDn(double);
		void SetDnc(double);
		void SetDt(double);
		void SetContact(double,double,double,double);
		void Output(void);
		void Describe(void);
        
		void SetNumber(int);
		int GetNumber(void);
         
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
        long numberSegments;
        double xmin,xmax,ymin,ymax;
		int fixedCrack,number;
		double crackFriction,crackDn,crackDnc,crackDt;
		bool customContact,hasTractionLaws;
		Vector initialDirection[2];
		bool allowAlternate[2];
        
};

extern CrackHeader *firstCrack;
extern int JGridSize,JContourType,JTerms;

#endif
