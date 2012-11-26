/********************************************************************************
    ElementBase.hpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _ELEMENTBASE_

#define _ELEMENTBASE_

#ifdef MPM_CODE
class MPMBase;

// Undefine to try axisymmetric GIMP using planar GIMP shape functions
//#define NONRADIAL_GIMP_AS

// Comment out to try uGIMP with no truncation AND lCPDI with no shrinkage
#define TRUNCATE

#endif

// element types
#define CS_TRIANGLE 1
#define FOUR_NODE_ISO 2
#define EIGHT_NODE_ISO 3
#define ISO_TRIANGLE 4
#define LINEAR_INTERFACE 5
#define QUAD_INTERFACE 6
#define EIGHT_NODE_ISO_BRICK 7
#define NINE_NODE_LAGRANGE 8

// other constants
#define TOLERANCE_RATIO 0.1

#define BMATRIX 2
// max nodes in an element + 1
#define MaxElFaces 4
#define MxFree 2

// neighbors
#define UNKNOWN_NEIGHBOR -2
#define NO_NEIGHBOR -1
        
class ElementBase : public LinkedObject
{
    public:
        int num;						// element number (1 based)
        int nodes[MaxElNd];				// 1 based node numbers in nodes[0] to nodes[NumberNodes()-1]
        double xmin,xmax,ymin,ymax;		// element extent
        int filled;						// flags for loaded material points
#ifdef MPM_CODE
        int *neighbors;					// Elements next to faces
#else
        int material;					// FEA material ID
        double strainEnergy;			// FEA strain energy
		int numGauss;					// number of gaussian points
		int gaussSet;					// which set of points to use
#endif
        
		// static variables
        static double gridTolerance;	// TOLERANCE_RATIO * minimum grid extent
#ifdef MPM_CODE
		static int useGimp;             // Code for GIMP method (0 is classic MPM special case of GIMP)
		static int analysisGimp;		// store GIMP option in case need to disable for a while
        static int numCPDINodes;        // number of nodes used by CPDI in particle domain
#endif
		
        // constructors and destructors
#ifdef MPM_CODE
        ElementBase(int,int *);
#else
        ElementBase(int,int *,int,double);
#endif
		virtual ~ElementBase();
        
        // prototypes for abstract methods (must override)
        virtual short ElementName(void) = 0;
        virtual int NumberNodes(void) = 0;
        virtual double GetArea(void) = 0;
        virtual double GetVolume(void) = 0;
		virtual int FaceNodes(void) = 0;
        virtual short PtInElement(Vector &) = 0;
        virtual void ShapeFunction(Vector *,int,double *,double *,double *,
                                    Vector *,double *,double *,double *) = 0;
#ifdef MPM_CODE
		virtual void ShapeFunction(Vector *,int,double *,double *,double *,double *) = 0;
		virtual void GetGimpNodes(int *,int *,int *,Vector *);
		virtual void GimpShapeFunction(Vector *,int,int *,int,double *,double *,double *,double *);
		virtual void GimpShapeFunctionAS(Vector *,int,int *,int,double *,double *,double *,double *);
#endif
									
        // prototypes of methods defined in ElementBase class (but may override)
		virtual void PrintElement(ostream &,int);
		virtual double GetThickness(void);
        virtual void SetThickness(double);
		virtual void GetXYZCentroid(Vector *);
        virtual double GetCenterX(void);
		virtual double GetDeltaX(void);
		virtual double GetDeltaY(void);
		virtual double GetDeltaZ(void);
		virtual bool IntersectsBox(double,double,double,double,double);
        virtual int NumberSides(void);
		int NodeIndex(int);
        virtual void FindExtent(void);
		virtual void FindCentroid(Vector *);
#ifdef MPM_CODE
		virtual void GetShapeFunctions(int *,double *,int *,Vector *,Vector *,MPMBase *);
		virtual void GetShapeFunctions(int *,double *,int *,Vector *,MPMBase *);
		virtual void GetShapeGradients(int *,double *,int *,Vector *,double *,double *,double *,MPMBase *);
		virtual void GetShapeFunctionsAndGradients(int *,double *,int *,Vector *,Vector *,double *,double *,double *,MPMBase *);
		virtual void GimpCompact(int *,int *,double *,double *,double *,double *);
		virtual void GetCPDIFunctions(int,CPDIDomain **,int *,int *,double *,double *,double *,double *);
        virtual void GridShapeFunctions(int *,int *,Vector *,double *);
		virtual bool OnTheEdge(void);
		virtual void GetListOfNeighbors(int *);
		virtual int NextNode(int);
        virtual int FindEdge(int,int);
        virtual int Neighbor(int);
		virtual void AllocateNeighborsArray(void);
		virtual int Orthogonal(double *,double *,double *);
        virtual int NearestNode(double,double,int *);
        virtual void MPMPoints(short,Vector *);
		virtual void GetXiPos(Vector *,Vector *);
#else
		virtual bool HasNode(int);
		virtual void DecrementNodeNums(int);
        virtual void MaxMinNode(int *,int *);
		virtual void MapNodes(int *);
		virtual void CalcEdgeLoads(double *,int,int,double *,int);
        virtual void Stiffness(int);
        virtual void ForceStress(double *,int,int);
        virtual void GetProperties(int);
        void IsoparametricStiffness(int);
		void ZeroUpperHalfStiffness(void);
		void FillLowerHalfStiffness(void);
        void IsoparametricForceStress(double *,int,int);
		virtual void ExtrapolateGaussStressToNodes(double [][5]);
        int WantElement(char,const vector< int > &);
		void LinearEdgeLoad(int,int,int,double *,double *,int);
		void QuadEdgeLoad(int,int,int,int,double *,double *,int);
		virtual bool BulkElement(void);
		virtual void MakeQuarterPointNodes(int,vector<int> &);
        virtual void SetAngle(double);
		void SetAngleInDegrees(double eAng);
		double GetAngleInDegrees(void);
#endif
		// class methods
#ifdef MPM_CODE
		static void ChangeGimp(int);
		static void RestoreGimp(void);
		static void AllocateNeighbors(void);
        static void InitializeCPDI(bool);
#else
		static void MoveCrackTipNodes(int);
#endif
		static double GetMinimumCellSize(void);

	protected:
#ifdef MPM_CODE
		bool pgElement;					// is it a parallelogram element?
		double pgTerm[6];				// precalculated terms for GetXiPos() speed
#else
        double angle;					// FEA material angle
#endif

#ifdef MPM_CODE
        virtual void GetCentroid(Vector *);
        virtual void GetCoordinates(Vector *,int,int *);
		virtual void GetNodes(int *,int *);
#endif

};

// List of elements stored as theElements[0] to theElements[nelems-1]
extern ElementBase **theElements;
extern int nelems;
#ifdef FEA_CODE
	// temporary globals used in FEA element calculations
    extern Vector ce[MaxElNd];
    extern double re[MxFree*MaxElNd];
    extern double se[MxFree*MaxElNd][MxFree*MaxElNd];
	extern double te[MaxElNd];
#endif

#endif

