/*********************************************************************
    MeshInfo.hpp
    Nairn Research Group MPM Code
    
    Created by John Nairn on 5/16/06.
    Copyright (c) 2006, All rights reserved.
	
	Dependencies
		none
*********************************************************************/

#ifndef _MESHINFO_

#define _MESHINFO_

class GridPatch;
class NodalPoint;
class MPMBase;

// grid type - all 3D ones above the marker, all known grid types > 0
//		UNKNOWN_GRID only while reading, it is set if not known before things start
// Square, Rectangular, Cubic, and Orthogonal all have equal element sizes
// Variable means rectangular (2D) or orthogonal (3D), but variable element sizes
// NOT_CARTENSIAN means unstrutured and elements not orthogonal. Code does not work in this
//		case without a lot of modification
enum {	UNKNOWN_GRID=-1,NOT_CARTESIAN=0,SQUARE_GRID,RECTANGULAR_GRID,VARIABLE_RECTANGULAR_GRID,
			BEGIN_3D_GRIDS,CUBIC_GRID,ORTHOGONAL_GRID,VARIABLE_ORTHOGONAL_GRID,NOT_CARTESIAN_3D };

class MeshInfo
{
    public:
		// properties of a regular grid
        double xmin,ymin,zmin;			// minimums (if from a grid)
		double xmax,ymax,zmax;			// maximums (if from a grid)
        double positionCutoff;          // cut off for normal contact when using positiong instead of displacements
		int xplane,yplane,zplane;		// node spacings in each plane
		Vector minParticleSize;
	
		// constructors
		MeshInfo(void);
		
		// methods
		void Output(int,bool);
        void OutputContactByDisplacements(bool);
		bool EdgeElement2D(int);
		bool EdgeElement3D(int);
		bool EdgeNode(int,char);
		void ListOfNeighbors2D(int,int *);
		void ListOfNeighbors3D(int,int *);
		GridPatch **CreatePatches(int,int);
		GridPatch **CreateOnePatch(int);
	
		// Accessors
		int GetPatchForElement(int);
		void SetCartesian(int,double,double,double);
		void SetElements(int,int,int,double,double,double,double,double,double);
		int GetCartesian(void);
		bool IsStructuredGrid(void) const;
        bool IsStructuredEqualElementsGrid(void) const;
        bool Is3DGrid(void);
		void GetGridPoints(int *,int *,int *);
		double GetDefaultThickness();
        bool GetContactByDisplacements(void);
        void SetContactByDisplacements(bool);
		double GetMinCellDimension(void);
		void TrackMinParticleSize(Vector);
		Vector GetGlobalMinParticleSize(void) const;
		double GetGlobalMinParticleLength(void) const;
		double GetThickness(void);                          // 2D only
		int FindElementFromPoint(const Vector *,MPMBase *);
		int BinarySearchForElement(int,double,int,int);
		void FindElementCoordinatesFromPoint(Vector *,int &,int &,int &);
		double GetCellVolume(NodalPoint *);
		void GetLocalCellSizes(NodalPoint *,double &,double &,double &,double &,double &,double &);
		double GetCellRatio(NodalPoint *,int);
		double GetAverageCellSize(MPMBase *);
		bool FindMeshLineForBCs(int,double,int,int &,double &,double &);
		Vector GetPerpendicularDistance(Vector *,NodalPoint *);
		void ScaleTangent(double,double,double,double &,double &);
   
        // use of these implies equal element sizes
        Vector GetCellSize(void);
			
	private:
		int cartesian;					// non-zero (NOT_CARTESIAN=0) is a regular grid
		bool equalElementSizes;			// true if all elements the same size
		int totalElems;					// total number of elements
		int horiz,vert,depth;			// number of elements in that direction (if from a grid)
		double cellMinSize;				// minimum cell length
        bool contactByDisplacements;    // true is using displacements, false if need to adjust normal COD
		int xpnum,ypnum,zpnum;			// patch grid size
		int xPatchSize,yPatchSize,zPatchSize;		// patch sizes in elements (last may differ)

		double cellVolume;				// cell volume when equal element sizes
		double avgCellSize;             // average cell size when equal element sizes
        Vector grid;                    // cell size when equal element sizes
};

extern MeshInfo mpmgrid;

#endif

