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

// grid type - all 3D ones above the marker, all known grid types > 0
// Square, Rect, Cubic, and orthogonal all have equal element sizes
// Variable means rectangular (2D) or orthogonal (3D), but variable element sizes
enum {	UNKNOWN_GRID=-1,NOT_CARTESIAN=0,SQUARE_GRID,RECTANGULAR_GRID,VARIABLE_RECTANGULAR_GRID,
			BEGIN_3D_GRIDS,CUBIC_GRID,ORTHOGONAL_GRID,VARIABLE_ORTHOGONAL_GRID };

class MeshInfo
{
    public:
		// properties of a regular grid
        double xmin,ymin,zmin;			// minimums (if from a grid)
		double xmax,ymax,zmax;			// maximums (if from a grid)
        double positionCutoff;          // cut off for normal contact when using positiong instead of displacements
		int xplane,yplane,zplane;		// node spacings in each plane
		
		// constructors
		MeshInfo(void);
		
		// methods
		void Output(int,bool);
        void OutputContactByDisplacements(void);
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
		double GetThickness(void);                          // 2D only
    
        // use of these implies equal element sizes
		int FindElementFromPoint(Vector *);
        double GetAverageCellSize(void);
        void SetParticleLength(int);
        double GetParametersForBCs(int axis,double *,double *);
        double GetCellVolume(void);
        double GetPerpendicularDistance(Vector *,Vector *,double);
		double GetNormalCODAdjust(Vector *,Vector *,double);
        Vector GetCellSize(void);
        double GetCellXSize(void);
		double GetCellYSize(void);
		double GetCellZSize(void);
        Vector GetParticleSize(void);
        double GetParticleXSize(void);
        double GetParticleYSize(void);
        double GetParticleZSize(void);
        double GetParticleSemiLength(void);
		
	private:
		int cartesian;					// non-zero (NOT_CARTESIAN=0) is a regular grid
		bool equalElementSizes;			// true is all elements the same size
		int totalElems;					// total number of elements
		int horiz,vert,depth;			// number of elements in that direction (if from a grid)
		double cellVolume;				// cell volume
        double avgCellSize;             // average cell size
		double cellMinSize;				// minimum cell length
        bool contactByDisplacements;    // TRUE is using displacements, false if need to adjust normal COD
		int xpnum,ypnum,zpnum;			// patch grid size
		int xPatchSize,yPatchSize,zPatchSize;		// patch sizes in elements (last may differ)

        double lp;                      // semilength of particle (lp in natural coordinates)
        Vector grid;                    // cell size when equal element sizes
        Vector part;                    // particle size in x, y, and z direction when equal element sizes
        Vector diag;                    // cell diagonal
};

extern MeshInfo mpmgrid;

#endif

