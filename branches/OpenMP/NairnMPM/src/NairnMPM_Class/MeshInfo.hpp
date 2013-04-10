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

// grid type - all 3D ones above the marker
enum {	UNKNOWN_GRID=-1,NOT_CARTESIAN=0,SQUARE_GRID,RECTANGULAR_GRID,
			BEGIN_3D_GRIDS,CUBIC_GRID,ORTHOGONAL_GRID };

class MeshInfo
{
    public:
		// properties of a regular grid
		double gridx,gridy,gridz;		// cell size
		double partx,party,partz,lp;	// semilength of particle (lp in natural coordinates)
		double diagx,diagy,diagz;		// cell diagonal
		int xplane,yplane,zplane;		// node spacings in each plane
		double xmin,ymin,zmin;			// minimums (if from a grid)
        double positionCutoff;             // cut off for normal contact when using positiong instead of displacements
		
		// constructors
		MeshInfo(void);
		
		// methods
		void Output(int,bool);
        void OutputContactByDisplacements(void);
		bool EdgeElement2D(int);
		bool EdgeElement3D(int);
		void ListOfNeighbors2D(int,int *);
		void ListOfNeighbors3D(int,int *);
        int FindElementFromPoint(Vector *);
		GridPatch **CreatePatches(int,int);
		
		// Accessors
		void SetCartesian(int,double,double,double);
		void SetElements(int,int,int,double,double,double);
		void SetParticleLength(int);
		int GetCartesian(void);
		double GetMinCellDimension(void);
		int CanDoGIMP(void);
		bool IsStructuredGrid(void);
        bool Is3DGrid(void);
		void GetGridPoints(int *,int *,int *);
		double GetCellVolume(void);
        double GetAverageCellSize(void);
		double GetThickness(void);
		double GetDefaultThickness();
        double GetNormalCODAdjust(Vector *,Vector *,double);
        double GetPerpendicularDistance(Vector *,Vector *,double);
        bool GetContactByDisplacements(void);
        void SetContactByDisplacements(bool);
		
	private:
		int cartesian;					// non-zero (NOT_CARTESIAN) is a regular grid
		int totalElems;					// total number of elements
		int horiz,vert,depth;			// number of elements in that direction (if from a grid)
		double cellVolume;				// cell volume
        double avgCellSize;             // average cell size
        bool contactByDisplacements;    // TRUE is using displacements, false if need to adjust normal COD
		int xpnum,ypnum,zpnum;			// patch grid size

};

extern MeshInfo mpmgrid;

#endif

