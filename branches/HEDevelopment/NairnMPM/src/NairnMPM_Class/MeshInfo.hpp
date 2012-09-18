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

// grid type - all 3D ones above the marker
enum {	UNKNOWN_GRID=-1,NOT_CARTESIAN=0,SQUARE_GRID,RECTANGULAR_GRID,
			BEGIN_3D_GRIDS,CUBIC_GRID,ORTHOGONAL_GRID };

class MeshInfo
{
    public:
		double gridx,gridy,gridz;
		double partx,party,partz;
		double diagx,diagy,diagz;
		int xplane,yplane,zplane;
		double xmin,ymin,zmin;			// minimums (if from a grid)
		
		// constructors
		MeshInfo(void);
		
		// methods
		void Output(int);
		bool EdgeElement2D(int);
		bool EdgeElement3D(int);
		void ListOfNeighbors2D(int,int *);
		void ListOfNeighbors3D(int,int *);
        int FindElementFromPoint(Vector *);
		
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
		double GetThickness(void);
		double GetDefaultThickness();
		
	private:
		int cartesian;
		int totalElems;
		int horiz,vert,depth;			// number of elements in that direction (if from a grid)
		double cellVolume;
};

extern MeshInfo mpmgrid;

#endif

