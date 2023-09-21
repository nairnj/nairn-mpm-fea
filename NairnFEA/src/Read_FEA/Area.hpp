/********************************************************************************
    Area.hpp
    NairnFEA
    
    Created by John Nairn on 9/28/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _JAN_AREA_

#define _JAN_AREA_

class Path;
class EightNodeIsoparam;

class Area
{
    public:
		int mat;
		double thick,angle;
		char *angleExpr;
	
        //  Constructors and Destructor
		Area(int,char *,double);
		~Area();
		
		// methods
		int AddPath(char *);
		const char *MeshElements(void);
		const char *MeshArea(void);
        const char *AreaError(const char *);
		const char *MeshInterface(void);
		int CheckIntervals(void);
		int PathsAvailable(void);
		double SignedArea(void);
		void MeshLine(double,double,double,double,double,int,int);
		void FindPtOnLine(double,double,double,double,double,int,int,double *,double *,
							double *,double *,double *);
		void LineLocate(double,double,double,double,double,double *,double *);
	
	private:
		int numPaths;
		Path *edges[4];
		Vector *areaNode;
		EightNodeIsoparam *areaElem;
};

extern Area *theArea;

#endif
