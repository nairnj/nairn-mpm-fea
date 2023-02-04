/********************************************************************************
 AreaOfInterest.hpp
 NairnMPM
 
 Created by John Nairn on 2/19/2016.
 Copyright (c) 2016 John A. Nairn, All rights reserved.
 
 Dependencies
	LinkedObject.hpp
 ********************************************************************************/

#ifndef _AREAOFINTEREST_
#define _AREAOFINTEREST_

#define GEOMETRIC_STYLE 0
#define LINEAR_STYLE 1

class AreaOfInterest
{
	public:
		double x1,x2,y1,y2,z1,z2;
	
		// constructors and destructors
		AreaOfInterest();
		AreaOfInterest(double,double,int,double,double,int,double,double,int);
	
		// methods
		bool IsBefore(AreaOfInterest *,int);
		AreaOfInterest *GetCellCounts(int,double &,double,double,int *,int);
		double GetAxisCounts(int,double &,double &,double &,int &,int &,int &,double,double,double,int *,int);
		double GetNextEdges(int,double &,double &);
		void Collapse(int,double,double);
		AreaOfInterest *MeshAreaOfInterest(int,double &,double,double,int *,double *,int);
		double MeshAxis(int,double,double,double,int,int,int,double,double,double,int *,double *,int);
	
		// links
		AreaOfInterest *GetNextAOI(int);
		void SetNextAOI(int,AreaOfInterest *);
	
	protected:
		int nx,ny,nz;
		double xcell,ycell,zcell;
		int nxBefore,nxAfter,nyBefore,nyAfter,nzBefore,nzAfter;
		AreaOfInterest *nextAOI[3];
};

#endif
