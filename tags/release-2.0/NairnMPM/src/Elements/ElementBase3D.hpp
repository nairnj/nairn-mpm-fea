/********************************************************************************
    ElementBase3D.hpp
    NairnMPM
    
    Created by John Nairn on John Nairn on 7/20/06.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		ElementBase.hpp
********************************************************************************/

#ifndef _ELEMENTBASE3D_

#define _ELEMENTBASE3D_

#include "Elements/ElementBase.hpp"

class ElementBase3D : public ElementBase
{
	public:
		double zmin,zmax;		// z extent
		
		// constructors
		ElementBase3D(long,long *);
			
		// methods
		virtual void FindExtent(void);
		virtual int FaceNodes(void);
		virtual void GetXYZCentroid(Vector *);
		virtual double GetDeltaZ(void);
		virtual bool IntersectsBox(double,double,double,double,double);
		virtual bool OnTheEdge(void);
		virtual void GetListOfNeighbors(int *);
};

#endif
