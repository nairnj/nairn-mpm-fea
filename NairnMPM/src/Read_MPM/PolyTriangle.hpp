/********************************************************************************
	PolyTriangle.hpp
	nairn-mpm-fea
 
	Created by John Nairn on 1/7/11.
	Copyright (c) 2011 John A. Nairn, All rights reserved.
 
	Dependencies
		none
********************************************************************************/

#ifndef _POLYTRIANGLE_

#define _POLYTRIANGLE_
#include "System/MPMPrefix.hpp"


enum { USE_NX=1,USE_NY,USE_NZ };

class PolyTriangle
{
	public:
		
		// constructors
		PolyTriangle(Vector,Vector,Vector);
		
		// methods
		int PointCrossesFace(Vector *,unsigned *);
		void GetExtents(Vector *,Vector *);
	
	private:
		Vector v0,v1,v2,n;
		double m11,m12,m21,m22;
		int style;
		Vector fmin,fmax;
};

#endif

