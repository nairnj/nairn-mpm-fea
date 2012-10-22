/********************************************************************************
	PolyTrianglechpp
	NairnMPM

	Created by John Nairn on 1/7/11.
	Copyright (c) 2011 John A. Nairn, All rights reserved.
********************************************************************************/

#include "PolyTriangle.hpp"

/********************************************************************************
	PolyTriangle: constructors and destructors
********************************************************************************/

// these are a triangle face ccw when looking from outside the polyhedron
PolyTriangle::PolyTriangle(Vector vp0,Vector vp1,Vector vp2)
{
	v0=vp0;
	v1=vp1;
	v2=vp2;
	
	// normal is (v1-v0) X (v2-v0)
	Vector a,b;
	SubVector(CopyVector(&a,&v1),&v0);
	SubVector(CopyVector(&b,&v2),&v0);
	CrossProduct(&n,&a,&b);
	
	if(!DbleEqual(n.x,0))
	{	m11 = b.z/n.x;
		m12 = -b.y/n.x;
		m21 = -a.z/n.x;
		m22 = a.y/n.x;
		style = USE_NX;
	}
	else if(!DbleEqual(n.z,0))
	{	m11 = b.y/n.z;
		m12 = -b.x/n.z;
		m21 = -a.y/n.z;
		m22 = a.x/n.z;
		style = USE_NZ;
	}
	else
	{	m11 = -b.z/n.y;
		m12 = b.x/n.y;
		m21 = a.z/n.y;
		m22 = -a.x/n.y;
		style = USE_NY;
	}
	
	fmin.x = fmin(fmin(v0.x,v1.x),v2.x);
	fmin.y = fmin(fmin(v0.y,v1.y),v2.y);
	fmin.z = fmin(fmin(v0.z,v1.z),v2.z);
	fmax.x = fmax(fmax(v0.x,v1.x),v2.x);
	fmax.y = fmax(fmax(v0.y,v1.y),v2.y);
	fmax.z = fmax(fmax(v0.z,v1.z),v2.z);
}

/********************************************************************************
	PolyTriangle: methods
********************************************************************************/

// return -1 if does not cross, 0 if point on the face, 1 if does cross
// Increment edges if point on edge of face (for nx!=0), but for nx==0
//    increment it if the ray passes through the face
int PolyTriangle::PointCrossesFace(Vector *p,unsigned *edges)
{
	double u,v;
	double d = (v0.x-p->x)*n.x + (v0.y-p->y)*n.y + (v0.z-p->z)*n.z;
	
	switch(style)
	{	case USE_NX:
			u = m11*(p->y-v0.y) + m12*(p->z-v0.z);
			if(u<0. || u>1.) return -1;
			v = m21*(p->y-v0.y) + m22*(p->z-v0.z);
			if(v<0 || v>1. || u+v>1.) return -1;
			// if d=0 and here, then point is on the face and consider as insize
			if(d==0.) break;
			// If signs are correct then ray passes through the faces, otherwise no
			if((d>0. && n.x>0.) || (d<0. && n.x<0.))
			{	if(u==0. || v==0. || u+v==1.) *edges++;
				return 1;
			}
			return -1;
		case USE_NZ:
			// nx=0, not cross if not on plane (i.e. if d!=0)
			if(d!=0.) return -1;
			// increment edges if ray passes through the face
			if(p->y>=fmin.y && p->y<=fmax.y && p->z>=fmin.z && p->z<=fmax.z && p->x<=fmin.x) *edges++;
			// now check if point is on this face
			u = m11*(p->x-v0.x) + m12*(p->y-v0.y);
			if(u<0. || u>1.) return -1;
			v = m21*(p->x-v0.x) + m22*(p->y-v0.y);
			if(v<0 || v>1. || u+v>1.) return -1;
			break;
		case USE_NY:
			// nx=nz=0, not cross if not on plane (i.e. if d!=0)
			if(d!=0.) return -1;
			// increment edges if ray passes through the face (no need to check y since all must be equal)
			if(p->z>=fmin.z && p->z<=fmax.z && p->x<=fmin.x) *edges++;
			// now check if point is on this face
			u = m11*(p->x-v0.x) + m12*(p->z-v0.z);
			if(u<0. || u>1.) return -1;
			v = m21*(p->x-v0.x) + m22*(p->z-v0.z);
			if(v<0 || v>1. || u+v>1.) return -1;
			break;
	}
	return 0;
}

// get extents in x direction
void PolyTriangle::GetExtents(Vector *amin,Vector *amax)
{	*amin = fmin;
	*amax = fmax;
}
