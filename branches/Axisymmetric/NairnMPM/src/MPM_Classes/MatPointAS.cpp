/********************************************************************************
	MatPointAS.hpp
	NairnMPM

	Created by John Nairn on 10/23/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
 
	When MPM is axisymmetric, the x axis becomes radial (or R) direction, the
		y axis because axial direction (or Z direction), and the z direction
		become the tangential direction (or theta direction)
 
	The thickness (in thick) is set to the initial radial position of the
		material point
********************************************************************************/

#include "MPM_Classes/MatPointAS.hpp"

#pragma mark MatPointAS::Constructors and Destructors

// Constructors
MatPointAS::MatPointAS() {}

// Constructors
MatPointAS::MatPointAS(int inElemNum,int theMatl,double angin,double thickin) : MatPoint2D(inElemNum,theMatl,angin,thickin)
{	
}

#pragma mark MatPoint2D::Accessors

// set original position (2D) (in mm)
void MatPointAS::SetOrigin(Vector *pt)
{	origpos.x=pt->x;
    origpos.y=pt->y;
	origpos.z=0.;
	thick = pt->x;
}

