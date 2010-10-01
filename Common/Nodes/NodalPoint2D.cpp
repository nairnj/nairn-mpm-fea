/********************************************************************************
    NodalPoint2D.cpp
    NairnMPM
    
    Created by John Nairn on 8/30/07.
    Copyright (c) 200y John A. Nairn, All rights reserved.
********************************************************************************/

#include "NodalPoint2D.hpp"

/********************************************************************************
	NodalPoint: Constructors and Destructor
********************************************************************************/

NodalPoint2D::NodalPoint2D(long nodeNum,double xPt,double yPt)
		 : NodalPoint(nodeNum)
{
    x=xPt;
    y=yPt;
	z=0.;
}

/********************************************************************************
	NodalPoint2D: Methods
********************************************************************************/

// write node to output file
void NodalPoint2D::PrintNodalPoint(ostream &os)
{
	char nline[200];
	sprintf(nline,"%5ld %15.7e %15.7e",num,x,y);
    os << nline << endl;
}