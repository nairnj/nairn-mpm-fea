/********************************************************************************
    NodalPoint3D.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 8/30/07.
    Copyright (c) 200y John A. Nairn, All rights reserved.
********************************************************************************/

#include "NodalPoint3D.hpp"

#pragma mark NodalPoin3D: Constructors and Destructor

NodalPoint3D::NodalPoint3D(int nodeNum,double xPt,double yPt,double zPt)
		 : NodalPoint(nodeNum)
{
    x=xPt;
    y=yPt;
	z=zPt;
}

#ifdef MPM_CODE
NodalPoint3D::NodalPoint3D(NodalPoint *real) : NodalPoint(real) {}
#endif

#pragma mark NodalPoint3D: Methods

// write node to output file
void NodalPoint3D::PrintNodalPoint(ostream &os)
{
    char nline[200];
	sprintf(nline,"%5d %15.7e %15.7e %15.7e",num,x,y,z);
    os << nline << endl;
}
