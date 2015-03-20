/********************************************************************************
    MoreNodalPoint.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Mar 15 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Nodes/NodalPoint.hpp"

/********************************************************************************
    More NodalPoint methods
********************************************************************************/

// allocate and initialize force field
void NodalPoint::InitForceField(void)
{
    fs=new ForceField;
    
    fs->stress.xx=0.;
    fs->stress.yy=0.;
    fs->stress.zz=0.;
    fs->stress.xy=0.;
    fs->force.x=0.;
    fs->force.y=0.;
    fs->numElems=0;
}

// allocate and initialize force field
void NodalPoint::PrintAvgStress(void)
{
    char fline[200];
    double factor=1./(double)fs->numElems;
    
    sprintf(fline,"%5d  %15.7e  %15.7e  %15.7e  %15.7e",num,
            fs->stress.xx*factor,fs->stress.yy*factor,
            fs->stress.zz*factor,fs->stress.xy*factor);
    cout << fline << endl;
}


