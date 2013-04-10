/********************************************************************************
    GridPatch.cpp
    NairnMPM

    Created by John Nairn on 4/10/13.
    Copyright (c) 2013 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Patches/GridPatch.hpp"


// globals
GridPatch **patches;            // list of patches
int GridPatch::ghostRows = 1;   // number of ghost rows. If needed, increase for higher strain limits

#pragma mark INITIALIZATION

// Constructors
GridPatch::GridPatch(int x1,int x2,int y1,int y2,int z1,int z2)
{	
    cout << "(" << x1 << "-" << x2 << "),(" << y1 << "-" << y2 << "),("  << z1 << "-" << z2 << ")" << endl;
    xmin = x1;
    xmax = x2;
    ymin = y1;
    ymax = y2;
    zmin = z1;
    zmax = z2;              // = 0 if 2D patch
}
