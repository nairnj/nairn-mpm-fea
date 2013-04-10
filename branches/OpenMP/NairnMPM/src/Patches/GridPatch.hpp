/********************************************************************************
    GridPatch.hpp
    NairnMPM

    Created by John Nairn on 4/10/13.
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        none
********************************************************************************/

#ifndef _GRIDPATCH_

#define _GRIDPATCH_

class MPMBase;

class GridPatch
{
    public:
        static int ghostRows;
    
        // constructors and destructors
        GridPatch(int,int,int,int,int,int);
    
    private:
        int xmin,xmax,ymin,ymax,zmin,zmax;      // element ranges
        MPMBase *firstNR;                       // first non-rigid particle
        MPMBase *firstRC;                       // first rigid contact particle
        MPMBase *firstRBC;                      // first rigid BC particle
};

extern GridPatch **patches;

#endif
