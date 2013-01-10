/********************************************************************************
    CrackLeaf.hpp
    NairnMPM

    Created by John Nairn on Wed Aug 22 2012
    Copyright (c) 2012 John A. Nairn, All rights reserved.

    Dependencies
        none
********************************************************************************/

#ifndef _CRACKLEAF_

#define _CRACKLEAF_

class CrackSegment;

// Crack Leaf Class
class CrackLeaf
{
    public:
        CrackLeaf *nextLeaf;
        CrackLeaf *parent;
        double cnear[4],cfar[4];
    
        // constructors
        CrackLeaf();
        CrackLeaf(CrackSegment *,CrackSegment *);
        CrackLeaf(CrackLeaf *,CrackLeaf *);
    
        // methods
        void GetLeafExtents(void);
        void ShiftLeft(CrackSegment *);
        CrackLeaf *AddSegment(CrackSegment *);
        CrackLeaf *AddLeaf(CrackLeaf *);
    
        // accessors
        void Describe(int level);
        void DescribeSegments(int);
        bool ChildrenAreSegments(void);
        void GetChildLeaves(CrackLeaf **,CrackLeaf **);
        void GetChildSegments(CrackSegment **,CrackSegment **);
    
    private:
        char *child1,*child2;
        bool terminalLeaf;
};


#endif