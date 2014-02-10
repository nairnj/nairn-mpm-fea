/********************************************************************************
    CrackLeaf.cpp
    nairn-mpm-fea

    Created by John Nairn on Wed Aug 22 2012
    Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Cracks/CrackLeaf.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Cracks/CrackHeader.hpp"

#pragma mark CrackLeaf: Constructors and Destructors

// Constructors
CrackLeaf::CrackLeaf()
{
}

// create Leaf for one (if scrk2->nextSeg==NULL) or two segments (otherwise)
//     of the crack.
// This is only called during initilization or when crack is propagating
CrackLeaf::CrackLeaf(CrackSegment *scrk1,CrackSegment *scrk2)
{
#ifdef HIERARCHICAL_CRACKS
    // next leaf in this level
    nextLeaf = NULL;
	parent = NULL;
    
    // child segment
    child1 = (char *)scrk1;
    scrk1->parent = this;
    
    // means the children are crack segments
    terminalLeaf = TRUE;
   
    int i;
    if(scrk2->nextSeg!=NULL)
    {   for(i=0;i<4;i++)
        {   cnear[i] = fmin(scrk1->cnear[i],scrk2->cnear[i]);
            cfar[i] = fmax(scrk1->cfar[i],scrk2->cfar[i]);
        }
        child2 = (char *)scrk2;
        scrk2->parent = this;
    }
    else
    {   for(i=0;i<4;i++)
        {   cnear[i] = scrk1->cnear[i];
            cfar[i] = scrk1->cfar[i];
        }
        child2 = NULL;
    }
#endif
}

// Create parent Leaf for two existing leaves.
// This is only called during initilization or when crack is propagating
CrackLeaf::CrackLeaf(CrackLeaf *leaf1,CrackLeaf *leaf2)
{
    // next leaf in this level
    nextLeaf = NULL;
	parent = NULL;
    
    // child leaf
    child1 = (char *)leaf1;
    leaf1->parent = this;
    
    // means the children are leaves too (not segments)
    terminalLeaf = FALSE;
   
    int i;
    if(leaf2!=NULL)
    {   for(i=0;i<4;i++)
        {   cnear[i] = fmin(leaf1->cnear[i],leaf2->cnear[i]);
            cfar[i] = fmax(leaf1->cfar[i],leaf2->cfar[i]);
        }
        child2 = (char *)leaf2;
        leaf2->parent = this;
    }
    else
    {   for(i=0;i<4;i++)
        {   cnear[i] = leaf1->cnear[i];
            cfar[i] = leaf1->cfar[i];
        }
        child2 = NULL;
    }
}

#pragma mark METHODS

// Recalculate extents when crack moves
void CrackLeaf::GetLeafExtents(void)
{
#ifdef HIERARCHICAL_CRACKS
    int i;
    if(terminalLeaf)
    {   CrackSegment *scrk1 = (CrackSegment *)child1;
        if(child2!=NULL)
        {   CrackSegment *scrk2 = (CrackSegment *)child2;
            for(i=0;i<4;i++)
            {   cnear[i] = fmin(scrk1->cnear[i],scrk2->cnear[i]);
                cfar[i] = fmax(scrk1->cfar[i],scrk2->cfar[i]);
            }
        }
        else
        {   for(i=0;i<4;i++)
            {   cnear[i] = scrk1->cnear[i];
                cfar[i] = scrk1->cfar[i];
            }
        }
    }
    else
    {   CrackLeaf *leaf1 = (CrackLeaf *)child1;
        if(child2!=NULL)
        {   CrackLeaf *leaf2 = (CrackLeaf *)child2;
            for(i=0;i<4;i++)
            {   cnear[i] = fmin(leaf1->cnear[i],leaf2->cnear[i]);
                cfar[i] = fmax(leaf1->cfar[i],leaf2->cfar[i]);
            }
        }
        else
        {   for(i=0;i<4;i++)
            {   cnear[i] = leaf1->cnear[i];
                cfar[i] = leaf1->cfar[i];
            }
        }
    }
#endif
}

// Shift to left to use provided segment and previous child segment
// Will need to redo extents in shifted structures
void CrackLeaf::ShiftLeft(CrackSegment *scrk1)
{
    // put scrk1 in child1 and move child1 to child2
    child2 = child1;
    child1 = (char *)scrk1;
#ifdef HIERARCHICAL_CRACKS
    scrk1->parent = this;
#endif
}

// Add segment to this leaf. If one cannot be added, create and return a new leaf
// If it can be added, return NULL
CrackLeaf *CrackLeaf::AddSegment(CrackSegment *scrk2)
{
    if(child2==NULL)
    {   child2 = (char *)scrk2;
#ifdef HIERARCHICAL_CRACKS
        scrk2->parent = this;
#endif
        return NULL;
    }
    else
        return new CrackLeaf(scrk2,scrk2->nextSeg);
}

// Add leaf to this leaf. If one cannot be added, create and return a new leaf
// If it can be added, return NULL
CrackLeaf *CrackLeaf::AddLeaf(CrackLeaf *leaf2)
{
    if(child2==NULL)
    {   child2 = (char *)leaf2;
        leaf2->parent = this;
        return NULL;
    }
    else
        return new CrackLeaf(leaf2,NULL);
}

#pragma mark ACCESSORS

// Describe hierarchical leaf structure with recursive calls
void CrackLeaf::Describe(int level)
{
    cout << "# Level = " << level << " box = ((" << cnear[0] << "," << cfar[0] << "),(" << cnear[1] << "," << cfar[1] << "))";
    if(terminalLeaf)
        cout << " (last leaf level)" << endl;
    else
    {   cout << endl;
        ((CrackLeaf *)child1)->Describe(level+1);
        if(child2!=NULL) ((CrackLeaf *)child2)->Describe(level+1);
    }
}

// Describe hierarchical leaf structure with recursive calls
void CrackLeaf::DescribeSegments(int level)
{
    if(level == 0) cout << "# ";
    
    if(terminalLeaf)
    {   CrackSegment *scrk1 = (CrackSegment *)child1;
        cout << "(" << scrk1->x << "," << scrk1->y << ")" ;
        if(child2!=NULL)
        {   scrk1 = (CrackSegment *)child2;
            cout << "(" << scrk1->x << "," << scrk1->y << ")" ;
        }
    }
    else
    {   ((CrackLeaf *)child1)->DescribeSegments(level+1);
        if(child2!=NULL) ((CrackLeaf *)child2)->DescribeSegments(level+1);
    }
    
    if(level == 0) cout << endl;

}


// return TRUE is childen are segments
bool CrackLeaf::ChildrenAreSegments(void) { return terminalLeaf; }

// return children cast as CrackLeaf objects
void CrackLeaf::GetChildLeaves(CrackLeaf **leaf1,CrackLeaf **leaf2)
{   *leaf1 = (CrackLeaf *)child1;
    *leaf2 = (CrackLeaf *)child2;
}

// return children cast as CrackSegment objects
// this always returns two segments even if child2 is NULL (which happens when seg2 is lastSeg in the crack)
void CrackLeaf::GetChildSegments(CrackSegment **seg1,CrackSegment **seg2)
{   *seg1 = (CrackSegment *)child1;
    *seg2 = (*seg1)->nextSeg;
}




