/********************************************************************************
    MaterialInterfaceNode.hpp
    NairnMPM
 
    Created by John Nairn on 5/3/12.
    Copyright (c) 2012 John A. Nairn, All rights reserved.
 
    Dependencies
        NodalPoint.hpp
********************************************************************************/

#ifndef _MATERIALINTERFACENODE_

#define _MATERIALINTERFACENODE_

#include "Nodes/NodalPoint.hpp"

class MaterialInterfaceNode
{
    public:
        // variables (changed in MPM time step)
        static MaterialInterfaceNode *currentNode;
    
        // constructors and destructors
        MaterialInterfaceNode(NodalPoint *,int,int,int,Vector *,double);
    
        // methods
        MaterialInterfaceNode *InterfaceForce(void);
        void SetPrevBC(MaterialInterfaceNode *);
        MaterialInterfaceNode *GetPrevBC(void);
        double GetInterfaceTraction(Vector *);
        void GetFieldInfo(int *,int *,int *);
    
        // class methods
        static void RemoveInterfaceNodes(void);
        static void InterfaceOnKnownNodes(void);
    
    private:
        NodalPoint *theNode;
        int vfld,mati,matipaired;
        MaterialInterfaceNode *prevBC;
        Vector traction;
        double energy;
};

#endif

