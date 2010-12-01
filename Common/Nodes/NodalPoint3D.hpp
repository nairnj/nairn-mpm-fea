/********************************************************************************
    NodalPoint3D.hpp
    NairnMPM
    
    Created by John Nairn on 8/30/07.
    Copyright (c) 200y John A. Nairn, All rights reserved.
	
	Dependencies
		NodalPoint.hpp
********************************************************************************/

#include "Nodes/NodalPoint.hpp"

#ifndef _NODALPOINT3D_

#define _NODALPOINT3D_

class NodalPoint3D : public NodalPoint
{
    public:
       
        // constructors and destructors
		NodalPoint3D(int,double,double,double);
		
		// methods
        virtual void PrintNodalPoint(ostream &);
		
};

#endif
