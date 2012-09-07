/********************************************************************************
    NodalPoint2D.hpp
    NairnMPM
    
    Created by John Nairn on 8/30/07.
    Copyright (c) 200y John A. Nairn, All rights reserved.
	
	Dependencies
		NodalPoint.hpp
********************************************************************************/

#include "Nodes/NodalPoint.hpp"

#ifndef _NODALPOINT2D_

#define _NODALPOINT2D_

class NodalPoint2D : public NodalPoint
{
    public:
       
        // constructors and destructors
		NodalPoint2D(int,double,double);

		// methods
        virtual void PrintNodalPoint(ostream &);
};

#endif
