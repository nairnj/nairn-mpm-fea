/********************************************************************************
    ContourPoint.hpp
    NairnMPM
    
    Created by John Nairn on Feb 13, 2003.
    Copyright (c) 2002 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _CONTOURPOINT_

#define _CONTOURPOINT_

class NodalPoint;

// line orientations
enum { ANGLED=-1,HORIZONTAL,VERTICAL};

class ContourPoint
{
    public:
        NodalPoint *node;
        ContourPoint *nextPoint;
        int orient;
        Vector norm;
        double ds;
        
        // constructors and destructors
        ContourPoint();
        ContourPoint(NodalPoint *);
        
        // methods
        int SetNextPoint(ContourPoint *);
        double Fraction(Vector &);

    private:
};

#endif
