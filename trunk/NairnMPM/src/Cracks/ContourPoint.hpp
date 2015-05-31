/********************************************************************************
    ContourPoint.hpp
    nairn-mpm-fea
    
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
		bool phantomNode;
    
        // constructors and destructors
        ContourPoint();
        ContourPoint(NodalPoint *);
		virtual ~ContourPoint();
	
        // methods
        int SetNextPoint(ContourPoint *);
        double Fraction(Vector &);
		void SetPhantomNode(bool);

    private:
};

#endif
