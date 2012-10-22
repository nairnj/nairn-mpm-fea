/********************************************************************************
    Keypoint.hpp
    NairnFEA
    
    Created by John Nairn on 6/23/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _KEYPOINT_

#define _KEYPOINT_

class Keypoint : public LinkedObject
{
    public:
		double x,y,temperature;
		char keyID[33];
		int node;
	
        //  Constructors and Destructor
		Keypoint(const char *,double,double);
		
		// methods
		void AssignToNode(void);
		int SamePtAs(Keypoint *);
		
		// accessors
		void GetKeypointXY(Vector *);
};

#endif
