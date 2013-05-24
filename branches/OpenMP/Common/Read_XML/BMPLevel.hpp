/********************************************************************************
    BMPLevel.hpp
    nairn-mpm-fea

    Created by John Nairn on Thu Dec 8 2004.
    Copyright (c) 2004 RSAC Software. All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _BMPLEVEL_

#define _BMPLEVEL_

class BMPLevel : public LinkedObject
{
    public:
		int mat,imin,imax;
		double angle,thickness,temperature,concentration;
		Vector vel;
        
        // constructors and destructors
        BMPLevel();
		BMPLevel(int,int,int);
        
		int Material(unsigned char);
		int Material(unsigned char,double);
		BMPLevel *ClearWeight(void);
		BMPLevel *MaximumWeight(double&,BMPLevel **);
		
    private:
		double weight;
};

extern BMPLevel *firstLevel,*currentLevel;

#endif
