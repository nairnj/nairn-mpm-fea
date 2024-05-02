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

#include "System/MPMPrefix.hpp"

class BMPLevel : public LinkedObject
{
    public:
		double angle,thickness,temperature,concentration;
		Vector vel;
		double currentDeltaT;
        
        // constructors and destructors
        BMPLevel();
		BMPLevel(int,int,int);
		void SetDefaults(void);
		BMPLevel *ClearWeight(void);
    
        // methods and accessors
		int Material(void) const;
		int Material(unsigned char) const;
		int Material(unsigned char,double);
		BMPLevel *MaximumWeight(double&,BMPLevel **);
		void *GetContextInfo(void);
		void SetContextInfo(void *info);
		void OutputLevel(void);
	
    private:
		int mat,imin,imax;
		double weight;
		void *contextInfo;
};

extern BMPLevel *firstLevel,*currentLevel;

#endif
