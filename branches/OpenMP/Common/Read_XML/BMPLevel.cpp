/********************************************************************************
    BMPLevel.cpp
    nairn-mpm-fea

    Created by John Nairn on Thu Dec 8 2004.
    Copyright (c) 2004 RSAC Software. All rights reserved.
********************************************************************************/

#include "Read_XML/BMPLevel.hpp"

BMPLevel *firstLevel=NULL;
BMPLevel *currentLevel=NULL;

/*******************************************************************
	BMPLevel: Constructors and Destructors
*******************************************************************/

// Constructors
BMPLevel::BMPLevel()
{
	imin=1;
	imax=0;
}

BMPLevel::BMPLevel(int matnum,int setmin,int setmax)
{
	mat=matnum;
	imin=setmin;
	imax=setmax;
	
	angle=0.;
	thickness=1.;
	ZeroVector(&vel);
	temperature=0.;
	concentration=0.;
}

/*******************************************************************
	BMPLevel: Methods
*******************************************************************/

// see if matches level
int BMPLevel::Material(unsigned char intensity)
{	return (intensity>=imin && intensity<=imax) ? mat : -1;	}

int BMPLevel::Material(unsigned char intensity,double pweight)
{
	if(intensity>=imin && intensity<=imax)
	{	weight+=pweight;
		return mat;
	}
	return -1;
}

// level weight
BMPLevel *BMPLevel::ClearWeight(void)
{	weight=0;
	return (BMPLevel *)nextObject;
}

// check for maximum weight
BMPLevel *BMPLevel::MaximumWeight(double& maxweight,BMPLevel **maxLevel)
{
	if(weight>maxweight)
	{	maxweight=weight;
		*maxLevel=this;
	}
	return (BMPLevel *)nextObject;
}
