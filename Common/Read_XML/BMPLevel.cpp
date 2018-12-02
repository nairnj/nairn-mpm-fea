/********************************************************************************
    BMPLevel.cpp
    nairn-mpm-fea

    Created by John Nairn on Thu Dec 8 2004.
    Copyright (c) 2004 RSAC Software. All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/BMPLevel.hpp"

BMPLevel *firstLevel=NULL;
BMPLevel *currentLevel=NULL;

#pragma mark BMPLevel::Constructors and Destructors

// Constructors
BMPLevel::BMPLevel()
{
	mat=0;
	imin=1;
	imax=0;
	SetDefaults();
}

BMPLevel::BMPLevel(int matnum,int setmin,int setmax)
{
	mat=matnum;
	imin=setmin;
	imax=setmax;
	SetDefaults();
}

// set all default values
void BMPLevel::SetDefaults(void)
{	angle=0.;
	thickness=1.;
	ZeroVector(&vel);
	temperature=0.;
	currentDeltaT=0.;
	concentration=0.;
	contextInfo = NULL;
}

// level weight
BMPLevel *BMPLevel::ClearWeight(void)
{	weight=0;
	return (BMPLevel *)nextObject;
}

#pragma mark BMPLevel::Methods and Accessors

// see if matches level and return mat if does, otherwise return -1
int BMPLevel::Material(void) const { return mat; }

// see if matches level and return mat if does, otherwise return -1
int BMPLevel::Material(unsigned char intensity) const
{	return (intensity>=imin && intensity<=imax) ? mat : -1;	}

// if intensity is range for this level, add to weight and return mat
// otherwise return -1
int BMPLevel::Material(unsigned char intensity,double pweight)
{	if(intensity>=imin && intensity<=imax)
	{	weight+=pweight;
		return mat;
	}
	return -1;
}

// check for maximum weight
BMPLevel *BMPLevel::MaximumWeight(double& maxweight,BMPLevel **maxLevel)
{	if(weight>maxweight)
	{	// Accept if not void (i.e. mat!=0), but if void, all accept if
		// void really has a higher weight
		if(mat!=0 || fabs(weight-maxweight)>1.e-10)
		{	maxweight=weight;
			*maxLevel=this;
		}
	}
	return (BMPLevel *)nextObject;
}

// context info pointer
void *BMPLevel::GetContextInfo(void) { return contextInfo; }
void BMPLevel::SetContextInfo(void *info) { contextInfo = info; }

// print level for thermal ramp
void BMPLevel::OutputLevel(void)
{	cout << "      " << imin << " to " << imax << " has dT = " << temperature << endl;
}

