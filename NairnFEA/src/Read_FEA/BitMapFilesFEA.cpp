/********************************************************************************
    BitMapFilesFEA.cpp - extra code for FEAReadHandler.cpp to generate elements
		from gray scale BMP files
    NairnFEA

    Created by John Nairn on Feb 5 2008.
    Copyright (c) 2008 RSAC Software. All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_FEA/FEAReadHandler.hpp"
#include "System/CommonArchiveData.hpp"
#include "Read_XML/RectController.hpp"
#include "Elements/ElementBase.hpp"
#include "System/FEAArchiveData.hpp"
#include "Read_XML/BMPLevel.hpp"

//-----------------------------------------------------------
// Check for bmp element, return false if not
//-----------------------------------------------------------
short FEAReadHandler::BMPFileInput(char *xName,const Attributes& attrs)
{
	// check for common commands
	if(BMPFileCommonInput(xName,attrs,MUST_BE_NO_BLOCK,false)) return true;
	
	return FALSE;
}

//-----------------------------------------------------------
// Subroutine to translate BMP file into material point
// throws std::bad_alloc, SAXException()
//-----------------------------------------------------------
void FEAReadHandler::TranslateBMPFiles(void)
{
	// file info and data
	unsigned char **rows,**angleRows=NULL;
	XYInfoHeader info,angleInfo;
	
	// read image file
	char *bmpFullPath=archiver->ExpandOutputPath(bmpFileName);
	ReadBMPFile(bmpFullPath,info,&rows);
	delete [] bmpFullPath;
	
	// angle file name
	bool setAngles = false;
	if(bmpAngleFileName[0][0]>0)
	{	setAngles=true;
		char *bmpFullAnglePath=archiver->ExpandOutputPath(bmpAngleFileName[0]);
		ReadBMPFile(bmpFullAnglePath,angleInfo,&angleRows);
		if(info.height!=angleInfo.height || info.width!=angleInfo.width)
			throw SAXException(BMPError("The image file and angle file sizes do not match.",bmpFileName));
		delete [] bmpFullAnglePath;
	}
	
	// get final image width, height, and size per pixel
	Vector pw;
	const char *msg = CommonReadHandler::DecodeBMPWidthAndHeight(info,bwidth,bheight,orig.z,pw,false);
	if(msg != NULL) throw SAXException(BMPError(msg,bmpFileName));
	pw.z = yflipped ? -1. : 1. ;
	
	// scan all elements - fill those that are in the area and have no material yet
	DomainMap map;
	ElementBase *elem;
	Vector center,del;
	for(int i=0;i<nelems;i++)
	{	elem=theElements[i]; 
        
		// skip if element already has a material
        if(elem->material!=NO_MATERIAL) continue;
		
		// get center (will skip unless the center of extent of the element is in the image)
		elem->GetXYZCentroid(&center);
		
		// half the extent of the element (z means 2D)
		del = MakeVector((elem->GetDeltaX())/2.,(elem->GetDeltaY())/2.,-1.);
		
		// Get range of rows and columns and their weights (or skip if not in the image)
		if(!MapDomainToImage(info,center,orig,del,pw,bwidth,bheight,map)) continue;
		
		// find maximum level and its material ID or none (a hole)
		int matID=-1;
		BMPLevel *nextLevel = FindBMPLevel(firstLevel,map,rows);
		if(nextLevel!=NULL) matID = nextLevel->Material();
		
		// set material ID if found a match
		if(matID>0)
		{	// set it
			elem->material=matID;
			
			// is there an angle image too?
			if(setAngles)
			{	double totalIntensity = FindAverageValue(map,angleRows);
				if(totalIntensity>0.)
				{	double matAngle=minAngle[0]+(totalIntensity-minIntensity[0])*angleScale[0];
					elem->SetAngleInDegrees(matAngle);
				}
			}
		}
	}

	// clean up
	for(int row=0;row<info.height;row++)
	{	delete [] rows[row];
		if(setAngles) delete [] angleRows[row];
	}
	delete [] rows;
	if(setAngles) delete [] angleRows;
}

