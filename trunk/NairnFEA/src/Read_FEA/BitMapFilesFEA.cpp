/********************************************************************************
    BitMapFilesFEA.cpp - extra code for FEAReadHandler.cpp to generate elements
		from gray scale BMP files
    NairnFEA

    Created by John Nairn on Feb 5 2008.
    Copyright (c) 2008 RSAC Software. All rights reserved.
********************************************************************************/

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
	if(BMPFileCommonInput(xName,attrs,MUST_BE_NO_BLOCK)) return TRUE;
	
	return FALSE;
}

//-----------------------------------------------------------
// Subroutine to translate BMP file into material point
//-----------------------------------------------------------
void FEAReadHandler::TranslateBMPFiles(void)
{
	unsigned char **rows,**angleRows;
	BMPInfoHeader info,angleInfo;
	bool setAngles=FALSE;
	
	// read image file
	char *bmpFullPath=archiver->ExpandOutputPath(bmpFileName);
	ReadBMPFile(bmpFullPath,info,&rows);
	delete [] bmpFullPath;
	
	// angle file name
	if(bmpAngleFileName[0]>0)
	{	setAngles=TRUE;
		char *bmpFullAnglePath=archiver->ExpandOutputPath(bmpAngleFileName);
		ReadBMPFile(bmpFullAnglePath,angleInfo,&angleRows);
		if(info.height!=angleInfo.height || info.width!=angleInfo.width)
			throw SAXException(BMPError("The image file and angle file sizes do not match.",bmpFileName));
		delete [] bmpFullAnglePath;
	}
	
	// provided mm per pixel
	double xpw=-1.,ypw=-1.;
	if(bheight<0. && bheight>-1.e8)
	{	ypw = -bheight;
		bheight = ypw*(double)info.height;
	}
	if(bwidth<0. && bwidth>-1.e8)
	{	xpw = -bwidth;
		bwidth = xpw*(double)info.width;
	}
	
	// total dimensions (if only one is known, find the other, never have both unknown)
	if(bheight<0) bheight=bwidth*(double)info.height/(double)info.width;
	if(bwidth<0) bwidth=bheight*(double)info.width/(double)info.height;
	
	// final mm per pixel (if needed)
	if(xpw<0.) xpw=bwidth/(double)info.width;
	if(ypw<0.) ypw=bheight/(double)info.height;
	
	// create a shape
	RectController *imageRect=new RectController(NO_BLOCK);
	imageRect->SetProperty("xmin",xorig);
	imageRect->SetProperty("xmax",xorig+bwidth);
	imageRect->SetProperty("ymin",yorig);
	imageRect->SetProperty("ymax",yorig+bheight);
	
	// scan all elements - fill those that are in the area and have no material yet
	int matID;
	BMPLevel *nextLevel;
		BMPLevel *maxLevel;
	double deltax,deltay;
	double rmin,rmax,wtr1,wtr2,rweight;
	double cmin,cmax,wtc1,wtc2,weight;
	int r1,r2,c1,c2;

	int i,row,col;
	ElementBase *elem;
	Vector center;
	for(i=0;i<nelems;i++)
	{	elem=theElements[i]; 
        
		// skip if element already has a material
        if(elem->material!=NO_MATERIAL) continue;
		
		// skip unless the center of extent of the element is in the image
		elem->GetXYZCentroid(&center);
		if(!imageRect->ContainsPoint(center)) continue;
		
		// half the extent of the volume of a particle
		deltax=(elem->GetDeltaX())/2.;
		deltay=(elem->GetDeltaY())/2.;
		
		// find range of rows and cols for pixels over this element extent
		rmin=(center.y-deltay-yorig)/ypw;
		rmax=(center.y+deltay-yorig)/ypw;
		r1=BMPIndex(rmin,info.height);
		r2=BMPIndex(rmax,info.height);
		if(r2==r1)
		{	wtr1=wtr2=1.;
		}
		else
		{	// fractional weight for first and last row in case not all within the extent
			wtr1=(double)(r1+1)-rmin;
			wtr2=rmax-(double)r2;
		}
		cmin=(center.x-deltax-xorig)/xpw;
		cmax=(center.x+deltax-xorig)/xpw;
		c1=BMPIndex(cmin,info.width);
		c2=BMPIndex(cmax,info.width);
		if(c2==c1)
		{	wtc1=wtc2=1.;
		}
		else
		{	// fractional weight for first and last col in case not all within the extent
			wtc1=(double)(c1+1)-cmin;
			wtc2=cmax-(double)c2;
		}
		
		// find material ID or none (a hole)
		// clear material weights
		nextLevel=firstLevel;
		while(nextLevel!=NULL) nextLevel=nextLevel->ClearWeight();

		// add weights
		Vector scanPt;
		scanPt.z=0.;
		for(row=r1;row<=r2;row++)
		{	scanPt.y=ypw*row+yorig;		// actual point
			if(row==r1)
				rweight=wtr1;
			else if(row==r2)
				rweight=wtr2;
			else
				rweight=1.;
			for(col=c1;col<=c2;col++)
			{	// skip if point is not in the element
				scanPt.x=xpw*col+xorig;
				if(!elem->PtInElement(scanPt)) continue;
				if(col==c1)
					weight=rweight*wtc1;
				else if(col==c2)
					weight=rweight*wtc2;
				else
					weight=rweight;
				
				// find material at this level
				// (note: last level with matID=0 catches empty space)
				nextLevel=firstLevel;
				while(nextLevel!=NULL)
				{	matID=nextLevel->Material(rows[row][col],weight);
					if(matID>=0) break;
					nextLevel=(BMPLevel *)nextLevel->GetNextObject();
				}
			}
		}
		
		// find maximum weight (matID=0 means max is empty space)
		weight=0.;
		maxLevel=NULL;
		nextLevel=firstLevel;
		while(nextLevel!=NULL) nextLevel=nextLevel->MaximumWeight(weight,&maxLevel);
		if(maxLevel!=NULL)
		{	matID=maxLevel->mat;
			nextLevel=maxLevel;
		}
		else
			matID=-1;

		// set material ID if found a match
		if(matID>0)
		{	// set it
			elem->material=matID;
			
			// is there an angle image?
			if(setAngles)
			{	// weight average of scanned pixels
				double totalWeight=0.;
				double totalIntensity=0.;
				for(row=r1;row<=r2;row++)
				{	scanPt.y=ypw*row+yorig;		// actual point
					if(row==r1)
						rweight=wtr1;
					else if(row==r2)
						rweight=wtr2;
					else
						rweight=1.;
					for(col=c1;col<=c2;col++)
					{	// skip if point if not in the element
						scanPt.x=xpw*col+xorig;
						if(!elem->PtInElement(scanPt)) continue;
						if(col==c1)
							weight=rweight*wtc1;
						else if(col==c2)
							weight=rweight*wtc2;
						else
							weight=rweight;
						
						totalIntensity+=weight*angleRows[row][col];
						totalWeight+=weight;
					}
				}
				if(totalWeight>0.)
				{	totalIntensity/=totalWeight;
					double matAngle=minAngle+(totalIntensity-minIntensity)*angleScale;
					elem->SetAngleInDegrees(matAngle);
				}
			}
		}
	}

	// done with shape
	delete imageRect;
	
	// clean up
	for(row=0;row<info.height;row++)
	{	delete [] rows[row];
		if(setAngles) delete [] angleRows[row];
	}
	delete [] rows;
	if(setAngles) delete [] angleRows;
}

