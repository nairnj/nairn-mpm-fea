/********************************************************************************
    BitMapFilesCommon.cpp - extra code for MPMReadHandler.cpp to generate particles
		from gray scale BMP files
    nairn-mpm-fea

    Created by John Nairn on Thu Dec 8 2004.
    Copyright (c) 2004 RSAC Software. All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include <fstream>

#include "Read_XML/CommonReadHandler.hpp"
#include "Read_XML/BMPLevel.hpp"
#include "Read_XML/MaterialController.hpp"
#include "Read_XML/XYBMPImporter.hpp"
#include "Read_XML/XYTXTImporter.hpp"
#include "Exceptions/CommonException.hpp"

#ifdef MPM_CODE
extern char rotationAxes[4];
extern char angleAxes[4];
#endif

#define WEIGHT_TOL 1.e-10

//-----------------------------------------------------------
// Check for bmp element, return false if not
// throws std::bad_alloc, SAXException()
//-----------------------------------------------------------
short CommonReadHandler::BMPFileCommonInput(char *xName,const Attributes& attrs,int expectedBlock,bool is3D)
{
    char *aName,*value;
    int i,numAttr;
	double aScaling;

    //-----------------------------------------------------------
    // Read BMP file name and resolution
    //-----------------------------------------------------------
    if(strcmp(xName,"BMP")==0)
	{	ValidateCommand(xName,expectedBlock,ANY_DIM);
        block=BMPBLOCK;
		bwidth=bheight=-1.e9;		// < -1.e8 means dimension was not specified
		bmpFileName[0]=0;
		bmpAngleFileName[0][0]=0;
		bmpAngleFileName[1][0]=0;
		bmpAngleFileName[2][0]=0;
		numAngles=0;
		orig = MakeVector(0.,0.,-1.e9);		// //  < -1.e8 means zlevel was not specified
        yflipped=false;
		aScaling=ReadUnits(attrs,LENGTH_UNITS);
        numAttr=(int)attrs.getLength();
#ifdef MPM_CODE
		angleAxes[0]=0;					// no angle files yet
		rotationAxes[0]=0;				// no rotations yet
        int maxRotations = is3D ? 3 : 1 ;
        bmpCustomPtsPerElement = 0;     // BMP commands can change it from using default
#endif
        for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"width")==0)
                sscanf(value,"%lf",&bwidth);
            else if(strcmp(aName,"height")==0)
                sscanf(value,"%lf",&bheight);
            else if(strcmp(aName,"name")==0)
				strcpy(bmpFileName,value);
#ifdef MPM_CODE
			// MPM allows up to three angle files
            else if(strcmp(aName,"angles")==0 || strcmp(aName,"anglesZ")==0)
			{	int rotNum = (int)strlen(angleAxes);
				if(rotNum==maxRotations)
					throw SAXException("Maximum of three angle files (for 3D) or one (for 2D) allowed in single BMP region.");
				strcpy(bmpAngleFileName[rotNum],value);
				strcat(angleAxes,"Z");
			}
            else if(strcmp(aName,"anglesX")==0 || strcmp(aName,"anglesY")==0)
			{	if(!is3D)
					throw SAXException("X and Y angle files only allowed for BMP regions in 3D simulations.");
				int rotNum = (int)strlen(angleAxes);
				if(rotNum==maxRotations)
					throw SAXException("Maximum of three angle files allowed in single BMP region.");
				strcpy(bmpAngleFileName[rotNum],value);
				angleAxes[rotNum]=aName[6];
				angleAxes[rotNum+1]=0;
			}
#else
			// FEA is 2D and can only rotate about z
            else if(strcmp(aName,"angles")==0 || strcmp(aName,"anglesZ")==0)
				strcpy(bmpAngleFileName[0],value);
#endif
            delete [] aName;
            delete [] value;
        }
		bwidth*=aScaling;
		bheight*=aScaling;
		if(bmpFileName[0]==0)
            throw SAXException("Bit-mapped file must specify the file in a name attribute.");
	}
	
    //-----------------------------------------------------------
    // Set origin on input data
    //-----------------------------------------------------------
    else if(strcmp(xName,"Origin")==0)
	{	ValidateCommand(xName,BMPBLOCK,ANY_DIM);
		aScaling=ReadUnits(attrs,LENGTH_UNITS);
        numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"x")==0)
			{	sscanf(value,"%lf",&orig.x);
				orig.x*=aScaling;
			}
            else if(strcmp(aName,"y")==0)
			{	sscanf(value,"%lf",&orig.y);
				orig.y*=aScaling;
			}
            else if(strcmp(aName,"z")==0)
			{	sscanf(value,"%lf",&orig.z);
				orig.z*=aScaling;
			}
            else if(strcmp(aName,"flipped")==0)
			{	if(strcmp(value,"yes")==0 || strcmp(value,"Yes")==0 || strcmp(value,"YES")==0 || strcmp(value,"1")==0 )
                    yflipped=true;
                else
                    yflipped=false;
			}
            delete [] aName;
            delete [] value;
        }
	}
	
    //-----------------------------------------------------------
    // Assign intensity to some material properties
    //-----------------------------------------------------------
    else if(strcmp(xName,"Intensity")==0)
	{	ValidateCommand(xName,BMPBLOCK,ANY_DIM);
        numAttr=(int)attrs.getLength();
		int mat=-1;
		int imin=-1;
		int imax=-1;
		double thisMinAngle=0.,maxAngle=0.;
		char matname[200];
		matname[0]=0;
        for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"mat")==0)
				sscanf(value,"%d",&mat);
			else if(strcmp(aName,"matname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(matname,value);
			}
            else if(strcmp(aName,"imin")==0)
                sscanf(value,"%d",&imin);
            else if(strcmp(aName,"imax")==0)
                sscanf(value,"%d",&imax);
            else if(strcmp(aName,"minAngle")==0)
                sscanf(value,"%lf",&thisMinAngle);
            else if(strcmp(aName,"maxAngle")==0)
                sscanf(value,"%lf",&maxAngle);
            delete [] aName;
            delete [] value;
        }
		// if gave a matname, it takes precedence over mat number
		if(strlen(matname)>0)
			mat = matCtrl->GetIDFromNewName(matname);
		if(imin<0 || imax<0)
			throw SAXException(XYFileError("<Intensity> has incomplete set of attributes.",bmpFileName));
		if(imin>imax)
			throw SAXException(XYFileError("<Intensity> range is not valid.",bmpFileName));
		if(mat>=0)
		{	BMPLevel *newLevel=new BMPLevel(mat,imin,imax);
			if(currentLevel==NULL)
				firstLevel=newLevel;
			else
				currentLevel->SetNextObject(newLevel);
			currentLevel=newLevel;
			block=INTENSITYBLOCK;
		}
		else
		{	if(numAngles==2)
				throw SAXException(XYFileError("Too many <Intensity> commands to set angle mappings.",bmpFileName));
            minAngle[numAngles]=thisMinAngle;
			angleScale[numAngles]=(maxAngle-thisMinAngle)/((double)imax-(double)imin);
			minIntensity[numAngles]=(double)imin;
			numAngles++;
		}
	}
	
	//-----------------------------------------------------------
    // Intensity properties for both MPM and FEA
    //-----------------------------------------------------------
	
	// thickness
    else if(strcmp(xName,"Thickness")==0)
	{	ValidateCommand(xName,INTENSITYBLOCK,MUST_BE_2D);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&currentLevel->thickness;
        gScaling=ReadUnits(attrs,LENGTH_UNITS);
    }
	
	// angle
    else if(strcmp(xName,"Angle")==0)
	{	ValidateCommand(xName,INTENSITYBLOCK,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&currentLevel->angle;
    }
	
	// temperature
    else if(strcmp(xName,"Temperature")==0)
	{
#ifdef FEA_CODE
		if(block==THERMAL) return false;
#endif
		ValidateCommand(xName,INTENSITYBLOCK,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&currentLevel->temperature;
    }
	
	else
		return false;
	
	return true;
}

//-----------------------------------------------------------
// Subroutine to set block on element end
// throws std::bad_alloc, SAXException()
//-----------------------------------------------------------
short CommonReadHandler::EndBMPInput(char *xName,int exitBlock)
{
    if(strcmp(xName,"BMP")==0)
	{	// make last intensity hole material
		if(currentLevel==NULL)
			throw SAXException(XYFileError("No intensity levels for materials were defined for <BMP> file",bmpFileName));
		
		// add level to grab left over intensities
		BMPLevel *newLevel=new BMPLevel(0,0,255);
		currentLevel->SetNextObject(newLevel);

		TranslateBMPFiles();		// scan file for material points or elements
        block=exitBlock;			// Must have been in MATERIALPOINTS
		
		// delete defined levels
		while(firstLevel!=NULL)
		{	currentLevel=(BMPLevel *)firstLevel->GetNextObject();
			delete firstLevel;
			firstLevel=currentLevel;
		}
		currentLevel=NULL;
    }

    else if(strcmp(xName,"Intensity")==0)
	{	block=BMPBLOCK;
	}
		
    else if(strcmp(xName,"Origin")==0 || strcmp(xName,"Thickness")==0 || strcmp(xName,"Angle")==0
				|| strcmp(xName,"Concentration")==0 || strcmp(xName,"Temperature")==0)
	{	// nothing to do
	}
	
	// not a BMP block element
	else
		return false;
		
	return true;
}

//-------------------------------------------------------------------
// Subroutine to translate BMP file into material points or elements
//-------------------------------------------------------------------
void CommonReadHandler::TranslateBMPFiles(void)
{
}

// set index into the BMP file, but make sure within bounds 0 to indexMax-1
int CommonReadHandler::BMPIndex(double value,int indexMax)
{
	int index=(int)value;
	if(index<0)
		index=0;
	else if(index>=indexMax)
		index=indexMax-1;
	return index;
}

//-----------------------------------------------------------
// Subroutine to read BMP file
// throws std::bad_alloc, SAXException()
//-----------------------------------------------------------
void *CommonReadHandler::ReadXYFile(char *bmpFullPath,XYInfoHeader &info,bool requireBytes,bool isReadingXML)
{
	// get file extension
	char ext[11];
	int dataType = BYTE_DATA;
	GetFileExtension(bmpFullPath,ext,10);
	
	// get file type from extension (case insensitive)
	XYFileImporter *xyFile = NULL;
	if(CIstrcmp(ext,"bmp")==0)
		xyFile = (XYFileImporter *)new XYBMPImporter(bmpFullPath,isReadingXML,dataType);
	else if(CIstrcmp(ext,"txt")==0 && !requireBytes)
	{	dataType = TEXT_DELIMITED;
		xyFile = (XYFileImporter *)new XYTXTImporter(bmpFullPath,isReadingXML,dataType);
	}
	else
	{	char *errMsg = XYFileError("The bit mapped file type is not recognized.",bmpFullPath);
		if(isReadingXML)
			throw SAXException(errMsg);
		else
			throw CommonException(errMsg,"CommonReadHandler::ReadXYFile");
	}
	
	// Read the file header
	xyFile->GetXYFileHeader(info);

	// buffer to read the file
	if(dataType==BYTE_DATA)
	{	unsigned char **rows=new (nothrow) unsigned char *[info.height];
		if(rows==NULL)
			xyFile->XYFileImportError("Out of memory reading bit mapped file.");

		for(int row=0;row<info.height;row++)
		{	rows[row]=new (nothrow) unsigned char[info.width];
			if(rows[row]==NULL)
				xyFile->XYFileImportError("Out of memory reading bit mapped file.");
		}
		
		// read the data
		xyFile->ReadXYFileData(rows,info);
	
		// close the file
		delete xyFile;
	
		// return the data pointer
		return rows;
	}
	
	else
	{	// read to floats (to save space)
		float **rows = new (nothrow) float *[info.height];
		if(rows==NULL)
			xyFile->XYFileImportError("Out of memory reading text file.");

		for(int row=0;row<info.height;row++)
		{	rows[row]=new (nothrow) float[info.width];
			if(rows[row]==NULL)
				xyFile->XYFileImportError("Out of memory reading bit text file.");
		}
		
		// read the data
		xyFile->ReadXYFileData((unsigned char **)rows,info);
	
		// close the file
		delete xyFile;
	
		// return the data pointer
		return rows;
	}
	
	return NULL;
}

// create error message with bmp file name
char *CommonReadHandler::XYFileError(const char *msg,const char *fileName)
{
	char *error=new (nothrow) char[strlen(msg)+strlen(fileName)+10];
	if(error != NULL)
	{	strcpy(error,msg);
		strcat(error," (file: ");
		strcat(error,fileName);
		strcat(error,")");
		return error;
	}
	else
		return (char *)msg;
}

// Give input width, height, zlevel and file info, decode into final image width, height, and zlevel and
// x and y length per pixel (note - dimensions are not scaled from mm)
// If error, return string
const char *CommonReadHandler::DecodeBMPWidthAndHeight(XYInfoHeader info,double &width,double &height,double &zlevel,Vector &pw,bool is3D)
{
	// height and width from input (<-1.e8 means not provided)
	//   other negative value means it is mm/pixel, positive is total size
	// if not provided, file might have it, but if provided, file not used
	pw.x=-1.;
	pw.y=-1.;
	if(width<-1.e8 && height<-1.e8)
	{	if(!info.knowsCellSize)
			return "Image file settings must specify width and/or height as size or pixels per mm.";
		width = info.width*info.xcell;
		height = info.height*info.ycell;
		pw.x = info.xcell;
		pw.y = info.ycell;
	}
	else
	{	// here means one or both was provided, but may be negative
        // if negative give -(length unts per pixel) (or data cell in text files)
		if(height<0. && height>-1.e8)
		{	pw.y = -height;
			height = pw.y*(double)info.height;
		}
		if(width<0. && width>-1.e8)
		{	pw.x = -width;
			width = pw.x*(double)info.width;
		}
	}
	
	// total dimensions (if only one is known, find the other, never have both unknown)
	if(height<0) height = width*(double)info.height/(double)info.width;
	if(width<0) width = height*(double)info.width/(double)info.height;
	
	// final length per pixel (if needed)
	if(pw.x<0.) pw.x = width/(double)info.width;
	if(pw.y<0.) pw.y = height/(double)info.height;
	
	// set zlevel (but only for 3D)
	if(is3D)
	{	// if not set, try to find it in the file info
		if(zlevel<-1.e8)
		{	if(!info.knowsCellSize)
				zlevel = 0.;
			else
				zlevel = info.zlevel;
		}
	}
	else
		zlevel = 0.;
	
	return NULL;
}

// Input: info is from the bit mapped file reading
//		Need: height, width
// spot is position to inspect
// orig = (xorigin,yorigin,zlevel)
// del = (deltax,deltay,deltaz) (deltaz>0 for 3D)
// pw = (xpw,ypw,+/-1) - pixel width, last is -1 to flip
// width and height of image
// return row and column ranges and weights of first and last one
bool CommonReadHandler::MapDomainToImage(XYInfoHeader info,Vector spot,Vector orig,Vector del,Vector pw,
										 double width,double height,DomainMap &map)
{
	// if point in the view area, then check it
	if((spot.x>=orig.x && spot.x<orig.x+width) && (spot.y>=orig.y && spot.y<orig.y+height))
	{	// for 3D, check z as well
		if(del.z>0. && (spot.z<orig.z-del.z || spot.z>=orig.z+del.z)) return false;
		
		// find range of rows and cols for pixels over this material point
		double rmin,rmax;
		if(pw.z<0.)
		{   rmin = (orig.y+height-spot.y-del.y)/pw.y;
			rmax = (orig.y+height-spot.y-del.y)/pw.y;
		}
		else
		{   rmin = (spot.y-del.y-orig.y)/pw.y;
			rmax = (spot.y+del.y-orig.y)/pw.y;
		}
		map.r1 = BMPIndex(rmin,info.height);
		map.r2 = BMPIndex(rmax,info.height);
		if(map.r2==map.r1)
		{	// a single row
			map.wtr1 = map.wtr2 = 1.;
		}
		else
		{	// fractional weight for first and last row in case not all within the extent
			map.wtr1 = (double)(map.r1+1)-rmin;
			map.wtr2 = rmax-(double)map.r2;
		}
		
		// ignore excellent alignment
		if(map.wtr1<WEIGHT_TOL)
		{	map.r1++;
			map.wtr1 = 1.;
		}
		if(map.wtr2<WEIGHT_TOL)
		{	map.r2--;
			map.wtr2 = 1.;
		}
		
		
		double cmin=(spot.x-del.x-orig.x)/pw.x;
		double cmax=(spot.x+del.x-orig.x)/pw.y;
		map.c1 = BMPIndex(cmin,info.width);
		map.c2 = BMPIndex(cmax,info.width);
		if(map.c2==map.c1)
		{	// a single column
			map.wtc1 = map.wtc2 = 1.;
		}
		else
		{	// fractional weight for first and last col in case not all within the extent
			map.wtc1 = (double)(map.c1+1)-cmin;
			map.wtc2 = cmax-(double)map.c2;
		}
		
		// ignore excellent alignment
		if(map.wtc1<WEIGHT_TOL)
		{	map.c1++;
			map.wtc1 = 1.;
		}
		if(map.wtc2<WEIGHT_TOL)
		{	map.c2--;
			map.wtc2 = 1.;
		}
		
		return true;
	}
	
	// not in the image
	return false;
}

// Find the most prominent level withing a domain on top of and bmp grid
BMPLevel *CommonReadHandler::FindBMPLevel(BMPLevel *startLevel,DomainMap map,unsigned char **rows)
{
	// clear level weights
	BMPLevel *nextLevel = startLevel;
	while(nextLevel!=NULL) nextLevel = nextLevel->ClearWeight();
		
	// add weights
	double rweight,weight;
	for(int row=map.r1;row<=map.r2;row++)
	{	// get weight for this row (first and last may differ)
		if(row==map.r1)
			rweight=map.wtr1;
		else if(row==map.r2)
			rweight=map.wtr2;
		else
			rweight=1.;
		for(int col=map.c1;col<=map.c2;col++)
		{	// get weight for this column (first and last may differ)
			if(col==map.c1)
				weight=rweight*map.wtc1;
			else if(col==map.c2)
				weight=rweight*map.wtc2;
			else
				weight=rweight;
				
			// find ID (i.e., material) at this level
			// (note: last level with ID=0 catches empty space when doing materials)
			nextLevel=startLevel;
			while(nextLevel!=NULL)
			{	int levID=nextLevel->Material(rows[row][col],weight);
				if(levID>=0) break;
				nextLevel=(BMPLevel *)nextLevel->GetNextObject();
			}
		}
	}
	
	// find level with maximum weight (may be empty space when doing materials)
	// it will be NULL if all weights are zero
	weight=0.;
	BMPLevel *maxLevel = NULL;
	nextLevel = startLevel;
	while(nextLevel!=NULL) nextLevel=nextLevel->MaximumWeight(weight,&maxLevel);
	return maxLevel;
}

// Search domain and find weight average of pixel values in the domain
// If find any value calculate average and return true, otherwise return false
double CommonReadHandler::FindAverageValue(DomainMap map,unsigned char **rows,int dataType,bool &hasWeight)
{
	// initialize
	double totalIntensity=0.;
	double totalWeight=0.;
	
	// scan the domain
	double rweight,weight;
	for(int row=map.r1;row<=map.r2;row++)
	{	// get weight for row (first and last may differ
		if(row==map.r1)
			rweight=map.wtr1;
		else if(row==map.r2)
			rweight=map.wtr2;
		else
			rweight=1.;
		
		for(int col=map.c1;col<=map.c2;col++)
		{	// get weigth for column (first and last may differ
			if(col==map.c1)
				weight=rweight*map.wtc1;
			else if(col==map.c2)
				weight=rweight*map.wtc2;
			else
				weight=rweight;
			
			// add pixel value
			if(dataType==BYTE_DATA)
				totalIntensity += weight*rows[row][col];
			else
				totalIntensity += weight*((float **)rows)[row][col];
			totalWeight += weight;
		}
	}
	
	// get weight average value
	if(totalWeight>0.)
	{	hasWeight = true;
		totalIntensity /= totalWeight;
	}
	return totalIntensity;
}


