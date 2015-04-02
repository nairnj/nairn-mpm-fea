/********************************************************************************
    BitMapFilesCommon.cpp - extra code for MPMReadHandler.cpp to generate particles
		from gray scale BMP files
    nairn-mpm-fea

    Created by John Nairn on Thu Dec 8 2004.
    Copyright (c) 2004 RSAC Software. All rights reserved.
********************************************************************************/

#include <fstream>

#include "Read_XML/CommonReadHandler.hpp"
#include "Read_XML/BMPLevel.hpp"
#include "Read_XML/MaterialController.hpp"

#ifdef MPM_CODE
extern char rotationAxes[4];
#endif

// file type allowed in BMP command
enum { BMP_INPUT_FILE=0 };

//-----------------------------------------------------------
// Check for bmp element, return false if not
//-----------------------------------------------------------
short CommonReadHandler::BMPFileCommonInput(char *xName,const Attributes& attrs,int expectedBlock)
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
		bmpAngleFileName[0]=0;
		xorig=yorig=0.;
        yflipped=false;
		zslice=0.;
		aScaling=ReadUnits(attrs,LENGTH_UNITS);
        numAttr=attrs.getLength();
#ifdef MPM_CODE
		rotationAxes[0]=0;				// no rotations yet
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
            else if(strcmp(aName,"angles")==0)
				strcpy(bmpAngleFileName,value);
            delete [] aName;
            delete [] value;
        }
		bwidth*=aScaling;
		bheight*=aScaling;
		if(bmpFileName[0]==0)
            throw SAXException("<BMP> must specify the file in a name attribute.");
	}
	
    //-----------------------------------------------------------
    // Set origin on input data
    //-----------------------------------------------------------
    else if(strcmp(xName,"Origin")==0)
	{	ValidateCommand(xName,BMPBLOCK,ANY_DIM);
		aScaling=ReadUnits(attrs,LENGTH_UNITS);
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"x")==0)
			{	sscanf(value,"%lf",&xorig);
				xorig*=aScaling;
			}
            else if(strcmp(aName,"y")==0)
			{	sscanf(value,"%lf",&yorig);
				yorig*=aScaling;
			}
            else if(strcmp(aName,"z")==0)
			{	sscanf(value,"%lf",&zslice);
				zslice*=aScaling;
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
        numAttr=attrs.getLength();
		int mat=-1;
		int imin=-1;
		int imax=-1;
		minAngle=0.;
		double maxAngle=0.;
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
                sscanf(value,"%lf",&minAngle);
            else if(strcmp(aName,"maxAngle")==0)
                sscanf(value,"%lf",&maxAngle);
            delete [] aName;
            delete [] value;
        }
		// if gave a matname, it takes precedence over mat number
		if(strlen(matname)>0)
			mat = matCtrl->GetIDFromNewName(matname);
		if(imin<0 || imax<0)
			throw SAXException(BMPError("<Intensity> has incomplete set of attributes.",bmpFileName));
		if(imin>imax)
			throw SAXException(BMPError("<Intensity> range is not valid.",bmpFileName));
		if(mat>=0)
		{	BMPLevel *newLevel=new BMPLevel(mat,imin,imax);
			if(newLevel==NULL)
				throw SAXException(BMPError("<Intensity> failed due to memory error.",bmpFileName));
			if(currentLevel==NULL)
				firstLevel=newLevel;
			else
				currentLevel->SetNextObject(newLevel);
			currentLevel=newLevel;
			block=INTENSITYBLOCK;
		}
		else
		{	angleScale=(maxAngle-minAngle)/((double)imax-(double)imin);
			minIntensity=(double)imin;
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
//-----------------------------------------------------------
short CommonReadHandler::EndBMPInput(char *xName,int exitBlock)
{
    if(strcmp(xName,"BMP")==0)
	{	// make last intensity hole material
		if(currentLevel==NULL)
			throw SAXException(BMPError("No intensity levels for materials were defined for <BMP> file",bmpFileName));
		
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
	{	if(currentLevel->concentration<0. || currentLevel->concentration>1.)
			throw SAXException(BMPError("Particle concentration potential for an intensity level must be from 0 to 1",bmpFileName));
        block=BMPBLOCK;
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
	int index=(int)floor(value);
	if(index<0)
		index=0;
	else if(index>=indexMax)
		index=indexMax-1;
	return index;
}

//-----------------------------------------------------------
// Subroutine to read BMP file
//-----------------------------------------------------------
void CommonReadHandler::ReadBMPFile(char *bmpFullPath,BMPInfoHeader &info,unsigned char ***theData)
{
	short mustReverse=false,status=true;
	FILE *fp;
	unsigned int i,ncolors;
	int row,col,rowBytes;
	
	// decide what computer I am on
    // byte order marker
    int test=1;
    char *testPtr=(char *)&test;
    if(*testPtr==0) mustReverse=true;
	
	// open the file
	if((fp=fopen(bmpFullPath,"r"))==NULL)
		throw SAXException(BMPError("The <BMP> file could not be opened.",bmpFullPath));
	
	// get file type from extension (all upper or all lower case)
	// If no extension, or unknown extension, assume to be BMP file
	int fileType = BMP_INPUT_FILE;
	int byteSize = 1;
	int endLoc = strlen(bmpFullPath)-1;
	int dotLoc = endLoc-1;
	while(dotLoc>0 && bmpFullPath[dotLoc]!='.') dotLoc--;
	if(dotLoc>0)
	{	char ext[11];
		dotLoc++;
		int ei = 0;
		while(ei<11 && dotLoc<=endLoc) ext[ei++]=bmpFullPath[dotLoc++];
		
		// look for supported extensions besides BMP
	}
	
	// Read the file header (depending on file type)
	
	if(fileType==BMP_INPUT_FILE)
	{	// read and check the header (individual reads due to Windows alignment issues)
		BMPHeader header;
		if(fread(header.type,2,1,fp)<1) status=false;
		if(fread(&header.size,4,1,fp)<1) status=false;
		if(fread(&header.reserved1,2,1,fp)<1) status=false;
		if(fread(&header.reserved2,2,1,fp)<1) status=false;
		if(fread(&header.offset,4,1,fp)<1) status=false;
		if(!status)
		{	fclose(fp);
			throw SAXException(BMPError("Error reading the specified <BMP> file.",bmpFileName));
		}
		if(mustReverse)
		{	Reverse((char *)&header.size,sizeof(int));
			Reverse((char *)&header.offset,sizeof(int));
		}
		
		// precalculate number of colors from header data
		ncolors=(header.offset-14-40)>>2;
		
		// exit if does not look like BMP file
		if(header.type[0]!='B' || header.type[1]!='M')
		{	fclose(fp);
			throw SAXException(BMPError("<BMP> file is not a valid bit map file.",bmpFileName));
		}
		
		// read information
		if(fread(&info.size,4,1,fp)<1) status=false;
		if(fread(&info.width,4,1,fp)<1) status=false;
		if(fread(&info.height,4,1,fp)<1) status=false;
		if(fread(&info.planes,2,1,fp)<1) status=false;
		if(fread(&info.bits,2,1,fp)<1) status=false;
		if(fread(&info.compression,4,1,fp)<1) status=false;
		if(fread(&info.imagesize,4,1,fp)<1) status=false;
		if(fread(&info.xresolution,4,1,fp)<1) status=false;
		if(fread(&info.yresolution,4,1,fp)<1) status=false;
		if(fread(&info.ncolors,4,1,fp)<1) status=false;
		if(fread(&info.importantcolors,4,1,fp)<1) status=false;
		if(!status)
		{	fclose(fp);
			throw SAXException(BMPError("Error reading the specified <BMP> file.",bmpFileName));
		}
		if(mustReverse)
		{	Reverse((char *)&info.size,sizeof(int));
			Reverse((char *)&info.width,sizeof(int));
			Reverse((char *)&info.height,sizeof(int));
			Reverse((char *)&info.planes,2);
			Reverse((char *)&info.bits,2);
			Reverse((char *)&info.compression,sizeof(int));
			Reverse((char *)&info.imagesize,sizeof(int));
			Reverse((char *)&info.xresolution,sizeof(int));
			Reverse((char *)&info.yresolution,sizeof(int));
			Reverse((char *)&info.ncolors,sizeof(int));
			Reverse((char *)&info.importantcolors,sizeof(int));
		}
		
		// correct for Photoshop lack of information
		if(info.imagesize==0)
			info.imagesize=header.size-header.offset;
		if(info.ncolors==0)
			info.ncolors=ncolors;
		
		// can I read this BMP file
		if(info.planes!=1)
		{	fclose(fp);
			throw SAXException(BMPError("Cannot read <BMP> files with more than 1 color plane.",bmpFileName));
		}
		if(info.compression!=0)
		{	fclose(fp);
			throw SAXException(BMPError("Cannot read compressed <BMP> files.",bmpFileName));
		}
		if(info.ncolors!=ncolors || (ncolors!=2 && ncolors!=16 && ncolors!=256))
		{	fclose(fp);
			throw SAXException(BMPError("<BMP> files does not appear to be 1, 4, or 8 bit grayscale image.",bmpFileName));
		}
		
		// read the 2, 16, or 256 colors
		unsigned char blueByte,greenByte,redByte,junkByte;
		for(i=0;i<ncolors;i++)
		{	if(fread(&blueByte,1,1,fp)<1) status=false;
			if(fread(&greenByte,1,1,fp)<1) status=false;
			if(fread(&redByte,1,1,fp)<1) status=false;
			if(fread(&junkByte,1,1,fp)<1) status=false;
			if(!status)
			{	fclose(fp);
				throw SAXException(BMPError("<BMP> color table is corrupted.",bmpFileName));
			}
			intensity[i]=((int)blueByte+(int)greenByte+(int)blueByte)/3;
		}
		
		// read one byte at a time (which might have multiple points) but each
		//	row is multiple of 4 bytes
		int colorsPerByte=1;
		if(ncolors==16)
			colorsPerByte=2;
		else if(ncolors==2)
			colorsPerByte=8;
		rowBytes=info.width/colorsPerByte;
		if(rowBytes%4!=0)
			rowBytes+=4-(rowBytes % 4);
		
		// BMP files do not know actual size
		info.knowsCellSize = false;
	}
	
	else
	{	// Read header of any other defined binary file types
	}
	
	// buffer to read the file
	unsigned char **rows=new unsigned char *[info.height];
	*theData=rows;
	if(rows==NULL)
	{	fclose(fp);
		throw SAXException(BMPError("Out of memory reading <BMP> file.",bmpFileName));
	}
	for(row=0;row<info.height;row++)
	{	rows[row]=new unsigned char[info.width];
		if(rows[row]==NULL)
		{	fclose(fp);
			throw SAXException(BMPError("Out of memory reading <BMP> file.",bmpFileName));
		}
	}
	
	// check size
	unsigned int expectSize=(unsigned int)(info.height*rowBytes);
	if(info.imagesize<expectSize)
	{	fclose(fp);
		throw SAXException(BMPError("<BMP> file size does not match expected data size.",bmpFileName));
	}
	
	// scan each row (depending on file type)
	if(fileType == BMP_INPUT_FILE)
	{	unsigned char dataByte,mask;
		int colByte;
		
		// loop over rows
		for(row=0;row<info.height;row++)
		{	// each row has chunks of rwoBytes bytes
			col=0;
			for(colByte=0;colByte<rowBytes;colByte++)
			{	// read next byte
				if(fread(&dataByte,byteSize,1,fp)<1)
				{	fclose(fp);
					throw SAXException(BMPError("Error reading image date of <BMP> file.",bmpFileName));
				}
				if(col>=info.width) continue;
				switch(ncolors)
				{	case 256:		// each byte is an index
						rows[row][col]=intensity[dataByte];
						col++;
						break;
					case 16:		// each byte has two indices
						rows[row][col]=intensity[(dataByte&0xF0)>>4];
						col++;
						if(col>=info.width) break;
						rows[row][col]=intensity[dataByte&0x0F];
						col++;
						break;
					case 2:
						mask=0x80;
						for(i=1;i<=8;i++)
						{	rows[row][col]=intensity[(dataByte&mask)>>(8-i)];
							col++;
							if(col>=info.width) break;
							mask>>=1;
						}
						break;
					default:
						break;
				}
			}
		}
	}
	
	else
	{	// Read any other supported binary type
	}
	
	// close the file
	fclose(fp);
}

// create error message with bmp file name
char *CommonReadHandler::BMPError(const char *msg,const char *fileName)
{
	char *error=new char[strlen(msg)+strlen(fileName)+10];
	strcpy(error,msg);
	strcat(error," (file: ");
	strcat(error,fileName);
	strcat(error,")");
	return error;
}


