/********************************************************************************
    BitMapFiles.cpp - extra code for MPMReadHandler.cpp to generate particles
		from gray scale BMP files
    nairn-mpm-fea

    Created by John Nairn on Thu Dec 8 2004.
    Copyright (c) 2004 RSAC Software. All rights reserved.
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Read_MPM/MPMReadHandler.hpp"
#include "Read_XML/BMPLevel.hpp"
#include "MPM_Classes/MatPoint2D.hpp"
#include "MPM_Classes/MatPointAS.hpp"
#include "MPM_Classes/MatPoint3D.hpp"
#include "Elements/ElementBase.hpp"
#include "Read_MPM/MpsController.hpp"
#include "System/ArchiveData.hpp"
#include "Read_XML/Expression.hpp"
#include "System/MPMPrefix.hpp"

extern char *angleExpr[3];
extern char rotationAxes[4];
extern char angleAxes[4];

//-----------------------------------------------------------
// Check for bmp element, return false if not
//-----------------------------------------------------------
short MPMReadHandler::BMPFileInput(char *xName,const Attributes& attrs)
{
	// check for common commands
	if(BMPFileCommonInput(xName,attrs,POINTSBLOCK,fmobj->IsThreeD())) return true;
	
	//-----------------------------------------------------------
    // Intensity properties for MPM only
	//		vel is handled in main MPMReadHandler
    //-----------------------------------------------------------
	
	// concentration
    if(strcmp(xName,"Concentration")==0)
	{	ValidateCommand(xName,INTENSITYBLOCK,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&currentLevel->concentration;
    }
	
	// not a BMP file element
	else
		return FALSE;

	return TRUE;
}

//-----------------------------------------------------------
// Subroutine to translate BMP file into material points
// throws std::bad_alloc, SAXException()
//-----------------------------------------------------------
void MPMReadHandler::TranslateBMPFiles(void)
{
	// file info and data
	unsigned char **rows,**angleRows = NULL,**angle2Rows = NULL,**angle3Rows = NULL;
	XYInfoHeader info,angleInfo,angleInfo2,angleInfo3;
	
	// read image file
	char *bmpFullPath=archiver->ExpandOutputPath(bmpFileName);
    rows = (unsigned char **)ReadXYFile(bmpFullPath,info,false,true);
	delete [] bmpFullPath;
	
	// angle file name (overrides other angle settings)
	bool setAngles = false;
	int numRotations=(int)strlen(rotationAxes);
	int fileRotations=(int)strlen(angleAxes);
	if(bmpAngleFileName[0][0]>0)
	{	setAngles = true;
		
		// first file always there
		char *bmpFullAnglePath=archiver->ExpandOutputPath(bmpAngleFileName[0]);
		angleRows = (unsigned char **)ReadXYFile(bmpFullAnglePath,angleInfo,true,true);
		if(info.height!=angleInfo.height || info.width!=angleInfo.width)
			throw SAXException(XYFileError("The image file and first angle file sizes do not match.",bmpFileName));
		delete [] bmpFullAnglePath;
		
		// was angle mapping set
		if(numAngles==0)
			throw SAXException(XYFileError("No mapping of pixels to angles for angle file were provided.",bmpFileName));
		
		if(fileRotations>1)
		{	bmpFullAnglePath=archiver->ExpandOutputPath(bmpAngleFileName[1]);
			angle2Rows = (unsigned char **)ReadXYFile(bmpFullAnglePath,angleInfo2,true,true);
			if(info.height!=angleInfo2.height || info.width!=angleInfo2.width)
				throw SAXException(XYFileError("The image file and second angle file sizes do not match.",bmpFileName));
			delete [] bmpFullAnglePath;
			if(numAngles<2)
			{	minAngle[1] = minAngle[0];
				minIntensity[1] = minIntensity[0];
				angleScale[1] = angleScale[0];
				numAngles++;
			}
		}
		
		if(fileRotations>2)
		{	bmpFullAnglePath=archiver->ExpandOutputPath(bmpAngleFileName[2]);
			angle3Rows = (unsigned char **)ReadXYFile(bmpFullAnglePath,angleInfo3,true,true);
			if(info.height!=angleInfo3.height || info.width!=angleInfo3.width)
				throw SAXException(XYFileError("The image file and second angle file sizes do not match.",bmpFileName));
			delete [] bmpFullAnglePath;
			if(numAngles<3)
			{	minAngle[2] = minAngle[1];
				minIntensity[2] = minIntensity[1];
				angleScale[2] = angleScale[1];
				numAngles++;
			}
		}
	}
	else if(numRotations>0)
	{	int i;
		for(i=0;i<numRotations;i++)
		{	if(!Expression::CreateFunction(angleExpr[i],i+1))
				throw SAXException("Invalid angle expression was provided in <RotateX(YZ)> command.");
		}
	}
	
	// get final width an height
	Vector pw;
	const char *msg = CommonReadHandler::DecodeBMPWidthAndHeight(info,bwidth,bheight,orig.z,pw,fmobj->IsThreeD());
	if(msg != NULL)
		throw SAXException(XYFileError("<BMP> command must specify width and/or height as size or pixels per mm.",bmpFileName));
    if(yflipped)
        pw.z = info.topDown==0 ? -1. : 1. ;
    else
        pw.z = info.topDown==0 ? 1. : -1. ;

    // if custom points, must be valid
    int usePtsPerElement = fmobj->ptsPerElement;
    int usePtsPerSide = fmobj->ptsPerSide;
    if(bmpCustomPtsPerElement>0)
    {   // switch
        usePtsPerElement = bmpCustomPtsPerElement;
        usePtsPerSide = bmpCustomPtsPerSide;
    }
    int numFound;
    
    // Length/semiscale is half particle with
    //    (semiscale=4 for 2D w 4 pts per element or 3D with 8 pts per element,
    //        or 2 (2D and 3D) if 1 particle in the element)
    double semiscale = 2*(double)usePtsPerSide;
    
	// variables
	Vector *mpos = new Vector[usePtsPerElement];
	DomainMap map;
	
    /* Parallelizing the following loop will speed up check meshes on Mac
        1. Will need copy of levels for each block (or pass copy of weight to BMPLevel methods)
            see BMPLevel methods: ClearWeight(),Material(double,double), MaximumWeight(double)
        2. mpCtrl->AddMaterialPoint() will need critical block
        3. FunctionMValue() when setting angles uses globals they are problem for parallel
    */
	
	// scan mesh and assign material points or angles
    int ptFlag=0;
    try
    {   for(int ii=1;ii<=nelems;ii++)
        {	// skip if image not in extent of element box
            ElementBase *elem=theElements[ii-1];
            if(!elem->IntersectsBox(orig,bwidth,bheight))
                continue;

            // load point coordinates
            elem->MPMPoints(usePtsPerSide,mpos,numFound,NULL,NULL);
        
            // particle radius (dimensioned) within volume of the element
			Vector del;
            del.x=(elem->GetDeltaX())/semiscale;
            del.y=(elem->GetDeltaY())/semiscale;
            if(fmobj->IsThreeD())
                del.z=elem->GetDeltaZ()/semiscale;
            else
                del.z=-1.;
 			
			// fill all points
            for(int k=0;k<numFound;k++)
            {   // default locations (skip if already filled)
                if(usePtsPerElement==fmobj->ptsPerElement)
                {   ptFlag=1<<k;
                    if(elem->filled&ptFlag) continue;
                }

                // if point in the view area, then check it
				if(MapDomainToImage(info,mpos[k],orig,del,pw,bwidth,bheight,map))
				{	// find maximum level and its material ID or none (a hole)
					int matID=-1;
					BMPLevel *nextLevel = FindBMPLevel(firstLevel,map,rows,info.dataType);
                    if(nextLevel!=NULL) matID = nextLevel->Material();
 					
                    // create a material point if one at this spot using matID and nextLevel
                    // Note that empty spaced is not marked as filled which allows superposition of
                    // images with different materials. If want to forcefully create a hole that
                    // cannot be filled by subsequent image, will need to define a new material
                    // type that can have matID for a hole. It will not create a point, but will mark
                    // the location as filled
					MPMBase *newMpt;
                    if(matID>0)
                    {	if(fmobj->IsThreeD())
                            newMpt=new MatPoint3D(ii,matID,nextLevel->angle);
                        else if(fmobj->IsAxisymmetric())
                            newMpt=new MatPointAS(ii,matID,nextLevel->angle,mpos[k].x);
                        else
                            newMpt=new MatPoint2D(ii,matID,nextLevel->angle,nextLevel->thickness);
                        newMpt->SetPosition(&mpos[k]);
                        newMpt->SetOrigin(&mpos[k]);
                        newMpt->SetVelocity(&nextLevel->vel);
						newMpt->SetDimensionlessByPts(usePtsPerSide);
                        mpCtrl->AddMaterialPoint(newMpt,nextLevel->concentration,nextLevel->temperature);
                        
						// is there an angle image too?
						if(setAngles)
						{	double matAngle[3];
							bool hasWeight = false;
							double totalIntensity = FindAverageValue(map,angleRows,angleInfo.dataType,hasWeight);
							matAngle[0] = minAngle[0]+(totalIntensity-minIntensity[0])*angleScale[0];
 							if(fileRotations>1)
							{	totalIntensity = FindAverageValue(map,angle2Rows,angleInfo2.dataType,hasWeight);
								matAngle[1] = minAngle[1]+(totalIntensity-minIntensity[1])*angleScale[1];
							}
							if(fileRotations>2)
							{	totalIntensity = FindAverageValue(map,angle3Rows,angleInfo3.dataType,hasWeight);
								matAngle[2] = minAngle[2]+(totalIntensity-minIntensity[2])*angleScale[2];
							}
							SetMptAnglesFromFunctions(angleAxes,matAngle,&mpos[k],newMpt);
						}
						else
						{	// If had Rotate commands then use them
							SetMptAnglesFromFunctions(rotationAxes,NULL,&mpos[k],newMpt);
						}
						
						// fill the spot
                        elem->filled|=ptFlag;
                    }
                }
            }
        }
    }
    catch(const char *msg)
    {   throw SAXException(msg);
    }
	
	// clean up
	for(int row=0;row<info.height;row++)
	{	delete [] rows[row];
		if(setAngles)
		{	delete [] angleRows[row];
			if(fileRotations>1) delete [] angle2Rows[row];
			if(fileRotations>2) delete [] angle3Rows[row];
		}
	}
	delete [] rows;
	if(setAngles)
	{	delete [] angleRows;
		if(fileRotations>1) delete [] angle2Rows;
		if(fileRotations>2) delete [] angle3Rows;
	}
	
	// angles if allocated
	for(int ii=0;ii<numRotations;ii++)
	{	delete [] angleExpr[ii];
	}
	Expression::DeleteFunction(-1);
	
	// vector array
	delete [] mpos;
		
}

// set current intensity velocity
void MPMReadHandler::SetLevelVelocity(double vx,double vy,double vz)
{
	currentLevel->vel.x=vx;
	currentLevel->vel.y=vy;
	currentLevel->vel.z=vz;
}

