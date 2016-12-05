/********************************************************************************
    BitMapFiles.cpp - extra code for MPMReadHandler.cpp to generate particles
		from gray scale BMP files
    nairn-mpm-fea

    Created by John Nairn on Thu Dec 8 2004.
    Copyright (c) 2004 RSAC Software. All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Read_MPM/MPMReadHandler.hpp"
#include "Read_XML/BMPLevel.hpp"
#include "MPM_Classes/MatPoint2D.hpp"
#include "MPM_Classes/MatPointAS.hpp"
#include "MPM_Classes/MatPoint3D.hpp"
#include "Elements/ElementBase.hpp"
#include "Read_MPM/MpsController.hpp"
#include "System/ArchiveData.hpp"
#include "Read_XML/mathexpr.hpp"

extern char *angleExpr[3];
extern char rotationAxes[4];

//-----------------------------------------------------------
// Check for bmp element, return false if not
//-----------------------------------------------------------
short MPMReadHandler::BMPFileInput(char *xName,const Attributes& attrs)
{
	// check for common commands
	if(BMPFileCommonInput(xName,attrs,POINTSBLOCK)) return TRUE;
	
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
// Subroutine to translate BMP file into material point
// throws std::bad_alloc, SAXException()
//-----------------------------------------------------------
void MPMReadHandler::TranslateBMPFiles(void)
{
	// file info and data
	unsigned char **rows,**angleRows = NULL;
	XYInfoHeader info,angleInfo;
	
	// read image file
	char *bmpFullPath=archiver->ExpandOutputPath(bmpFileName);
	ReadBMPFile(bmpFullPath,info,&rows);
	delete [] bmpFullPath;
	
	// angle file name (overrides other angle settings)
	bool setAngles = false;
	int numRotations=(int)strlen(rotationAxes);
	if(bmpAngleFileName[0]>0)
	{	setAngles = true;
		char *bmpFullAnglePath=archiver->ExpandOutputPath(bmpAngleFileName);
		ReadBMPFile(bmpFullAnglePath,angleInfo,&angleRows);
		if(info.height!=angleInfo.height || info.width!=angleInfo.width)
			throw SAXException(BMPError("The image file and angle file sizes do not match.",bmpFileName));
		delete [] bmpFullAnglePath;
	}
	else if(numRotations>0)
	{	int i;
		for(i=0;i<numRotations;i++)
		{	char *expr=new char[strlen(angleExpr[i])+1];
			strcpy(expr,angleExpr[i]);
			if(!CreateFunction(expr,i+1))
				throw SAXException("Invalid angle expression was provided in <RotateX(YZ)> command.");
			delete [] expr;
		}
	}
	
	// get final width an height
	Vector pw;
	const char *msg = CommonReadHandler::DecodeBMPWidthAndHeight(info,bwidth,bheight,orig.z,pw,fmobj->IsThreeD());
	if(msg != NULL)
		throw SAXException(BMPError("<BMP> command must specify width and/or height as size or pixels per mm.",bmpFileName));
	pw.z = yflipped ? -1. : 1. ;
	
	// Length/semiscale is half particle with
	//	(semiscale=4 for 2D w 4 pts per element or 3D with 8 pts per element,
	//		or 2 (2D and 3D) if 1 particle in the element)
	double semiscale;
	if(fmobj->IsThreeD())
		semiscale=2.*pow((double)fmobj->ptsPerElement,1./3.);
	else
		semiscale=2.*sqrt((double)fmobj->ptsPerElement);
    
	// variables
	Vector mpos[MaxElParticles];
	DomainMap map;
	
    /* Parallelizing the following loop will speed up check meshes on Mac
        1. Will need copy of levels for each block (or pass copy of weight to BMPLevel methods)
            see BMPLevel methods: ClearWeight(),Material(double,double), MaximumWeight(double)
        2. mpCtrl->AddMaterialPoint() will need critical block
        3. FunctionValue() when setting angles uses globals they are problem for parallel
    */
	
	// scan mesh and assign material points or angles
    try
    {   for(int ii=1;ii<=nelems;ii++)
        {	// skip if image not in extent of element box
            ElementBase *elem=theElements[ii-1];
            if(!elem->IntersectsBox(orig,bwidth,bheight))
                continue;
            
            // load point coordinates
            elem->MPMPoints(fmobj->ptsPerElement,mpos);
        
            // particle size withing volume of the element
			Vector del;
            del.x=(elem->GetDeltaX())/semiscale;
            del.y=(elem->GetDeltaY())/semiscale;
            if(fmobj->IsThreeD())
                del.z=elem->GetDeltaZ()/semiscale;
            else
                del.z=-1.;
			
            for(int k=0;k<fmobj->ptsPerElement;k++)
            {	int ptFlag=1<<k;
            
                // skip if already filled
                if(elem->filled&ptFlag) continue;

                // if point in the view area, then check it
				if(MapDomainToImage(info,mpos[k],orig,del,pw,bwidth,bheight,map))
				{	// find maximum level and its material ID or none (a hole)
					int matID=-1;
					BMPLevel *nextLevel = FindBMPLevel(firstLevel,map,rows);
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
						newMpt->SetDimensionlessByPts(fmobj->ptsPerElement);
                        mpCtrl->AddMaterialPoint(newMpt,nextLevel->concentration,nextLevel->temperature);
                        
						// is there an angle image too?
						if(setAngles)
						{	double totalIntensity = FindAverageValue(map,rows);
							if(totalIntensity>0.)
							{	double matAngle=minAngle+(totalIntensity-minIntensity)*angleScale;
								newMpt->SetAnglez0InDegrees(matAngle);
							}
						}
						else
						{	// If had Rotate commands then use them
							SetMptAnglesFromFunctions(numRotations,&mpos[k],newMpt);
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
		if(setAngles) delete [] angleRows[row];
	}
	delete [] rows;
	if(setAngles) delete [] angleRows;
	
	// angles if allocated
	for(int ii=0;ii<numRotations;ii++)
	{	delete [] angleExpr[ii];
	}
	DeleteFunction(-1);
		
}

// set current intensity velocity
void MPMReadHandler::SetLevelVelocity(double vx,double vy,double vz)
{
	currentLevel->vel.x=vx;
	currentLevel->vel.y=vy;
	currentLevel->vel.z=vz;
}

