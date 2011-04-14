/********************************************************************************
    BitMapFiles.cpp - extra code for MPMReadHandler.cpp to generate particles
		from gray scale BMP files
    NairnMPM

    Created by John Nairn on Thu Dec 8 2004.
    Copyright (c) 2004 RSAC Software. All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Read_MPM/MPMReadHandler.hpp"
#include "Read_XML/BMPLevel.hpp"
#include "MPM_Classes/MatPoint2D.hpp"
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
//-----------------------------------------------------------
void MPMReadHandler::TranslateBMPFiles(void)
{
	unsigned char **rows,**angleRows;
	BMPInfoHeader info,angleInfo;
	bool setAngles=FALSE;
	int numRotations=strlen(rotationAxes);
	
	// read image file
	char *bmpFullPath=archiver->ExpandInputPath(bmpFileName);
	ReadBMPFile(bmpFullPath,info,&rows);
	delete [] bmpFullPath;
	
	// angle file name (overrides other angle settings)
	if(bmpAngleFileName[0]>0)
	{	setAngles=TRUE;
		char *bmpFullAnglePath=archiver->ExpandInputPath(bmpAngleFileName);
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
	
	// resolutions
	if(bheight<0) bheight=bwidth*(double)info.height/(double)info.width;
	if(bwidth<0) bwidth=bheight*(double)info.width/(double)info.height;
	double xpw=bwidth/(double)info.width;
	double ypw=bheight/(double)info.height;
	
	// variables for scanning BMP file
	int matID;
	MPMBase *newMpt;
	int ii,k;
	Vector mpos[MaxElParticles];
	int ptFlag;
	BMPLevel *nextLevel;
	double deltax,deltay,deltaz,semiscale;
	double rmin,rmax,wtr1,wtr2,rweight;
	double cmin,cmax,wtc1,wtc2,weight;
	int r1,r2,c1,c2;
	BMPLevel *maxLevel;
	int row,col;
	ElementBase *elem;
	
	// Length/semiscale is half particle with
	//	(semiscale=4 for 2D w 4 pts per element or 3D with 8 pts per element,
	//		or 2 if 1 particle in the element)
	if(fmobj->IsThreeD())
		semiscale=2.*pow((double)fmobj->ptsPerElement,1./3.);
	else
		semiscale=2.*sqrt((double)fmobj->ptsPerElement);
	
	// scan mesh and assign material points or angles
    for(ii=1;ii<=nelems;ii++)
	{	// skip if image not in extend of element box
		elem=theElements[ii-1];
		if(!elem->IntersectsBox(xorig,yorig,bwidth,bheight,zslice))
			continue;
		
		// load point coordinates
		elem->MPMPoints(fmobj->ptsPerElement,mpos);
	
		// half the extent of the volume of a particle
		deltax=(elem->GetDeltaX())/semiscale;
		deltay=(elem->GetDeltaY())/semiscale;
		if(fmobj->IsThreeD())
			deltaz=elem->GetDeltaZ()/semiscale;
		else
			deltaz=1.;
		if(deltaz<0.) deltaz=1.;
		
        for(k=0;k<fmobj->ptsPerElement;k++)
		{	ptFlag=1<<k;
		
			// skip if alreadyu filled
            if(elem->filled&ptFlag) continue;
			
			// if point in the view area, then check it
			if((mpos[k].x>=xorig && mpos[k].x<xorig+bwidth)
					&& (mpos[k].y>=yorig && mpos[k].y<yorig+bheight)
						 && (mpos[k].z>=zslice-deltaz && mpos[k].z<zslice+deltaz))
			{	// find range of rows and cols for pixels over this material point
				rmin=(mpos[k].y-deltay-yorig)/ypw;
				rmax=(mpos[k].y+deltay-yorig)/ypw;
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
				cmin=(mpos[k].x-deltax-xorig)/xpw;
				cmax=(mpos[k].x+deltax-xorig)/xpw;
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
				for(row=r1;row<=r2;row++)
				{	if(row==r1)
						rweight=wtr1;
					else if(row==r2)
						rweight=wtr2;
					else
						rweight=1.;
					for(col=c1;col<=c2;col++)
					{	if(col==c1)
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

				// create a material point if one at this spot using matID and nextLevel
				// Note that empty spaced is not marked as filled which allows superposition of
				// images with different materials. If want to forcefully create a hole that
				// cannot be filled by subsequent image, will need to define a new material
				// type that can have matID for a hole. It will not create a point, but will mark
				// the location as filled
				if(matID>0)
				{	if(fmobj->IsThreeD())
						newMpt=new MatPoint3D(ii,matID,nextLevel->angle);
					else
						newMpt=new MatPoint2D(ii,matID,nextLevel->angle,nextLevel->thickness);
					newMpt->SetPosition(&mpos[k]);
					newMpt->SetOrigin(&mpos[k]);
					newMpt->SetVelocity(&nextLevel->vel);
					mpCtrl->AddMaterialPoint(newMpt,nextLevel->concentration,nextLevel->temperature);
					
					// is there an angle image?
					if(setAngles)
					{	// weight average of scanned pixels
						double totalWeight=0.;
						double totalIntensity=0.;
						for(row=r1;row<=r2;row++)
						{	if(row==r1)
								rweight=wtr1;
							else if(row==r2)
								rweight=wtr2;
							else
								rweight=1.;
							for(col=c1;col<=c2;col++)
							{	if(col==c1)
									weight=rweight*wtc1;
								else if(col==c2)
									weight=rweight*wtc2;
								else
									weight=rweight;
								
								totalIntensity+=weight*angleRows[row][col];
								totalWeight+=weight;
							}
						}
						totalIntensity/=totalWeight;
						double matAngle=minAngle+(totalIntensity-minIntensity)*angleScale;
						newMpt->SetAnglez0InDegrees(matAngle);
					}
					else
					{	// If had Rotate commands then use them
						SetMptAnglesFromFunctions(numRotations,&mpos[k],newMpt);
					}
					elem->filled|=ptFlag;
				}
            }
        }
    }
	
	// clean up
	for(row=0;row<info.height;row++)
	{	delete [] rows[row];
		if(setAngles) delete [] angleRows[row];
	}
	delete [] rows;
	if(setAngles) delete [] angleRows;
	
	// angles if allocated
	for(ii=0;ii<numRotations;ii++)
	{	delete [] angleExpr[ii];
		DeleteFunction(ii+1);
	}
		
}

// set current intensity velocity
void MPMReadHandler::SetLevelVelocity(double vx,double vy,double vz)
{
	currentLevel->vel.x=vx;
	currentLevel->vel.y=vy;
	currentLevel->vel.z=vz;
}

