/********************************************************************************
    LineController.cpp
    NairnFEA
    
    Created by John Nairn on 6/62/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/ElementBase.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Read_XML/LineController.hpp"

/********************************************************************************
	LineController: Constructors and Destructor
********************************************************************************/

LineController::LineController(int block) : ShapeController(block)
{
}

LineController::LineController(int block,double x1,double x2,double y1,double y2,double tolerate)
		: ShapeController(block,x1,x2,y1,y2,tolerance)
{
	// globals for length and tolerance
	distanceSq=(xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin);
	dtolerance=sqrt(distanceSq)*tolerate;
}

// called after initialization is done
void LineController::FinishSetup(void)
{
	// globals for length and tolerance
	distanceSq=(xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin);
	dtolerance=sqrt(distanceSq)*tolerance;
}

/********************************************************************************
	LineController: methods
********************************************************************************/

// Deterime if point is sufficiently close to the line from
//   (xmin,ymin) to (xmax,ymax). Must be within rectangle a
//   distance tolerance from the line in all directions 
bool LineController::PtOnShape(Vector v)
{
    double ypr=(xmax-xmin)*(v.y-ymin)-(ymax-ymin)*(v.x-xmin);
    
    if(ypr>dtolerance) return FALSE;
    if(ypr<-dtolerance) return FALSE;
    
    double xpr=(xmax-xmin)*(v.x-xmin)+(ymax-ymin)*(v.y-ymin);
    
    if(xpr<-dtolerance) return FALSE;
    if(xpr>distanceSq+dtolerance) return FALSE;
    
    return TRUE;
}

// set a property
void LineController::SetProperty(char *aName,char *value,CommonReadHandler *reader)
{
	if(strcmp(aName,"tolerance")==0)
	{	if(value[0]=='*')
		{	sscanf(&value[1],"%lf",&tolerance);
			tolerance*=ElementBase::GetMinimumCellSize();
		}
		else
		{	sscanf(value,"%lf",&tolerance);
			tolerance*=distScaling;
		}
	}
	else
		ShapeController::SetProperty(aName,value,reader);
}

// set tolerance (no need of scaling)
void LineController::SetTolerance(double value) { tolerance=value; }

