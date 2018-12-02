/********************************************************************************
    MatRegionFEA.cpp - extra code for FEAReadHandler.cpp to generate elements
        from shape commands mapped to a mesh of elements
    NairnFEA

    Created by John Nairn on 4/20/12.
    Copyright (c) 2012 RSAC Software. All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_FEA/FEAReadHandler.hpp"
#include "Read_XML/MaterialController.hpp"
#include "Read_XML/RectController.hpp"
#include "Read_XML/OvalController.hpp"
#include "Read_XML/PolygonController.hpp"
#include "Read_XML/LineController.hpp"
#include "Read_XML/ArcController.hpp"
#include "Elements/ElementBase.hpp"
#include "Read_XML/Expression.hpp"

static int regionMatNum;
static double regionThick;
static char *regionAngleExpr = NULL;

// Check for material region commands
// throws std::bad_alloc, SAXException()
short FEAReadHandler::MatRegionInput(char *xName,const Attributes& attrs)
{
    int i,numAttr;
    char *value,*aName;
    
    // Begin material shape regions
    if(strcmp(xName,"Body")==0 || strcmp(xName,"Hole")==0)
    {   ValidateCommand(xName,MUST_BE_NO_BLOCK,ANY_DIM);
        block=MATREGIONBLOCK;
        
		regionMatNum = strcmp(xName,"Hole")==0 ? HOLE_MATERIAL : NO_MATERIAL;
		regionThick = 1.0;           // optional
		char matname[200];
		matname[0]=0;
        
    	numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"mat")==0)
				sscanf(value,"%d",&regionMatNum);
			else if(strcmp(aName,"matname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(matname,value);
			}
            else if(strcmp(aName,"angle")==0)
			{	if(regionAngleExpr!=NULL) delete [] regionAngleExpr;
				regionAngleExpr=new char[strlen(value)+1];
				strcpy(regionAngleExpr,value);
			}
			else if(strcmp(aName,"thick")==0)
				sscanf(value,"%lf",&regionThick);
            delete [] aName;
            delete [] value;
        }
		
		// if gave a matname, it takes precedence over mat number
		if(strlen(matname)>0)
			regionMatNum = matCtrl->GetIDFromNewName(matname);
		if(regionMatNum<HOLE_MATERIAL)
			throw SAXException("Missing or invalid material number for Area.");
        
        // create function if being used
        if(regionAngleExpr!=NULL)
		{	if(!Expression::CreateFunction(regionAngleExpr,1))
                throw SAXException("The expression for material angle is not a valid function");
        }
    }
	
    // basic material region shapes
    // Rect and Oval done now, Polygon finished later
    else if(strcmp(xName,"Oval")==0 || strcmp(xName,"Rect")==0 || strcmp(xName,"Polygon")==0
			|| strcmp(xName,"Line")==0 || strcmp(xName,"Arc")==0)
	{	// only allowed in BODYPART or within a parent BODY_SHAPE
		if(block!=MATREGIONBLOCK && block!=BODY_SHAPE)
			ThrowCompoundErrorMessage(xName,"command found at invalid location","");
		
		// only allowed if BODY_SHAPE has active shape
		if(block==BODY_SHAPE && theShape==NULL)
			throw SAXException("Subordinate shape found without a parent shape.");
		
		// create shape
		ShapeController *newShape;
		if(strcmp(xName,"Rect")==0)
        {	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
            newShape = new RectController(MATREGIONBLOCK);
        }
        else if(strcmp(xName,"Oval")==0)
        {	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
            newShape = new OvalController(MATREGIONBLOCK);
        }
        else if(strcmp(xName,"Polygon")==0)
        {	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
            newShape = new PolygonController(MATREGIONBLOCK);
        }
        else if(strcmp(xName,"Line")==0)
        {	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
            newShape = new LineController(MATREGIONBLOCK,true);
        }
        else if(strcmp(xName,"Arc")==0)
        {	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
            newShape = new ArcController(MATREGIONBLOCK);
        }
		else
		{	// shape that contains no points
			newShape = new ShapeController(BODYPART);
		}
		
		// attributes
		newShape->SetScaling(ReadUnits(attrs,LENGTH_UNITS));
        numAttr = (int)attrs.getLength();
        for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
			newShape->SetProperty(aName,value,this);
            delete [] aName;
            delete [] value;
        }
		newShape->FinishSetup();
		
		// If in BODY_SHAPE then add this as child of current shape
		if(block==BODY_SHAPE)
			newShape->SetParentShape(theShape);
		block=BODY_SHAPE;
		theShape = newShape;
	}
    
	// add to polygon body object
    else if(strcmp(xName,"pt")==0)
	{	ValidateCommand(xName,BODY_SHAPE,MUST_BE_2D);
		if(theShape == NULL)
			throw SAXException("Body object <Ppt> command occurred without an active 2D body shape.");
		theShape->SetScaling(ReadUnits(attrs,LENGTH_UNITS));
		numAttr = (int)attrs.getLength();
		for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
			theShape->SetParameter(aName,value);
			delete [] aName;
			delete [] value;
		}
		theShape->FinishParameter();
	}

	// add to arec to oval body object
	else if(strcmp(xName,"arc")==0)
	{	ValidateCommand(xName,BODY_SHAPE,MUST_BE_2D);
		if(theShape == NULL)
			throw SAXException("Body object <arc> command occurred without an active 2D body shape.");
		numAttr=(int)attrs.getLength();
		for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
			theShape->SetParameter(aName,value);
			delete [] aName;
			delete [] value;
		}
		if(!theShape->FinishParameter())
			throw SAXException("<arc> command has invalid start and/or end angle.");
	}

	// no a material region command
    else
        return false;
    
    // was handled
    return true;

}

// Subroutine to finish up mat region commands
// throws SAXException()
short FEAReadHandler::EndMatRegionInput(char *xName,int exitBlock)
{
    // finish the matregion block
    if(strcmp(xName,"Body")==0 || strcmp(xName,"Hole")==0)
    {   block=exitBlock;
        if(regionAngleExpr!=NULL)
        {   delete [] regionAngleExpr;
            regionAngleExpr = NULL;
        }
    }
    
    // mesh shape now
    else if(strcmp(xName,"Oval")==0 || strcmp(xName,"Rect")==0 || strcmp(xName,"Polygon")==0
				|| strcmp(xName,"Line")==0 || strcmp(xName,"Arc")==0)
	{	if(strcmp(xName,"Polygon")==0)
		{	if(!theShape->HasAllParameters())
				throw SAXException("<Polygon> must have at least 3 subordinate <Ppt> commands.");
		}
		
		// If parent shape then done
		ShapeController *parentShape = theShape->GetParentShape();
		if(parentShape==NULL)
		{	SetRegionElements();
			delete theShape;
			theShape = NULL;
			block=MATREGIONBLOCK;
		}
		
		// if subordinate shape, add to parent shape
		else
		{	parentShape->AddCutoutShape(theShape);
			// return to parent
			theShape = parentShape;
			// stay in BODY_SHAPE block
		}
	}
    
	// not a MatRegion block element
	else
		return FALSE;
    
	return TRUE;
}

// Subroutine to finish up mat region commands
void FEAReadHandler::SetRegionElements(void)
{   
    int i;
    double angle=0;
    
    theShape->resetElementEnumerator();
    ElementBase *elem;
    while((i=theShape->nextElement())>=0)
    {   // get material ID
        elem = theElements[i];
        int matID = elem->material;

        // skip if not empty
        if(matID!=NO_MATERIAL) continue;
        
        // calculate angle
        if(regionAngleExpr!=NULL)
        {   Vector midPt;
            elem->GetXYZCentroid(&midPt);
			angle=Expression::FunctionValue(1,midPt.x,midPt.y,0.,0.,0.,0.);
        }
        
        // assign material, thickness, an angle
        elem->material = regionMatNum;
        elem->SetThickness(regionThick);
        elem->SetAngleInDegrees(angle);
    }
    
    
}
