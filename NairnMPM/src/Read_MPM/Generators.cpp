/********************************************************************************
    Generators.cpp - extra code for MPMReadHandler.cpp to generate mesh
        from xml commands (instead of from point by point commands)
    NairnMPM

    Created by John Nairn on Thu Apr 11 2002.
    Copyright (c) 2002 RSAC Software. All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Read_MPM/MPMReadHandler.hpp"
#include "Elements/ElementBase.hpp"
#include "MPM_Classes/MatPoint2D.hpp"
#include "MPM_Classes/MatPointAS.hpp"
#include "MPM_Classes/MatPoint3D.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Exceptions/StrX.hpp"
#include "Read_XML/ArcController.hpp"
#include "Read_XML/RectController.hpp"
#include "Read_XML/OvalController.hpp"
#include "Read_XML/BoxController.hpp"
#include "Read_MPM/MpsController.hpp"
#include "Elements/FourNodeIsoparam.hpp"
#include "Elements/EightNodeIsoparamBrick.hpp"
#include "Nodes/NodalPoint2D.hpp"
#include "Nodes/NodalPoint3D.hpp"
#include "Read_MPM/CrackController.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Read_MPM/SphereController.hpp"
#include "Read_XML/PolygonController.hpp"
#include "Read_MPM/PolyhedronController.hpp"
#include "Read_XML/mathexpr.hpp"
#include "Read_XML/ElementsController.hpp"
#include "Read_XML/MaterialController.hpp"


// Global variables for Generator.cpp (first letter all capitalized)
double Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,Rhoriz=1.,Rvert=1.,Rdepth=1.,Z2DThickness;
double pConc,pTemp,Angle,Thick;
int Nhoriz=0,Nvert=0,Ndepth=0,MatID;
double cellHoriz=-1.,cellVert=-1.,cellDepth=-1.;
Vector Vel;
char *angleExpr[3];
char rotationAxes[4];

// to allow new elements, add attribute or command to change this for grid
// Create those element in ElementsController::MeshElement()
int mpm2DElement = FOUR_NODE_ISO;
int mpm3DElement = EIGHT_NODE_ISO_BRICK;

enum { LINE_SHAPE=0,CIRCLE_SHAPE };		// shape for generating crack segments, Liping Xue

// Function prototypes
void swap(double&,double&,double&,double&,double&,double&);

short MPMReadHandler::GenerateInput(char *xName,const Attributes& attrs)
{
    char *aName,*value;
    int i,numAttr;
	double aScaling;

    //-----------------------------------------------------------
    // Read into grid parameters (Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)
    //-----------------------------------------------------------
    if(strcmp(xName,"Grid")==0)
	{	ValidateCommand(xName,MESHBLOCK,ANY_DIM);
		if(meshType!=UNKNOWN_MESH)
            throw SAXException("<Grid> cannot be mixed with an explicit mesh.");
		aScaling=ReadUnits(attrs,LENGTH_UNITS);
        block=GRIDBLOCK;
        numAttr=attrs.getLength();
		Xmin=Xmax=Ymin=Ymax=Zmin=Zmax=0.;
		Z2DThickness=1.0;
        for(i=0;i<numAttr;i++) {
            aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"xmin")==0)
                sscanf(value,"%lf",&Xmin);
            else if(strcmp(aName,"xmax")==0)
                sscanf(value,"%lf",&Xmax);
            else if(strcmp(aName,"ymin")==0)
                sscanf(value,"%lf",&Ymin);
            else if(strcmp(aName,"ymax")==0)
                sscanf(value,"%lf",&Ymax);
            else if(strcmp(aName,"zmin")==0)
                sscanf(value,"%lf",&Zmin);
            else if(strcmp(aName,"zmax")==0)
                sscanf(value,"%lf",&Zmax);
            else if(strcmp(aName,"thickness")==0)
                sscanf(value,"%lf",&Z2DThickness);
			else if(strcmp(aName,"type")==0)
			{	int mpmElement;
				sscanf(value,"%d",&mpmElement);
				if(fmobj->IsThreeD())
				{	// only one type is allowed in 3D
					if(mpmElement!=EIGHT_NODE_ISO_BRICK)
						throw SAXException("<Grid type='#'> attribute is invalid number for 3D MPM calculations.");
				}
				else
				{	// only one type is allowed, but can switch to allow two once that element does gimp
					//if(mpmElement!=FOUR_NODE_ISO && mpmElement!=NINE_NODE_LAGRANGE)
					if(mpmElement!=FOUR_NODE_ISO)
						throw SAXException("<Grid type='#'> attribute is invalid number for 2D MPM calculations.");
					mpm2DElement=mpmElement;
				}
			}
            delete [] aName;
            delete [] value;
        }
		Xmin*=aScaling;
		Xmax*=aScaling;
		Ymin*=aScaling;
		Ymax*=aScaling;
		Zmin*=aScaling;
		Zmax*=aScaling;
        swap(Xmin,Xmax,Ymin,Ymax,Zmin,Zmax);
    }
    
    //-----------------------------------------------------------
    // Read into cell parameters (nx,rx,ny,ry,nz,rz)
    //-----------------------------------------------------------
    else if(strcmp(xName,"Horiz")==0)
	{	ValidateCommand(xName,GRIDBLOCK,ANY_DIM);
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++) {
            aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"nx")==0)
                sscanf(value,"%d",&Nhoriz);
            else if(strcmp(aName,"rx")==0)
                sscanf(value,"%lf",&Rhoriz);
            else if(strcmp(aName,"cellsize")==0)
                sscanf(value,"%lf",&cellHoriz);
            delete [] aName;
            delete [] value;
        }
    }

    else if(strcmp(xName,"Vert")==0)
	{	ValidateCommand(xName,GRIDBLOCK,ANY_DIM);
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++) {
            aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"ny")==0)
                sscanf(value,"%d",&Nvert);
            else if(strcmp(aName,"ry")==0)
                sscanf(value,"%lf",&Rvert);
 			else if(strcmp(aName,"cellsize")==0)
                sscanf(value,"%lf",&cellVert);
           	delete [] aName;
            delete [] value;
        }
    }

    else if(strcmp(xName,"Depth")==0)
	{	ValidateCommand(xName,GRIDBLOCK,MUST_BE_3D);
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++) {
            aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"nz")==0)
                sscanf(value,"%d",&Ndepth);
            else if(strcmp(aName,"rz")==0)
                sscanf(value,"%lf",&Rdepth);
            else if(strcmp(aName,"cellsize")==0)
                sscanf(value,"%lf",&cellDepth);
           	delete [] aName;
            delete [] value;
        }
    }

    //--------------------------------------------------------------
    // Read into body properties (materID, angle, thick, velocities)
    //--------------------------------------------------------------
    else if(strcmp(xName,"Body")==0 || strcmp(xName,"Hole")==0)
	{	ValidateCommand(xName,POINTSBLOCK,ANY_DIM);
        block=BODYPART;
        MatID=0;
		char matname[200];
		matname[0]=0;
        if(strcmp(xName,"Body")==0)
        {   Thick=mpmgrid.GetDefaultThickness();
            Angle=Vel.x=Vel.y=Vel.z=0.0;
			pConc=0.0;
			pTemp=thermal.reference;
            numAttr=attrs.getLength();
			rotationAxes[0]=0;			// no rotations yet
            for(i=0;i<numAttr;i++)
			{	aName=XMLString::transcode(attrs.getLocalName(i));
                value=XMLString::transcode(attrs.getValue(i));
                if(strcmp(aName,"mat")==0)
                    sscanf(value,"%d",&MatID);
                else if(strcmp(aName,"matname")==0)
				{	if(strlen(value)>199) value[200]=0;
                    strcpy(matname,value);
				}
                else if(strcmp(aName,"angle")==0)
                    sscanf(value,"%lf",&Angle);
                else if(strcmp(aName,"thick")==0)
                    sscanf(value,"%lf",&Thick);
                else if(strcmp(aName,"vx")==0)
                    sscanf(value,"%lf",&Vel.x);
                else if(strcmp(aName,"vy")==0)
                    sscanf(value,"%lf",&Vel.y);
                else if(strcmp(aName,"vz")==0)
                    sscanf(value,"%lf",&Vel.z);
                else if(strcmp(aName,"conc")==0)
				{	sscanf(value,"%lf",&pConc);
					if(pConc<0. || pConc>1.)
						throw SAXException("Material point concentration potential must be from 0 to 1");
				}
				else if(strcmp(aName,"wtconc")==0)
				{	sscanf(value,"%lf",&pConc);
					pConc=-pConc;
					if(pConc>0.)
						throw SAXException("Material point weight fraction concentration must be >= 0");
				}
                else if(strcmp(aName,"temp")==0)
                    sscanf(value,"%lf",&pTemp);
                delete [] aName;
                delete [] value;
            }
			// if gave a matname, it takes precedence over mat number
			if(strlen(matname)>0)
				MatID = matCtrl->GetIDFromNewName(matname);
            if(MatID<=0)
                throw SAXException("Positive material ID needed in <Body> element.");
            if(Thick<=0. && !fmobj->IsThreeD() && !fmobj->IsAxisymmetric())
                throw SAXException("Positive material thickness needed in <Body> element for planar 2D analyses.");
        }
    }

    //-----------------------------------------------------------
    // Read into geometry parameters for Body shape objects
    //-----------------------------------------------------------
	else if(strcmp(xName,"RotateZ")==0)
	{	if(block!=BODYPART && block!=BMPBLOCK)
            ValidateCommand(xName,BAD_BLOCK,ANY_DIM);
		int rotationNum=strlen(rotationAxes);
		if(rotationNum>=3)
			throw SAXException("Maximum of three rotations allowed within a single <Body>.");
		strcat(rotationAxes,"Z");
		input=TEXT_BLOCK;
        inputID=CHAR_ARRAY;
		inputPtr=NULL;
	}
	
	else if(strcmp(xName,"RotateY")==0 || strcmp(xName,"RotateX")==0)
	{	ValidateCommand(xName,NO_BLOCK,MUST_BE_3D);
		if(block!=BODYPART && block!=BMPBLOCK)
            ValidateCommand(xName,BAD_BLOCK,ANY_DIM);
		int rotationNum=strlen(rotationAxes);
		if(rotationNum>=3)
			throw SAXException("Maximum of three rotations allowed within a single <Body>.");
		rotationAxes[rotationNum]=xName[6];
		rotationAxes[rotationNum+1]=0;
		input=TEXT_BLOCK;
        inputID=CHAR_ARRAY;
		inputPtr=NULL;
	}
	
	else if(strcmp(xName,"Unrotate")==0)
	{	ValidateCommand(xName,BODYPART,ANY_DIM);
		int numRotations=strlen(rotationAxes);
		for(i=0;i<numRotations;i++) delete [] angleExpr[i];
		rotationAxes[0]=0;
	}
		
    //-----------------------------------------------------------
    // Read into geometry parameters for Body shape objects
    //-----------------------------------------------------------
    else if(strcmp(xName,"Oval")==0 || strcmp(xName,"Rect")==0 || strcmp(xName,"Polygon")==0
				|| strcmp(xName,"Sphere")==0 || strcmp(xName,"Box")==0 || strcmp(xName,"Cylinder")==0
				|| strcmp(xName,"Polyhedron")==0 )
	{	if(strcmp(xName,"Rect")==0)
		{	ValidateCommand(xName,BODYPART,MUST_BE_2D);
            theShape = new RectController(BODYPART);
 		}
		else if(strcmp(xName,"Oval")==0)
		{	ValidateCommand(xName,BODYPART,MUST_BE_2D);
            theShape = new OvalController(BODYPART);
		}
		else if(strcmp(xName,"Box")==0 || strcmp(xName,"Cylinder")==0)
		{	ValidateCommand(xName,BODYPART,MUST_BE_3D);
			theShape = new BoxController(BODYPART);
		}
 		else if(strcmp(xName,"Sphere")==0)
		{	ValidateCommand(xName,BODYPART,MUST_BE_3D);
			theShape = new SphereController(BODYPART);
		}
		else if(strcmp(xName,"Polygon")==0)
		{	ValidateCommand(xName,BODYPART,MUST_BE_2D);
			theShape = new PolygonController(BODYPART);
		}
		else if(strcmp(xName,"Polyhedron")==0)
		{	ValidateCommand(xName,BODYPART,MUST_BE_3D);
			theShape = new PolyhedronController(BODYPART);
		}
		theShape->SetScaling(ReadUnits(attrs,LENGTH_UNITS));
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
			theShape->SetProperty(aName,value,this);
            delete [] aName;
            delete [] value;
        }
		
		// finish up and if body is done, generate points now
		if(theShape->FinishSetup())
		{	MPMPts();
			delete theShape;
			theShape=NULL;
		}
		else
			block=BODY_SHAPE;
	}

    //-----------------------------------------------------------
    // Read into geometry parameters for crack segments
    //       Line (xmin,ymin,xmax,ymax)
	//       Circle (xmin,ymin,xmax,ymax)
	// Liping Xue
    //-----------------------------------------------------------
    else if(strcmp(xName,"Line")==0 | strcmp(xName,"Circle")==0 )
	{	ValidateCommand(xName,CRACKLIST,MUST_BE_2D);
		int crackShape=LINE_SHAPE;
		if(strcmp(xName,"Circle")==0)
			crackShape=CIRCLE_SHAPE;
        numAttr=attrs.getLength();
		aScaling=ReadUnits(attrs,LENGTH_UNITS);
		Xmin=Xmax=Ymin=Ymax=0.;
		int resolution=1;
		double start_angle=0.,end_angle=360.;
		int startTip=-1,endTip=-1;
		MatID=0;
		char lawname[200];
		lawname[0]=0;
        for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"xmin")==0)
                sscanf(value,"%lf",&Xmin);
            else if(strcmp(aName,"xmax")==0)
                sscanf(value,"%lf",&Xmax);
            else if(strcmp(aName,"ymin")==0)
                sscanf(value,"%lf",&Ymin);
            else if(strcmp(aName,"ymax")==0)
                sscanf(value,"%lf",&Ymax);
            else if(strcmp(aName,"resolution")==0)
			{	sscanf(value,"%d",&resolution);
				if(resolution<1) resolution=1;
			}
            else if(strcmp(aName,"start_angle")==0)
                sscanf(value,"%lf",&start_angle);
            else if(strcmp(aName,"end_angle")==0)
                sscanf(value,"%lf",&end_angle);
            else if(strcmp(aName,"start_tip")==0)
			{	sscanf(value,"%d",&startTip);
				if(startTip<0 && startTip!=-2) startTip=-1;
			}
            else if(strcmp(aName,"end_tip")==0)
			{	sscanf(value,"%d",&endTip);
				if(endTip<0 && endTip!=-2) endTip=-1;
			}
			else if(strcmp(aName,"mat")==0)
				sscanf(value,"%d",&MatID);
			else if(strcmp(aName,"matname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(lawname,value);
			}
			
            delete [] aName;
            delete [] value;
        }
		// if gave a lawname, it takes precedence over mat number
		if(strlen(lawname)>0)
			MatID = matCtrl->GetIDFromNewName(lawname);
		Xmin*=aScaling;
		Xmax*=aScaling;
		Ymin*=aScaling;
		Ymax*=aScaling;
		// create crack segment
		MPMCracks(crackShape,resolution,start_angle,end_angle,startTip,endTip);
    }

	// add to polygon body object
    else if(strcmp(xName,"pt")==0)
	{	ValidateCommand(xName,BODY_SHAPE,MUST_BE_2D);
		if(theShape == NULL)
			throw SAXException("Body object <pt> command occurred without an active 2D body shape.");
		theShape->SetScaling(ReadUnits(attrs,LENGTH_UNITS));
		numAttr=attrs.getLength();
		for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
			theShape->SetParameter(aName,value);
			delete [] aName;
			delete [] value;
		}
		theShape->FinishParameter();
	}
	
	// triclinic option
	else if(strcmp(xName,"faces")==0)
	{	ValidateCommand(xName,BODY_SHAPE,MUST_BE_3D);
		if(theShape == NULL)
			throw SAXException("Body object <faces> command occurred without an active 3D body shape.");
		theShape->SetScaling(ReadUnits(attrs,LENGTH_UNITS));
		theShape->SetParameter("style","");
		numAttr=attrs.getLength();
		for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
			theShape->SetParameter(aName,value);
			delete [] aName;
			delete [] value;
		}
		if(!theShape->FinishParameter())
			throw SAXException("Body object <faces> command did not specify the style for the data.");
		input=POLYHEDRON_BLOCK;
	}
	
    // Store a line and tolerance into the current BC shape
    else if(strcmp(xName,"BCLine")==0)
	{	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
    	if(block!=GRIDBCHEADER && block!=PARTICLEBCHEADER)
            ValidateCommand(xName,BAD_BLOCK,ANY_DIM);
		LineController *theLine=new LineController(block);
        numAttr=attrs.getLength();
		theLine->SetScaling(ReadUnits(attrs,LENGTH_UNITS));
		theLine->SetTolerance(ElementBase::gridTolerance);
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
			theLine->SetProperty(aName,value,this);
            delete [] aName;
            delete [] value;
        }
		theLine->FinishSetup();
		theShape=theLine;
        block=BCSHAPE;
    }

    // Store an arc and tolerance into the current BC shape
    else if(strcmp(xName,"BCArc")==0)
	{	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
    	if(block!=GRIDBCHEADER && block!=PARTICLEBCHEADER)
            ValidateCommand(xName,BAD_BLOCK,ANY_DIM);
		ArcController *theArc=new ArcController(block);
        numAttr=attrs.getLength();
		theArc->SetScaling(ReadUnits(attrs,LENGTH_UNITS));
		theArc->SetTolerance(ElementBase::gridTolerance);
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
			theArc->SetProperty(aName,value,this);
            delete [] aName;
            delete [] value;
        }
		theArc->FinishSetup();
		theShape=theArc;
        block=BCSHAPE;
    }

    // Store a box into the current BC shape
    else if(strcmp(xName,"BCBox")==0)
	{	ValidateCommand(xName,NO_BLOCK,MUST_BE_3D);
    	if(block!=GRIDBCHEADER && block!=PARTICLEBCHEADER)
            ValidateCommand(xName,BAD_BLOCK,ANY_DIM);
		BoxController *theBox=new BoxController(block);
        numAttr=attrs.getLength();
		theBox->SetScaling(ReadUnits(attrs,LENGTH_UNITS));
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
			theBox->SetProperty(aName,value,this);
            delete [] aName;
            delete [] value;
        }
        theBox->FinishSetup();
		theShape=theBox;
        block=BCSHAPE;
    }
	
    // Read into velocity boundary conditions on current line (direction, value)
    else if(strcmp(xName,"DisBC")==0)
	{	ValidateCommand(xName,BCSHAPE,ANY_DIM);
		if(!theShape->RequiredBlock(GRIDBCHEADER))
			ValidateCommand(xName,BAD_BLOCK,ANY_DIM);
        double dispvel=0.0,ftime=0.0;
		int dof=0,style=CONSTANT_VALUE;
		char *function=NULL;
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"dir")==0)
                sscanf(value,"%d",&dof);
            else if(strcmp(aName,"disp")==0 || strcmp(aName,"vel")==0)
                sscanf(value,"%lf",&dispvel);
            else if(strcmp(aName,"time")==0 || strcmp(aName,"ftime")==0)
                sscanf(value,"%lf",&ftime);
            else if(strcmp(aName,"style")==0)
                sscanf(value,"%d",&style);
            else if(strcmp(aName,"function")==0)
			{	if(function!=NULL) delete [] function;
				function=new char[strlen(value)+1];
				strcpy(function,value);
			}
            delete [] aName;
            delete [] value;
        }
		if(fmobj->IsThreeD())
		{	if(dof<1 || dof>3)
				throw SAXException("'dir' in DisBC element must be 1, 2, or 3 for 3D analyses.");
		}
        else if(dof>2 || dof<0)
            throw SAXException("'dir' in DisBC element must be 0, 1, or 2 for 2D analyses.");
        
        // check all nodes
		theShape->resetNodeEnumerator();
		while((i=theShape->nextNode()))
		{	NodalVelBC *newVelBC=new NodalVelBC(nd[i]->num,dof,style,dispvel,ftime);
			newVelBC->SetFunction(function);
			velocityBCs->AddObject(newVelBC);
		}
		if(function!=NULL) delete [] function;
    }

	// Allow to specify total BC on particles
	else if(strcmp(xName,"net")==0)
	{	ValidateCommand(xName,BCSHAPE,ANY_DIM);
		theShape->setNetBC(true);
	}

	// Allow to specify BC on per particle basis
	else if(strcmp(xName,"perParticle")==0)
	{	ValidateCommand(xName,BCSHAPE,ANY_DIM);
		theShape->setNetBC(false);
	}
	
    // Read into concentration or temperature boundary conditions
	//  on current shape
    else if(strcmp(xName,"ConcBC")==0 || strcmp(xName,"TempBC")==0)
	{	ValidateCommand(xName,BCSHAPE,ANY_DIM);
		if(!theShape->RequiredBlock(GRIDBCHEADER))
			ValidateCommand(xName,BAD_BLOCK,ANY_DIM);
		int tempBC=strcmp(xName,"TempBC")==0 ? TRUE : FALSE;
        double bcvalue=0.0,ftime=0.0;
        int style=CONSTANT_VALUE;
		char *function=NULL;
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"value")==0)
                sscanf(value,"%lf",&bcvalue);
            else if(strcmp(aName,"time")==0)
                sscanf(value,"%lf",&ftime);
            else if(strcmp(aName,"style")==0)
                sscanf(value,"%d",&style);
            else if(strcmp(aName,"function")==0)
			{	if(function!=NULL) delete [] function;
				function=new char[strlen(value)+1];
				strcpy(function,value);
			}
            delete [] aName;
            delete [] value;
        }
        
        // check all nodes
		if(tempBC)
		{	theShape->resetNodeEnumerator();
			while((i=theShape->nextNode()))
			{	NodalTempBC *newTempBC=new NodalTempBC(nd[i]->num,style,bcvalue,ftime);
				newTempBC->SetFunction(function);
				tempBCs->AddObject(newTempBC);
			}
		}
		else
		{	theShape->resetNodeEnumerator();
			while((i=theShape->nextNode()))
			{	NodalConcBC *newConcBC=new NodalConcBC(nd[i]->num,style,bcvalue,ftime);
				newConcBC->SetFunction(function);
				concBCs->AddObject(newConcBC);
			}
		}
		if(function!=NULL) delete [] function;
    }

    //-----------------------------------------------------------
    // Read a rectangular shape for boundary conditions
    //-----------------------------------------------------------
    else if(strcmp(xName,"LdRect")==0)
	{	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
    	if(block!=GRIDBCHEADER && block!=PARTICLEBCHEADER)
            ValidateCommand(xName,BAD_BLOCK,ANY_DIM);
		RectController *theRect=new RectController(block);
		theRect->SetScaling(ReadUnits(attrs,LENGTH_UNITS));
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
			theRect->SetProperty(aName,value,this);
            delete [] aName;
            delete [] value;
        }
		theRect->FinishSetup();
		theShape=theRect;
        block=BCSHAPE;
    }

    // Read into load/traction boundary conditions (direction,style,values)
	//  or concentration flux BC on current shape
    else if(strcmp(xName,"LoadBC")==0 || strcmp(xName,"ConcFluxBC")==0 || strcmp(xName,"TractionBC")==0)
	{	if(strcmp(xName,"LoadBC")==0 || strcmp(xName,"TractionBC")==0)
			ValidateCommand(xName,BCSHAPE,ANY_DIM);
		else
			ValidateCommand(xName,BCSHAPE,ANY_DIM);
		if(!theShape->RequiredBlock(PARTICLEBCHEADER))
			ValidateCommand(xName,BAD_BLOCK,ANY_DIM);
        int dof=0,style=1,face=1;
        double ftime=0.0,load=0.0;
		char *function=NULL;
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"dir")==0)
                sscanf(value,"%d",&dof);
            else if(strcmp(aName,"style")==0)
                sscanf(value,"%d",&style);
            else if(strcmp(aName,"face")==0)
                sscanf(value,"%d",&face);
            else if(strcmp(aName,"load")==0 || strcmp(aName,"value")==0 || strcmp(aName,"stress")==0)
                sscanf(value,"%lf",&load);
            else if(strcmp(aName,"ftime")==0 || strcmp(aName,"time")==0 || strcmp(aName,"bath")==0)
                sscanf(value,"%lf",&ftime);
            else if(strcmp(aName,"function")==0)
			{	if(function!=NULL) delete [] function;
				function=new char[strlen(value)+1];
				strcpy(function,value);
			}
            delete [] aName;
            delete [] value;
        }
		if(fmobj->IsThreeD())
        {   // checks for 3D
			if(dof<1 || dof>3)
			{	throw SAXException("'dir' in <LoadBC>, <TractionBC>, or <ConcFluxBC> must be 1-3 for 3D analyses.");
			}
            if(strcmp(xName,"TractionBC")==0)
            {	if(face<1 || face>6)
                    throw SAXException("'face' in <TractionBC> element must be 1 to 6 for 3D analyses.");
            }
		}
        else
        {   // checks for 2D
            if(dof>2 || dof<1)
			{	if(dof<11 || dof>12 || strcmp(xName,"TractionBC")!=0)
				{	throw SAXException("'dir' in <LoadBC> or <ConcFluxBC> must be 1 or 2 and in <TractionBC> must be 1, 2, 11, or 12 for 2D analyses.");
				}
			}
            if(strcmp(xName,"TractionBC")==0)
            {	if(face<1 || face>4)
                    throw SAXException("'face' in <TractionBC> element must be 1 to 4 for 3D analyses.");
            }
        }

		
		// separate search depending on type
		if(strcmp(xName,"LoadBC")==0)
		{   // check each material point
			MatPtLoadBC *newLoadBC;
			theShape->resetParticleEnumerator();
			double normalize=theShape->particleCount();
			while((i=theShape->nextParticle())>=0)
			{   newLoadBC=new MatPtLoadBC(i+1,dof,style);
				newLoadBC->ftime=ftime;
				newLoadBC->SetFunction(function);
				if(function==NULL)
					newLoadBC->value=load/normalize;
				else
					newLoadBC->value=1./normalize;
				mpLoadCtrl->AddObject(newLoadBC);
			}
		}
		else if(strcmp(xName,"TractionBC")==0)
        {   if(style==SILENT)
                throw SAXException("Tractions boundary conditions cannot use silent style.");
            
		    // check each material point
			MatPtTractionBC *newTractionBC;
			theShape->resetParticleEnumerator();
			while((i=theShape->nextParticle())>=0)
			{   newTractionBC=new MatPtTractionBC(i+1,dof,face,style);
				newTractionBC->ftime=ftime;
				newTractionBC->SetFunction(function);
				newTractionBC->value=load;
				mpTractionCtrl->AddObject(newTractionBC);
			}
		}
		else if(strcmp(xName,"ConcFluxBC")==0)
        {	if(dof==2 && style!=FUNCTION_VALUE)
				throw SAXException("Coupled concentration flux condition must use function style");
			else if(style==SILENT && fmobj->IsThreeD())
				throw SAXException("Silent <ConcFluxBC> is currently only available in 2D analyses.");
				
			// check each material point
			MatPtFluxBC *newFluxBC;
			theShape->resetParticleEnumerator();
			while((i=theShape->nextParticle())>=0)
			{   newFluxBC=new MatPtFluxBC(i+1,dof,style,face);
				newFluxBC->value=load;
				newFluxBC->ftime=ftime;
				newFluxBC->SetFunction(function);
				mpConcFluxCtrl->AddObject(newFluxBC);
			}
		}
		if(function!=NULL) delete [] function;
    }
	
	// Not recognized as a generator element
	else
		return FALSE;
    
    return TRUE;
    
} //End of Generators.cpp

//-----------------------------------------------------------
// Subroutine to set block on element end
//-----------------------------------------------------------
short MPMReadHandler::EndGenerator(char *xName)
{
    if(strcmp(xName,"Grid")==0)
    {	grid();                     // Generate grid (1 per analysis)
		SetGIMPBorderAsHoles();     // if GIMP, mark edge elements as filled (i.e., no material pts allowed)
        block=MESHBLOCK;            // Must have been in MESHBLOCK
    }

    else if(strcmp(xName,"Body")==0)
	{	int numRotations=strlen(rotationAxes);
		int i;
		for(i=0;i<numRotations;i++)
		{	delete [] angleExpr[i];
			DeleteFunction(i+1);
		}
    	block=POINTSBLOCK;
	}
    
	else if(strcmp(xName,"RotateZ")==0 || strcmp(xName,"RotateY")==0 || strcmp(xName,"RotateX")==0)
	{	if(inputPtr==NULL)
			throw SAXException("No rotation angle or expression was provided in <RotateX(YZ)> command.");
		int rotationNum=strlen(rotationAxes);
		angleExpr[rotationNum-1]=inputPtr;
	}
	
    else if(strcmp(xName,"Hole")==0)
    	block=POINTSBLOCK;
    
    else if(strcmp(xName,"Polygon")==0)
	{	if(!theShape->HasAllParameters())
			throw SAXException("<Polygon> must have at least 3 subordinate <pt> commands.");
		MPMPts();
		delete theShape;
		theShape = NULL;
		block=BODYPART;
	}
    
    else if(strcmp(xName,"Polyhedron")==0)
	{	if(!theShape->HasAllParameters())
            throw SAXException("<Polyhedron> must have at least 4 faces.");
		MPMPts();
		delete theShape;
		theShape = NULL;
		block=BODYPART;
	}

    else if(strcmp(xName,"BCLine")==0 || strcmp(xName,"LdRect")==0 || strcmp(xName,"BCLine")==0 || strcmp(xName,"BCBox")==0)
	{	block=theShape->GetSourceBlock();
		delete theShape;
		theShape=NULL;
    }
    
	// Not recognized as a generator element
	else
		return FALSE;
		
	return TRUE;
}

//-----------------------------------------------------------
// Subroutine for creating mpm object (mpm)
//-----------------------------------------------------------
void MPMReadHandler::MPMPts(void)
{
	Vector ppos[MaxElParticles];
    MPMBase *newMpt;
    int i,k,ptFlag;
	int numRotations=strlen(rotationAxes);
	
	if(numRotations>0)
	{	Angle=0.;			// rotations override angle settings
		for(i=0;i<numRotations;i++)
		{	char *expr=new char[strlen(angleExpr[i])+1];
			strcpy(expr,angleExpr[i]);
			if(!CreateFunction(expr,i+1))
				throw SAXException("Invalid angle expression was provided in <RotateX(YZ)> command.");
			delete [] expr;
		}
	}

    for(i=1;i<=nelems;i++)
	{	theElements[i-1]->MPMPoints(fmobj->ptsPerElement,ppos);
        for(k=0;k<fmobj->ptsPerElement;k++)
		{	ptFlag=1<<k;
            if(theElements[i-1]->filled&ptFlag) continue;
            if(theShape->ContainsPoint(ppos[k]))
			{	if(MatID>0)
				{	if(fmobj->IsThreeD())
						newMpt=new MatPoint3D(i,MatID,Angle);
					else if(fmobj->IsAxisymmetric())
						newMpt=new MatPointAS(i,MatID,Angle,ppos[k].x);
					else
						newMpt=new MatPoint2D(i,MatID,Angle,Thick);
                    newMpt->SetPosition(&ppos[k]);
                    newMpt->SetOrigin(&ppos[k]);
					newMpt->SetVelocity(&Vel);
					SetMptAnglesFromFunctions(numRotations,&ppos[k],newMpt);
					mpCtrl->AddMaterialPoint(newMpt,pConc,pTemp);
                }
                theElements[i-1]->filled|=ptFlag;
            }
        }
    }
}

//------------------------------------------------------------------
// If just created a GIMP grid, make all border elements as holes
//------------------------------------------------------------------
void MPMReadHandler::SetGIMPBorderAsHoles(void)
{
	// not needed (or wanted) if not using GIMP
	if(ElementBase::useGimp == POINT_GIMP) return;
	
    int i,k,ptFlag;
    for(i=1;i<=nelems;i++)
	{	if(!theElements[i-1]->OnTheEdge()) continue;
        for(k=0;k<fmobj->ptsPerElement;k++)
		{	ptFlag=1<<k;
			theElements[i-1]->filled|=ptFlag;
        }
    }
}

//-----------------------------------------------------------
// Subroutine for creating crack segments
//-----------------------------------------------------------
void MPMReadHandler::MPMCracks(int crackShape,int resolution,double start_angle,double end_angle,int startTip,int endTip)
{
	int tipMatnum;
	double x, y, dx, dy, a, b, x0, y0, startSita, deltaSita;
	
	if(crackShape==LINE_SHAPE){
		dx=(Xmax-Xmin)/(double)resolution;
		dy=(Ymax-Ymin)/(double)resolution;
		tipMatnum=startTip;
		for(int i=0;i<=resolution;i++){
			x=Xmin+dx*i;
			y=Ymin+dy*i;
			if(!crackCtrl->AddSegment(new CrackSegment(x,y,tipMatnum,MatID)))
				throw SAXException("Crack not in the mesh or out of memory adding a crack segment.");
			tipMatnum = i!=resolution-1 ? -1 : endTip;
		}
	}
	else if(crackShape==CIRCLE_SHAPE){
		x0=(Xmin+Xmax)/2.0;
		y0=(Ymin+Ymax)/2.0;
		a=Xmax-x0;
		b=Ymax-y0;
		deltaSita=(end_angle-start_angle)/180.*PI_CONSTANT/(double)resolution;
		startSita=start_angle/180.*PI_CONSTANT;
		tipMatnum=startTip;
		for(int i=0;i<=resolution;i++){
			x=x0+a*cos(i*deltaSita+startSita);
			y=y0+b*sin(i*deltaSita+startSita);
			if(!crackCtrl->AddSegment(new CrackSegment(x,y,tipMatnum,MatID)))
				throw SAXException("Crack not in the mesh or out of memory adding a crack segment.");
			tipMatnum = i!=resolution-1 ? -1 : endTip;
		}
	}
}

//------------------------------------------------------------
// Find last loaded point
//-------------------------------------------------------------
char *MPMReadHandler::LastBC(char *firstBC)
{
	if(firstBC==NULL) return (char *)NULL;
	
    BoundaryCondition *lastBC=(BoundaryCondition *)firstBC;
	while(TRUE)
	{   if(lastBC->GetNextObject()==NULL) break;
		lastBC=(BoundaryCondition *)lastBC->GetNextObject();
	}
	return (char *)lastBC;
}

//-----------------------------------------------------------
// Subroutine for creating grid objects (nd, theElements)
//-----------------------------------------------------------
void MPMReadHandler::grid()
{
    int i,j,k,curPt=0,curEl=0;
	int node,element,enode[MaxElNd]={0};
    double xpt,ypt,zpt;
    
    // use cell sizes instead
    if(cellHoriz>0.)
    {	Nhoriz=(int)((Xmax-Xmin)/cellHoriz);
    	double newMax=Xmin+Nhoriz*cellHoriz;
    	if(newMax<Xmax)
    	{	Xmax=newMax+cellHoriz;
    		Nhoriz++;
    	}
    }
    if(cellVert>0.)
    {	Nvert=(int)((Ymax-Ymin)/cellVert);
    	double newMax=Ymin+Nvert*cellVert;
    	if(newMax<Ymax)
    	{	Ymax=newMax+cellVert;
    		Nvert++;
    	}
    }
    if(cellDepth>0. && fmobj->IsThreeD())
    {	Ndepth=(int)((Zmax-Zmin)/cellDepth);
    	double newMax=Zmin+Ndepth*cellDepth;
    	if(newMax<Zmax)
    	{	Zmax=newMax+cellDepth;
    		Ndepth++;
    	}
    }

	// check has grid settingsd
    if(Nhoriz<1 || Nvert<1 || (Ndepth<1 && fmobj->IsThreeD()))
        throw SAXException("Number of grid elements in all direction must be >= 1.");
	
	// save the limits before extra GIMP layer
	mxmin=Xmin;
	mxmax=Xmax;
	mymin=Ymin;
	mymax=Ymax;
	mzmin=Zmin;
	mzmax=Zmax;
		
	// allow for GIMP
	if(ElementBase::useGimp != POINT_GIMP)
	{	double cell=(Xmax-Xmin)/(double)Nhoriz;
		Nhoriz+=2;
		Xmin-=cell;
		Xmax+=cell;
		cell=(Ymax-Ymin)/(double)Nvert;
		Ymin-=cell;
		Ymax+=cell;
		Nvert+=2;
		if(fmobj->IsThreeD())
		{	cell=(Zmax-Zmin)/(double)Ndepth;
			Zmin-=cell;
			Zmax+=cell;
			Ndepth+=2;
		}
	}
    
	if(DbleEqual(Rhoriz,1.) && DbleEqual(Rvert,1.) && (!fmobj->IsThreeD() || DbleEqual(Rdepth,1.)))
	{	double zparam,gridz = 0.;
		if(fmobj->IsThreeD())
		{	zparam = Zmin;
			gridz = (Zmax-Zmin)/(double)Ndepth;
		}
		else if(fmobj->IsAxisymmetric())
			zparam = 1.0;
		else
			zparam = Z2DThickness;
		mpmgrid.SetCartesian(TRUE,(Xmax-Xmin)/(double)Nhoriz,(Ymax-Ymin)/(double)Nvert,gridz);
		mpmgrid.SetElements(Nhoriz,Nvert,Ndepth,Xmin,Ymin,zparam);
	}
	else
		mpmgrid.SetCartesian(FALSE,0.,0.,0.);
	
    // number of nodes and elements
	if(fmobj->IsThreeD())
    {	nnodes=(Nhoriz+1)*(Nvert+1)*(Ndepth+1);
		nelems=Nhoriz*Nvert*Ndepth;
	}
	else if(mpm2DElement==FOUR_NODE_ISO)
    {	nnodes=(Nhoriz+1)*(Nvert+1);
		nelems=Nhoriz*Nvert;
	}
	else
	{	// 9 node Lagranging and side nodes and internal nodes adding one per element in each direction
    	nnodes=(Nhoriz+1+Nhoriz)*(Nvert+1+Nvert);
		nelems=Nhoriz*Nvert;
	}
	
     // space for nodal points (1 based)
    curPt=1;
    nd=(NodalPoint **)malloc(sizeof(NodalPoint *)*(nnodes+1));
    if(nd==NULL) throw SAXException("Out of memory allocating space for nodes.");

    // space for elements (0 based)
    curEl=0;
    theElements=(ElementBase **)malloc(sizeof(ElementBase *)*nelems);
    if(theElements==NULL) throw SAXException("Out of memory allocating space for elements.");

    // create the elements
	if(fmobj->IsThreeD())
	{	// create the nodes: vary x first, then y, last z
		for(k=0;k<=Ndepth;k++)
		{	zpt=Zmin+(double)k*(Zmax-Zmin)/(double)Ndepth;
			for(j=0;j<=Nvert;j++)
			{	ypt=Ymin+(double)j*(Ymax-Ymin)/(double)Nvert;
				for(i=0;i<=Nhoriz;i++)
				{	xpt=Xmin+(double)i*(Xmax-Xmin)/(double)Nhoriz;
					zpt=Zmin+(double)k*(Zmax-Zmin)/(double)Ndepth;
					node=k*(Nhoriz+1)*(Nvert+1)+j*(Nhoriz+1)+(i+1);
					nd[curPt]=new NodalPoint3D(node,xpt,ypt,zpt);
					curPt++;
				}
			}
		}
		// create elements by z planes and then blocks at constant y
		int zplane=(Nhoriz+1)*(Nvert+1);
		for(k=1;k<=Ndepth;k++)
		{	for(j=1;j<=Nvert;j++)
			{	for(i=1;i<=Nhoriz;i++)
				{	element=(k-1)*Nhoriz*Nvert+(j-1)*Nhoriz+i;
					enode[0]=(k-1)*zplane+(j-1)*(Nhoriz+1)+i;
					enode[1]=enode[0]+1;
					enode[3]=enode[0]+(Nhoriz+1);
					enode[2]=enode[3]+1;
					enode[4]=enode[0]+zplane;
					enode[5]=enode[1]+zplane;
					enode[6]=enode[2]+zplane;
					enode[7]=enode[3]+zplane;
					theElements[curEl]=ElementsController::MeshElement(mpm3DElement,element,enode);
					theElements[curEl]->FindExtent();
					curEl++;
				}
			}
		}
	}
	else if(mpm2DElement==FOUR_NODE_ISO)
	{	// create the nodes vary x first, y second
		for(j=0;j<=Nvert;j++)
		{	ypt=Ymin+(double)j*(Ymax-Ymin)/(double)Nvert;
			for(i=0;i<=Nhoriz;i++)
			{	xpt=Xmin+(double)i*(Xmax-Xmin)/(double)Nhoriz;
				node=j*(Nhoriz+1)+(i+1);
				nd[curPt]=new NodalPoint2D(node,xpt,ypt);
				curPt++;
			}
		}
		// create the elements in blocks at constant y
		for(j=1;j<=Nvert;j++)
		{	for(i=1;i<=Nhoriz;i++)
			{	element=(j-1)*Nhoriz+i;
				enode[0]=(j-1)*(Nhoriz+1)+i;
				enode[1]=enode[0]+1;
				enode[2]=enode[0]+(Nhoriz+2);
				enode[3]=enode[2]-1;
				theElements[curEl]=ElementsController::MeshElement(mpm2DElement,element,enode);
				theElements[curEl]->FindExtent();
				curEl++;
			}
		}
	}
	else
	{	// 2D Lagrange adds side and internal nodes
		// create the nodes vary x first, y second
		int NvertHalfCells=2*Nvert,NhorizHalfCells=2*Nhoriz;
		for(j=0;j<=NvertHalfCells;j++)
		{	ypt=Ymin+(double)j*(Ymax-Ymin)/(double)NvertHalfCells;
			for(i=0;i<=NhorizHalfCells;i++)
			{	xpt=Xmin+(double)i*(Xmax-Xmin)/(double)NhorizHalfCells;
				node=j*(NhorizHalfCells+1)+(i+1);
				nd[curPt]=new NodalPoint2D(node,xpt,ypt);
				curPt++;
			}
		}
		// create the elements in blocks at constant y
		for(j=1;j<=Nvert;j++)
		{	for(i=1;i<=Nhoriz;i++)
			{	element=(j-1)*Nhoriz+i;
				enode[0]=(j-1)*2*(NhorizHalfCells+1)+2*i-1;
				enode[1]=enode[0]+2;
				enode[2]=enode[1]+2*(NhorizHalfCells+1);
				enode[3]=enode[2]-2;
				enode[4]=enode[0]+1;
				enode[5]=enode[1]+(NhorizHalfCells+1);
				enode[6]=enode[3]+1;
				enode[7]=enode[5]-2;
				enode[8]=enode[7]+1;
				theElements[curEl]=ElementsController::MeshElement(mpm2DElement,element,enode);
				theElements[curEl]->FindExtent();
				curEl++;
			}
		}
	}
}

//------------------------------------------------------------------
// If just created a GIMP grid, make all border elements as holes
//------------------------------------------------------------------
void SetMptAnglesFromFunctions(int numRotations,Vector *mpos,MPMBase *newMpt)
{
	if(numRotations==0) return;
	
	// evaluate each angle
	int i;
	double rotAngle[3];
	for(i=0;i<numRotations;i++)
	{	rotAngle[i]=FunctionValue(i+1,mpos->x,mpos->y,mpos->z,0.,0.,0.);
	}
	
	// supported rotation schemes
	
	// Single angle is trivial
	if(strcmp(rotationAxes,"Z")==0)
	{	newMpt->SetAnglez0InDegrees(rotAngle[0]);
	}
	
	else if(strcmp(rotationAxes,"Y")==0)
	{	newMpt->SetAngley0InDegrees(rotAngle[0]);
	}
	
	else if(strcmp(rotationAxes,"X")==0)
	{	newMpt->SetAnglex0InDegrees(rotAngle[0]);
	}
	
	// all double angle options supported
	
	else if(strcmp(rotationAxes,"XY")==0)
	{	double xx=rotAngle[0]*PI_CONSTANT/180.;
		double yy=rotAngle[1]*PI_CONSTANT/180.;
		double R11=cos(yy);
		double R12=sin(xx)*sin(yy);
		double R13=cos(xx)*sin(yy);
		double R21=0.;
		double R22=cos(xx);
		double R23=-sin(xx);
		double R31=-sin(yy);
		double R33=cos(xx)*cos(yy);
		ConvertToZYX(newMpt,R11,R12,R13,R21,R22,R23,R31,R33);
	}
	
	else if(strcmp(rotationAxes,"XZ")==0)
	{	double xx=rotAngle[0]*PI_CONSTANT/180.;
		double zz=rotAngle[1]*PI_CONSTANT/180.;
		double R11=cos(zz);
		double R12=-cos(xx)*sin(zz);
		double R13=sin(xx)*sin(zz);
		double R21=sin(zz);
		double R22=cos(xx)*cos(zz);
		double R23=-cos(zz)*sin(xx);
		double R31=0.;
		double R33=cos(xx);
		ConvertToZYX(newMpt,R11,R12,R13,R21,R22,R23,R31,R33);
	}
	
	else if(strcmp(rotationAxes,"YX")==0)
	{	newMpt->SetAngley0InDegrees(rotAngle[0]);
		newMpt->SetAnglex0InDegrees(rotAngle[1]);
	}
	
	else if(strcmp(rotationAxes,"YZ")==0)
	{	double yy=rotAngle[0]*PI_CONSTANT/180.;
		double zz=rotAngle[1]*PI_CONSTANT/180.;
		double R11=cos(yy)*cos(zz);
		double R12=-sin(zz);
		double R13=cos(zz)*sin(yy);
		double R21=cos(yy)*sin(zz);
		double R22=cos(zz);
		double R23=sin(zz)*sin(yy);
		double R31=-sin(yy);
		double R33=cos(yy);
		ConvertToZYX(newMpt,R11,R12,R13,R21,R22,R23,R31,R33);
	}
	
	else if(strcmp(rotationAxes,"ZX")==0)
	{	newMpt->SetAnglez0InDegrees(rotAngle[0]);
		newMpt->SetAnglex0InDegrees(rotAngle[1]);
	}
	
	else if(strcmp(rotationAxes,"ZY")==0)
	{	newMpt->SetAnglez0InDegrees(rotAngle[0]);
		newMpt->SetAngley0InDegrees(rotAngle[1]);
	}
	
	// Triple Angle options
	
	else if(strcmp(rotationAxes,"ZYX")==0)
	{	newMpt->SetAnglez0InDegrees(rotAngle[0]);
		newMpt->SetAngley0InDegrees(rotAngle[1]);
		newMpt->SetAnglex0InDegrees(rotAngle[2]);
	}
	
	else if(strcmp(rotationAxes,"ZYZ")==0)
	{	double z=rotAngle[0]*PI_CONSTANT/180.;
		double y=rotAngle[1]*PI_CONSTANT/180.;
		double x=rotAngle[2]*PI_CONSTANT/180.;		// second z
		double R11 = cos(x)*cos(y)*cos(z) - sin(x)*sin(z);
		double R12 = -cos(z)*sin(x) - cos(x)*cos(y)*sin(z);
		double R13 = cos(x)*sin(y);
		double R21 = cos(y)*cos(z)*sin(x) + cos(x)*sin(z);
		double R22 = cos(x)*cos(z) - cos(y)*sin(x)*sin(z);
		double R23 = sin(x)*sin(y);
		double R31 = -cos(z)*sin(y);
		double R33 = cos(y);
		ConvertToZYX(newMpt,R11,R12,R13,R21,R22,R23,R31,R33);
	}
	
	else
	{	char badAngles[200];
		strcpy(badAngles,"'");
		strcat(badAngles,rotationAxes);
		strcat(badAngles,"' is an unsupported rotation scheme for <RotateX(Y)(Z)> commands.");
		throw SAXException(badAngles);
	}
}

// Decompose matrix, input most elements (never need R32) but if R13 is 1, only need R12, R21, and R22 (may not be worth a check)
void ConvertToZYX(MPMBase *newMpt,double R11,double R12,double R13,double R21,double R22,double R23,double R31,double R33)
{
	// this as two roots (angle and PI-angle), but either one can be used (i.e., two solutions for
	//	xyz are x,y,x and -PI+x,PI-y,-PI+z. Here we take the first one
	double y=asin(R13);
	double x,z;
	
	// Special case for y=±pi/2
	if(DbleEqual(fabs(R13),1))
	{	// can only deterimine x+z, set x=0 and find z
		x=0.;
		z=asin(R21);
		if(!DbleEqual(R22,cos(z))) z=PI_CONSTANT-z;
	}
	
	// assume y is not ±pi/2
	else
	{	if(!DbleEqual(R11,0))
		{	z=atan(-R12/R11);
			if(!DbleEqual(R11,cos(y)*cos(z))) z-=PI_CONSTANT;
		}
		else
		{	// Cos[z]=0, Sin[z]=±1, Cos[y]­0
			if(DbleEqual(R12,-cos(y)))
				z=PI_CONSTANT/2.;		// Sin[z]=1
			else
				z=-PI_CONSTANT/2.;		// Sin[z]=-1
		}
		
		if(!DbleEqual(R33,0))
		{	x=atan(-R23/R33);
			if(!DbleEqual(R33,cos(y)*cos(x))) x-=PI_CONSTANT;
		}
		else
		{	// Cos[x]=0, Sin[x]=±1, Cos[y]­0
			if(DbleEqual(R23,-cos(y)))
				x=PI_CONSTANT/2.;		// Sin[x]=1
			else
				x=-PI_CONSTANT/2.;		// Sin[x]=-1
		}
	}
	
	// set decomposed angles
	newMpt->SetAnglez0(z);
	newMpt->SetAngley0(y);
	newMpt->SetAnglex0(x);
}

//------------------------------------------------------------
// swap Xmin and Xmax, Ymin and Ymax if Xmin>Xmax or Ymin>Ymax
// if 3D, do same for Zmin and Zmax
//-------------------------------------------------------------
void swap(double& Xmin, double& Xmax, double& Ymin, double& Ymax, double& Zmin, double& Zmax)
{  
    double temp;
    if(Xmin>Xmax) {
        temp=Xmax;
        Xmax=Xmin;
        Xmin=temp;
    }
    if(DbleEqual(Xmin,Xmax)) {
        throw SAXException("xmax can not equal xmin in input parameters.");
    }
    if(Ymin>Ymax) {
        temp=Ymax;
        Ymax=Ymin;
        Ymin=temp;
    }
    if(DbleEqual(Ymin,Ymax)) {
        throw SAXException("ymax can not equal ymin in input parameters.");
    }
	if(fmobj->np!=THREED_MPM) return;
	
    if(Zmin>Zmax) {
        temp=Zmax;
        Zmax=Zmin;
        Zmin=temp;
    }
    if(DbleEqual(Zmin,Zmax)) {
        throw SAXException("zmax can not equal zmin in input parameters.");
    }
}

