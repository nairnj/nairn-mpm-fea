/********************************************************************************
    FEAReadHandler.cpp
    NairnFEA
    
    Created by John Nairn on Feb 2 2003.
    Copyright (c) 2003 John A. Nairn. All rights reserved.
********************************************************************************/

#include "NairnFEA_Class/NairnFEA.hpp"
#include "Read_FEA/FEAReadHandler.hpp"
#include "Elements/ElementBase.hpp"
#include "Boundary_Conditions/NodalDispBC.hpp"
#include "Boundary_Conditions/NodalLoad.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Boundary_Conditions/EdgeBC.hpp"
#include "Read_FEA/Keypoint.hpp"
#include "Read_FEA/KeypointsController.hpp"
#include "Read_FEA/Path.hpp"
#include "Read_FEA/PathsController.hpp"
#include "Read_FEA/Area.hpp"
#include "Read_XML/NodesController.hpp"
#include "Read_XML/PointController.hpp"
#include "Read_FEA/PathBCController.hpp"
#include "Read_XML/ElementsController.hpp"
#include "Read_FEA/EdgeBCController.hpp"
#include "Read_XML/LineController.hpp"
#include "Read_XML/MaterialController.hpp"
#include "Read_FEA/NodalDispBCController.hpp"
#include "Read_FEA/NodalLoadController.hpp"
#include "Read_FEA/ConstraintController.hpp"
#include "Boundary_Conditions/Constraint.hpp"
#include <algorithm>

/********************************************************************************
	FEAReadHandler: Constructors and Destructor
********************************************************************************/

FEAReadHandler::FEAReadHandler()
{
	resequence=-1;
}

FEAReadHandler::~FEAReadHandler()
{
	// delete keypoints and paths now
	if(keyPts!=NULL) delete keyPts;
	if(paths!=NULL) delete paths;
	
	// remove unused elements and check material types
	RemoveEmptyElements();
	
	// find periodic nodes
	if(fmobj->periodic.dof>0) FindPeriodicNodes();
	
	// resequence if requested (but not yet supported when periodic)
	if(resequence>0) ResequenceNodes();
	delete theNodes;
	
	// sort and remove redundancy
	if(fmobj->selectedNodes.size()>0)
	{	std::sort(fmobj->selectedNodes.begin(),fmobj->selectedNodes.end());
		vector< int >::iterator iter=std::unique(fmobj->selectedNodes.begin(),fmobj->selectedNodes.end());
		// returned iter is beginning of left over data, which needs to be deleted
		if(iter!=fmobj->selectedNodes.end())
			fmobj->selectedNodes.erase(iter);
		// remove a 0 if there which corresponds to removed nodes
		iter=std::remove(fmobj->selectedNodes.begin(),fmobj->selectedNodes.end(),0);
		// returned iter is beginning of left over data, which needs to be deleted
		if(iter!=fmobj->selectedNodes.end())
			fmobj->selectedNodes.erase(iter);
	}
}

/********************************************************************************
	FEAReadHandler: Methods
********************************************************************************/

// Custom FEA element start
bool FEAReadHandler::myStartElement(char *xName,const Attributes& attrs)
{
    char *value,*aName;
    char matName[100];
    int i,numAttr;
    double x,y,temp;
	
    //-------------------------------------------------------
    // <Header> data for FEA
	
	// select output options
    if(strcmp(xName,"Output")==0)
	{	ValidateCommand(xName,HEADER,MUST_BE_2D);
         input=OUTFLAGS_BLOCK;
    }
    
	// Select a node for output when explicitly listed
    else if(strcmp(xName,"Select")==0)
	{	ValidateCommand(xName,HEADER,MUST_BE_2D);
		value=ReadTagValue("node",attrs);
		if(value!=NULL)
		{	int nodeNum;
			sscanf(value,"%d",&nodeNum);
			fmobj->selectedNodes.push_back(nodeNum);
			delete [] value;
		}
	}
    
    //-------------------------------------------------------
	// <Mesh> section
	
	// begin mesh section
	else if(strcmp(xName,"Mesh")==0)
	{	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
		block=MESHBLOCK;
		theNodes=new NodesController();			// deleted when Mesh is done
		theElems=new ElementsController();		// deleted when Mesh is done
	}
	
	// Keypoints list - implies generated mesh
	else if(strcmp(xName,"Keypoints")==0)
	{	ValidateCommand(xName,MESHBLOCK,MUST_BE_2D);
    	if(meshType!=UNKNOWN_MESH)
            throw SAXException("Invalid use of <Keypoints> block.");
		keyPts=new KeypointsController();		// deleted when XML parsing is done
	    block=KEYPOINTBLOCK;
		meshType=GENERATED_MESH;
	}
	
	// define a path
	else if(strcmp(xName,"Path")==0)
	{	ValidateCommand(xName,MESHBLOCK,MUST_BE_2D);
    	if(meshType!=GENERATED_MESH)
            throw SAXException("Invalid use of <Path> block.");
		if(paths==NULL) paths=new PathsController();		// delete when XML parsing is done
		block=PATHBLOCK;
		matName[0]=0;
		double ratio=1.;
		int intervals=1;
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"intervals")==0)
				sscanf(value,"%d",&intervals);
            else if(strcmp(aName,"ratio")==0)
				sscanf(value,"%lf",&ratio);
			else if(strcmp(aName,"id")==0)
				strcpy(matName,value);
            delete [] aName;
            delete [] value;
        }
		
		// check and add new path
		if(!paths->ValidName(matName))
			ThrowCatErrorMessage("Duplicate or invalid Path ID",matName);
		paths->AddObject(new Path(matName,intervals,ratio));	// Path delete when paths is deleted
	}
	
	// define an area
	else if(strcmp(xName,"Area")==0)
	{	ValidateCommand(xName,MESHBLOCK,MUST_BE_2D);
    	if(meshType!=GENERATED_MESH || paths==NULL)
            throw SAXException("Invalid use of <Area> block.");
		block=AREABLOCK;
		int matnum=0;			// must set number or name
		char matname[200];
		matname[0]=0;
		char *angleExpr=NULL;		// optional
		double thick=1.0;		// optional
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"mat")==0)
				sscanf(value,"%d",&matnum);
			else if(strcmp(aName,"matname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(matname,value);
			}
            else if(strcmp(aName,"angle")==0)
			{	if(angleExpr!=NULL) delete [] angleExpr;
				angleExpr=new char[strlen(value)+1];
				strcpy(angleExpr,value);
			}
			else if(strcmp(aName,"thick")==0)
				sscanf(value,"%lf",&thick);
            else if(strcmp(aName,"type")==0)
			{	if(!theElems->SetElemIDStr(value))
					throw SAXException("Invalid or incompatible element type.");
			}
			else if(strcmp(aName,"flip")==0)
				theElems->SetFlipTriangles(value);
            delete [] aName;
            delete [] value;
        }
		
		// if gave a matname, it takes precedence over mat number
		if(strlen(matname)>0)
			matnum = matCtrl->GetIDFromNewName(matname);
		if(matnum<0)
			throw SAXException("Missing or invalid material number for Area.");
		theArea=new Area(matnum,angleExpr,thick);		// deleted when Area is done
	}
	
	// point in NodeList or Keypoints list
    else if(strcmp(xName,"pt")==0)
	{	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
        x=y=temp=0.;
		matName[0]=0;
    	numAttr=attrs.getLength();
		
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"x")==0)
				sscanf(value,"%lf",&x);
            else if(strcmp(aName,"y")==0)
				sscanf(value,"%lf",&y);
			else if(strcmp(aName,"temp")==0)
				sscanf(value,"%lf",&temp);
			else if(strcmp(aName,"id")==0)
				strcpy(matName,value);
            delete [] aName;
            delete [] value;
        }

		// make node or keypoint
        if(block==NODELIST)
			theNodes->AddNode(x,y,(double)0.,temp);
		else if(block==KEYPOINTBLOCK)
		{	if(!keyPts->ValidName(matName))
				ThrowCatErrorMessage("Duplicate or invalid Keypoint ID",matName);
			keyPts->AddObject(new Keypoint(matName,x,y));		// deleted when keyPts is deleted
		}
        else
			ValidateCommand(xName,BAD_BLOCK,MUST_BE_2D);
    }
	
	// Keypoint in path definition
    else if(strcmp(xName,"keypt")==0)
	{	ValidateCommand(xName,PATHBLOCK,MUST_BE_2D);
		value=ReadTagValue("id",attrs);
		if(value!=NULL)
		{	if(!paths->AddKeypoint(value))
				ThrowCatErrorMessage("Invalid keypoint or too many keypoints for a path",value);
			delete [] value;
		}
		else
			throw SAXException("<keypt> name in an <Path> not provided.");
	}
    
	// Path in area definition
    else if(strcmp(xName,"path")==0)
	{	ValidateCommand(xName,AREABLOCK,MUST_BE_2D);
    	if(theArea==NULL)
            throw SAXException("<path> must be within an <Area> element.");
		value=ReadTagValue("id",attrs);
		if(value!=NULL)
		{	if(!theArea->AddPath(value))
				ThrowCatErrorMessage("Invalid or disconnected path added to an area",value);
			delete [] value;
		}
		else
			throw SAXException("<path> name in an <Area> not provided.");
	}

	// Explicit element definition
    else if(strcmp(xName,"elem")==0)
	{	ValidateCommand(xName,ELEMENTLIST,MUST_BE_2D);
    	input=NODE_BLOCK;
        elemMat=1;
        elemAngle=0.;
        elemThick=1.;
		int elemTypeSet=FALSE;
    	numAttr=attrs.getLength();
		char elmat[200];
		elmat[0]=0;
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"type")==0)
			{	if(!theElems->SetElemIDStr(value))
					throw SAXException("Invalid or incompatible element type.");
				elemTypeSet=TRUE;
			}
            else if(strcmp(aName,"matl")==0 || strcmp(aName,"mat")==0)
                sscanf(value,"%d",&elemMat);
			else if(strcmp(aName,"matname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(elmat,value);
			}
            else if(strcmp(aName,"angle")==0)
                sscanf(value,"%lf",&elemAngle);
            else if(strcmp(aName,"thick")==0)
                sscanf(value,"%lf",&elemThick);
            delete [] aName;
            delete [] value;
        }
		// if gave a matname, it takes precedence over mat number
		if(strlen(elmat)>0)
			elemMat = matCtrl->GetIDFromNewName(elmat);
		if(!elemTypeSet)
            throw SAXException("<elem> does not specify element type.");
    }

    //-----------------------------------------------------------
    // <GridBCs> section
	
	// begin GridBCs section
    else if(strcmp(xName,"GridBCs")==0)
	{	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
		if(firstDispBC!=NULL || firstLoadBC!=NULL || firstEdgeBC!=NULL)
			throw SAXException("More than one <GridBCs> block");
	    block=GRIDBCHEADER;
		dispBCCtrl=new NodalDispBCController();		// deleted when GridBCs ends
		loadBCCtrl=new NodalLoadController();		// deleted when GridBCs ends
		edgeBCCtrl=new EdgeBCController();			// deleted when GridBCs ends
    }
	
	// LoadBCs section
    else if(strcmp(xName,"LoadBCs")==0)
	{	ValidateCommand(xName,GRIDBCHEADER,MUST_BE_2D);
	    block=LOADEDNODES;
    }
        
	// EdgeBCs section (stress on edges)
    else if(strcmp(xName,"EdgeBCs")==0)
	{	ValidateCommand(xName,GRIDBCHEADER,MUST_BE_2D);
    	block=LOADEDFACES;
    }
        
    // Store a line and tolerance into the current line (in theShape)
    else if(strcmp(xName,"BCLine")==0)
	{	ValidateCommand(xName,GRIDBCHEADER,MUST_BE_2D);
        numAttr=attrs.getLength();
        double tolerance=ElementBase::gridTolerance;
		double xmin=0.,xmax=0.,ymin=0.,ymax=0.;
		bool select=FALSE;
		Path *bcPath=NULL;
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"x1")==0)
                sscanf(value,"%lf",&xmin);
            else if(strcmp(aName,"y1")==0)
                sscanf(value,"%lf",&ymin);
            else if(strcmp(aName,"x2")==0)
                sscanf(value,"%lf",&xmax);
            else if(strcmp(aName,"y2")==0)
                sscanf(value,"%lf",&ymax);
            else if(strcmp(aName,"tolerance")==0)
                sscanf(value,"%lf",&tolerance);
            else if(strcmp(aName,"path")==0)
			{	if(paths!=NULL) bcPath=paths->FindPath(value);
				if(bcPath==NULL)
					throw SAXException("Invalid path in <BCLine>");
			}
            else if(strcmp(aName,"select")==0)
                select=TRUE;
            delete [] aName;
            delete [] value;
        }
		
		if(bcPath==NULL)
			theShape=new LineController(block,xmin,xmax,ymin,ymax,tolerance);
		else
			theShape=new PathBCController(block,bcPath);
        block=BCSHAPE;
		if(select)
		{	theShape->resetNodeEnumerator();
			while((i=theShape->nextNode()))
				fmobj->selectedNodes.push_back(i);
		}
    }

    // Store a line and tolerance into the current line (in theShape)
    else if(strcmp(xName,"BCPt")==0 || strcmp(xName,"Resequence")==0 || strcmp(xName,"Cracktip")==0)
	{	ValidateCommand(xName,GRIDBCHEADER,MUST_BE_2D);
		if(strcmp(xName,"Cracktip")==0)
		{	if(dispBCCtrl->numObjects>0 || loadBCCtrl->numObjects>0 || edgeBCCtrl->numObjects>0)
				throw SAXException("<Cracktip> must be first command in <GridBCs> section");
		}
        numAttr=attrs.getLength();
		double xpt=0.,ypt=0.;
		int nearestNode=-1;
		bool select=FALSE;
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"x")==0)
                sscanf(value,"%lf",&xpt);
            else if(strcmp(aName,"y")==0)
                sscanf(value,"%lf",&ypt);
            else if(strcmp(aName,"keypt")==0)
			{	if(keyPts!=NULL)
				{	Keypoint *bcKey=keyPts->FindKeypoint(value);
					if(bcKey!=NULL)
						nearestNode=bcKey->node;
				}
				if(nearestNode<=0)
					throw SAXException("Invalid keypoint in <BCPt> or <Resequence>");
			}
            else if(strcmp(aName,"select")==0)
                select=TRUE;
            delete [] aName;
            delete [] value;
        }
		
		// find nearest node (if needed)
		if(nearestNode<=0)
		{	double dist,minDist=1e30;
			nearestNode=0;
			for(i=1;i<=nnodes;i++)
			{	dist=(xpt-nd[i]->x)*(xpt-nd[i]->x)+(ypt-nd[i]->y)*(ypt-nd[i]->y);
				if(dist<minDist)
				{	minDist=dist;
					nearestNode=i;
				}
			}
		}

		if(strcmp(xName,"BCPt")==0)
		{	// create special case line to reference the node (deleted when BCLine done)
			theShape=new PointController(block,nearestNode);
			block=BCSHAPE;
			if(select) fmobj->selectedNodes.push_back(nearestNode);
		}
		else if(strcmp(xName,"Resequence")==0)
			resequence=nearestNode;
		else
			ElementBase::MoveCrackTipNodes(nearestNode);
    }

	// Set periodic in one direction
    else if(strcmp(xName,"Periodic")==0)
	{	ValidateCommand(xName,GRIDBCHEADER,MUST_BE_2D);
		int dof=0;
		bool fixDu=FALSE,fixDudy=FALSE;
		double du=0.,dudy=0.;
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"dof")==0)
				dof=GetDOFAttribute(value);
			else if(strcmp(aName,"delta")==0)
			{	fixDu=TRUE;
				sscanf(value,"%lf",&du);
			}
			else if(strcmp(aName,"slope")==0 || strcmp(aName,"shear")==0)
			{	fixDudy=TRUE;
				sscanf(value,"%lf",&dudy);
			}
            delete [] aName;
            delete [] value;
        }
		
		// x direction
		if(dof==1)
		{	if(fmobj->periodic.dof&0x01)
				throw SAXException("<Periodic> in x direction was already set.");
			if(fmobj->IsAxisymmetric())
				throw SAXException("<Periodic> in r direction for axisymmetric problems is not allowed.");
			fmobj->periodic.dof|=dof;
			fmobj->periodic.fixDu=fixDu;
			fmobj->periodic.du=du;
			fmobj->periodic.fixDudy=fixDudy;
			fmobj->periodic.dudy=dudy;
		}
		
		// y direction
		else if(dof==2 || dof==3)
		{	if(fmobj->periodic.dof&0x02)
				throw SAXException("<Periodic> in y or z (axisymmetric) direction was already set.");
			if(dof==3)
			{	if(fmobj->IsAxisymmetric())
					dof=2;
				else
					throw SAXException("<Periodic> has invalid dof attribute.");
			}
			fmobj->periodic.dof|=dof;
			fmobj->periodic.fixDv=fixDu;
			fmobj->periodic.dv=du;
			if(fmobj->IsAxisymmetric())
			{	if(fixDudy && !DbleEqual(dudy,0.))
					throw SAXException("Axisymmetric problems cannot set non-=zero shear or slope periodicity.");
				fmobj->periodic.fixDvdx=TRUE;
				fmobj->periodic.dvdx=0.;
			}
			else
			{	fmobj->periodic.fixDvdx=fixDudy;
				fmobj->periodic.dvdx=dudy;
			}
		}
		
		// error
		else
            throw SAXException("<Periodic> has missing or invalid dof attribute.");
    }

	// Set fixed displacement at a node
    else if(strcmp(xName,"fix")==0)
	{	ValidateCommand(xName,FIXEDNODES,MUST_BE_2D);
		int node=0;
		int dof=0;
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"node")==0)
				sscanf(value,"%d",&node);
            else if(strcmp(aName,"dof")==0)
				dof=GetDOFAttribute(value);
            delete [] aName;
            delete [] value;
        }
        if(node==0 || dof==0)
            throw SAXException("<fix> missing required attribute.");
        
        // create object and get input
        NodalDispBC *newDispBC=new NodalDispBC(node,dof);
		if(dispBCCtrl->AddObject(newDispBC))
    	{	input=DOUBLE_NUM;
			inputPtr=(char *)newDispBC->GetValuePtr();
		}
    }
    
	// rotate coordinates at a node for skewed BCs
    else if(strcmp(xName,"rotate")==0)
	{	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
    	if(block!=FIXEDNODES && block!=BCSHAPE)
            ValidateCommand(xName,BAD_BLOCK,MUST_BE_2D);
		int node=0;
		int axis=0;
		double angle=0.;
		bool hasAngle=FALSE;
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"node")==0)
				sscanf(value,"%d",&node);
            else if(strcmp(aName,"axis")==0)
				axis=GetDOFAttribute(value);
            else if(strcmp(aName,"angle")==0)
			{	sscanf(value,"%lf",&angle);
				hasAngle=TRUE;
			}
            delete [] aName;
            delete [] value;
        }
        if((node==0 && block==FIXEDNODES) || axis==0 || (!hasAngle && block==BCSHAPE))
            throw SAXException("<rotate> missing required attribute.");
		if(axis!=3)
            throw SAXException("<rotate> axis can only be 'z' or '3' in 2D calculations.");
        
        // create object and get input
		if(block==FIXEDNODES)
		{	NodalDispBC *newDispBC=new NodalDispBC(-node,axis);
			if(dispBCCtrl->AddObject(newDispBC))
			{	if(hasAngle)
					newDispBC->angle=angle;
				else
				{	input=DOUBLE_NUM;
					inputPtr=(char *)&newDispBC->angle;
				}
			}
		}
		else
		{	// set one node or check all node
			NodalDispBC *newDispBC;
			const char *msg=theShape->startNodeEnumerator(ROTATE_PATH_NODES,0);
			if(msg!=NULL) throw SAXException(msg);
			while((i=theShape->nextNode()))
			{	newDispBC=new NodalDispBC(-nd[i]->num,axis);
				if(dispBCCtrl->AddObject(newDispBC))
					newDispBC->angle=angle;
			}
		}
    }
	
	// Apply a load to a node
    else if(strcmp(xName,"load")==0)
	{	ValidateCommand(xName,LOADEDNODES,MUST_BE_2D);
		int node=0;
		int dof=0;
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"node")==0)
				sscanf(value,"%d",&node);
            else if(strcmp(aName,"dof")==0)
				dof=GetDOFAttribute(value);
            delete [] aName;
            delete [] value;
        }
        if(node<=0 || dof<=0)
            throw SAXException("<load> missing required attribute.");
        
        // create object and get input
        NodalLoad *newLoadBC=new NodalLoad(node,dof);
		if(loadBCCtrl->AddObject(newLoadBC))
    	{	input=DOUBLE_NUM;
			inputPtr=(char *)newLoadBC->GetValuePtr();
		}
    }
    
	// Apply stress to face of an element
    else if(strcmp(xName,"stress")==0)
	{	ValidateCommand(xName,LOADEDFACES,MUST_BE_2D);
		int elem=0,face=0,dir=0;
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"elem")==0)
				sscanf(value,"%d",&elem);
			if(strcmp(aName,"face")==0)
 				sscanf(value,"%d",&face);
			else if(strcmp(aName,"dir")==0)
				dir=GetDOFAttribute(value);
            delete [] aName;
            delete [] value;
        }
        if(elem<=0 || face<=0 || dir<=0 || dir>2 || elem>nelems)
            throw SAXException("<stress> has missing or invalid attribute.");
        
        // create object
		if(edgeBCCtrl->AddObject(new EdgeBC(elem,face,dir)))
			input=STRESS_LIST;
    }

    // Read into displacements boundary conditions on current line (direction, value or function)
    else if(strcmp(xName,"DisBC")==0)
	{	ValidateCommand(xName,BCSHAPE,MUST_BE_2D);
        double disp=0.;
		int dof=0;
		char *function=NULL;
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"dof")==0)
				dof=GetDOFAttribute(value);
            else if(strcmp(aName,"disp")==0)
                sscanf(value,"%lf",&disp);
            else if(strcmp(aName,"function")==0)
			{	if(function!=NULL) delete [] function;
				function=new char[strlen(value)+1];
				strcpy(function,value);
			}
            delete [] aName;
            delete [] value;
        }
        if(dof>2 || dof<1)
            throw SAXException("'dof' must be either 1 or 2 in DisBC element.");
        
		// set all nodes on the line
		NodalDispBC *newDispBC;
		const char *msg=theShape->startNodeEnumerator(FIX_PATH_NODES,dof);
		if(msg!=NULL) throw SAXException(msg);
		while((i=theShape->nextNode()))
		{	newDispBC=new NodalDispBC(nd[i]->num,dof);
			if(dispBCCtrl->AddObject(newDispBC))
				newDispBC->SetValue(disp,function);
		}
		if(function!=NULL) delete [] function;
    }

    // Read into load boundary conditions (dirction, value or function
    else if(strcmp(xName,"LoadBC")==0)
	{	ValidateCommand(xName,BCSHAPE,MUST_BE_2D);
        int dof=0;
        double load=0.0;
		char *function=NULL;
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"dof")==0)
				dof=GetDOFAttribute(value);
            else if(strcmp(aName,"load")==0)
                sscanf(value,"%lf",&load);
            else if(strcmp(aName,"function")==0)
			{	if(function!=NULL) delete [] function;
				function=new char[strlen(value)+1];
				strcpy(function,value);
			}
            delete [] aName;
            delete [] value;
        }
        if(dof>2 || dof<1)
            throw SAXException("'dof' must be either 1 or 2 in <LoadBC> element.");
		
		// set one node or check all node
		NodalLoad *newLoadBC;
		const char *msg=theShape->startNodeEnumerator(LOAD_PATH_NODES,dof);
		if(msg!=NULL) throw SAXException(msg);
		while((i=theShape->nextNode()))
		{	newLoadBC=new NodalLoad(nd[i]->num,dof);
			if(loadBCCtrl->AddObject(newLoadBC))
				newLoadBC->SetValue(load,function);
		}
		if(function!=NULL) delete [] function;
    }

    // Read into stress boundary conditions (direction,style,values)
    else if(strcmp(xName,"StressBC")==0)
	{	ValidateCommand(xName,BCSHAPE,MUST_BE_2D);
		Path *thePath=(Path *)(theShape->GetContextInfo());
		if(thePath==NULL)
            throw SAXException("<StressBC> can only be applied to paths.");
        int dir=0;
		int nloads=0;
        double load[3]={0.,0.,0.};
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"dir")==0)
				dir=GetDOFAttribute(value);
            else if(strcmp(aName,"stress")==0)
			{	int j=0;
				nloads=1;
				while(value[j]!=0)
				{	if(value[j]==',') nloads++;
					j++;
				}
				if(nloads==1)
					sscanf(value,"%lf",&load[0]);
				else if(nloads==2)
					sscanf(value,"%lf,%lf",&load[0],&load[1]);
				else if(nloads==3)
					sscanf(value,"%lf,%lf,%lf",&load[0],&load[1],&load[2]);
				else
					throw SAXException("Too many loads in a <StressBC> element.");
			}
            delete [] aName;
            delete [] value;
        }
        if(dir>2 || dir<1)
            throw SAXException("'dir' must be either 'n' or 't' in <StressBC> element.");
		if(nloads==0)
			throw SAXException("No stress given for <StressBC> element.");
		
		// set one node or check all node
		const char *msg=thePath->VerifyBCOption(LOAD_PATH_EDGES,dir);
		if(msg!=NULL) throw SAXException(msg);
		thePath->AddEdgeBCsToPath(dir,nloads,load);
    }

    //-------------------------------------------------------------------
    // BMP file generator commands (must come before Temperature below)
    else if(BMPFileInput(xName,attrs))
    {
    }

    //-------------------------------------------------------
    // <Thermal> section
    
	// set temperature expression
    else if(strcmp(xName,"Temperature")==0)
	{	ValidateCommand(xName,THERMAL,MUST_BE_2D);
    	input=TEXT_BLOCK;
        inputID=TEMPERATURE_EXPR;
    }
    
	// set stress free temperature
    else if(strcmp(xName,"StressFreeTemp")==0)
	{	ValidateCommand(xName,THERMAL,MUST_BE_2D);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&fmobj->stressFreeTemperature;
    }

	//-------------------------------------------------------
	// Unknown element
	else
		return FALSE;
    
    return TRUE;
}

// End an element
void FEAReadHandler::myEndElement(char *xName)
{
	if(strcmp(xName,"GridBCs")==0)
	{	// set r==0 nodes if axisymmetric
		if(fmobj->np==AXI_SYM)
		{	int i;
			for(i=1;i<=nnodes;i++)
			{	if(DbleEqual(0.,nd[i]->x))
					dispBCCtrl->AddObject(new NodalDispBC(nd[i]->num,1));
			}
 		}
		firstDispBC=(NodalDispBC *)dispBCCtrl->firstObject;
		firstLoadBC=(NodalLoad *)loadBCCtrl->firstObject;
		firstEdgeBC=(EdgeBC *)edgeBCCtrl->firstObject;
		delete dispBCCtrl;
		delete loadBCCtrl;
		delete edgeBCCtrl;
	}
	
	else if(strcmp(xName,"LoadBCs")==0 || strcmp(xName,"EdgeBCs")==0)
	{   block=GRIDBCHEADER;
	}
	
	else if(strcmp(xName,"Keypoints")==0)
	{	block=MESHBLOCK;
	}
	
	else if(strcmp(xName,"Path")==0)
	{	if(!paths->ValidPath())
			throw SAXException("A path has insufficient keypoints.");
		block=MESHBLOCK;
	}
	
	else if(strcmp(xName,"Area")==0)
	{	const char *errMsg=theArea->MeshElements();
		if(errMsg!=NULL)
			throw SAXException(errMsg);
		delete theArea;
		theArea=NULL;
		block=MESHBLOCK;
	}
	
    else if(strcmp(xName,"BCLine")==0 || strcmp(xName,"BCPt")==0)
	{	block=theShape->GetSourceBlock();
		delete theShape;
		theShape=NULL;
    }
	
    else if(EndBMPInput(xName,MESHBLOCK))
    {   
    }
}

// Decode block of characters for FEA input
void FEAReadHandler::myCharacters(char *xData,const unsigned int length)
{
    switch(input)
	{	case TEXT_BLOCK:
			switch(inputID)
			{	case TEMPERATURE_EXPR:
					if(fmobj->temperatureExpr!=NULL) delete [] fmobj->temperatureExpr;
					fmobj->temperatureExpr=new char[strlen(xData)+1];
					strcpy(fmobj->temperatureExpr,xData);
					break;
				default:
					break;
			}
			break;
        
        case NODE_BLOCK:
			if(!theElems->CreateElement(xData,elemMat,elemAngle,elemThick))
				throw SAXException("Unknown element type was found in <elem> command.");
            break;
        
		case STRESS_LIST:
			if(!edgeBCCtrl->SetStress(xData))
				throw SAXException("Unknown number of stresses found in <stress> command.");
			break;
		
		case OUTFLAGS_BLOCK:
			int i;
			for(i=0;i<NUMBER_OUT;i++)
			{	if(xData[i]==0) break;
				fmobj->outFlags[i]=xData[i];
			}
			break;
        
        default:
            break;
    }
}

/********************************************************************************
	 Translate dof attribute into integer. It can be 1 (or n or x),
	 2 (or t or y), or 3 (or z). Anything else returns zero
********************************************************************************/

int FEAReadHandler::GetDOFAttribute(char *value)
{	
	int i=0,dof=0;
	while(value[i]!=0)
	{	if(value[i]==' ' || value[i]=='\t') continue;
		if(value[i]=='1' || value[i]=='x' || value[i]=='n' || value[i]=='X' || value[i]=='N')
			dof=1;
		else if(value[i]=='2' || value[i]=='y' || value[i]=='t' || value[i]=='Y' || value[i]=='T')
			dof=2;
		else if(value[i]=='3' || value[i]=='z' || value[i]=='Z')
			dof=3;
		break;
	}
	return dof;
}

/********************************************************************************
	 Resequence nodes starting at node resequence
********************************************************************************/

#define MAX_CONNECTIVITY 41

typedef struct {
	int degree;
	int cons[MAX_CONNECTIVITY];
} ConnectRec;

void FEAReadHandler::ResequenceNodes(void)
{
	int i,j,k,l,numnds;
	
	// set up data structures
	ConnectRec *nList=new ConnectRec[nnodes];	// Nodal connectivities
	int *theLevel=new int[nnodes];				// list of nodes in a level
	int *lastLevel=new int[nnodes];				// list of nodes in previous level
	bool *levelFlags=new bool[nnodes];			// flags to remember nodes in level
	bool *mapFlags=new bool[nnodes];			// bits to remember nodes that have already been mapped
	int *nodeMap=new int[nnodes];				// map of resequenced nodes
	if(nList==NULL || theLevel==NULL || lastLevel==NULL || levelFlags==NULL
				|| mapFlags==NULL || nodeMap==NULL)
		throw SAXException("Out of memory resequencing the nodes.");
	
	for(i=0;i<nnodes;i++)
	{	nList[i].degree=0;			// zero the degrees
		mapFlags[i]=FALSE;			// clear all flags
	}

	/* Invert element list to nodal connectivity list (note list has node #'s-1)
		nList[0 to nnode-1] says that node is connected to .degree nodes
												 listed in .cons[0] to .degree-1
	*/
	int node1,node2,degree;
	bool addNode;
	for(i=0;i<nelems;i++)
	{	numnds=theElements[i]->NumberNodes();
		for(j=1;j<=numnds;j++)
		{	node1=theElements[i]->NodeIndex(j);
			degree=nList[node1].degree;
			for(k=1;k<=numnds;k++)
			{	if(k==j) continue;
				node2=theElements[i]->NodeIndex(k);
				addNode=TRUE;
				for(l=0;l<degree;l++)
				{	if(node2==nList[node1].cons[l])
					{	addNode=FALSE;
						break;
					}
				}
				if(!addNode) continue;
				if(degree>=MAX_CONNECTIVITY)
				{	char msg[100];
					sprintf(msg,"Mesh too highly connected for resequencing at node %d.",node1+1);
					throw SAXException(msg);
				}
				nList[node1].cons[degree++]=node2;
			}
			nList[node1].degree=degree;
		}
	}
	
	// variables
	int mapped,inLastLevel,inLevel,lognb2;
	
	// Start node map using requested resequence node (-1 from 1 based to zero based here)
	nodeMap[0]=resequence-1;
	mapped=1;
	mapFlags[nodeMap[0]]=TRUE;
	
	// Initialize level structure
	inLastLevel=1;
	lastLevel[0]=nodeMap[0];
	
	// Loop through levels until map is done
	while(TRUE)
	{	// Calculate new level from last level
		for(i=0;i<nnodes;i++) levelFlags[i]=FALSE;
		inLevel=0;
		for(i=0;i<inLastLevel;i++)
		{	node1=lastLevel[i];
			for(j=0;j<nList[node1].degree;j++)
			{	node2=nList[node1].cons[j];
				if(!mapFlags[node2] && !levelFlags[node2])
				{	theLevel[inLevel++]=node2;
					levelFlags[node2]=TRUE;
				}
			}
		}
		
		// If inLevel is zero then all done
		if(inLevel==0) break;
		
		// Sort by degree - shell sort - Numerical Recipes in C, pg 244
		lognb2=(int)(log((double)inLevel)*1.442695022+1.0e-5);	// log base 2
		k=inLevel;
		for(l=1;l<=lognb2;l++)
		{	k>>=1;		// divide by 2
			for(j=k;j<inLevel;j++)
			{	i=j-k;
				degree=nList[theLevel[j]].degree;
				node1=theLevel[j];
				
				// Back up until find insertion point
				while(i>=0 && nList[theLevel[i]].degree>degree)
				{	theLevel[i+k]=theLevel[i];
					i-=k;
				}
				
				// Insert point
				theLevel[i+k]=node1;
			}
		}
		
		// Assign maps to nodes
		for(i=0;i<inLevel;i++)
		{	nodeMap[mapped++]=theLevel[i];
			mapFlags[theLevel[i]]=TRUE;
		}
		
		// Copy level to last level
		inLastLevel=inLevel;
		for(i=0;i<inLevel;i++)
			lastLevel[i]=theLevel[i];
	}
	
	// delete data structures (except nodeMap)
	delete [] nList;
	delete [] theLevel;
	delete [] lastLevel;
	delete [] mapFlags;
	delete [] levelFlags;
	
	// verify found all nodes
	if(mapped!=nnodes)
		throw SAXException("Resequencing failed because of disconnected nodes. Turn resequencing off and plot the mesh to find them,\n  watching out for areas that touch but do not share nodes.");
	
	// reverse the node map
	inLevel=(nnodes-2)/2;
	for(i=0;i<=inLevel;i++)
	{	k=nodeMap[i];
		j=nnodes-1-i;
		nodeMap[i]=nodeMap[j];
		nodeMap[j]=k;
	}
	
	// move original nodes to resequenced nodes
	// Node: revMap 1 based and node numbers in the map are 1 based too
	int *revMap=new int[nnodes+1];
	for(i=1;i<=nnodes;i++)
		revMap[nodeMap[i-1]+1]=i;
	delete [] nodeMap;
	
	// remap nodes
	free(nd);
	theNodes->SetNodeArray(revMap);

	// map nodes in all data structures
	//  ... Elements, Nodal Displacements, Nodal Loads, selected nodes, periodic nodes
	for(i=0;i<nelems;i++) theElements[i]->MapNodes(revMap);
	
	NodalDispBC *nextBC=firstDispBC;
	while(nextBC!=NULL)
		nextBC=nextBC->MapNodes(revMap);
	
	NodalLoad *nextLoad=firstLoadBC;
	while(nextLoad!=NULL)
		nextLoad=nextLoad->MapNodes(revMap);
		
	Constraint *nextConstraint=firstConstraint;
	while(nextConstraint!=NULL)
		nextConstraint=nextConstraint->MapNodes(revMap);
		
	for(i=0;i<(int)fmobj->selectedNodes.size();i++)
		fmobj->selectedNodes[i]=revMap[fmobj->selectedNodes[i]];

	// all done
	delete [] revMap;
}

/********************************************************************************
	 Find periodics nodes if requested
********************************************************************************/

void FEAReadHandler::FindPeriodicNodes(void)
{
	vector< int > xPeriodicNodes;	// pairs of nodes periodic in x
	vector< int > yPeriodicNodes;	// pairs of nodes periodic in y
	
	// should have elements by now
	if(nelems<0) return;
	
	// find range of all nodes
	double xmin=theElements[0]->xmin;
	double xmax=theElements[0]->xmax;
	double ymin=theElements[0]->ymin;
	double ymax=theElements[0]->ymax;
	int i;
	for(i=1;i<nelems;i++)
	{	xmin=min(xmin,theElements[i]->xmin);
		xmax=max(xmax,theElements[i]->xmax);
		ymin=min(ymin,theElements[i]->ymin);
		ymax=max(ymax,theElements[i]->ymax);
	}
	fmobj->periodic.xmin=xmin;
	fmobj->periodic.xmax=xmax;
	fmobj->periodic.ymin=ymin;
	fmobj->periodic.ymax=ymax;
	
	// find all nodes at xmin and/or ymin
	vector< int > xMinNodes;
	vector< int > yMinNodes;
	for(i=1;i<=nnodes;i++)
	{	if(fmobj->periodic.dof&0x01)
		{	if(DbleEqual(xmin,nd[i]->x))
				xMinNodes.push_back(i);
		}
		if(fmobj->periodic.dof&0x02)
		{	if(DbleEqual(ymin,nd[i]->y))
				yMinNodes.push_back(i);
		}
	}
	
	// match to all nodes at xmax and/or ymax
	// Stored as (min1,max1,min2,max2, ...)
	vector< int >::iterator iter;
	for(i=1;i<=nnodes;i++)
	{	if(fmobj->periodic.dof&0x01)
		{	if(DbleEqual(xmax,nd[i]->x))
			{	unsigned int numOrig=xMinNodes.size();
				for(iter=xMinNodes.begin(); iter!=xMinNodes.end(); ++iter)
				{	if(DbleEqual(nd[i]->y,nd[*iter]->y))
					{	xPeriodicNodes.push_back(*iter);
						xPeriodicNodes.push_back(i);
						xMinNodes.erase(iter);
						break;
					}
				}
				if(xMinNodes.size()==numOrig)
					throw SAXException("Found unmatched periodic nodes in x on the right.");
			}
		}
		if(fmobj->periodic.dof&0x02)
		{	if(DbleEqual(ymax,nd[i]->y))
			{	unsigned int numOrig=yMinNodes.size();
				for(iter=yMinNodes.begin(); iter!=yMinNodes.end(); ++iter)
				{	if(DbleEqual(nd[i]->x,nd[*iter]->x))
					{	yPeriodicNodes.push_back(*iter);
						yPeriodicNodes.push_back(i);
						yMinNodes.erase(iter);
						break;
					}
				}
				if(yMinNodes.size()==numOrig)
					throw SAXException("Found unmatched periodic nodes in y (or z) on the top.");
			}
		}
	}
	
	// error if not all were matched
	if(xMinNodes.size()>0)
		throw SAXException("Found unmatched periodic nodes in x on the left.");
	if(yMinNodes.size()>0)
		throw SAXException("Found unmatched periodic nodes in y on the bottom.");
	
	// error if too few (need at least two pairs)
	if(fmobj->periodic.dof&0x01 && xPeriodicNodes.size()<4)
		throw SAXException("Found too few periodic nodes in x (need at least two pairs).");
	if(fmobj->periodic.dof&0x02 && yPeriodicNodes.size()<4)
		throw SAXException("Found too few periodic nodes in y (or z) (need at least two pairs).");
	
	// move min and max in x direction to positions 1 and 2
	if(xPeriodicNodes.size()>0)
	{	double rymin=nd[xPeriodicNodes[0]]->y,rymax=nd[xPeriodicNodes[0]]->y;
		for(i=0;i<(int)xPeriodicNodes.size();i+=2)
		{	rymin=fmin(rymin,nd[xPeriodicNodes[i]]->y);
			rymax=fmax(rymax,nd[xPeriodicNodes[i]]->y);
		}
		
		for(i=0;i<(int)xPeriodicNodes.size();i+=2)
		{	if(DbleEqual(rymin,nd[xPeriodicNodes[i]]->y))
			{	if(i!=0)
				{	int left1=xPeriodicNodes[0];
					int right1=xPeriodicNodes[1];
					xPeriodicNodes[0]=xPeriodicNodes[i];
					xPeriodicNodes[1]=xPeriodicNodes[i+1];
					xPeriodicNodes[i]=left1;
					xPeriodicNodes[i+1]=right1;
				}
			}
			else if(DbleEqual(rymax,nd[xPeriodicNodes[i]]->y))
			{	if(i!=2)
				{	int left2=xPeriodicNodes[2];
					int right2=xPeriodicNodes[3];
					xPeriodicNodes[2]=xPeriodicNodes[i];
					xPeriodicNodes[3]=xPeriodicNodes[i+1];
					xPeriodicNodes[i]=left2;
					xPeriodicNodes[i+1]=right2;
				}
			}
		}
	}
	
	// move min and max in y direction to positions 1 and 2
	if(yPeriodicNodes.size()>0)
	{	double rxmin=nd[yPeriodicNodes[0]]->x,rxmax=nd[yPeriodicNodes[0]]->x;
		for(i=0;i<(int)yPeriodicNodes.size();i+=2)
		{	rxmin=fmin(rxmin,nd[yPeriodicNodes[i]]->x);
			rxmax=fmax(rxmax,nd[yPeriodicNodes[i]]->x);
		}
		
		for(i=0;i<(int)yPeriodicNodes.size();i+=2)
		{	if(DbleEqual(rxmin,nd[yPeriodicNodes[i]]->x))
			{	if(i!=0)
				{	int bottom1=yPeriodicNodes[0];
					int top1=yPeriodicNodes[1];
					yPeriodicNodes[0]=yPeriodicNodes[i];
					yPeriodicNodes[1]=yPeriodicNodes[i+1];
					yPeriodicNodes[i]=bottom1;
					yPeriodicNodes[i+1]=top1;
				}
			}
			else if(DbleEqual(rxmax,nd[yPeriodicNodes[i]]->x))
			{	if(i!=2)
				{	int bottom2=yPeriodicNodes[2];
					int top2=yPeriodicNodes[3];
					yPeriodicNodes[2]=yPeriodicNodes[i];
					yPeriodicNodes[3]=yPeriodicNodes[i+1];
					yPeriodicNodes[i]=bottom2;
					yPeriodicNodes[i+1]=top2;
				}
			}
		}
	}

	// set up constraints
	constraintCtrl=new ConstraintController();
	Constraint *newConstraint;
	
	// periodic in x only (not axisymmetric)
	if(fmobj->periodic.dof==1)
	{	// v dof Ri : vRi = vLi + vR1 - vL1 for i=2 to m
		for(i=3;i<(int)xPeriodicNodes.size();i+=2)
		{	newConstraint=new Constraint(xPeriodicNodes[i],2);
			newConstraint->AddNode(xPeriodicNodes[i-1],(double)1.);
			newConstraint->AddNode(xPeriodicNodes[1],(double)1.);
			newConstraint->AddNode(xPeriodicNodes[0],(double)-1.);
			constraintCtrl->AddObject(newConstraint);
		}
		
		if(fmobj->periodic.fixDu && fmobj->periodic.fixDudy)
		{	// constraints at R1, R2, R3, ... i=1 to m
			double yR1=nd[xPeriodicNodes[1]]->y;
			double ybar=(nd[xPeriodicNodes[3]]->y+yR1)/2.;		// (yR2+yR1)/2
			for(i=1;i<(int)xPeriodicNodes.size();i+=2)
			{	// u dof at Ri : uRi = uLi + du + dudy*(yRi-ybar) 
				newConstraint=new Constraint(xPeriodicNodes[i],1);
				newConstraint->AddNode(xPeriodicNodes[i-1],(double)1.);
				double q=fmobj->periodic.du + fmobj->periodic.dudy*(nd[xPeriodicNodes[i]]->y-ybar);
				newConstraint->AddConstant(q/1000.);
				constraintCtrl->AddObject(newConstraint);
			}
		}
		
		else if(fmobj->periodic.fixDu && !fmobj->periodic.fixDudy)
		{	// constraints at R2, R3, R4, ... i=2 to m
			double yR1=nd[xPeriodicNodes[1]]->y;
			double ybar=(nd[xPeriodicNodes[3]]->y+yR1)/2.;		// (yR2+yR1)/2
			for(i=3;i<(int)xPeriodicNodes.size();i+=2)
			{	// u dof at Ri : uRi = uLi + (1-gi)*(uR1-uL1) + gi*du 
				double gi=(nd[xPeriodicNodes[i]]->y-yR1)/(ybar-yR1);
				newConstraint=new Constraint(xPeriodicNodes[i],1);
				newConstraint->AddNode(xPeriodicNodes[i-1],(double)1.);
				newConstraint->AddNode(xPeriodicNodes[1],(1.-gi));
				newConstraint->AddNode(xPeriodicNodes[0],-(1.-gi));
				newConstraint->AddConstant(gi*fmobj->periodic.du/1000.);
				constraintCtrl->AddObject(newConstraint);
			}
		}
		
		else if(!fmobj->periodic.fixDu && fmobj->periodic.fixDudy)
		{	// constraints at R2, R3, R4, ... i=2 to m
			double yR1=nd[xPeriodicNodes[1]]->y;
			for(i=3;i<(int)xPeriodicNodes.size();i+=2)
			{	// u dof at Ri : uRi = uLi + (uR1-uL1) + dudy*(yRi-yR1)
				newConstraint=new Constraint(xPeriodicNodes[i],1);
				newConstraint->AddNode(xPeriodicNodes[i-1],(double)1.);
				newConstraint->AddNode(xPeriodicNodes[1],(double)1.);
				newConstraint->AddNode(xPeriodicNodes[0],(double)-1.);
				double q=fmobj->periodic.dudy*(nd[xPeriodicNodes[i]]->y-yR1);
				newConstraint->AddConstant(q/1000.);
				constraintCtrl->AddObject(newConstraint);
			}
		}
		
		else
		{	// contraints at R3, R4, ... i=3 to m
			double yR1=nd[xPeriodicNodes[1]]->y;
			double ydiff=nd[xPeriodicNodes[3]]->y-yR1;		// yR2-yR1
			for(i=5;i<(int)xPeriodicNodes.size();i+=2)
			{	// x dof uRi = uLi + (1-fi)*(uR1-uL1) + fi*(uR2-uL2) 
				double fi=(nd[xPeriodicNodes[i]]->y-yR1)/ydiff;
				newConstraint=new Constraint(xPeriodicNodes[i],1);
				newConstraint->AddNode(xPeriodicNodes[i-1],(double)1.);
				newConstraint->AddNode(xPeriodicNodes[1],(1.-fi));
				newConstraint->AddNode(xPeriodicNodes[0],-(1.-fi));
				newConstraint->AddNode(xPeriodicNodes[3],fi);
				newConstraint->AddNode(xPeriodicNodes[2],-fi);
				constraintCtrl->AddObject(newConstraint);
			}
		}
	}
	
	// periodic in y only for planar or z only for axisymmetric problems
	else if(fmobj->periodic.dof==2)
	{	// x dof Ti : vTi = vBi + vT1 - vB1 for i=2 to m
		//            vTi = vBi got i=1 to m if axisymmetric (or 2 to m if rmin on the origin - need to avoid r=0 BC there)
		int ibeg=3;
		if(fmobj->IsAxisymmetric())
		{	if(!DbleEqual(nd[yPeriodicNodes[0]]->x,0.)) ibeg=1;
			fmobj->periodic.fixDvdx=TRUE;
			fmobj->periodic.dvdx=0.;
		}
		for(i=ibeg;i<(int)yPeriodicNodes.size();i+=2)
		{	newConstraint=new Constraint(yPeriodicNodes[i],1);
			newConstraint->AddNode(yPeriodicNodes[i-1],(double)1.);
			if(!fmobj->IsAxisymmetric())
			{	newConstraint->AddNode(yPeriodicNodes[1],(double)1.);
				newConstraint->AddNode(yPeriodicNodes[0],(double)-1.);
			}
			constraintCtrl->AddObject(newConstraint);
		}
		
		if(fmobj->periodic.fixDv && fmobj->periodic.fixDvdx)
		{	// constraints at T1, T2, T3, ... i=1 to m
			double xT1=nd[yPeriodicNodes[1]]->x;
			double xbar=(nd[yPeriodicNodes[3]]->x+xT1)/2.;		// (xT2+xT1)/2
			for(i=1;i<(int)yPeriodicNodes.size();i+=2)
			{	// v dof at Ti : vTi = vBi + dv + dvdx*(xTi-xbar) 
				newConstraint=new Constraint(yPeriodicNodes[i],2);
				newConstraint->AddNode(yPeriodicNodes[i-1],(double)1.);
				double q=fmobj->periodic.dv + fmobj->periodic.dvdx*(nd[yPeriodicNodes[i]]->x-xbar);
				newConstraint->AddConstant(q/1000.);
				constraintCtrl->AddObject(newConstraint);
			}
		}
		
		else if(fmobj->periodic.fixDv && !fmobj->periodic.fixDvdx)
		{	// constraints at T2, T3, T4, ... i=2 to m
			double xT1=nd[yPeriodicNodes[1]]->x;
			double xbar=(nd[yPeriodicNodes[3]]->x+xT1)/2.;		// (xT2+xT1)/2
			for(i=3;i<(int)yPeriodicNodes.size();i+=2)
			{	// v dof at Ti : vTi = vBi + (1-gi)*(vT1-uB1) + gi*dv 
				double gi=(nd[yPeriodicNodes[i]]->x-xT1)/(xbar-xT1);
				newConstraint=new Constraint(yPeriodicNodes[i],2);
				newConstraint->AddNode(yPeriodicNodes[i-1],(double)1.);
				newConstraint->AddNode(yPeriodicNodes[1],(1.-gi));
				newConstraint->AddNode(yPeriodicNodes[0],-(1.-gi));
				newConstraint->AddConstant(gi*fmobj->periodic.dv/1000.);
				constraintCtrl->AddObject(newConstraint);
			}
		}
		
		else if(!fmobj->periodic.fixDv && fmobj->periodic.fixDvdx)
		{	// constraints at T2, T3, T4, ... i=2 to m
			double xT1=nd[yPeriodicNodes[1]]->x;
			for(i=3;i<(int)yPeriodicNodes.size();i+=2)
			{	// v dof at Ti : vTi = vBi + (vT1-vB1) + dvdx*(xTi-xT1)
				newConstraint=new Constraint(yPeriodicNodes[i],2);
				newConstraint->AddNode(yPeriodicNodes[i-1],(double)1.);
				newConstraint->AddNode(yPeriodicNodes[1],(double)1.);
				newConstraint->AddNode(yPeriodicNodes[0],(double)-1.);
				double q=fmobj->periodic.dvdx*(nd[yPeriodicNodes[i]]->x-xT1);
				newConstraint->AddConstant(q/1000.);
				constraintCtrl->AddObject(newConstraint);
			}
		}
		
		else
		{	// contraints at T3, RT4, ... i=3 to m
			double xT1=nd[yPeriodicNodes[1]]->x;
			double xdiff=nd[yPeriodicNodes[3]]->x-xT1;		// xT2-xT1
			for(i=5;i<(int)yPeriodicNodes.size();i+=2)
			{	// y dof vTi = vBi + (1-fi)*(vT1-vB1) + fi*(vT2-vB2) 
				double fi=(nd[yPeriodicNodes[i]]->x-xT1)/xdiff;
				newConstraint=new Constraint(yPeriodicNodes[i],2);
				newConstraint->AddNode(yPeriodicNodes[i-1],(double)1.);
				newConstraint->AddNode(yPeriodicNodes[1],(1.-fi));
				newConstraint->AddNode(yPeriodicNodes[0],-(1.-fi));
				newConstraint->AddNode(yPeriodicNodes[3],fi);
				newConstraint->AddNode(yPeriodicNodes[2],-fi);
				constraintCtrl->AddObject(newConstraint);
			}
		}
	}
	
	// periodic in both x and y (not axisymmetric)
	else
	{	// u dof Ri : uRi = uLi + uR1 - uL1 for i=2 to m
		//          : uRi = uLi + du        for i=1 to m if fixed du
		int ibeg = fmobj->periodic.fixDu ? 1 : 3;
		for(i=ibeg;i<(int)xPeriodicNodes.size();i+=2)
		{	newConstraint=new Constraint(xPeriodicNodes[i],1);
			newConstraint->AddNode(xPeriodicNodes[i-1],(double)1.);
			if(fmobj->periodic.fixDu)
				newConstraint->AddConstant(fmobj->periodic.du/1000.);
			else
			{	newConstraint->AddNode(xPeriodicNodes[1],(double)1.);
				newConstraint->AddNode(xPeriodicNodes[0],(double)-1.);
			}
			constraintCtrl->AddObject(newConstraint);
		}
	
		// v dof Ri : vRi = vLi + vR1 - vL1 for i=2 to m
		//          : vRi = vLi + dvdx      for i=1 to m if fixed dvdx
		ibeg = fmobj->periodic.fixDvdx ? 1 : 3;
		for(i=ibeg;i<(int)xPeriodicNodes.size();i+=2)
		{	newConstraint=new Constraint(xPeriodicNodes[i],2);
			newConstraint->AddNode(xPeriodicNodes[i-1],(double)1.);
			if(fmobj->periodic.fixDvdx)
				newConstraint->AddConstant(fmobj->periodic.dvdx/1000.);
			else
			{	newConstraint->AddNode(xPeriodicNodes[1],(double)1.);
				newConstraint->AddNode(xPeriodicNodes[0],(double)-1.);
			}
			constraintCtrl->AddObject(newConstraint);
		}
		
		// u dof Ti : uTi = uBi + uT1 - uB1 for i=2 to m
		//          : uTi = uBi + dudy      for i=1 to m if fixed dudy
		ibeg = fmobj->periodic.fixDudy ? 1 : 3;
		for(i=ibeg;i<(int)yPeriodicNodes.size();i+=2)
		{	newConstraint=new Constraint(yPeriodicNodes[i],1);
			newConstraint->AddNode(yPeriodicNodes[i-1],(double)1.);
			if(fmobj->periodic.fixDudy)
				newConstraint->AddConstant(fmobj->periodic.dudy/1000.);
			else
			{	newConstraint->AddNode(yPeriodicNodes[1],(double)1.);
				newConstraint->AddNode(yPeriodicNodes[0],(double)-1.);
			}
			constraintCtrl->AddObject(newConstraint);
		}
	
		// v dof Ti : vTi = vBi + vT1 - vB1 for i=2 to m
		//          : vTi = vBi + dv        for i=1 to m if fixed dv
		ibeg = fmobj->periodic.fixDv ? 1 : 3;
		for(i=ibeg;i<(int)yPeriodicNodes.size();i+=2)
		{	newConstraint=new Constraint(yPeriodicNodes[i],2);
			newConstraint->AddNode(yPeriodicNodes[i-1],(double)1.);
			if(fmobj->periodic.fixDv)
				newConstraint->AddConstant(fmobj->periodic.dv/1000.);
			else
			{	newConstraint->AddNode(yPeriodicNodes[1],(double)1.);
				newConstraint->AddNode(yPeriodicNodes[0],(double)-1.);
			}
			constraintCtrl->AddObject(newConstraint);
		}
	}
	
	firstConstraint=(Constraint *)constraintCtrl->firstObject;
	delete constraintCtrl;
}

/********************************************************************************
	 Remove empty elements (material=0) amnd their nodes
	 Also verify other elements have valid material type
********************************************************************************/

void FEAReadHandler::RemoveEmptyElements(void)
{
	int nmatl=matCtrl->numObjects;
	
	int j,i=0;
	bool emptyElements=FALSE;
	while(i<nelems)
	{	int matID=theElements[i]->material;
	
		// invalid materials
		if(matID<0 || matID>nmatl)
			throw SAXException("An element has an unknown material type.");
		
		// if valid material continue
		if(matID!=0)
		{	i++;
			continue;
		}
		
		// found empty element that needs to be removed, which may trigger removing nodes
		// and boundary conditions
		emptyElements=TRUE;
		delete theElements[i];
		
		// compact element list
		nelems--;
		for(j=i;j<nelems;j++)
		{	theElements[j]=theElements[j+1];
			theElements[j]->num--;
		}
	}
		
	// see if any nodes are now not in any element
	if(emptyElements)
	{	int nodeNum;
		i=1;
		while(i<=nnodes)
		{	nodeNum=nd[i]->num;
			
			// search remaining elements, if found then continue
			for(j=0;j<nelems;j++)
			{	if(theElements[j]->HasNode(nodeNum)) break;
			}
			if(j<nelems)
			{	i++;
				continue;
			}
			
			// found dangling node
			delete nd[i];
			nnodes--;
			for(j=i;j<=nnodes;j++)
			{	nd[j]=nd[j+1];
				nd[j]->num--;
			}
			
			// theNodes nodes controller may still need linkages for resequencing
			if(i==1)
				theNodes->firstObject=nd[1];
			else if(i>nnodes)
			{	theNodes->lastObject=nd[nnodes];
				nd[nnodes]->SetNextObject(NULL);
			}
			else
				nd[i-1]->SetNextObject(nd[i]);
			
			// renumber nodes in elements
			for(j=0;j<nelems;j++)
			{	theElements[j]->DecrementNodeNums(nodeNum);
			}
			
			// nodal displacement BCs
			NodalDispBC *nextBC=firstDispBC;
			NodalDispBC *prevBC=NULL;
			NodalDispBC *tempBC=NULL;
			bool deleted;
			while(nextBC!=NULL)
			{	tempBC=(NodalDispBC *)nextBC->DecrementNodeNum(nodeNum,prevBC,&deleted);
				if(deleted)
				{	if(prevBC==NULL)
					{	// the first one was deleted, reset it
						firstDispBC=tempBC;
					}
				}
				else
					prevBC=nextBC;
				nextBC=tempBC;
			}
			
			// nodal load BCs
			NodalLoad *nextLoadBC=firstLoadBC;
			NodalLoad *prevLoadBC=NULL;
			NodalLoad *tempLoadBC=NULL;
			while(nextLoadBC!=NULL)
			{	tempLoadBC=(NodalLoad *)nextLoadBC->DecrementNodeNum(nodeNum,prevLoadBC,&deleted);
				if(deleted)
				{	if(prevLoadBC==NULL)
					{	// the first one was deleted, reset it
						firstLoadBC=tempLoadBC;
					}
				}
				else
					prevLoadBC=nextLoadBC;
				nextLoadBC=tempLoadBC;
			}
			
			// selected nodes
			for(j=0;j<(int)fmobj->selectedNodes.size();j++)
			{	if(fmobj->selectedNodes[j]>nodeNum)
					fmobj->selectedNodes[j]--;
				else if(fmobj->selectedNodes[j]==nodeNum)
					fmobj->selectedNodes[j]=0;
			}
					
		}
		
		// verify the mesh is well connected
	}

}
