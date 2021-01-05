/********************************************************************************
    Area.cpp
    NairnFEA
    
    Created by John Nairn on 9/28/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_FEA/Area.hpp"
#include "Read_FEA/PathsController.hpp"
#include "Read_XML/ElementsController.hpp"
#include "Read_XML/NodesController.hpp"
#include "Read_FEA/KeypointsController.hpp"
#include "Read_FEA/Path.hpp"
#include "Read_FEA/Keypoint.hpp"
#include "Elements/ElementBase.hpp"
#include "Elements/EightNodeIsoparam.hpp"
#include "Read_XML/Expression.hpp"

Area *theArea=NULL;

/********************************************************************************
	Area: Constructors and Destructor
********************************************************************************/

// Initialize, but assumes name already 32 characters or less and unique
Area::Area(int matNum,char *matAngle,double thickness)
{
	mat=matNum;
	angleExpr=matAngle;
	angle=0.;
	thick=thickness;
	numPaths=0;
	areaNode=NULL;
	areaElem=NULL;
}

// Destructor
Area::~Area()
{
	if(areaNode!=NULL)
		delete [] areaNode;
	if(areaElem!=NULL)
		delete areaElem;
	if(angleExpr!=NULL)
		delete [] angleExpr;
}

/********************************************************************************
	Area: methods
********************************************************************************/

// add key point to the last path
int Area::AddPath(char *pathName)
{   
	// someday may allow more general areas
	if(numPaths>=4) return FALSE;
	
	// find the path object
	edges[numPaths]=paths->FindPath(pathName);
	if(edges[numPaths]==NULL) return FALSE;
	
	// make sure path not already in the area
	int i;
	for(i=0;i<numPaths;i++)
	{	if(edges[i]==edges[numPaths])
			return FALSE;
	}
	
	// make sure it connects to previous path and reorient if needed
	if(numPaths>0)
	{	if(edges[numPaths-1]->LastKeypoint()!=edges[numPaths]->FirstKeypoint())
		{	if(edges[numPaths-1]->LastKeypoint()==edges[numPaths]->LastKeypoint())
				edges[numPaths]->ReorientPath();
			else if(numPaths==1)
			{	if(edges[0]->FirstKeypoint()==edges[1]->FirstKeypoint())
					edges[0]->ReorientPath();
				else if(edges[0]->FirstKeypoint()==edges[1]->LastKeypoint())
				{	edges[0]->ReorientPath();
					edges[1]->ReorientPath();
				}
				else
				{	/* could be second path of interface elements
						If there are eventually valid, they will have already
						been meshed into areas and thus should be oriented
						in opposite directions (first as last keys should be the
						same) and should be reoriented
					*/
					if(!theElems->InterfaceElements()) return FALSE;
					if(edges[0]->FirstKeypoint()->SamePtAs(edges[1]->LastKeypoint()))
					{	edges[0]->ReorientPath();
						edges[1]->ReorientPath();
					}
					else
						return FALSE;
				}
			}
			else
				return FALSE;
		}
	}
	
	// all done
	numPaths++;
	return TRUE;
}

// Mesh the area into nodes and elements or return error message if fails
// assumes edges already oriented correctly
const char *Area::MeshElements(void)
{
	// exactly 2 paths will be for an interface
	if(numPaths==2)
		return MeshInterface();
		
	// Must be 4-sided area. In future should support more
	if(numPaths!=4)
		return "Code can only mesh quadrilateral areas (4 paths) or interfaces (2 paths).";
	
	// cannot be interface elements
	if(theElems->InterfaceElements())
		return "Areas cannot be meshed into interface elements.";
	
	// must be connected
	if(edges[numPaths-1]->LastKeypoint()!=edges[0]->FirstKeypoint())
		return "Area does not define enclosed area.";
	
	// is element provided
	if(theElems->CurrentElemID()<0)
		return "No element type defined for area.";
	
	// Check intervals
	if(!CheckIntervals())
		return "Number of intervals on sides of the area are not valid";
	
	// Check path usage
	if(!PathsAvailable())
		return "Paths cannot be meshed into more than 2 areas";
	
	// Check ccw
	if(SignedArea()<0.)
		return "The paths must circumnavigate the area in a counter-clockwise direction";

	// create function if being used
	if(angleExpr!=NULL)
	{	if(!Expression::CreateFunction(angleExpr,1))
			return "The expression for material angle is not a valid function";
	}
	
	// mesh the error
	return MeshArea();
}

// mesh two overlapping paths into interface elements
const char *Area::MeshInterface(void)
{
	// must be interface elements
	if(!theElems->InterfaceElements())
		return "Interfaces can only be meshed into interface elements.";
		
	// must mesh into area before an interface
	if(edges[0]->face==0 || edges[1]->face==0)
		return "An interface elements path has not yet be meshed into elements.";
	
	// can not interface an interior path or one already with an interface
	if(edges[0]->face<0 || edges[1]->face<0)
		return "An interface elements is not on a boundary (i.e., it is an interior path).";

	// paths are already oriented in same direction (not ccw around interface) - verify more
	
	// same intervals
	if(edges[0]->intervals!=edges[1]->intervals)
		return "Two interface paths have different numbers of nodes.";
	
	// same ratio
	if(!DbleEqual(edges[0]->ratio,1./edges[1]->ratio))
		return "The nodes on two interface paths are at different locations.";
		
	// set faces to -face  to indicate path used twice - Area and an Interface
	edges[0]->face=-edges[0]->face;
	edges[1]->face=-edges[1]->face;
				
	// add the interface elements - 1 for each path interval
	int i,numInt=edges[0]->intervals,rInt;
	int eNode[MaxElNd];
	int lnameEl=theElems->CurrentElemID();
	for(i=1;i<=numInt;i++)
	{	// nodes 1 and 2 or 1, 2, and 3 on first path
		if(lnameEl==LINEAR_INTERFACE)
		{	eNode[1]=edges[0]->nodeAtIndex(i);
			eNode[2]=edges[0]->nodeAtIndex(i+1);
		}
		else
		{	eNode[1]=edges[0]->nodeAtIndex(2*i-1);
			eNode[2]=edges[0]->nodeAtIndex(2*i);
			eNode[3]=edges[0]->nodeAtIndex(2*i+1);
		}

		// nodes 3 and 4 or 4, 5, and 6 on second path (in reverse order) */
		rInt=numInt-i+1;
		if(lnameEl==LINEAR_INTERFACE)
		{	eNode[3]=edges[1]->nodeAtIndex(rInt);
			eNode[4]=edges[1]->nodeAtIndex(rInt+1);
		}
		else
		{	eNode[4]=edges[1]->nodeAtIndex(2*rInt-1);
			eNode[5]=edges[1]->nodeAtIndex(2*rInt);
			eNode[6]=edges[1]->nodeAtIndex(2*rInt+1);
		}
		
		// add the element
		if(!theElems->MeshElement(eNode,mat,(double)0.,thick))
			return "Unable to create elements in the interface area";
	}
	
	return NULL;
}

// Mesh an area that has been checked already as valid and is quadrilateral
// throws std::bad_alloc
const char *Area::MeshArea(void)
{
	int i;
	int eNode[MaxElNd];
	
	// space for area nodes and element object
	areaNode=new Vector[2*numPaths];
	for(i=0;i<8;i++) eNode[i]=i+1;
	areaElem=new EightNodeIsoparam(1,eNode,1,(double)0.0,(double)0.0);
	
	// put nodes at the corners
	for(i=0;i<numPaths;i++)
	{	edges[i]->FirstKeypointToNode();
		edges[i]->GetFirstAndMidXY(&areaNode[i],&areaNode[i+numPaths]);
	}
		
	// place nodes along the edges
	for(i=0;i<numPaths;i++)
	{	if(!edges[i]->IsMeshed())
		{	// Get line in dimensionless coordinates
			double x1=1.,y1=1.,x2=1.,y2=1.;
			switch(i)
			{	case 0:
					x1=y1=y2=-1;
					break;
				case 1:
					y1=-1;
					break;
				case 2:
					x2=-1;
					break;
				case 3:	
					x1=x2=y2=-1;
					break;
				default:
					break;
			}
			
			// Do main nodes (note: 1 interval is special case)
			if(edges[i]->intervals==1)
			{	// Only 1 interval - no need to get nodes
				edges[i]->firstMainNode=0;
			}
			else
			{	edges[i]->firstMainNode=theNodes->NextNodeNumber();
				MeshLine(x1,y1,x2,y2,edges[i]->ratio,edges[i]->intervals,FALSE);
			}
			
			/* Do mid side nodes */
			if(theElems->HasMidsideNodes())
			{	edges[i]->firstMidsideNode=theNodes->NextNodeNumber();
				MeshLine(x1,y1,x2,y2,edges[i]->ratio,edges[i]->intervals,TRUE);
			}
		}
	}
	
	/*-------------------------------------------------------------------
	 	Use the following section when the numbers of elements are not
	 		that same on opposite sides
		Divide element as follows:
		 
			   _____na__________
			  |\               |
			  |  \+5           |+3      na, nb, nc are numbers of intervals
			  |    \ nc        5 nc                
			  |      \         |        1. create new keypoints (#'s on lines)
		  nb  |        6___+6__3        2. create new paths (+numbers
			  |  1st    |  na  |+2      3. Call MeshArea recursively for the
			  |   nb  +4|      4 nb          newly created areas with equal
			  |   +0    |   +1 |             interval numbers on opposites sides
			  -----1----0---2--
				  nc         na
	-------------------------------------------------------------------*/
	if(edges[0]->intervals!=edges[2]->intervals)
	{	// area intervals (have been ordered such that n1>n3, which implies n2>n4)
		int n1=edges[0]->intervals;
		int n3=edges[2]->intervals;
		int n4=edges[3]->intervals;
		int nc=n1-n3;
		int na=n3;
		int nb=n4;
		
		// keypoint at junction of two new paths along path 1
		double targ,t1arg,x1arg,x2arg,y1arg,y2arg,xkey,ykey;
		Keypoint *addedKeys[7];
		FindPtOnLine((double)-1.,(double)-1.,(double)1.,(double)-1.,edges[0]->ratio,
				edges[0]->intervals,nc,&targ,&x1arg,&x2arg,&y1arg,&y2arg);
		t1arg=targ;
		LineLocate(targ,x1arg,x2arg,y1arg,y2arg,&xkey,&ykey);
		addedKeys[0]=new Keypoint("",xkey,ykey);
		addedKeys[0]->node=edges[0]->firstMainNode+(nc-1)*edges[0]->nodeIncrement;
		
		// Create keypoint in middle of added path 0
		LineLocate((targ-1.)/2.,x1arg,x2arg,y1arg,y2arg,&xkey,&ykey);
		addedKeys[1]=new Keypoint("",xkey,ykey);
		
		// Create keypoint in middle of added path 1 */
		LineLocate((1.+targ)/2.,x1arg,x2arg,y1arg,y2arg,&xkey,&ykey);
		addedKeys[2]=new Keypoint("",xkey,ykey);
		
		// keypoint at junction of two new paths along path 2
		FindPtOnLine((double)1.,(double)-1.,(double)1.,(double)1.,edges[1]->ratio,
				edges[1]->intervals,nb,&targ,&x1arg,&x2arg,&y1arg,&y2arg);
		LineLocate(targ,x1arg,x2arg,y1arg,y2arg,&xkey,&ykey);
		addedKeys[3]=new Keypoint("",xkey,ykey);
		addedKeys[3]->node=edges[1]->firstMainNode+(nb-1)*edges[1]->nodeIncrement;
		
		// Create keypoint in middle of added path 2
		LineLocate((targ-1.)/2.,x1arg,x2arg,y1arg,y2arg,&xkey,&ykey);
		addedKeys[4]=new Keypoint("",xkey,ykey);
		
		// Create keypoint in middle of added path 3
		LineLocate((1.+targ)/2.,x1arg,x2arg,y1arg,y2arg,&xkey,&ykey);
		addedKeys[5]=new Keypoint("",xkey,ykey);
		
		// keypoint at middle of area
		FindPtOnLine(t1arg,(double)-1.,t1arg,(double)1.,edges[1]->ratio,
				edges[1]->intervals,nb,&targ,&x1arg,&x2arg,&y1arg,&y2arg);
		LineLocate(targ,x1arg,x2arg,y1arg,y2arg,&xkey,&ykey);
		addedKeys[6]=new Keypoint("",xkey,ykey);
		
		// First new path along path 1
		Path *addedPaths[7];
		double r,a1,a4;
		GetLineParameters(edges[0]->intervals,edges[0]->ratio,&r,&a1);
		addedPaths[0]=new Path("",nc,GetRatio(nc,r));
		addedPaths[0]->SetKeys(edges[0]->FirstKeypoint(),addedKeys[1],addedKeys[0]);
		addedPaths[0]->firstMainNode=edges[0]->firstMainNode;
		addedPaths[0]->firstMidsideNode=edges[0]->firstMidsideNode;
		addedPaths[0]->nodeIncrement=edges[0]->nodeIncrement;
		edges[0]->subPath[0]=addedPaths[0];
		
		// Second new path along path 1
		addedPaths[1]=new Path("",na,GetRatio(na,r));
		addedPaths[1]->SetKeys(addedKeys[0],addedKeys[2],edges[0]->LastKeypoint());
		addedPaths[1]->firstMainNode=edges[0]->firstMainNode+nc*edges[0]->nodeIncrement;
		addedPaths[1]->firstMidsideNode=edges[0]->firstMidsideNode+nc*edges[0]->nodeIncrement;
		addedPaths[1]->nodeIncrement=edges[0]->nodeIncrement;
		edges[0]->subPath[1]=addedPaths[1];
		
		// First new path along path 2
		GetLineParameters(edges[1]->intervals,edges[1]->ratio,&r,&a4);
		addedPaths[2]=new Path("",nb,GetRatio(nb,r));
		addedPaths[2]->SetKeys(edges[1]->FirstKeypoint(),addedKeys[4],addedKeys[3]);
		addedPaths[2]->firstMainNode=edges[1]->firstMainNode;
		addedPaths[2]->firstMidsideNode=edges[1]->firstMidsideNode;
		addedPaths[2]->nodeIncrement=edges[1]->nodeIncrement;
		edges[1]->subPath[0]=addedPaths[2];

		// Second new path along path 2
		addedPaths[3]=new Path("",nc,GetRatio(nc,r));
		addedPaths[3]->SetKeys(addedKeys[3],addedKeys[5],edges[1]->LastKeypoint());
		addedPaths[3]->firstMainNode=edges[1]->firstMainNode+nb*edges[1]->nodeIncrement;
		addedPaths[3]->firstMidsideNode=edges[1]->firstMidsideNode+nb*edges[1]->nodeIncrement;
		addedPaths[3]->nodeIncrement=edges[1]->nodeIncrement;
		edges[1]->subPath[1]=addedPaths[3];
		
		// New Path from path 1 to central point
		addedPaths[4]=new Path("",nb,addedPaths[2]->ratio);
		addedPaths[4]->SetKeys(addedKeys[0],addedKeys[6],NULL);
		
		// New Path from central point to upper-left corner
		GetLineParameters(addedPaths[0]->intervals,addedPaths[0]->ratio,&r,&a1);
		GetLineParameters(addedPaths[3]->intervals,addedPaths[3]->ratio,&r,&a4);
		x2arg=addedKeys[0]->x-(edges[0]->FirstKeypoint())->x;
		y2arg=addedKeys[0]->y-(edges[0]->FirstKeypoint())->y;
		double p1d1=a1*a1*(x2arg*x2arg+y2arg*y2arg);
		double p1d2=p1d1/(addedPaths[0]->ratio*addedPaths[0]->ratio);
		x2arg=(edges[1]->LastKeypoint())->x-addedKeys[3]->x;
		y2arg=(edges[1]->LastKeypoint())->y-addedKeys[3]->y;
		double p4d1=a4*a4*(x2arg*x2arg+y2arg*y2arg);
		double p4d2=p4d1/(addedPaths[3]->ratio*addedPaths[3]->ratio);
		double ratio=sqrt(p1d2 + p4d1)/sqrt(p1d1 + p4d2);
		addedPaths[5]=new Path("",nc,ratio);
		addedPaths[5]->SetKeys(addedKeys[6],edges[3]->FirstKeypoint(),NULL);

		// New Path from path 2 to central point
		addedPaths[6]=new Path("",na,1./addedPaths[1]->ratio);
		addedPaths[6]->SetKeys(addedKeys[3],addedKeys[6],NULL);
		
		// clear and store items
		delete [] areaNode;
		delete areaElem;
		areaNode=NULL;
		areaElem=NULL;
		Path *oldEdges[4];
		for(i=0;i<=3;i++) oldEdges[i]=edges[i];
		
		// Call MeshArea three times (do not flip 1st or 3rd areas)
		bool oldFlip=theElems->FlipTriangles();
		theElems->SetFlipTriangles("0");
		edges[0]=addedPaths[0];
		edges[1]=addedPaths[4];
		edges[2]=addedPaths[5];
		edges[3]=oldEdges[3];
		const char *msg=MeshArea();
		if(msg!=NULL) return msg;
		
		theElems->SetFlipTriangles(oldFlip);
		edges[0]=addedPaths[1];
		edges[1]=addedPaths[2];
		edges[2]=addedPaths[6];
		edges[3]=addedPaths[4];
		edges[3]->ReorientPath();
		msg=MeshArea();
		if(msg!=NULL) return msg;

		theElems->SetFlipTriangles("0");
		edges[0]=addedPaths[6];
		edges[1]=addedPaths[3];
		edges[2]=oldEdges[2];
		edges[3]=addedPaths[5];
		edges[0]->ReorientPath();
		edges[3]->ReorientPath();
		msg=MeshArea();
		if(msg!=NULL) return msg;
		
		theElems->SetFlipTriangles(oldFlip);		// restore
		// add to list so they get deleted when no longer needed
		for(i=0;i<=6;i++)
		{	keyPts->AddObject(addedKeys[i]);
			paths->AddObject(addedPaths[i]);
		}
		return NULL;
	}

	/* Mesh nodes that are interior to the area
	*/
	double rLeft,aLeft,rRight,aRight;
	int firstInterior=0,firstInteriorMidside=0;
	double yLeft,yRight,yLLast,yRLast;
	double ratio,botFraction,yLeftLast,yRightLast,yLeftMid,yRightMid;
	
	// Linear interpolation of ratios
	double ratioTop=1./edges[2]->ratio;
	double ratioBot=edges[0]->ratio;
	
	// Define interval spacing along left and right edges
	GetLineParameters(edges[3]->intervals,1/edges[3]->ratio,&rLeft,&aLeft);
	GetLineParameters(edges[1]->intervals,edges[1]->ratio,&rRight,&aRight);
	
	// Do main nodes interior to area (if any are there) */
	if(edges[3]->intervals>1 && edges[0]->intervals>1)
	{	firstInterior=theNodes->NextNodeNumber();
		
		// Loop over intervals
		yLeft=yRight=-1.;
		yLLast=aLeft;
		yRLast=aRight;
		for(i=1;i<edges[3]->intervals;i++)
		{	yLeft+=yLLast;
			yLLast*=rLeft;
			yRight+=yRLast;
			yRLast*=rRight;
			botFraction=(1.-(yLeft+yRight)/2.)/2.;
			ratio=botFraction*ratioBot+(1-botFraction)*ratioTop;
			MeshLine((double)-1.,yLeft,(double)1.,yRight,ratio,edges[0]->intervals,FALSE);
		}
	}

	// Do mid side nodes (if needed)
	if(theElems->HasMidsideNodes() && (edges[0]->intervals>1 || edges[3]->intervals>1))
	{	firstInteriorMidside=theNodes->NextNodeNumber();
		
		// Loop over intervals
		yLeft=yRight=yLeftLast=yRightLast=-1.;
		yLLast=aLeft;
		yRLast=aRight;
		for(i=1;i<=edges[3]->intervals;i++)
		{	yLeft+=yLLast;
			yLLast*=rLeft;
			yRight+=yRLast;
			yRLast*=rRight;
			
			// Do mid side nodes on right edges of elements
			if(edges[0]->intervals>1)
			{	yLeftMid=(yLeftLast+yLeft)/2.;
				yRightMid=(yRightLast+yRight)/2.;
				botFraction=(1.-(yLeftMid+yRightMid)/2.)/2.;
				ratio=botFraction*ratioBot+(1-botFraction)*ratioTop;
				MeshLine((double)-1.,yLeftMid,(double)1.,yRightMid,ratio,edges[0]->intervals,FALSE);
			}
			
			// Do mid side nodes to tops of elements
			if(i!=edges[3]->intervals)
			{	botFraction=(1.-(yLeft+yRight)/2.)/2.;
				ratio=botFraction*ratioBot+(1-botFraction)*ratioTop;
				MeshLine((double)-1.,yLeft,(double)1.,yRight,ratio,edges[0]->intervals,TRUE);
			}
			
			// Save current point
			yLeftLast=yLeft;
			yRightLast=yRight;
		}
	}
	
	/* Mesh the area or get nodes in each element
	*/
	int row,rows=edges[1]->intervals;
	int col,cols=edges[0]->intervals;
	
	/* Calculate element numbers along surronding paths -
		First element means first one on boundary */
	int elemNum=theElems->numObjects;
	if(theElems->ElementSides()==3)
	{	edges[0]->firstElem=elemNum+1;
		edges[0]->elemIncrement=2;
		if(theElems->FlipTriangles())
			edges[1]->firstElem=elemNum+2*cols;
		else
			edges[1]->firstElem=elemNum+2*cols-1;
		edges[1]->elemIncrement=2*cols;
		edges[2]->firstElem=elemNum+2*rows*cols;
		edges[2]->elemIncrement=-2;
		if(theElems->FlipTriangles())
			edges[3]->firstElem=elemNum+2*(rows-1)*cols+1;
		else
			edges[3]->firstElem=elemNum+2*(rows-1)*cols+2;
		edges[3]->elemIncrement=-2*cols;
	}
	else if(theElems->ElementSides()==4)
	{	edges[0]->firstElem=elemNum+1;
		edges[0]->elemIncrement=1;
		edges[1]->firstElem=elemNum+cols;
		edges[1]->elemIncrement=cols;
		edges[2]->firstElem=elemNum+rows*cols;
		edges[2]->elemIncrement=-1;
		edges[3]->firstElem=elemNum+(rows-1)*cols+1;
		edges[3]->elemIncrement=-cols;
	}
	
	// Calculate the elements
	for(row=1;row<=rows;row++)
	{	for(col=1;col<=cols;col++)
		{	/* Assume all nodes are interior to area and get
				the nodes - wrong ones will be corrected latter */
			eNode[3]=firstInterior+(row-1)*(cols-1)+(col-1);
			eNode[4]=eNode[3]-1;
			eNode[1]=eNode[4]-(cols-1);
			eNode[2]=eNode[1]+1;
			if(theElems->HasMidsideNodes())
			{	eNode[6]=firstInteriorMidside+(row-1)*(2*cols-1)+(col-1);
				eNode[7]=eNode[6]+cols-1;
				eNode[8]=eNode[6]-1;
				eNode[5]=eNode[6]-cols;
			}
				
			// Correct nodes that are actually on an edge for first row
			if(row==1)
			{	// First assume rows and cols>1
				if(col==1)
				{	eNode[1]=(edges[0]->FirstKeypoint())->node;
					eNode[2]=edges[0]->firstMainNode;
					eNode[4]=edges[3]->firstMainNode+edges[3]->nodeIncrement*(rows-2);
					if(theElems->HasMidsideNodes())
					{	eNode[5]=edges[0]->firstMidsideNode;
						eNode[8]=edges[3]->firstMidsideNode+edges[3]->nodeIncrement*(rows-1);
					}
				}
				else if(col==cols)
				{	eNode[1]=edges[0]->firstMainNode+edges[0]->nodeIncrement*(cols-2);
					eNode[2]=(edges[0]->LastKeypoint())->node;
					eNode[3]=edges[1]->firstMainNode;
					if(theElems->HasMidsideNodes())
					{	eNode[5]=edges[0]->firstMidsideNode+edges[0]->nodeIncrement*(cols-1);
						eNode[6]=edges[1]->firstMidsideNode;
					}
				}
				else
				{	eNode[1]=edges[0]->firstMainNode+edges[0]->nodeIncrement*(col-2);
					eNode[2]=eNode[1]+edges[0]->nodeIncrement;
					if(theElems->HasMidsideNodes())
						eNode[5]=edges[0]->firstMidsideNode+edges[0]->nodeIncrement*(col-1);
				}
				
				/* Final correction for rows or cols equal to 1
					(Note: both will not be 1 (trapped earlier)) */
				if(rows==1)
				{	if(col==1)	// Note: cols won't be 1
					{	eNode[3]=edges[2]->firstMainNode+edges[2]->nodeIncrement*(cols-2);
						eNode[4]=(edges[3]->FirstKeypoint())->node;
						if(theElems->HasMidsideNodes())
							eNode[7]=edges[2]->firstMidsideNode+edges[2]->nodeIncrement*(cols-1);
					}
					else if(col==cols)
					{	eNode[3]=(edges[2]->FirstKeypoint())->node;
						eNode[4]=edges[2]->firstMainNode;
						if(theElems->HasMidsideNodes())
							eNode[7]=edges[2]->firstMidsideNode;
					}
					else
					{	eNode[3]=edges[2]->firstMainNode+edges[2]->nodeIncrement*(cols-col-1);
						eNode[4]=eNode[3]+edges[2]->nodeIncrement;
						if(theElems->HasMidsideNodes())
							eNode[7]=edges[2]->firstMidsideNode+edges[2]->nodeIncrement*(cols-col);
					}
				}
				else if(cols==1)	// Note rows can not also be 1
				{	eNode[2]=(edges[1]->FirstKeypoint())->node;
					eNode[3]=edges[1]->firstMainNode;
					if(theElems->HasMidsideNodes())
						eNode[6]=edges[1]->firstMidsideNode;
				}
			}
			
			// Correct nodes that are actually on an edge for last row
			else if(row==rows)
			{	if(col==1)
				{	eNode[1]=edges[3]->firstMainNode;
					eNode[3]=edges[2]->firstMainNode+edges[2]->nodeIncrement*(cols-2);
					eNode[4]=(edges[2]->LastKeypoint())->node;
					if(theElems->HasMidsideNodes())
					{	eNode[7]=edges[2]->firstMidsideNode+edges[2]->nodeIncrement*(cols-1);
						eNode[8]=edges[3]->firstMidsideNode;
					}
				}
				else if(col==cols)
				{	eNode[2]=edges[1]->firstMainNode+edges[1]->nodeIncrement*(rows-2);
					eNode[3]=(edges[1]->LastKeypoint())->node;
					eNode[4]=edges[2]->firstMainNode;
					if(theElems->HasMidsideNodes())
					{	eNode[6]=edges[1]->firstMidsideNode+edges[1]->nodeIncrement*(rows-1);
						eNode[7]=edges[2]->firstMidsideNode;
					}
				}
				else
				{	eNode[3]=edges[2]->firstMainNode+edges[2]->nodeIncrement*(cols-col-1);
					eNode[4]=eNode[3]+edges[2]->nodeIncrement;
					if(theElems->HasMidsideNodes())
						eNode[7]=edges[2]->firstMidsideNode+edges[2]->nodeIncrement*(cols-col);
				}
				
				/* Final correction cols equal to 1
					(Note: rows won't be 1 - it will come through in
						first in if(row==1) section */
				if(cols==1)
				{	eNode[3]=(edges[2]->FirstKeypoint())->node;
					eNode[2]=edges[1]->firstMainNode+edges[1]->nodeIncrement*(rows-2);
					if(theElems->HasMidsideNodes())
						eNode[6]=edges[1]->firstMidsideNode+edges[1]->nodeIncrement*(rows-1);
				}
			}
			
			// Correct nodes that are actually on an edge for interior rows
			else
			{	if(col==1)
				{	eNode[4]=edges[3]->firstMainNode+edges[3]->nodeIncrement*(rows-row-1);
					eNode[1]=eNode[4]+edges[3]->nodeIncrement;
					if(theElems->HasMidsideNodes())
						eNode[8]=edges[3]->firstMidsideNode+edges[3]->nodeIncrement*(rows-row);
				}
				
				if(col==cols)
				{	eNode[2]=edges[1]->firstMainNode+edges[1]->nodeIncrement*(row-2);
					eNode[3]=eNode[2]+edges[1]->nodeIncrement;
					if(theElems->HasMidsideNodes())
						eNode[6]=edges[1]->firstMidsideNode+edges[1]->nodeIncrement*(row-1);
				}
			}
			
			// Special case for 1 element area
			if(rows==1 && cols==1)
			{	eNode[1]=(edges[0]->FirstKeypoint())->node;
				eNode[2]=(edges[1]->FirstKeypoint())->node;
				eNode[3]=(edges[2]->FirstKeypoint())->node;
				eNode[4]=(edges[3]->FirstKeypoint())->node;
				if(theElems->HasMidsideNodes())
				{	eNode[5]=edges[0]->firstMidsideNode;
					eNode[6]=edges[1]->firstMidsideNode;
					eNode[7]=edges[2]->firstMidsideNode;
					eNode[8]=edges[3]->firstMidsideNode;
				}
			}
			
			// calculate angle
			if(angleExpr!=NULL)
			{	int imax= theElems->HasMidsideNodes() ? 4 : 8 ;
				Vector midPt;
				theNodes->MidPoint(&eNode[1],imax,&midPt);
				angle=Expression::FunctionValue(1,midPt.x,midPt.y,0.,0.,0.,0.);
			}
			
			if(!theElems->MeshElement(eNode,mat,angle,thick))
				return "Unable to create elements in the meshing area";
			
		}
	}

	/* adjust path faces for some elements */
	if(theElems->ElementSides()==3)
	{	if(!theElems->FlipTriangles())
		{	if(edges[2]->face>0) edges[2]->face=1;
			if(edges[3]->face>0) edges[3]->face=2;
		}
		else
		{	if(edges[3]->face>0) edges[3]->face=3;
		}
	}
	
	// all done
	return NULL;
}

/* Mesh a line in dimensionless coordinates from (x1,y1) to (x2,y2)
	ratio is ratio of first to last interval size
	numIntervals is number of intervals (>1)
*/
void Area::MeshLine(double x1,double y1,double x2,double y2,double ratio,int numIntervals,int doingMidNodes)
{
	double a,r,t,fxn[8],xNode,yNode,last;
	double x1arg,x2arg,y1arg,y2arg,xlast,ylast;
	int numMesh,i,j;
	Vector xi;
	
	GetLineParameters(numIntervals,ratio,&r,&a);
	x1arg=(x2+x1)/2.;
	x2arg=(x2-x1)/2.;
	y1arg=(y2+y1)/2.;
	y2arg=(y2-y1)/2.;
	
	numMesh=numIntervals-1;
	if(doingMidNodes) numMesh++;
	
	t=-1.+a;
	last=a;
	xlast=x1;
	ylast=y1;
	for(i=1;i<=numMesh;i++)
	{	// vector xi has dimensionless coordinates of next node
		xi.x=x1arg+t*x2arg;
		xi.y=y1arg+t*y2arg;
		if(doingMidNodes)
		{	xNode=(xi.x+xlast)/2.;
			yNode=(xi.y+ylast)/2.;
			xlast=xi.x;
			ylast=xi.y;
			xi.x=xNode;
			xi.y=yNode;
		}
		
		// Convert to real coordinates using shape functions
		areaElem->ShapeFunction(&xi,FALSE,fxn,NULL,NULL,NULL,NULL,NULL,NULL);
		xNode=0.;
		yNode=0.;
		for(j=0;j<areaElem->NumberNodes();j++)
		{	xNode+=fxn[j]*areaNode[j].x;
			yNode+=fxn[j]*areaNode[j].y;
		}
		
		// add node to list
		theNodes->AddNode(xNode,yNode,(double)0.,(double)0.0);
		
		last*=r;
		t+=last;
	}
}

/* Copy algorithm of Mesh Line to find a point along a path
	(x1,y1) to (x2,y2) line in dimensionless coordinates
	ratio is element size ratio
	numIntervals in number of intervals on full line
	numFind is which node point to find
	returns t (relative position along line) and parameters
		for subsequent call to LineLocate()
*/
void Area::FindPtOnLine(double x1,double y1,double x2,double y2,double ratio,
		int numIntervals,int numFind,double *targ,double *x1arg,
		double *x2arg,double *y1arg,double *y2arg)
{
	double a,r;
	GetLineParameters(numIntervals,ratio,&r,&a);
	*x1arg=(x2+x1)/2.;
	*x2arg=(x2-x1)/2.;
	*y1arg=(y2+y1)/2.;
	*y2arg=(y2-y1)/2.;
	
	// find dimensionless position of point numFind along the path
	if(fabs(ratio-1.)>0.001)
		*targ=-1. + a*(1.-pow(r,(double)numFind))/(1.-r);
	else
		*targ=-1. + a*(double)numFind;
}

// Call with output from Find Pt on Line
void Area::LineLocate(double t,double x1arg,double x2arg,double y1arg,double y2arg,double *xkey,double *ykey)
{
	double fxn[8],xNode,yNode;
	int j;
	Vector xi;
	
	xi.x=x1arg+t*x2arg;
	xi.y=y1arg+t*y2arg;
	areaElem->ShapeFunction(&xi,FALSE,fxn,NULL,NULL,NULL,NULL,NULL,NULL);
	xNode=0.;
	yNode=0.;
	for(j=0;j<8;j++)
	{	xNode+=fxn[j]*areaNode[j].x;
		yNode+=fxn[j]*areaNode[j].y;
	}
	*xkey=xNode;
	*ykey=yNode;
}

/*	Check that opposite side intervals match. They must have
	   1. Opposite faces with same number
			(n1+n4)==(n2+n3) and n1==n3
	   2. Or sum at opposite corners the same
			(n1+n4)==(n2+n3) where n1!=n3
			or (n1+n2)==(n3+n4) where n1!=n3
		  For latter, reorder edges to make (n1+n4)==(n2+n3) 
		  Then make n1>n3, reorder if not
		  Note: when done will always have n2>n4 too
*/
int Area::CheckIntervals(void)
{
	int n1=edges[0]->intervals;
	int n2=edges[1]->intervals;
	int n3=edges[2]->intervals;
	int n4=edges[3]->intervals;
	if((n1+n4)!=(n2+n3))
	{	if((n1+n2)==(n3+n4))
		{	// OK, but reorder faces
			Path *ptemp=edges[0];
			edges[0]=edges[1];
			edges[1]=edges[2];
			edges[2]=edges[3];
			edges[3]=ptemp;
			n1=edges[0]->intervals;
			n2=edges[1]->intervals;
			n3=edges[2]->intervals;
			n4=edges[3]->intervals;
		}
		else
			return FALSE;
	}
	
	// Valid intervals, but may want to reorder
	if(n1<n3)
	{	Path *ptemp=edges[0];
		edges[0]=edges[2];
		edges[2]=ptemp;
		ptemp=edges[1];
		edges[1]=edges[3];
		edges[3]=ptemp;
	}
	
	return TRUE;
}

/* Make sure paths can be used in solid Area meshing
	face setting depends on previous path uses
		never used = 0
		Area = element face #
		Area Area = INTERNAL_PATH
		Area Interface = -(element face #)
	No other combinations are legal
*/
int Area::PathsAvailable(void)
{
	int i;
	for(i=0;i<numPaths;i++)
	{	if(edges[i]->face==0)
			edges[i]->face=i+1;					// first use of path
		else if(edges[i]->face>0)
			edges[i]->face=INTERNAL_PATH;		// was single Area use
		else
			return FALSE;						// used twice (Area Area) or (Area Interface)
	}
	return TRUE;
}

// signed area of the enclosed area (based only on corner points)
double Area::SignedArea(void)
{
	int i;
	double ccwarea=0.;
	Keypoint *key1,*key2;
	
	for(i=0;i<numPaths;i++)
	{	key1=edges[i]->FirstKeypoint();
		key2=edges[i]->LastKeypoint();
		ccwarea+=key1->x*key2->y - key1->y*key2->x;
	}
	
	return 0.5*ccwarea;
}
