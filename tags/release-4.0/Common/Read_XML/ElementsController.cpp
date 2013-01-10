/********************************************************************************
    ElementsController.cpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_XML/ElementsController.hpp"
#include "Elements/FourNodeIsoparam.hpp"
#include "Elements/Lagrange2D.hpp"
#ifdef FEA_CODE
	#include "Elements/EightNodeIsoparam.hpp"
	#include "Elements/SixNodeTriangle.hpp"
	#include "Elements/CSTriangle.hpp"
	#include "Elements/LinearInterface.hpp"
	#include "Elements/QuadInterface.hpp"
#else
	#include "NairnMPM_Class/NairnMPM.hpp"
	#include "Elements/EightNodeIsoparamBrick.hpp"
#endif
#include "Read_XML/NodesController.hpp"
#include "Read_XML/CommonReadHandler.hpp"

ElementsController *theElems=NULL;

/********************************************************************************
	ElementsController: Constructur and Destructor
********************************************************************************/

ElementsController::ElementsController(void)  : ParseController()
{
	currentElemID=-1;
#ifdef FEA_CODE
	flipTriangles=FALSE;
#endif
}

/********************************************************************************
	ElementsController: methods
********************************************************************************/

void ElementsController::AddElement(ElementBase *newElem)
{
	AddObject(newElem);
	
	// set some properties
	newElem->num=numObjects;
}

// set array of elements when done
int ElementsController::SetElementArray(void)
{
	theElements=(ElementBase **)MakeObjectArray(0);
	if(theElements==NULL) return FALSE;
	
	// fill the array (zero based)
	ElementBase *anElem=(ElementBase *)firstObject;
	nelems=0;
	while(anElem!=NULL)
	{	theElements[nelems]=anElem;
		theElements[nelems]->FindExtent();
		nelems++;
		anElem=(ElementBase *)anElem->GetNextObject();
	}
	
	// done
	return TRUE;
}

// create specfic element from node list in character string
#ifdef MPM_CODE
int ElementsController::CreateElement(char *xData)
#else
int ElementsController::CreateElement(char *xData,int elemMat,double elemAngle,double elemThick)
#endif
{
	vector<double> pts;
	
	// read the data into a vector of doubles
	if(!CommonReadHandler::GetFreeFormatNumbers(xData,pts,1.0)) return FALSE;
	int numnds = pts.size();
	if(numnds<=0 || numnds>MaxElNd) return FALSE;
	
	// convert to integers
	int i,eNode[MaxElNd];
	for(i=0;i<numnds;i++) eNode[i] = int(pts[i]+0.1);
	
	// create elements
	ElementBase *newElem=NULL;

	switch(currentElemID)
	{	
#ifdef MPM_CODE
		// All MPM Element types
		case FOUR_NODE_ISO:
			if(numnds!=4) break;
			newElem=new FourNodeIsoparam(1,eNode);
			break;
			
		case EIGHT_NODE_ISO_BRICK:
			if(numnds!=8) break;
			newElem=new EightNodeIsoparamBrick(1,eNode);
			break;
#else
		// All FEA Element types
		case FOUR_NODE_ISO:
			if(numnds!=4) break;
			newElem=new FourNodeIsoparam(1,eNode,elemMat,elemAngle,elemThick);
			break;
			
		case EIGHT_NODE_ISO:
			if(numnds!=8) break;
			newElem=new EightNodeIsoparam(1,eNode,elemMat,elemAngle,elemThick);
			break;

		case ISO_TRIANGLE:
			if(numnds!=6) break;
			newElem=new SixNodeTriangle(1,eNode,elemMat,elemAngle,elemThick);
			break;
			
		case CS_TRIANGLE:
			if(numnds!=3) break;
			newElem=new CSTriangle(1,eNode,elemMat,elemAngle,elemThick);
			break;
			
		case NINE_NODE_LAGRANGE:
			if(numnds!=9) break;
			newElem=new Lagrange2D(1,eNode,elemMat,elemAngle,elemThick);
			break;
			
		case LINEAR_INTERFACE:
			if(numnds!=4) break;
			newElem=new LinearInterface(1,eNode,elemMat,elemAngle,elemThick);
			break;
			
		case QUAD_INTERFACE:
			if(numnds!=6) break;
			newElem=new QuadInterface(1,eNode,elemMat,elemAngle,elemThick);
			break;
			
#endif
		default:
			break;
	}
	
	// add element and return the result
	if(newElem==NULL) return FALSE;
	AddElement(newElem);
	return TRUE;
}

#ifdef MPM_CODE
// Create MPM element(s) from node numbers calculated in MPM grid() method
ElementBase *ElementsController::MeshElement(int elemID,int element,int *enode)
{
	ElementBase *newElem=NULL;
	
	switch(elemID)
	{	case FOUR_NODE_ISO:
			newElem = new FourNodeIsoparam(element,enode);
			break;
		
		case NINE_NODE_LAGRANGE:
			newElem = new Lagrange2D(element,enode);
			break;
		
		case EIGHT_NODE_ISO_BRICK:
			newElem = new EightNodeIsoparamBrick(element,enode);
			break;
		
		default:
			break;
	}
	
	return newElem;
}
			
#else
// Create FEA element(s) from node numbers calculated in FEA meshing routine
int ElementsController::MeshElement(int *eNode,int elemMat,double elemAngle,double elemThick)
{
	ElementBase *newElem=NULL;
	
	switch(currentElemID)
	{	case FOUR_NODE_ISO:
			newElem=new FourNodeIsoparam(1,&eNode[1],elemMat,elemAngle,elemThick);
			if(newElem==NULL) return FALSE;
			AddElement(newElem);
			break;
			
		case EIGHT_NODE_ISO:
			newElem=new EightNodeIsoparam(1,&eNode[1],elemMat,elemAngle,elemThick);
			if(newElem==NULL) return FALSE;
			AddElement(newElem);
			break;
		
		case NINE_NODE_LAGRANGE:
			int lNode[9],i;
			Vector qmidPt;
			theNodes->MidPoint(&eNode[1],8,&qmidPt);
			theNodes->AddNode(qmidPt.x,qmidPt.y,(double)0.,(double)0.0);
			for(i=1;i<=8;i++) lNode[i-1]=eNode[i];
			lNode[8]=theNodes->numObjects;
			newElem=new Lagrange2D(1,lNode,elemMat,elemAngle,elemThick);
			if(newElem==NULL) return FALSE;
			AddElement(newElem);
			break;
		
		case ISO_TRIANGLE:
			int tNode[9];
			Vector midPt;
			theNodes->MidPoint(&eNode[1],8,&midPt);
			theNodes->AddNode(midPt.x,midPt.y,(double)0.,(double)0.0);
			if(!FlipTriangles())
			{	tNode[1]=eNode[1];
				tNode[2]=eNode[2];
				tNode[3]=eNode[3];
				tNode[4]=eNode[5];
				tNode[5]=eNode[6];
				tNode[6]=theNodes->numObjects;
				newElem=new SixNodeTriangle(1,&tNode[1],elemMat,elemAngle,elemThick);
				if(newElem==NULL) return FALSE;
				AddElement(newElem);
				
				tNode[1]=eNode[3];
				tNode[2]=eNode[4];
				tNode[3]=eNode[1];
				tNode[4]=eNode[7];
				tNode[5]=eNode[8];
				newElem=new SixNodeTriangle(1,&tNode[1],elemMat,elemAngle,elemThick);
				if(newElem==NULL) return FALSE;
				AddElement(newElem);
			}
			else
			{	tNode[1]=eNode[1];
				tNode[2]=eNode[2];
				tNode[3]=eNode[4];
				tNode[4]=eNode[5];
				tNode[5]=theNodes->numObjects;
				tNode[6]=eNode[8];
				newElem=new SixNodeTriangle(1,&tNode[1],elemMat,elemAngle,elemThick);
				if(newElem==NULL) return FALSE;
				AddElement(newElem);
				
				tNode[1]=eNode[4];
				tNode[2]=eNode[2];
				tNode[3]=eNode[3];
				tNode[4]=theNodes->numObjects;
				tNode[5]=eNode[6];
				tNode[6]=eNode[7];
				newElem=new SixNodeTriangle(1,&tNode[1],elemMat,elemAngle,elemThick);
				if(newElem==NULL) return FALSE;
				AddElement(newElem);
			}
			break;
		
		case CS_TRIANGLE:
			if(!FlipTriangles())
			{	newElem=new CSTriangle(1,&eNode[1],elemMat,elemAngle,elemThick);
				if(newElem==NULL) return FALSE;
				AddElement(newElem);
				
				eNode[5]=eNode[1];
				newElem=new CSTriangle(1,&eNode[3],elemMat,elemAngle,elemThick);
				if(newElem==NULL) return FALSE;
				AddElement(newElem);
			}
			else
			{	int holdNode=eNode[3];
				eNode[3]=eNode[4];
				newElem=new CSTriangle(1,&eNode[1],elemMat,elemAngle,elemThick);
				if(newElem==NULL) return FALSE;
				AddElement(newElem);
				
				eNode[1]=eNode[4];
				eNode[3]=holdNode;
				newElem=new CSTriangle(1,&eNode[1],elemMat,elemAngle,elemThick);
				if(newElem==NULL) return FALSE;
				AddElement(newElem);
			}
			break;
		
		case LINEAR_INTERFACE:
			newElem=new LinearInterface(1,&eNode[1],elemMat,elemAngle,elemThick);
			if(newElem==NULL) return FALSE;
			AddElement(newElem);
			break;
		
		case QUAD_INTERFACE:
			newElem=new QuadInterface(1,&eNode[1],elemMat,elemAngle,elemThick);
			if(newElem==NULL) return FALSE;
			AddElement(newElem);
			break;
		
		default:
			return FALSE;
	}
	
	return TRUE;
}
#endif

// verify new element type is valid and compatible with previous one
int ElementsController::ElementsCompatible(int previousType)
{
	// if previous was not used, then make sure it is valid
	if(previousType<0)
	{	switch(currentElemID)
		{	case FOUR_NODE_ISO:
			case EIGHT_NODE_ISO:
			case ISO_TRIANGLE:
			case CS_TRIANGLE:
			case LINEAR_INTERFACE:
			case QUAD_INTERFACE:
			case EIGHT_NODE_ISO_BRICK:
			case NINE_NODE_LAGRANGE:
				return TRUE;
				break;
			
			default:
				break;
		}
		return FALSE;
	}
	
	// if the same then always OK
	if(currentElemID==previousType) return TRUE;
	
	// if different, look for compatibility
	switch(currentElemID)
	{	case FOUR_NODE_ISO:
			if(previousType==CS_TRIANGLE) return TRUE;
			if(previousType==LINEAR_INTERFACE) return TRUE;
			break;
		
		case EIGHT_NODE_ISO:
			if(previousType==ISO_TRIANGLE) return TRUE;
			if(previousType==QUAD_INTERFACE) return TRUE;
			if(previousType==NINE_NODE_LAGRANGE) return TRUE;
			break;
		
		case ISO_TRIANGLE:
			if(previousType==EIGHT_NODE_ISO) return TRUE;
			if(previousType==QUAD_INTERFACE) return TRUE;
			if(previousType==NINE_NODE_LAGRANGE) return TRUE;
			break;
		
		case CS_TRIANGLE:
			if(previousType==FOUR_NODE_ISO) return TRUE;
			if(previousType==LINEAR_INTERFACE) return TRUE;
			break;
		
		case LINEAR_INTERFACE:
			if(previousType==CS_TRIANGLE) return TRUE;
			if(previousType==FOUR_NODE_ISO) return TRUE;
			break;
		
		case QUAD_INTERFACE:
			if(previousType==EIGHT_NODE_ISO) return TRUE;
			if(previousType==ISO_TRIANGLE) return TRUE;
			if(previousType==NINE_NODE_LAGRANGE) return TRUE;
			break;
		
		case NINE_NODE_LAGRANGE:
			if(previousType==ISO_TRIANGLE) return TRUE;
			if(previousType==QUAD_INTERFACE) return TRUE;
			if(previousType==EIGHT_NODE_ISO) return TRUE;
			break;

		default:
			break;
	}
	return FALSE;
}

#ifdef FEA_CODE
// does current element type have mid side nodes
int ElementsController::HasMidsideNodes(void)
{
	switch(currentElemID)
	{	case EIGHT_NODE_ISO:
		case ISO_TRIANGLE:
		case QUAD_INTERFACE:
		case NINE_NODE_LAGRANGE:
			return TRUE;
		
		default:
			break;
	}
	return FALSE;
}
#endif

#ifdef FEA_CODE
// number of sides for proposed element type
int ElementsController::ElementSides(void)
{
	switch(currentElemID)
	{	case EIGHT_NODE_ISO:
		case FOUR_NODE_ISO:
		case NINE_NODE_LAGRANGE:
			return 4;
		case ISO_TRIANGLE:
		case CS_TRIANGLE:
			return 3;
		case LINEAR_INTERFACE:
		case QUAD_INTERFACE:
			return 2;
		default:
			break;
	}
	return 0;
}
#endif

/********************************************************************************
	ElementsController: accessors
********************************************************************************/

int ElementsController::CurrentElemID(void) { return currentElemID; }
int ElementsController::SetElemIDStr(char *value)
{	int tempID;
	sscanf(value,"%d",&tempID);
	return SetCurrentElemID(tempID);
}
int ElementsController::SetCurrentElemID(int elemID)
{
	int oldElemID=currentElemID;
	currentElemID=elemID;
#ifdef MPM_CODE
	// limit elements allowed in MPM
	if(fmobj->IsThreeD())
	{	if(elemID!=EIGHT_NODE_ISO_BRICK) return FALSE;
		return TRUE;
	}
	else if(elemID!=FOUR_NODE_ISO && elemID!=NINE_NODE_LAGRANGE)
		return FALSE;
#endif
	return ElementsCompatible(oldElemID);
}

#ifdef FEA_CODE
bool ElementsController::FlipTriangles(void) { return flipTriangles; }
void ElementsController::SetFlipTriangles(const char *value)
{	int tempID;
	sscanf(value,"%d",&tempID);
	flipTriangles = (tempID==0) ? FALSE : TRUE ;
}
void ElementsController::SetFlipTriangles(bool newvalue) { flipTriangles=newvalue; }
int ElementsController::InterfaceElements(void)
{	return (currentElemID==LINEAR_INTERFACE || currentElemID==QUAD_INTERFACE);
}
#endif
