/********************************************************************************
    ElementsController.cpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_XML/ElementsController.hpp"
#include "Elements/FourNodeIsoparam.hpp"
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
	int eNode[MaxElNd];
	ElementBase *newElem=NULL;

	switch(currentElemID)
	{	case FOUR_NODE_ISO:
			sscanf(xData,"%d,%d,%d,%d",&eNode[0],&eNode[1],&eNode[2],&eNode[3]);
#ifdef MPM_CODE
			newElem=new FourNodeIsoparam(1,eNode);
#else
			newElem=new FourNodeIsoparam(1,eNode,elemMat,elemAngle,elemThick);
#endif
			break;
			
#ifdef MPM_CODE
		case EIGHT_NODE_ISO_BRICK:
			sscanf(xData,"%d,%d,%d,%d,%d,%d,%d,%d",&eNode[0],&eNode[1],
					&eNode[2],&eNode[3],&eNode[4],&eNode[5],&eNode[6],&eNode[7]);
			newElem=new EightNodeIsoparamBrick(1,eNode);
			break;
#endif

#ifdef FEA_CODE
		case EIGHT_NODE_ISO:
			sscanf(xData,"%d,%d,%d,%d,%d,%d,%d,%d",&eNode[0],&eNode[1],
					&eNode[2],&eNode[3],&eNode[4],&eNode[5],&eNode[6],&eNode[7]);
			newElem=new EightNodeIsoparam(1,eNode,elemMat,elemAngle,elemThick);
			break;

		case ISO_TRIANGLE:
			sscanf(xData,"%d,%d,%d,%d,%d,%d",&eNode[0],&eNode[1],
					&eNode[2],&eNode[3],&eNode[4],&eNode[5]);
			newElem=new SixNodeTriangle(1,eNode,elemMat,elemAngle,elemThick);
			break;
			
		case CS_TRIANGLE:
			sscanf(xData,"%d,%d,%d",&eNode[0],&eNode[1],&eNode[2]);
			newElem=new CSTriangle(1,eNode,elemMat,elemAngle,elemThick);
			break;
			
		case LINEAR_INTERFACE:
			sscanf(xData,"%d,%d,%d,%d",&eNode[0],&eNode[1],&eNode[2],&eNode[3]);
			newElem=new LinearInterface(1,eNode,elemMat,elemAngle,elemThick);
			break;
			
		case QUAD_INTERFACE:
			sscanf(xData,"%d,%d,%d,%d,%d,%d",&eNode[0],&eNode[1],
					&eNode[2],&eNode[3],&eNode[4],&eNode[5]);
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

#ifdef FEA_CODE
// Create element(s) from node numbers calculated in meshing routine
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

#ifdef FEA_CODE
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
			break;
		
		case ISO_TRIANGLE:
			if(previousType==EIGHT_NODE_ISO) return TRUE;
			if(previousType==QUAD_INTERFACE) return TRUE;
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
			break;
		
		default:
			break;
	}
	return FALSE;
}
#endif

#ifdef FEA_CODE
// does current element type have mid side nodes
int ElementsController::HasMidsideNodes(void)
{
	switch(currentElemID)
	{	case EIGHT_NODE_ISO:
		case ISO_TRIANGLE:
		case QUAD_INTERFACE:
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
#ifdef FEA_CODE
	int oldElemID=currentElemID;
	currentElemID=elemID;
	return ElementsCompatible(oldElemID);
#else
	if(fmobj->IsThreeD())
	{	if(elemID!=EIGHT_NODE_ISO_BRICK)
			return FALSE;
	}
	else if(elemID!=FOUR_NODE_ISO)
		return FALSE;
	currentElemID=elemID;
	return TRUE;
#endif
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
