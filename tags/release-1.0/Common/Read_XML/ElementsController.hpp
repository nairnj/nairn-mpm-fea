/********************************************************************************
    ElementsController.hpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		ParseController.hpp
********************************************************************************/

#ifndef _ELEMENTSCONTROLLER_

#define _ELEMENTSCONTROLLER_

#include "Read_XML/ParseController.hpp"
class ElementBase;

class ElementsController : public ParseController
{
    public:

        //  Constructors and Destructor
		ElementsController(void);
	
		// methods
		void AddElement(ElementBase *);
		int SetElementArray(void);
#ifdef MPM_CODE
		int CreateElement(char *);
#else
		int CreateElement(char *,int,double,double);
		int MeshElement(long *,int,double,double);
		int ElementsCompatible(int);
		int HasMidsideNodes(void);
		int ElementSides();
		bool FlipTriangles(void);
		void SetFlipTriangles(const char *);
		void SetFlipTriangles(bool);
		int InterfaceElements(void);
#endif		
		// accessors
		int SetElemIDStr(char *);
		int SetCurrentElemID(int);
		int CurrentElemID(void);
	
	private:
		int currentElemID;
#ifdef FEA_CODE
		bool flipTriangles;
#endif
};

extern ElementsController *theElems;

#endif
