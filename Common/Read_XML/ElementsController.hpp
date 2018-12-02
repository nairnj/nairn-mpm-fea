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
		ElementBase **SetElementArray(int &,bool);
		int ElementsCompatible(int);
#ifdef MPM_CODE
		int CreateElement(char *);
		static ElementBase *MeshElement(int,int,int *);
#else
		int CreateElement(char *,int,double,double);
		int MeshElement(int *,int,double,double);
		int HasMidsideNodes(void);
		int ElementSides();
		bool FlipTriangles(void);
		void SetFlipTriangles(const char *);
		void SetFlipTriangles(bool);
		int InterfaceElements(void);
#endif		
		// accessors
		bool SetElemIDStr(char *,int);
		bool SetCurrentElemID(int,int);
		int CurrentElemID(void);
	
	private:
		int currentElemID;
#ifdef FEA_CODE
		bool flipTriangles;
#endif
};

extern ElementsController *theElems;

#endif
