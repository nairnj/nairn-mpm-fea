/********************************************************************************
    FEAReadHandler.hpp
    NairnFEA
    
    Created by John Nairn on Feb 2 2002.
    Copyright (c) 2003 John A. Nairn. All rights reserved.    
	
	Dependencies
		CommonReadHandler.hpp
********************************************************************************/

#ifndef _FEAREADHANDLER_

#define _FEAREADHANDLER_

#define NO_MATERIAL 0
#define HOLE_MATERIAL -1

#include "Read_XML/CommonReadHandler.hpp"

class FEAReadHandler : public CommonReadHandler
{
    public:
        //  Constructors and Destructor
        FEAReadHandler();
        ~FEAReadHandler();
    
        // Handlers for processing FEA data
		virtual CommonAnalysis *GetCommonAnalysis(void);
		virtual bool myStartElement(char *,const Attributes&);
		virtual void myEndElement(char *);
		virtual void myCharacters(char *,const unsigned int);
		virtual void TranslateBMPFiles(void);
		virtual NodesController *DeleteTheNodes(void);
	
		// My Methods
		int GetDOFAttribute(char *);
		void ResequenceNodes(void);
		void FindPeriodicNodes(void) ;
		void RemoveEmptyElements(void);
		short BMPFileInput(char *,const Attributes&);
    
        short MatRegionInput(char *,const Attributes&);
        short EndMatRegionInput(char *,int);
        void SetRegionElements(void);
				         
    private:
        int elemMat,resequence;
        double elemAngle,elemThick;
};

#endif




