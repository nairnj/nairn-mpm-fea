/********************************************************************************
    PathBCController.hpp
    NairnFEA
    
    Created by John Nairn on 8/8/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		ShapeController.hpp
********************************************************************************/

#ifndef _PATHBCCONTROLLER_

#define _PATHBCCONTROLLER_

#include "Read_XML/ShapeController.hpp"

class Path;

class PathBCController : public ShapeController
{
    public:
	
		// contructors
		PathBCController(int,Path *);
        virtual bool FinishSetup(void);
	
		// methods
		virtual const char *startNodeEnumerator(int,int);
		virtual int nextNode(void);
    
        // accessors
        virtual const char *GetShapeName();
        virtual char *GetContextInfo(void);

	private:
		Path *thePath;
	
};

#endif

