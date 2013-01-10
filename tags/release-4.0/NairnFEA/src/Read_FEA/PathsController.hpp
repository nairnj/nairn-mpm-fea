/********************************************************************************
    PathsController.hpp
    NairnFEA
    
    Created by John Nairn on 6/23/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		ParseController.hpp
********************************************************************************/

#ifndef _PATHSCONTROLLER_

#define _PATHSCONTROLLER_

#include "Read_XML/ParseController.hpp"
class Path;

class PathsController : public ParseController
{
    public:
	
		~PathsController(void);
		
		// methods
		int AddKeypoint(char *);
		int ValidPath(void);
		int ValidName(char *);
		Path *FindPath(char *);

};

extern PathsController *paths;

#endif

