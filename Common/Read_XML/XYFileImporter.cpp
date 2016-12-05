/********************************************************************************
	XYFileImporter.cpp
	nairn-mpm-fea

	Created by John Nairn on Thu 9/24/2016.
	Copyright (c) 2016 RSAC Software. All rights reserved.
 
	Parent class to support reading x-y digital date from
	various file types
********************************************************************************/

#include "stdafx.h"
#include "XYFileImporter.hpp"

#pragma mark XYFileImporter: Constructors and Destructors

// Constructors
XYFileImporter::XYFileImporter()
{	fp = NULL;
	fullPath = NULL;
}

XYFileImporter::XYFileImporter(char *filePath)
{
	// save file name
	fullPath = new char[strlen(filePath)+1];
	strcpy(fullPath,filePath);
	
	// open the file
	if((fp=fopen(fullPath,"r"))==NULL)
		XYFileError("The bit mapped file could not be opened.");
	
	// decide what computer I am on
	// byte order marker
	mustReverse=false;
	int test=1;
	char *testPtr=(char *)&test;
	if(*testPtr==0) mustReverse=true;

}

// Destructor (and it is virtual)
XYFileImporter::~XYFileImporter()
{	if(fp!=NULL) fclose(fp);
	if(fullPath!=NULL) delete [] fullPath;
}

#pragma mark XYFileImporter: Methods

// create error message with bmp file name
// throws SAXException or std::bad_alloc
void XYFileImporter::XYFileError(const char *msg)
{
	if(fp!=NULL) fclose(fp);
	char *error = new (nothrow) char[strlen(msg)+strlen(fullPath)+10];
	if(error != NULL)
	{	strcpy(error,msg);
		strcat(error," (file: ");
		strcat(error,fullPath);
		strcat(error,")");
		throw SAXException(error);
	}
	else
		throw SAXException(msg);
}

