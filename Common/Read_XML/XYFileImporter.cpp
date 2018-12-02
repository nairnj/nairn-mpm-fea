/********************************************************************************
	XYFileImporter.cpp
	nairn-mpm-fea

	Created by John Nairn on Thu 9/24/2016.
	Copyright (c) 2016 RSAC Software. All rights reserved.
 
	Parent class to support reading x-y digital date from
	various file types
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/XYFileImporter.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark XYFileImporter: Constructors and Destructors

// Constructors
XYFileImporter::XYFileImporter()
{	fp = NULL;
	fullPath = NULL;
	readingXML = true;
}

XYFileImporter::XYFileImporter(char *filePath,bool isReadingXML)
{
	readingXML = isReadingXML;
	
	// save file name
	fullPath = new char[strlen(filePath)+1];
	strcpy(fullPath,filePath);
	
	// open the file
	if((fp=fopen(fullPath,"r"))==NULL)
		XYFileImportError("The bit mapped file could not be opened.");
	
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
// throws SAXException or std::bad_alloc or CommonException
void XYFileImporter::XYFileImportError(const char *msg)
{
	if(fp!=NULL) fclose(fp);
	char *error = new (nothrow) char[strlen(msg)+strlen(fullPath)+10];
	if(error != NULL)
	{	strcpy(error,msg);
		strcat(error," (file: ");
		strcat(error,fullPath);
		strcat(error,")");
		if(readingXML)
			throw SAXException(error);
		else
			throw CommonException(error,"XYFileImporter::XYFileImportError");
	}
	else
		if(readingXML)
			throw SAXException(msg);
		else
			throw CommonException(msg,"XYFileImporter::XYFileImportError");
}

