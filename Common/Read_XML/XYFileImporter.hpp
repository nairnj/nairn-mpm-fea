/********************************************************************************
	XYFileImporter.hpp
	nairn-mpm-fea

	Created by John Nairn on Thu 9/24/2016.
	Copyright (c) 2016 RSAC Software. All rights reserved.

	Dependencies
		CommonReadHandler.hpp
********************************************************************************/

#ifndef __NairnMPM__XYFileImporter__
#define __NairnMPM__XYFileImporter__

#include <fstream>
#include "Read_XML/CommonReadHandler.hpp"


class XYFileImporter
{
	public:
	
		// constructors and destructors
		XYFileImporter();
		XYFileImporter(char *,bool);
		virtual ~XYFileImporter();
		
		// virtual methods
		virtual void GetXYFileHeader(XYInfoHeader &) = 0;
		virtual void ReadXYFileData(unsigned char **,XYInfoHeader &) = 0;

		// methods
		void XYFileImportError(const char *msg);
		
	protected:
		char *fullPath;
		FILE *fp;
		int mustReverse;
		bool readingXML;
};

#endif
