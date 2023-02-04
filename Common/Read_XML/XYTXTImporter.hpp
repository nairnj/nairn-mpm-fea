/********************************************************************************
	XYTXTImporter.hpp
	nairn-mpm-fea

	Created by John Nairn on Thu 8/30/2021.
	Copyright (c) 2021 RSAC Software. All rights reserved.

	Dependencies
		XYFileImporter.hpp, CommonReadHandler.hpp
********************************************************************************/

#ifndef __NairnMPM__XYTXTImporter__
#define __NairnMPM__XYTXTImporter__

#include "Read_XML/XYFileImporter.hpp"

class XYTXTImporter : public XYFileImporter
{
	public:
	
		// constructors
		XYTXTImporter(char *,bool,int);
	
		// methods
		virtual void GetXYFileHeader(XYInfoHeader &);
		virtual void ReadXYFileData(unsigned char **,XYInfoHeader &);
	
	protected:
		long fileLength;
		unsigned char *buffer;
		char delim;
};

#endif

