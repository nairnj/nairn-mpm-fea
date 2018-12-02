/********************************************************************************
	XYBMPImporter.hpp
	nairn-mpm-fea

	Created by John Nairn on Thu 9/24/2016.
	Copyright (c) 2016 RSAC Software. All rights reserved.

	Dependencies
		XYFileImporter.hpp, CommonReadHandler.hpp
********************************************************************************/

#ifndef __NairnMPM__XYBMPImporter__
#define __NairnMPM__XYBMPImporter__

#include "Read_XML/XYFileImporter.hpp"

typedef struct
{	unsigned char type[2];					 // Should be 'BM'
	unsigned int size;                       // File size in Intel order bytes
	unsigned short int reserved1, reserved2;
	unsigned int offset;                     // Offset in bytes to image data
	unsigned short int planes;		// Number of colour planes
	unsigned int compression;		// Compression type
	unsigned short int bits;		// Bits per pixel
	unsigned int importantcolors;	// Important colors
	unsigned int ncolors;			// Number of colors
	int xresolution,yresolution;	// Pixels per meter
	unsigned int imagesize;			// Image size in bytes (not counting header)
} BMPHeader;

class XYBMPImporter : public XYFileImporter
{
	public:
	
		// constructors
		XYBMPImporter(char *,bool);

		// methods
		virtual void GetXYFileHeader(XYInfoHeader &);
		virtual void ReadXYFileData(unsigned char **,XYInfoHeader &);
	
	private:
		unsigned int intensity[256];
		BMPHeader header;

};

#endif
