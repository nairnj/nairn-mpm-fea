/********************************************************************************
	XYTXTImporter.cpp
	nairn-mpm-fea

	Created by John Nairn on Thu 8/30/2021.
	Copyright (c) 2021 RSAC Software. All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/XYTXTImporter.hpp"

#pragma mark XYTXTImporter: Constructors and Destructors

// create an instance
XYTXTImporter::XYTXTImporter(char *filePath,bool isReadingXML,int dataType) :
				XYFileImporter(filePath,isReadingXML,dataType)
{
}

#pragma mark XYFPGImporter: Methods

// Get the file header
// For text file, read to buffer and count rows and columns
// Assume all rows have the same number of columns
// throws SAXException or CommonException
void XYTXTImporter::GetXYFileHeader(XYInfoHeader &info)
{
	info.dataType = FLOAT_DATA;
	
	// get file length and read it
	if(fseek(fp,0L,SEEK_END)!=0)
		XYFileImportError("Error reading the specified text file.");
	fileLength=ftell(fp);
	rewind(fp);
	
	// create buffer
	buffer=(unsigned char *)malloc(fileLength);
	if(buffer==NULL)
		XYFileImportError("Memory error reading the specified text file.");
	
	// read the file
	size_t result = fread(buffer, 1, fileLength, fp);
	if(result != fileLength)
	{	if(ferror(fp)==0)
		{	// assume Windows has dropped CRs
			fileLength = result;
		}
		else
			XYFileImportError("Error reading all data from a text file.");
	}
	
	// count columns and rows
	info.width = 1;
	info.height = 0;
	
	for(int i=0;i<fileLength;i++)
	{	// if new line, increment height
		//   (but ignore ignore zero-character lines, including CR/LF endings)
		if(buffer[i]=='\n' || buffer[i]=='\r')
		{	if(i>0)
			{	if(buffer[i-1]!='\r' && buffer[i-1]!='\n')
					info.height++;
			}
		}
		else if(buffer[i]=='\t' || buffer[i]==',')
		{	// all tabs or commas, but convert to tabs
			buffer[i] = '\t';
			
			// count columns (first row only)
			if(info.height==0) info.width++;
		}
	}
	
	if(info.height==0 || info.width==1)
		XYFileImportError("Text file has no rows or no delimited columns");
	
	// first row of text file is top of the image
	info.topDown = 1;
	
	// text files do not know their cell size, user must supply
	info.knowsCellSize = false;
}

// Read all XY data into an array
// throws SAXException or CommonException
void XYTXTImporter::ReadXYFileData(unsigned char **rows,XYInfoHeader &info)
{
	// read into floats
	float **frows = (float **)rows;
	
	unsigned char *p = buffer;
	unsigned char *lastp = buffer+fileLength;
	for(int row=0;row<info.height;row++)
	{	// get all columns
		for(int col=0;col<info.width;col++)
		{	// get first float
			sscanf((char *)p,"%f",&frows[row][col]);

			// scan to next delimiter
			p++;
			while(*p!='\t' && *p!='\r' && *p!='\n' && p<lastp) p++;
			if((*p!='\t' || p>=lastp) && col<info.width-1)
			{	// row has too few columns
				XYFileImportError("Row of a text file has too few columns.");
			}
			
			p++;
		}
		
		// skip CR/LF or blank lines
		while((*p=='\r' || *p=='\n') && p<lastp) p++;
		
		// error if reached end of data before the last row
		if(p>=lastp && row<info.height-1)
		{	// row has too few columns
			XYFileImportError("Text files has too few rows.");
		}
	}
}
