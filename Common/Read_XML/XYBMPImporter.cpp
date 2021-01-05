/********************************************************************************
	XYBMPImporter.cpp
	nairn-mpm-fea

	Created by John Nairn on Thu 9/24/2016.
	Copyright (c) 2016 RSAC Software. All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/XYBMPImporter.hpp"

#pragma mark XYBMPImporter: Constructors and Destructors

// create an instance
XYBMPImporter::XYBMPImporter(char *filePath,bool isReadingXML)  : XYFileImporter(filePath,isReadingXML)
{
}

#pragma mark XYBMPImporter: Methods

// Get the file header
// throws SAXException or CommonException
void XYBMPImporter::GetXYFileHeader(XYInfoHeader &info)
{
	// read and check the header (individual reads due to Windows alignment issues)
	// 14 byte header
	short status=true;
	if(fread(header.type,2,1,fp)<1) status=false;
	if(fread(&header.size,4,1,fp)<1) status=false;
	if(fread(&header.reserved1,2,1,fp)<1) status=false;
	if(fread(&header.reserved2,2,1,fp)<1) status=false;
	if(fread(&header.offset,4,1,fp)<1) status=false;
	if(!status) XYFileImportError("Error reading the specified <BMP> file.");
	if(mustReverse)
	{	Reverse((char *)&header.size,sizeof(int));
		Reverse((char *)&header.offset,sizeof(int));
	}
	
	// exit if does not look like BMP file
	if(header.type[0]!='B' || header.type[1]!='M')
		XYFileImportError("<BMP> file is not a valid bit map file.");
	
	// read information (40 byte header)
	unsigned int header2Size;
	if(fread(&header2Size,4,1,fp)<1) status=false;
	if(fread(&info.width,4,1,fp)<1) status=false;
	if(fread(&info.height,4,1,fp)<1) status=false;
	if(fread(&header.planes,2,1,fp)<1) status=false;
	if(fread(&header.bits,2,1,fp)<1) status=false;
	if(fread(&header.compression,4,1,fp)<1) status=false;
	if(fread(&header.imagesize,4,1,fp)<1) status=false;
	if(fread(&header.xresolution,4,1,fp)<1) status=false;
	if(fread(&header.yresolution,4,1,fp)<1) status=false;
	if(fread(&header.ncolors,4,1,fp)<1) status=false;
	if(fread(&header.importantcolors,4,1,fp)<1) status=false;
	if(!status) XYFileImportError("Error reading the specified <BMP> file.");
	if(mustReverse)
	{	Reverse((char *)&header2Size,sizeof(int));
		Reverse((char *)&info.width,sizeof(int));
		Reverse((char *)&info.height,sizeof(int));
		Reverse((char *)&header.planes,2);
		Reverse((char *)&header.bits,2);
		Reverse((char *)&header.compression,sizeof(int));
		Reverse((char *)&header.imagesize,sizeof(int));
		Reverse((char *)&header.xresolution,sizeof(int));
		Reverse((char *)&header.yresolution,sizeof(int));
		Reverse((char *)&header.ncolors,sizeof(int));
		Reverse((char *)&header.importantcolors,sizeof(int));
	}
    
    // negative means image stored top to bottom, commads can flip if needed, here read all bottom to top
    info.topDown = 0;
    if(info.height<0)
    {   info.height = -info.height;
        info.topDown = 1;
    }
    
    // calculate number of colors from table in the header
    // header should have table and 14+header2Size byte header, subtract get table size in bytes
    // Each entry has 4 bytes (RGBA), so divide by 4 to get number of colors
    unsigned int ncolors = (header.offset-14-header2Size)>>2;

	// correct for Photoshop lack of information
	if(header.imagesize==0)
    {   // files size minus offset to bitmap dat}
		header.imagesize=header.size-header.offset;
    }
	if(header.ncolors==0)
        header.ncolors = ncolors;
    
	// can I read this BMP file
	if(header.planes!=1) XYFileImportError("Cannot read <BMP> files with more than 1 color plane.");
	if(header.compression!=0) XYFileImportError("Cannot read compressed <BMP> files.");
    
    // using a color table
    if(header.ncolors>0)
    {   if(header.ncolors!=ncolors || (ncolors!=2 && ncolors!=4 && ncolors!=16 && ncolors!=256))
            XYFileImportError("<BMP> files is neither a color table (<=8 bits) nor a 24- or 32-bit pixel map image.");
        
        // read the 2, 16, 64, or 256 colors in a color table
        // (i.e., bytes in image data have 8, 4, 2, or 1 pixels per byte)
        unsigned char blueByte,greenByte,redByte,junkByte;
        for(unsigned int i=0;i<header.ncolors;i++)
        {	if(fread(&blueByte,1,1,fp)<1) status=false;
            if(fread(&greenByte,1,1,fp)<1) status=false;
            if(fread(&redByte,1,1,fp)<1) status=false;
            if(fread(&junkByte,1,1,fp)<1) status=false;
            if(!status) XYFileImportError("<BMP> color table is corrupted.");
            intensity[i]=((int)redByte+(int)greenByte+(int)blueByte)/3;
        }
        
        // read one byte at a time (which might have multiple points) but each
        //	row must be padded to multiple of 4 bytes
        int colorsPerByte=1;
        if(header.ncolors==16)
            colorsPerByte=2;
        else if(header.ncolors==4)
            colorsPerByte=4;
        else if(header.ncolors==2)
            colorsPerByte=8;
        info.rowBytes = info.width/colorsPerByte;
        
        // if needed, padd row to multiple of 4 bytes
        if(info.rowBytes%4 != 0)
            info.rowBytes += 4-(info.rowBytes % 4);
        
        // BMP files do not know actual size
        // User must supply in the BMP command
        info.knowsCellSize = false;
        
        // check size
        unsigned int expectSize=(unsigned int)(info.height*info.rowBytes);
        if(header.imagesize<expectSize)
            XYFileImportError("BMP file size does not match expected data size.");
    }
    else
    {   // no colors means raw RGB data and no color table
        if(header.bits!=32 && header.bits!=24)
            XYFileImportError("<BMP> files is neither a color table (<=8 bits) nor a 24- or 32-bit pixel map image.");
        
        // check size
		info.rowBytes = header.bits == 32 ? info.width*4 : info.width*3;
        unsigned int expectSize=(unsigned int)(info.height*info.rowBytes);
        if(header.imagesize<expectSize)
            XYFileImportError("BMP file size does not match expected data size.");
    }
}

// Read all XY data into an array
// throws SAXException or CommonException
void XYBMPImporter::ReadXYFileData(unsigned char **rows,XYInfoHeader &info)
{   int colByte,row,col;
	
    if(header.ncolors>0)
    {   // color table image
        // Note color table stored with RGBA values if need to use in the future, here all reduced to gray scale
        unsigned char dataByte,mask;
        for(row=0;row<info.height;row++)
        {	// each row has chunks of rwoByte bytes
            col=0;
            for(colByte=0;colByte<info.rowBytes;colByte++)
            {	// read one byte at a time
                if(fread(&dataByte,1,1,fp)<1)
                    XYFileImportError("Error reading bit mapped file data.");
                
                // skip if width is not multiple of 4
                if(col>=info.width) continue;
                
                // unpack to intensity level
                switch(header.ncolors)
                {	case 256:		// each byte is an index
                        rows[row][col]=intensity[dataByte];
                        col++;
                        break;
                    case 16:		// each byte has two indices
                        rows[row][col]=intensity[(dataByte&0xF0)>>4];
                        col++;
                        if(col>=info.width) break;
                        rows[row][col]=intensity[dataByte&0x0F];
                        col++;
                        break;
                    case 4:        // each byte has 4 indices
                        rows[row][col]=intensity[(dataByte&0xC0)>>6];
                        col++;
                        if(col>=info.width) break;
                        rows[row][col]=intensity[(dataByte&0x30)>>4];
                        col++;
                        if(col>=info.width) break;
                        rows[row][col]=intensity[(dataByte&0x0C)>>2];
                        col++;
                        if(col>=info.width) break;
                        rows[row][col]=intensity[dataByte&0x03];
                        col++;
                        break;
                    case 2:			// each byte has eight indices
                        mask=0x80;
                        for(unsigned int i=1;i<=8;i++)
                        {	rows[row][col]=intensity[(dataByte&mask)>>(8-i)];
                            col++;
                            if(col>=info.width) break;
                            mask>>=1;
                        }
                        break;
                    default:
                        break;
                }
            }
        }
    }
    else
    {   // 32 bit pixel map, loop over rows get 4 bytes per entry
        unsigned char blueByte,greenByte,redByte,junkByte;
		int rgbSize = header.bits == 32 ? 4 : 3;
        bool status=true;
        for(row=0;row<info.height;row++)
        {    // each row rgbSize bytes per pixel
            col=0;
            for(colByte=0;colByte<info.rowBytes;colByte+=rgbSize)
            {   // read ARGB bytes in little endian order as BGRA or BGR if 24 bit
                // It is possible R and B or flipped here
                if(fread(&blueByte,1,1,fp)<1) status=false;
                if(fread(&greenByte,1,1,fp)<1) status=false;
                if(fread(&redByte,1,1,fp)<1) status=false;
				if(rgbSize == 4)
				{	if (fread(&junkByte, 1, 1, fp) < 1) status = false;
				}
                if(!status) XYFileImportError("Error reading bit mapped file data.");
                
                // for now convert to gray scale - future could make use of separate RGBA values
                unsigned int grayValue = ((int)redByte+(int)greenByte+(int)blueByte)/3;

                // assign to array
                rows[row][col]=(unsigned char)grayValue;
                col++;
            }
        }
    }
}
