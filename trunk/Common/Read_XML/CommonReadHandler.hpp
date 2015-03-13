/********************************************************************************
    CommonReadHandler.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Jan 13 2006.
    Copyright (c) 2006 John A. Nairn. All rights reserved.    
	
	Dependencies
		none
********************************************************************************/

#ifndef _COMMONREADHANDLER_

#define _COMMONREADHANDLER_

/* All Xerces includes are here. Any file that includes this header (directly or
	indirectly) might be affected by version changes in Xerces

	Current files are:
		CommonReadHandler.cpp, ShapeController.cpp, BitMapFilesCommon.cpp
	Header Files: MPMReadHandler.hpp, FEAReadHandler.hpp
	which adds files
		CommonAnalysis.cpp, Generators.cpp, BitMapFiles.cpp, BitMapFilesFEA.cpp
*/
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/sax/SAXException.hpp>

XERCES_CPP_NAMESPACE_USE		// for Xerces >=2.3.0

// Input blocks
#define BAD_BLOCK -1
enum { NO_BLOCK=0,HEADER,MPMHEADER,NODELIST,MESHBLOCK,POINTSBLOCK,
        ELEMENTLIST,MATLPTS,MATERIAL,GRAVITY,FIXEDNODES,LOADEDPOINTS,
        CRACKHEADER,MULTIMATERIAL,CRACKLIST,THERMAL,CUSTOMTASKS,TASKPARAMETERS,
		GRIDBCHEADER,PARTICLEBCHEADER,BMPBLOCK,INTENSITYBLOCK,
		LOADEDNODES,LOADEDFACES,KEYPOINTBLOCK,PATHBLOCK,AREABLOCK,MATREGIONBLOCK,
		GRIDBLOCK=1000,BODYPART,BCSHAPE,BODY_SHAPE,
        MUST_BE_NO_BLOCK=2000 };

// input IDs
enum { DESCRIPTION=0,TEMPERATURE_EXPR,ARCHIVEROOT_NAME,CHAR_ARRAY,UNIQUE_ARCHIVEROOT_NAME };

// mesh types
enum { UNKNOWN_MESH=0,EXPLICIT_MESH,GENERATED_MESH };

// unit types
enum { SEC_UNITS=0,LENGTH_UNITS,VELOCITY_UNITS,MASS_UNITS };

// dimension requirement
enum { ANY_DIM=0,MUST_BE_2D,MUST_BE_3D };

typedef struct
{	unsigned char type[2];					// Should be 'BM'
	unsigned int size;                       // File size in Intel order bytes
	unsigned short int reserved1, reserved2;
	unsigned int offset;                     // Offset in bytes to image data
} BMPHeader;

typedef struct
{	unsigned int size;				// Header size in bytes
	int width,height;				// Width and height of image
	unsigned short int planes;		// Number of colour planes
	unsigned short int bits;		// Bits per pixel
	unsigned int compression;		// Compression type
	unsigned int imagesize;			// Image size in bytes
	int xresolution,yresolution;	// Pixels per meter
	unsigned int ncolors;			// Number of colors
	unsigned int importantcolors;	// Important colors
} BMPInfoHeader;

class CommonReadHandler : public DefaultHandler
{
    public:
        int input,inputID;
        
        //  Constructors and Destructor
        CommonReadHandler();
        ~CommonReadHandler();
    
        //  Handlers for the SAX ContentHandler interface
        void startElement(const XMLCh* const,const XMLCh* const,const XMLCh* const,const Attributes&);
        void endElement(const XMLCh *const,const XMLCh *const,const XMLCh *const);
        void characters(const XMLCh* const,const XMLSize_t);

        // prototypes for abstract methods (must override)
        virtual bool myStartElement(char *,const Attributes&) = 0;
		virtual void myEndElement(char *) = 0;
		virtual void myCharacters(char *,const unsigned int) = 0;
		
		// My Methods
		char *ReadTagValue(const char *,const Attributes&);
        void ReadTagNumber(int *,const Attributes&);
		double ReadNumericAttribute(const char *,const Attributes&,double);
        double ReadUnits(const Attributes&,int type);
		short BMPFileCommonInput(char *,const Attributes&,int);
		short EndBMPInput(char *xName,int);
		virtual void TranslateBMPFiles(void);
		int BMPIndex(double,int);
		void ReadBMPFile(char *,BMPInfoHeader &,unsigned char ***);
		char *BMPError(const char *,const char *);
		
        //  Handlers for the SAX ErrorHandler interface
        void warning(const SAXParseException&);
        void error(const SAXParseException&);
        void fatalError(const SAXParseException&);
		void ThrowCatErrorMessage(const char *,const char *);
		void ThrowCompoundErrorMessage(const char *,const char *,const char *);
		void ValidateCommand(char *,int,int);
    
		// accessors
		double ReadX(char *,double);
		double ReadY(char *,double);
		double ReadZ(char *,double);
		double ReadGridPoint(char *,double,double,double);
	
		// class methods
		static bool GetFreeFormatNumbers(char *,vector<double> &,double);
	
    protected:
        int block,meshType;
        char *inputPtr;
        double gScaling;
        double mxmin,mxmax,mymin,mymax,mzmin,mzmax;
		
		// bmp file globals
		double bwidth,bheight,xorig,yorig,zslice;
        bool yflipped;
		double minAngle,minIntensity,angleScale;
		char bmpFileName[300],bmpAngleFileName[300];
		unsigned int intensity[256];
};

#endif

