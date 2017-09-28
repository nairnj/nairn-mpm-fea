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
	Header Files: MPMReadHandler.hpp, FEAReadHandler.hpp, XYFileImporter.hpp
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

typedef struct
{	int width;						// Width of image
	int height;						// Height of image
	bool knowsCellSize;				// true if file gives xcell, ycell, and zlevel
	double xcell;					// x pixels size (if file has it, ignored if XML gives it)
	double ycell;					// y pixels size (if file has it, ignored if XML gives it)
	double zlevel;					// z value (if in the file, ignored if XML gives it)
	int rowBytes;					// bytes in row (i.e., width * bytes per pixel)
	int version;					// file version
} XYInfoHeader;

typedef struct
{	int r1,r2;						// range of rows in a bit map
	double wtr1,wtr2;				// weight for first and last row
	int c1,c2;						// range of columns in a bit map
	double wtc1,wtc2;				// weigth for first and last column
} DomainMap;

class BMPLevel;

// Input blocks
#define BAD_BLOCK -1
enum { NO_BLOCK=0,HEADER,MPMHEADER,NODELIST,MESHBLOCK,POINTSBLOCK,
        ELEMENTLIST,MATLPTS,MATERIAL,GRAVITY,FIXEDNODES,LOADEDPOINTS,
        CRACKHEADER,MULTIMATERIAL,CRACKLIST,THERMAL,CUSTOMTASKS,TASKPARAMETERS,
		GRIDBCHEADER,PARTICLEBCHEADER,BMPBLOCK,INTENSITYBLOCK,
		LOADEDNODES,LOADEDFACES,KEYPOINTBLOCK,PATHBLOCK,AREABLOCK,MATREGIONBLOCK,
		GRIDBLOCK=1000,BODYPART,BCSHAPE,BODY_SHAPE,MEMBRANEPART,
        MUST_BE_NO_BLOCK=2000 };

// input IDs
enum { DESCRIPTION=0,TEMPERATURE_EXPR,ARCHIVEROOT_NAME,CHAR_ARRAY,UNIQUE_ARCHIVEROOT_NAME };

// mesh types
enum { UNKNOWN_MESH=0,EXPLICIT_MESH,GENERATED_MESH };

// unit types
enum { SEC_UNITS=0,LENGTH_UNITS,VELOCITY_UNITS,MASS_UNITS };

// dimension requirement
enum { ANY_DIM=0,MUST_BE_2D,MUST_BE_3D };

// data type for input files
// For image or other file input
enum { BYTE_DATA=0,SHORT_DATA,INT_DATA,FLOAT_DATA,DOUBLE_DATA };

class CommonReadHandler : public DefaultHandler
{
    public:
        int input,inputID;
        
        //  Constructors and Destructor
        CommonReadHandler();
        ~CommonReadHandler();
		void FinishUp(void);
    
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
		short BMPFileCommonInput(char *,const Attributes&,int,bool);
		short EndBMPInput(char *xName,int);
		virtual void TranslateBMPFiles(void);
	
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
		static void *ReadXYFile(char *,XYInfoHeader &,int);
		static char *XYFileError(const char *,const char *);
		static const char *DecodeBMPWidthAndHeight(XYInfoHeader,double &,double &,double &,Vector &,bool);
		static bool MapDomainToImage(XYInfoHeader,Vector,Vector,Vector,Vector,double,double,DomainMap &);
		static int BMPIndex(double,int);
		static BMPLevel *FindBMPLevel(BMPLevel *,DomainMap,unsigned char **);
		static double FindAverageValue(DomainMap,unsigned char **);
	
    protected:
        int block,meshType;
        char *inputPtr;
        double gScaling;
        double mxmin,mxmax,mymin,mymax,mzmin,mzmax;
		
		// bmp file globals
		Vector orig;
		double bwidth,bheight;
        bool yflipped;
		double minAngle[3],minIntensity[3],angleScale[3];
		int numAngles;
		char bmpFileName[300],bmpAngleFileName[3][300];
};

#endif

