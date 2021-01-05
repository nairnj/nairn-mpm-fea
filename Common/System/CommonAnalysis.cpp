/********************************************************************************
    CommonAnalysis.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Jan 13, 2006.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/CommonReadHandler.hpp"
#include "Exceptions/CommonException.hpp"
#include "Exceptions/StrX.hpp"
#include "Materials/MaterialBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Elements/ElementBase.hpp"
#include "System/UnitsController.hpp"

// reeading globals
static bool doSchema=TRUE;
static bool schemaFullChecking=FALSE;
extern const char *svninfo;

#pragma mark CommonAnalysis: Constructors and Destructor

CommonAnalysis::CommonAnalysis()
{
	validate=FALSE;				// validate with DTD file
	reverseBytes=FALSE;			// reverse bytes in archived files
	description=NULL;
	SetDescription("\nNo Description");
	
	np=-1;						// analysis method to be set
	nfree=2;					// 2D analysis
	
	// If non-zero, this flags can change the calculation
	int i;
	for(i=0;i<NUMBER_DEVELOPMENT_FLAGS;i++) dflag[i]=0;
}

CommonAnalysis::~CommonAnalysis()
{
	delete [] description;
}

#pragma mark CommonAnalysis: Running Analysis

// Start comuptational mechanics calculations
void CommonAnalysis::StartAnalysis(bool abort)
{
	// start output file
	StartResultsOutput();
	
	// finish with analysis specific items
	CMStartResultsOutput();
	
	// prepare for computational mechanics
	CMPreparations();
	
	// run the computational mechancis
	CMAnalysis(abort);
	
}

// start results file before proceeding to analysis
// throws CommonException()
void CommonAnalysis::StartResultsOutput(void)
{
    //--------------------------------------------------
    // Analysis Title
	PrintAnalysisTitle();
	cout << "Written by: Nairn Research Group, Oregon State University\n"
        << "Date: " << __DATE__ << "\n"
        << "Source: " << svninfo << "\n"
		<< "Units: ";
		UnitsController::OutputUnits();
		cout << "\n"
#ifdef _OPENMP
        << "Processors: " << numProcs << "\n"
#endif
        << endl;
    
    //--------------------------------------------------
    // Description
    PrintSection("ANALYSIS DESCRIPTION");
    cout << GetDescription() << endl;

	/*	NairnMPM: Current Flags in Use
	 *   dflag[3]==xxyyzz - custom patch method for parallel (ignored if xx*yy*zz != numProcs)
	 *
	 */
	int i;
	bool hasFlags=FALSE;
	for(i=0;i<NUMBER_DEVELOPMENT_FLAGS;i++)
	{	if(dflag[i]!=0)
		{	cout << "Development Flag #" << i << " = " << dflag[i] << endl;
			cout << "... WARNING: unrecognized developer flag specified" << endl;
			hasFlags=TRUE;
		}
	}
	if(hasFlags) cout << endl;
	
    //--------------------------------------------------
    // Analysis Type
	PrintAnalysisMethod();
	
    // start section (nodes and elements in grid)
	PrintSection(NodesAndElementsTitle());
    
    //---------------------------------------------------
    // Nodes
	char hline[200];
    if(nnodes==0 || nelems<0)
        throw CommonException("No nodes or elements were defined in the input file.","CommonAnalysis::StartResultsOutput");
    sprintf(hline,"Nodes: %d       Elements: %d\n",nnodes,nelems);
    cout << hline;
    sprintf(hline,"DOF per node: %d     ",nfree);
    cout << hline;
	GetAnalysisType(np,hline);
	cout << hline << endl << endl;
	
    //---------------------------------------------------
    // Nodes
	sprintf(hline,"NODAL POINT COORDINATES (in %s)",UnitsController::Label(CULENGTH_UNITS));
    PrintSection(hline);
	ArchiveNodalPoints(np);

    //---------------------------------------------------
    // Elements
    PrintSection("ELEMENT DEFINITIONS");
	ArchiveElements(np);
    
    //---------------------------------------------------
    // Defined Materials
	const char *err;
#ifdef MPM_CODE
	sprintf(hline,"DEFINED MATERIALS\n       (Note: moduli and stresses in %s, thermal exp. coeffs in ppm/K)\n       (      density in %s)",
			UnitsController::Label(PRESSURE_UNITS),UnitsController::Label(DENSITY_UNITS));
#else
	sprintf(hline,"DEFINED MATERIALS\n       (Note: moduli in %s, thermal exp. coeffs in ppm/K)",
			UnitsController::Label(PRESSURE_UNITS));
#endif
    PrintSection(hline);
    if(theMaterials==NULL)
        throw CommonException("No materials were defined.","CommonAnalysis::StartResultsOutput");
    for(i=0;i<nmat;i++)
	{	err = theMaterials[i]->VerifyAndLoadProperties(np);
    	theMaterials[i]->PrintMaterial(i+1);
		cout << endl;
        if(err!=NULL)
		{	cout << "Invalid material properties\n   " << err << endl;
            throw CommonException("Invalid material properties - see output file for details","CommonAnalysis::StartResultsOutput");
		}
   }
	
	//---------------------------------------------------
	// initialize timers
#ifdef _OPENMP
	startTime = omp_get_wtime();
#else
    time(&startTime);
#endif
	startCPU=clock();
}

// print analysis type (if needed)
void CommonAnalysis::PrintAnalysisMethod(void) {}

// Do an iniitial preparations before starting the analysis
// throws CommonException()
void CommonAnalysis::CMPreparations(void) {}

#pragma mark CommonAnalysis: Base Methods

// Main entry to read file and decode into objects
// xmlFile are command arg - relative or full path, dos path in windows
int CommonAnalysis::ReadFile(const char *xmlFile,bool useWorkingDir)
{
	// set directory of input file
	SetInputDirPath(xmlFile,useWorkingDir);
	
    // Initialize the XML4C2 system
	SAX2XMLReader* parser;
    try
    {	XMLPlatformUtils::Initialize();
	
		// Create a SAX parser object.
		parser=XMLReaderFactory::createXMLReader();
		
		// Validation (first means yes or no, second means to skip if no DTD, even if yes)
		parser->setFeature(XMLString::transcode("http://xml.org/sax/features/validation"),GetValidate());
		parser->setFeature(XMLString::transcode("http://apache.org/xml/features/validation/dynamic"),true);
		if(!GetValidate())
		{	parser->setFeature(XMLString::transcode("http://apache.org/xml/features/nonvalidating/load-external-dtd"),
									GetValidate());
		}

		// default settings
		parser->setFeature(XMLString::transcode("http://apache.org/xml/features/validation/schema"),doSchema);
		parser->setFeature(XMLString::transcode("http://apache.org/xml/features/validation/schema-full-checking"),schemaFullChecking);
    }

    catch(const XMLException& toCatch)
    {   cerr << "\nAn error occurred while initializing Xerces: "
            << StrX(toCatch.getMessage()) << endl;
        return XercesInitErr;
    }
	
	catch(std::bad_alloc&)
	{	cerr << "\nMemory error initializing Xerces" << endl;
		return XercesInitErr;
	}
	
    catch( ... )
    {	cerr << "\nUnexpected error initializing Xerces" << endl;
		return XercesInitErr;
    }

    // parce the file
	CommonReadHandler *handler = GetReadHandler();
    try
	{	// create the handler
		parser->setContentHandler(handler);
		parser->setErrorHandler(handler);
		parser->parse(xmlFile);
		
		// Reading done, now clean up and exit
		delete parser;
		handler->FinishUp();
		delete handler;
		XMLPlatformUtils::Terminate();
    }

    catch(const XMLException& toCatch)
	{	char *message = XMLString::transcode(toCatch.getMessage());
        cerr << "\nParcing error: " << message << endl;
        return ReadFileErr;
    }
    
    catch(const SAXException& err)
	{	char *message = XMLString::transcode(err.getMessage());
    	cerr << "\nInput file error: " << message << endl;
        return ReadFileErr;
    }
	
	catch(std::bad_alloc& ba)
	{	cerr << "\nMemory error: " << ba.what() << endl;
		return ReadFileErr;
	}
	
	catch(const char *errMsg)
	{	cerr << errMsg << endl;
		return ReadFileErr;
	}
    
    catch( ... )
    {	cerr << "Unexpected error occurred while reading XML file" << endl;
		return ReadFileErr;
    }
    
    return noErr;
}

// print conde version to standard output
void CommonAnalysis::CoutCodeVersion(void)
{
	cout << CodeName() << " " << version << "." << subversion << " build " << buildnumber << endl;
}

// elapsed actual time
double CommonAnalysis::ElapsedTime(void)
{
#ifdef _OPENMP
	return omp_get_wtime()-startTime;
#else
	time_t currentTime;
    time(&currentTime);
    return (double)difftime(currentTime,startTime);
#endif
}

// elapsed CPU time
double CommonAnalysis::CPUTime(void)
{
	return (double)(clock()-startCPU)/CLOCKS_PER_SEC;
}

#pragma mark CommonAnalysis: accessors

// validation option
void CommonAnalysis::SetValidate(bool setting) { validate=setting; }
bool CommonAnalysis::GetValidate(void) { return validate; }

// reverse bytes (MPM only)
void CommonAnalysis::SetReverseBytes(bool setting) { reverseBytes=setting; }
bool CommonAnalysis::GetReverseBytes(void) { return reverseBytes; }

// set description
// throws std::bad_alloc
void CommonAnalysis::SetDescription(const char *descrip)
{	if(description!=NULL) delete [] description;
    description=new char[strlen(descrip)+1];
    strcpy(description,descrip);
}
// point to second character because first is always new line character
char *CommonAnalysis::GetDescription(void) { return &description[1]; }

// is it 3D analysis
bool CommonAnalysis::IsThreeD(void) { return np==THREED_MPM; }

// is it axisymmetric (FEA or MPM)
bool CommonAnalysis::IsAxisymmetric(void) { return np==AXI_SYM || np==AXISYMMETRIC_MPM; }

// string about MPM additions (when printing out analysis type
const char *CommonAnalysis::MPMAugmentation(void) { return ""; }

// set number of processes (0 for serial code)
void CommonAnalysis::SetNumberOfProcessors(int npr) { numProcs = npr; }
int CommonAnalysis::GetNumberOfProcessors(void) { return numProcs; }
int CommonAnalysis::GetTotalNumberOfPatches(void) { return numProcs>1 ? numProcs : 1 ; }

#pragma mark Methods to keep Xerces contact in fewer files

// throws SAXException()
void ThrowSAXException(const char *msg)
{	throw SAXException(msg);
}

// throw an SAXException with format string and one string variable
// throws SAXException()
void ThrowSAXException(const char *frmt,const char *msg)
{	char eline[500];
	sprintf(eline,frmt,msg);
	throw SAXException(eline);
}





