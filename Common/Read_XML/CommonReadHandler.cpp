/********************************************************************************
    CommonReadHandler.cpp
    NairnFEA
    
    Created by John Nairn on Jan 13 2006.
    Copyright (c) 2006 John A. Nairn. All rights reserved.
********************************************************************************/

#ifdef MPM_CODE
	#include "NairnMPM_Class/NairnMPM.hpp"
#else
	#include "NairnFEA_Class/NairnFEA.hpp"
#endif
#include "Read_XML/CommonReadHandler.hpp"
#include "Read_XML/MaterialController.hpp"
#include "Read_XML/NodesController.hpp"
#include "Read_XML/ElementsController.hpp"
#include "Elements/ElementBase.hpp"
#include "Exceptions/StrX.hpp"

/********************************************************************************
	CommonReadHandler: Constructors and Destructor
********************************************************************************/

CommonReadHandler::CommonReadHandler()
{
	meshType=UNKNOWN_MESH;
    block=NO_BLOCK;
	matCtrl=new MaterialController();
    gScaling=1.;
    mxmin=mymin=mzmin=0.;
    mxmax=mymax=mzmax=-1.;
};

CommonReadHandler::~CommonReadHandler()
{
	// set material array
	const char *emsg = matCtrl->SetMaterialArray();
	if(emsg!=NULL)
		throw SAXException(emsg);
	delete matCtrl;
}

/********************************************************************************
	CommonReadHandler: Methods
********************************************************************************/

// Start a new element
void CommonReadHandler::startElement(const XMLCh* const uri,const XMLCh* const localname,
            const XMLCh* const qname,const Attributes& attrs)
{
	// decode tag name
    char *xName=XMLString::transcode(localname);
	input=NO_INPUT;
	
    // Handle all elements common between MPM and FEA analysis

    //-------------------------------------------------------
    // <Header> elements
    if(strcmp(xName,"Header")==0)
		block=HEADER;
        
    else if(strcmp(xName,"Description")==0)
    {	if(block!=HEADER)
            throw SAXException("<Description> must be within the <Header> element.");
        input=TEXT_BLOCK;
		inputID=DESCRIPTION;
    }
    
    else if(strcmp(xName,"Analysis")==0)
	{   if(block!=HEADER)
            throw SAXException("<Analysis> must be within the <Header> element.");
    	input=ANALYSIS_NUM;
    }
	
	// begin a material
    else if(strcmp(xName,"DevelFlag")==0)
    {	if(block!=HEADER)
			throw SAXException("<DevelFlag> must be within the <Header> element.");
		double flagNumDble=ReadNumericAttribute("Number",attrs,(double)0.0);
    	int flagNum=(int)(flagNumDble+0.5);
		if(!flagNum<0 || flagNum>=NUMBER_DEVELOPMENT_FLAGS)
			throw SAXException("The <DevelFlag> 'Number' must be from 0 to 4");
		input=INT_NUM;
		inputPtr=(char *)&fmobj->dflag[flagNum];
    }
    //-------------------------------------------------------
    // <Mesh> block
	
    // Node list
    else if(strcmp(xName,"NodeList")==0)
	{	ValidateCommand(xName,MESHBLOCK,ANY_DIM);
    	if(meshType!=UNKNOWN_MESH)
            throw SAXException("<NodeList> can not be used with a generated mesh.");
    	block=NODELIST;
		meshType=EXPLICIT_MESH;
		if(theNodes==NULL) theNodes=new NodesController();
    }
        
    // Element list
    else if(strcmp(xName,"ElementList")==0)
	{	ValidateCommand(xName,MESHBLOCK,ANY_DIM);
		if(meshType!=EXPLICIT_MESH)
            throw SAXException("<ElementList> cannot be used with a generated mesh.");
    	block=ELEMENTLIST;
		// MPM only uses when in ElementList (which is uncommon because does not suppport GIMP)
		if(theElems==NULL) theElems=new ElementsController();
	}
        
    //-----------------------------------------------------------
	// <GridBCs> section
	
    // DisplacementBC section
    else if(strcmp(xName,"DisplacementBCs")==0)
	{	ValidateCommand(xName,GRIDBCHEADER,ANY_DIM);
	    block=FIXEDNODES;
    }
        
    //-------------------------------------------------------
    // <Thermal> section
	
	// begin Thermal section
    else if(strcmp(xName,"Thermal")==0)
    {	block=THERMAL;
    }
    
    //-------------------------------------------------------
    // <Material> Definitions
	
	// begin a material
    else if(strcmp(xName,"Material")==0)
    {	block=MATERIAL;
        int matID=0;		// invalid ID unless it is set
		char matName[200],*aName,*value;
		matName[0]=0;		// to get the name of the material
    	int i,numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"Type")==0)
                sscanf(value,"%d",&matID);
            else if(strcmp(aName,"Name")==0)
			{	if(strlen(value)>199) value[200]=0;
                strcpy(matName,value);
			}
            delete [] aName;
            delete [] value;
        }
		if(strlen(matName)==0)
			throw SAXException("<Material> must be named using 'Name' atttribute.");
		if(!matCtrl->AddMaterial(matID,matName))
			ThrowCatErrorMessage("Invalid material: either undefined or not allowed for current analysis type",matName);
    }

	// begin a material
    else if(strcmp(xName,"color")==0)
    {	if(block!=MATERIAL)
			throw SAXException("<color> must be within a <Material> definition.");
		float redClr=-1.,greenClr=-1.,blueClr=-1.;
		char *aName,*value;
    	int i,numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"red")==0)
                sscanf(value,"%f",&redClr);
            else if(strcmp(aName,"green")==0)
                sscanf(value,"%f",&greenClr);
            else if(strcmp(aName,"blue")==0)
                sscanf(value,"%f",&blueClr);
            delete [] aName;
            delete [] value;
        }
		if(redClr>0.)
		{	if(redClr>1.) redClr=1.;
			if(greenClr<0.) greenClr=redClr;
			if(blueClr<0.) blueClr=redClr;
			if(greenClr>1.) greenClr=1.;
			if(blueClr>1.) blueClr=1.;
			matCtrl->SetMatColor(redClr,greenClr,blueClr);
		}
	}
	
    //-------------------------------------------------------
    // Analysis-specific elements
	else if(!myStartElement(xName,attrs))
	{	// look for material property
		if(block==MATERIAL)
		{	inputPtr=matCtrl->InputPointer(xName,input);
			if(inputPtr==NULL)
				ThrowCatErrorMessage("Unrecognized material property was found",xName);
		}
		
		// or an invlid tag
		else if(!strcmp(xName,"JANFEAInput")==0)
			ThrowCatErrorMessage("Unrecognized input element found",xName);
	}
	
	
	// delete tag string
	delete [] xName;
}

// End an element
void CommonReadHandler::endElement(const XMLCh *const uri,const XMLCh *const localname,const XMLCh *const qname)
{
    char *xName=XMLString::transcode(localname);
    
    if(strcmp(xName,"Material")==0 || strcmp(xName,"Thermal")==0 )
    {	block=NO_BLOCK;
    }
	
	else if(strcmp(xName,"Header")==0)
	{	block=NO_BLOCK;
		if(fmobj->np<0)
			throw SAXException("Header did not specify required <Analysis> type.");
	}
	
	else if(strcmp(xName,"Mesh")==0)
	{	if(theNodes!=NULL)
		{	if(!theNodes->SetNodeArray(&mxmin,&mxmax,&mymin,&mymax,&mzmin,&mzmax))
				throw SAXException("Memory error in <NodeList> or no nodes.");
			//delete theNodes;		// delete later, in case needed for limits or resequencing in FEA
		}
		if(theElems!=NULL)
		{	if(!theElems->SetElementArray())
				throw SAXException("Memory error in <ElementList> or no elements.");
			delete theElems;
		}
		block=NO_BLOCK;
	}
	
	else if(strcmp(xName,"ElementList")==0 || strcmp(xName,"NodeList")==0)
	{	block=MESHBLOCK;
	}
	
	else if(strcmp(xName,"DisplacementBCs")==0)
	{   block=GRIDBCHEADER;
	}
	
    else
		myEndElement(xName);
    
    delete [] xName;
    input=NO_INPUT;
    inputPtr=NULL;
}

// Decode block of characters if input!=NO_INPUT
void CommonReadHandler::characters(const XMLCh* const chars,const XMLSize_t length)
{
    if(input==NO_INPUT) return;
    char *xData=XMLString::transcode(chars);
    switch(input)
    {	case TEXT_BLOCK:
			switch(inputID)
			{	case DESCRIPTION:
					fmobj->SetDescription(xData);
					break;
				case CHAR_ARRAY:
					if(inputPtr!=NULL) delete [] inputPtr;
					inputPtr=new char[strlen(xData)+1];
					strcpy(inputPtr,xData);
					// the end element method should assign its pointer to inputPtr
					break;
				default:
					myCharacters(xData,length);
					break;
			}
			break;
        
        case INT_NUM:
            sscanf(xData,"%d",(int *)inputPtr);
            break;
            
        case DOUBLE_NUM:
            sscanf(xData,"%lf",(double *)inputPtr);
			*((double *)inputPtr)*=gScaling;
			gScaling=1.;
            break;
		
		case ANALYSIS_NUM:
            sscanf(xData,"%d",&fmobj->np);
			if(!fmobj->ValidAnalysisType())
				throw SAXException("Unknown <Analysis> type in the <Header>.");
			break;
			
		case NOT_NUM:
			break;
		
        default:
			myCharacters(xData,length);
            break;
    }
    delete [] xData;
}

// Read tag with name *theTag and return string (or NULL if not found)
// Caller must delete [] it
char *CommonReadHandler::ReadTagValue(const char *theTag,const Attributes& attrs)
{
	int i,numAttr=attrs.getLength();
    char *aName,*value;
	
	for(i=0;i<numAttr;i++)
	{   aName=XMLString::transcode(attrs.getLocalName(i));
		if(strcmp(aName,theTag)==0)
		{   value=XMLString::transcode(attrs.getValue(i));
			delete [] aName;
			return value;
		}
		delete [] aName;
	}
	return NULL;
}

// Read tag with name *theTag and return its value as a double
// If not found, return the supplied default value
// Main use is for command with a single numeric tag
double CommonReadHandler::ReadNumericAttribute(const char *theTag,const Attributes& attrs,double defaultValue)
{
	int i,numAttr=attrs.getLength();
    char *aName,*value;
	double attrValue=defaultValue;
	
	for(i=0;i<numAttr;i++)
	{   aName=XMLString::transcode(attrs.getLocalName(i));
		if(strcmp(aName,theTag)==0)
		{   value=XMLString::transcode(attrs.getValue(i));
			sscanf(value,"%lf",&attrValue);
			delete [] aName;
			delete [] value;
			break;
		}
		delete [] aName;
	}
	return attrValue;
}

// read the 'number' attribute in current tag, ignore other attributes
void CommonReadHandler::ReadTagNumber(int *myval,const Attributes& attrs)
{
    int i,aval;
    char *aName,*value;
    int numAttr=attrs.getLength();
	*myval=0;
	
    for(i=0;i<numAttr;i++)
    {   aName=XMLString::transcode(attrs.getLocalName(i));
        if(strcmp(aName,"number")==0)
        {   value=XMLString::transcode(attrs.getValue(i));
            sscanf(value,"%i",&aval);
            *myval=aval;
            delete [] value;
			delete [] aName;
			break;
        }
        delete [] aName;
    }
}

// Read 1st 'units' attribute and return a scaling factor an input quantity
//     to convert to standard analysis units
// The options are for time, length, and velocity units
// If no 'units' found or invalid one found, return 1.
double CommonReadHandler::ReadUnits(const Attributes& attrs,int type)
{
    int i,numAttr=attrs.getLength();
    char *aName,*value;
    double attrScale=1.;
    
    for(i=0;i<numAttr;i++)
    {	aName=XMLString::transcode(attrs.getLocalName(i));
        if(strcmp(aName,"units")==0)
        {   value=XMLString::transcode(attrs.getValue(i));
			switch(type)
			{	case TIME_UNITS:
					// convert to seconds
					if(strcmp(value,"msec")==0)
						attrScale=1.e-3;
					else if(strcmp(value,"ms")==0)
						attrScale=1.e-3;
					else if(strcmp(value,"microsec")==0)
						attrScale=1.e-6;
					break;
					
				case LENGTH_UNITS:
					// convert to mm
 					if(strcmp(value,"m")==0)
						attrScale=1.e3;
					else if(strcmp(value,"cm")==0)
						attrScale=10.;
					else if(strcmp(value,"microns")==0)
						attrScale=1.e-3;
					else if(strcmp(value,"in")==0)
						attrScale=25.4;
					else if(strcmp(value,"ft")==0)
						attrScale=12.*25.4;
					break;
					
				case VELOCITY_UNITS:
					// convert to mm/sec
					if(strcmp(value,"m/sec")==0)
						attrScale=1.e3;
					else if(strcmp(value,"mm/msec")==0)
						attrScale=1.e3;
					else if(strcmp(value,"cm/sec")==0)
						attrScale=10.;
					else if(strcmp(value,"in/sec")==0)
						attrScale=25.4;
					else if(strcmp(value,"ft/sec")==0)
						attrScale=12.*25.4;
					break;
				
				case MASS_UNITS:
					// convert to g
					if(strcmp(value,"kg")==0)
						attrScale=1.e3;
					else if(strcmp(value,"mg")==0)
						attrScale=1.e-3;
					else if(strcmp(value,"lbs")==0)
						attrScale=453.594;
					else if(strcmp(value,"oz")==0)
						attrScale=28.3495;
					break;
					
				default:
					break;
			}
			delete [] aName;
            delete [] value;
            break;
        }
        delete [] aName;
    }
    return attrScale;
}

/********************************************************************************
	 CommonReadHandler: Overrides of the SAX ErrorHandler interface
		and custom error messages
********************************************************************************/

// throw error message "msg (comment)"
void CommonReadHandler::ThrowCatErrorMessage(const char *msg,const char *comment)
{  
	char *errMsg=new char[strlen(msg)+strlen(comment)+5];
	strcpy(errMsg,msg);
	strcat(errMsg," (");
	strcat(errMsg,comment);
	strcat(errMsg,").");
	throw SAXException(errMsg);
}

// throw error message "msg extra (comment)" and extra and/or comment may be empty
void CommonReadHandler::ThrowCompoundErrorMessage(const char *msg,const char *extra,const char *comment)
{  
	int msgLen=strlen(msg)+1;
	if(strlen(extra)>0) msgLen+=strlen(extra)+1;
	char *errMsg=new char[msgLen];
	strcpy(errMsg,msg);
	if(strlen(extra)>0)
	{	strcat(errMsg," ");
		strcat(errMsg,extra);
	}
	if(strlen(comment)>0)
		ThrowCatErrorMessage(errMsg,comment);
	else
		throw SAXException(errMsg);
}

// error if not correct dimension
void CommonReadHandler::ValidateCommand(char *cmdName,int blockNeed,int dimNeed)
{
	// check block location
	if(blockNeed!=NO_BLOCK && block!=blockNeed)
		ThrowCompoundErrorMessage(cmdName,"command found at invalid location","");
		
	// only allowed in 3D calculations
	if(dimNeed==MUST_BE_3D)
	{	if(!fmobj->IsThreeD())
			ThrowCompoundErrorMessage(cmdName,"command only allowed in 3D analyses","");
	}
	
	// only allowed in 2D calculations
	else if(dimNeed==MUST_BE_2D)
	{	if(fmobj->IsThreeD())
			ThrowCompoundErrorMessage(cmdName,"command only allowed in 2D analyses","");
	}
	
}

// overrides SAX error Handler
void CommonReadHandler::error(const SAXParseException& e)
{	cerr << "\nError parsing the XML input file:" << endl;
	cerr << "\nError in file " << StrX(e.getSystemId())
		 << ", at line " << e.getLineNumber()
		 << ", char " << e.getColumnNumber()
         << ": " << StrX(e.getMessage()) << endl;
    throw "\nFix the input file and retry.";
}

// overrides SAX fatalError Handler
void CommonReadHandler::fatalError(const SAXParseException& e)
{	cerr << "\nFatal Error parsing the XML input file:" << endl;
    cerr << "\nFatal Error in file " << StrX(e.getSystemId())
		 << ", at line " << e.getLineNumber()
		 << ", char " << e.getColumnNumber()
         << ": " << StrX(e.getMessage()) << endl;
    throw "\nFix the input file and retry.";
}

// overrides SAX warning Handler
void CommonReadHandler::warning(const SAXParseException& e)
{
    cerr << "\nWarning in file " << StrX(e.getSystemId())
		 << ", at line " << e.getLineNumber()
		 << ", char " << e.getColumnNumber()
         << ": " << StrX(e.getMessage()) << endl;
}


#pragma mark ACCESSORS

/********************************************************************************
	 Read coordinate allowing max or min for current grid
	 After mesh is known
********************************************************************************/

// x value
double CommonReadHandler::ReadX(char *value,double distScaling)
{	return ReadGridPoint(value,distScaling,mxmin,mxmax);
}

// y value
double CommonReadHandler::ReadY(char *value,double distScaling)
{	return ReadGridPoint(value,distScaling,mymin,mymax);
}

// z value
double CommonReadHandler::ReadZ(char *value,double distScaling)
{	return ReadGridPoint(value,distScaling,mzmin,mzmax);
}

// generic value
double CommonReadHandler::ReadGridPoint(char *value,double distScaling,double axisMin,double axisMax)
{	double dval;
	
	// look for max, max+. max-, min, min+, min- optionally with a number
	if(strlen(value)>=3 && value[0]=='m')
	{	if(value[1]=='a' && value[2]=='x')
		{	sscanf(value,"%*3c%lf",&dval);
			if(DbleEqual(dval,0.))
			{	dval=0.;
				if(strlen(value)>3)
				{	if(value[3]=='-')
						dval=-1.;
					else if(value[3]=='+')
						dval=1.;
				}
			}
			return axisMax+dval*ElementBase::GetMinimumCellSize();
		}
		else if(value[1]=='i' && value[2]=='n')
		{	sscanf(value,"%*3c%lf",&dval);
			if(DbleEqual(dval,0.))
			{	dval=0.;
				if(strlen(value)>3)
				{	if(value[3]=='-')
						dval=-1.;
					else if(value[3]=='+')
						dval=1.;
				}
			}
			return axisMin+dval*ElementBase::GetMinimumCellSize();
		}
	}
	
	// just read a number
	sscanf(value,"%lf",&dval);
	return dval*distScaling;
}

#pragma mark ACCESSORS

// decode string delimited by white space, commands, or tabs and add intervening numbers
// to the input vector (which is cleared first). Any other characters trigger an error
bool CommonReadHandler::GetFreeFormatNumbers(char *nData,vector<double> &values,double scaling)
{
	int offset=0,numOffset=0;
	double dval;
	char numstr[200];
	
	// clear input array
	values.clear();
	
	while(nData[offset]!=0)
	{	// skip delimiters
		char nc = nData[offset];
		if(nc==' ' || nc == '\n' || nc == '\r' || nc == '\t' || nc == ':' || nc == ';' || nc == ',')
		{	// decode number now
			if(numOffset)
			{	numstr[numOffset]=0;
				sscanf(numstr,"%lf",&dval);
				values.push_back(dval*scaling);
				numOffset=0;
			}
			offset++;
			continue;
		}
		
		// valid number characters - exit on error
		if((nc<'0' || nc>'9') && nc!='+' && nc!='-' && nc!='e' && nc!='E' && nc!='.')
			return FALSE;
		
		// add to number string
		numstr[numOffset++] = nc;
				
		// go on
		offset++;
	}
	
	// decode last number now
	if(numOffset)
	{	numstr[numOffset]=0;
		sscanf(numstr,"%lf",&dval);
		values.push_back(dval);
	}
	
	return TRUE;
}
		
