/********************************************************************************
    CommonReadHandler.cpp
    NairnFEA
    
    Created by John Nairn on Jan 13 2006.
    Copyright (c) 2006 John A. Nairn. All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "System/UnitsController.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Read_XML/MaterialController.hpp"
#include "Read_XML/NodesController.hpp"
#include "Read_XML/ElementsController.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Exceptions/StrX.hpp"
#include "Read_XML/Expression.hpp"

#pragma mark CommonReadHandler: Constructors and Destructor

// main constuctor
CommonReadHandler::CommonReadHandler()
{
	meshType=UNKNOWN_MESH;
    block=NO_BLOCK;
	matCtrl=new MaterialController();
    gScaling=1.;
    mxmin=mymin=mzmin=0.;
    mxmax=mymax=mzmax=-1.;
};

// destructor
CommonReadHandler::~CommonReadHandler()
{
	delete matCtrl;
}

// set material array before deleting
void CommonReadHandler::FinishUp(void)
{	const char *emsg = matCtrl->SetMaterialArray();
	if(emsg!=NULL)
		throw SAXException(emsg);
}

#pragma mark CommonReadHandler: Methods

// Start a new element
// throws std::bad_alloc, SAXException()
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
	
	// developed flag
	else if(strcmp(xName,"RandSeed")==0)
	{	if(block!=HEADER)
			throw SAXException("<RandSeed> must be within the <Header> element.");
		input=INT_NUM;
		CommonAnalysis *obj = GetCommonAnalysis();
		inputPtr=(char *)&obj->randseed;
	}

	// developed flag
    else if(strcmp(xName,"DevelFlag")==0)
    {	if(block!=HEADER)
			throw SAXException("<DevelFlag> must be within the <Header> element.");
		double flagNumDble=ReadNumericAttribute("Number",attrs,(double)0.0);
    	int flagNum=(int)(flagNumDble+0.5);
		// must be from 0 to NUMBER_DEVELOPMENT_FLAGS-1
		if(flagNum<0 || flagNum>=NUMBER_DEVELOPMENT_FLAGS)
			throw SAXException("The <DevelFlag> 'Number' must be from 0 to 19");
		input=INT_NUM;
		CommonAnalysis *obj = GetCommonAnalysis();
		inputPtr=(char *)&obj->dflag[flagNum];
    }
	
	else if(strcmp(xName,"ConsistentUnits")==0)
	{   if(block!=HEADER)
			throw SAXException("<ConsistentUnits> must be within the <Header> element.");
		char length[10],mass[10],timeu[10];
		strcpy(length,"");
		strcpy(mass,"");
		strcpy(timeu,"");
		char *aName,*value;
    	int i,numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
			if(strlen(value)>9)
				throw SAXException("<ConsistentUnits> length, mass, or time attribute is invalid.");
            if(strcmp(aName,"length")==0)
                strcpy(length,value);
            else if(strcmp(aName,"mass")==0)
                strcpy(mass,value);
            else if(strcmp(aName,"time")==0)
                strcpy(timeu,value);
            delete [] aName;
            delete [] value;
        }
		if(strlen(length)==0 && strlen(mass)==0 && strlen(timeu)==0)
		{	strcpy(length,"L");
			strcpy(mass,"M");
			strcpy(timeu,"T");
		}
		if(!UnitsController::SetConsistentUnits(length, mass, timeu))
			throw SAXException("Duplicated <ConsistentUnits> command or one of the units is invalid.");
    }

	
    //-------------------------------------------------------
    // <Mesh> block
	
    // Node list
    else if(strcmp(xName,"NodeList")==0)
	{	// mesh or could be 3D crack in MPM
		if(block==MESHBLOCK)
		{	if(meshType!=UNKNOWN_MESH)
				throw SAXException("<NodeList> cannot be used with a generated mesh or used a second time for creating a mesh.");
			block=NODELIST;
			meshType=EXPLICIT_MESH;
			if(theNodes==NULL) theNodes=new NodesController();
		}
		else
			throw SAXException("<NodeList> must be within a <Mesh> element.");
    }
        
    // Element list
    else if(strcmp(xName,"ElementList")==0)
	{	// mesh or could be 3D crack in MPM
		if(block==MESHBLOCK)
		{	if(meshType!=EXPLICIT_MESH)
				throw SAXException("<ElementList> cannot be used with a generated mesh.");
			block=ELEMENTLIST;
			if(theElems==NULL) theElems=new ElementsController();
		}
		else
			throw SAXException("<ElementList> must be within a <Mesh> element.");
		
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
    	int i,numAttr=(int)attrs.getLength();
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
			ThrowCatErrorMessage("Invalid material: either undefined type or not allowed for current analysis type",matName);
    }

	// begin a material
    else if(strcmp(xName,"color")==0)
    {	if(block!=MATERIAL)
			throw SAXException("<color> must be within a <Material> definition.");
		float redClr=-1.,greenClr=-1.,blueClr=-1.,alpha=1.;
		char *aName,*value;
    	int i,numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"red")==0)
                sscanf(value,"%f",&redClr);
            else if(strcmp(aName,"green")==0)
                sscanf(value,"%f",&greenClr);
            else if(strcmp(aName,"blue")==0)
                sscanf(value,"%f",&blueClr);
            else if(strcmp(aName,"alpha")==0)
                sscanf(value,"%f",&alpha);
            delete [] aName;
            delete [] value;
        }
		if(redClr>=0.)
		{	if(redClr>1.) redClr=1.;
			if(greenClr<0.) greenClr=redClr;
			if(blueClr<0.) blueClr=redClr;
			if(greenClr>1.) greenClr=1.;
			if(blueClr>1.) blueClr=1.;
            if(alpha<0.) alpha=0.;
            if(alpha>1.) alpha=1.;
			matCtrl->SetMatColor(redClr,greenClr,blueClr,alpha);
		}
	}
	
    //-------------------------------------------------------
    // Analysis-specific elements
	else if(!myStartElement(xName,attrs))
	{	// look for material property
		if(block==MATERIAL)
		{	inputPtr=matCtrl->InputPointer(xName,input,gScaling);
			if(inputPtr==NULL)
			{	char msg[255];
				strcpy(msg,"Unrecognized ");
				strcat(msg,matCtrl->MaterialType());
				strcat(msg," material property was found");
				ThrowCatErrorMessage(msg,xName);
			}
		}
		
		// or an invalid tag
		else if(!(strcmp(xName,"JANFEAInput")==0))
			ThrowCatErrorMessage("Unrecognized input element found",xName);
	}
	
	
	// delete tag string
	delete [] xName;
}

// End an element
// throws SAXException()
void CommonReadHandler::endElement(const XMLCh *const uri,const XMLCh *const localname,const XMLCh *const qname)
{
    char *xName=XMLString::transcode(localname);
    
    if(strcmp(xName,"Material")==0 || strcmp(xName,"Thermal")==0 )
    {	block=NO_BLOCK;
    }
	
	else if(strcmp(xName,"Header")==0)
	{	block=NO_BLOCK;
		CommonAnalysis *obj = GetCommonAnalysis();
		if(obj->np<0)
			throw SAXException("Header did not specify required <Analysis> type.");
	}
	
	else if(strcmp(xName,"Mesh")==0)
	{	// Convert nodes and elements into arrays
		if(theNodes!=NULL)
		{	nd = theNodes->SetNodeArray(nnodes,&mxmin,&mxmax,&mymin,&mymax,&mzmin,&mzmax);;
			if(nd==NULL)
				throw SAXException("Memory error in <NodeList> or no nodes in the list.");
			theNodes = DeleteTheNodes();
		}
		if(theElems!=NULL)
		{	theElements = theElems->SetElementArray(nelems,block==MESHBLOCK);
			if(theElements==NULL)
				throw SAXException("Memory error in <ElementList> or no elements.");
			delete theElems;
			theElems = NULL;
		}
		
		block = NO_BLOCK;
	}
	
	else if(strcmp(xName,"NodeList")==0)
	{	// return to parent block
		block = MESHBLOCK;
	}
	
	else if(strcmp(xName,"ElementList")==0)
	{	// return to parent block
		block = MESHBLOCK;
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
// throws SAXException()
void CommonReadHandler::characters(const XMLCh* const chars,const XMLSize_t length)
{
    if(input==NO_INPUT) return;
    char *xData=XMLString::transcode(chars);
    switch(input)
    {	case TEXT_BLOCK:
			switch(inputID)
			{	case DESCRIPTION:
				{	CommonAnalysis *obj = GetCommonAnalysis();
					obj->SetDescription(xData);
					break;
				}
				case CHAR_ARRAY:
					if(inputPtr!=NULL) delete [] inputPtr;
					inputPtr=new char[strlen(xData)+1];
					strcpy(inputPtr,xData);
					// the end element method should assign its pointer to inputPtr
					break;
				default:
					myCharacters(xData,(int)length);
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
        
        case EXPRESSION_STR:
            // if needed, cannot be NULL, empty, or a duplicate
            if(strlen(xData)==0)
            {   ThrowSAXException("A material expression with zero-length string");
                return;
            }
            if(*((Expression **)inputPtr)!=NULL)
            {    ThrowSAXException("Duplicate expression found for the same material property");
                return;
            }
            *((Expression **)inputPtr) =  Expression::CreateExpression(xData,"Material expression is not a valid function");
            break;
		
		case ANALYSIS_NUM:
		{	CommonAnalysis *obj = GetCommonAnalysis();
			sscanf(xData,"%d",&obj->np);
			if(!obj->ValidAnalysisType())
				throw SAXException("Unknown <Analysis> type in the <Header>.");
			break;
		}

		case NOT_NUM:
			break;
		
        default:
			myCharacters(xData,(int)length);
            break;
    }
    delete [] xData;
}

// Read tag with name *theTag and return string (or NULL if not found)
// Caller must delete [] it
char *CommonReadHandler::ReadTagValue(const char *theTag,const Attributes& attrs)
{
	int i,numAttr=(int)attrs.getLength();
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
	int i,numAttr=(int)attrs.getLength();
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
    int numAttr=(int)attrs.getLength();
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

// Read first 'units' attribute and return a scaling factor an input quantity
//     to convert to standard analysis units
// The options are for time, length, and velocity units
// If no 'units' found or invalid one found, return 1.
double CommonReadHandler::ReadUnits(const Attributes& attrs,int type)
{
    int i,numAttr=(int)attrs.getLength();
    char *aName,*value;
    double attrScale=1.;
    
    for(i=0;i<numAttr;i++)
    {	aName=XMLString::transcode(attrs.getLocalName(i));
        if(strcmp(aName,"units")==0)
        {   value=XMLString::transcode(attrs.getValue(i));
			attrScale = UnitsController::UnitsAttribute(value,type);
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
// throws SAXException()
void CommonReadHandler::ThrowCatErrorMessage(const char *msg,const char *comment)
{  
	char *errMsg=new (nothrow) char[strlen(msg)+strlen(comment)+5];
	if(errMsg != NULL)
	{	strcpy(errMsg,msg);
		strcat(errMsg," (");
		strcat(errMsg,comment);
		strcat(errMsg,").");
		throw SAXException(errMsg);
	}
	else
		throw SAXException(msg);
}

// throw error message "msg extra (comment)" and extra and/or comment may be empty
// throws SAXException()
void CommonReadHandler::ThrowCompoundErrorMessage(const char *msg,const char *extra,const char *comment)
{  
	int msgLen=(int)strlen(msg)+1;
	if(strlen(extra)>0) msgLen+=(int)strlen(extra)+1;
	char *errMsg=new (nothrow) char[msgLen];
	if(errMsg != NULL)
	{	strcpy(errMsg,msg);
		if(strlen(extra)>0)
		{	strcat(errMsg," ");
			strcat(errMsg,extra);
		}
		if(strlen(comment)>0)
			ThrowCatErrorMessage(errMsg,comment);
		else
			throw SAXException(errMsg);
	}
	else
		throw SAXException(msg);
}

// error if not correct dimension
void CommonReadHandler::ValidateCommand(char *cmdName,int blockNeed,int dimNeed)
{
    // check those that require NO_BLOCK
    if(blockNeed==MUST_BE_NO_BLOCK)
    {   if(block!=NO_BLOCK)
            ThrowCompoundErrorMessage(cmdName,"command found at invalid location","");
    }
    
	// check block location
	else if(blockNeed!=NO_BLOCK && block!=blockNeed)
		ThrowCompoundErrorMessage(cmdName,"command found at invalid location","");
		
	// only allowed in 3D calculations
	if(dimNeed==MUST_BE_3D)
	{	CommonAnalysis *obj = GetCommonAnalysis();
		if(!obj->IsThreeD())
			ThrowCompoundErrorMessage(cmdName,"command only allowed in 3D analyses","");
	}
	
	// only allowed in 2D calculations
	else if(dimNeed==MUST_BE_2D)
	{	CommonAnalysis *obj = GetCommonAnalysis();
		if(obj->IsThreeD())
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

// delete nodes controller
NodesController *CommonReadHandler::DeleteTheNodes(void)
{
	delete theNodes;
	return NULL;
}

#pragma mark CommonReadHandler: Accessors

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
	{	dval = 0.;
		
		// max+, max-, max+f, or max-f, where f is an unsigned float
		if(value[1]=='a' && value[2]=='x')
		{	if(strlen(value)==4)
			{	// old style can only be + or -
				if(value[3]=='-')
					dval=-1.;
				else if(value[3]=='+')
					dval=1.;
				else if(value[3]>=0x30 && value[3]<=0x39)
				{	// allow max0 to max9 to be max+#
					dval=(float)(value[3]-0x30);
				}
			}
			else if(strlen(value)>4)
				sscanf(value,"%*3c%lf",&dval);
			return axisMax+dval*ElementBase::GetMinimumCellSize();
		}
		
		// min+, min-, min+f, or min-f, where f is an unsigned float
		else if(value[1]=='i' && value[2]=='n')
		{	if(strlen(value)==4)
			{	// old style can only be + or -
				if(value[3]=='-')
					dval=-1.;
				else if(value[3]=='+')
					dval=1.;
				else if(value[3]>=0x30 && value[3]<=0x39)
				{	// allow min0 to min9 to be min+#
					dval=(float)(value[3]-0x30);
				}
			}
			else  if(strlen(value)>4)
				sscanf(value,"%*3c%lf",&dval);
			return axisMin+dval*ElementBase::GetMinimumCellSize();
		}
	}
	
	// just read a number
	sscanf(value,"%lf",&dval);
	return dval*distScaling;
}

// decode string delimited by white space, puctuation (, : ;), or tabs and add intervening numbers
// to the input vector (which is cleared first). Any other characters trigger an error
// not thread safe due to push_back()
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
		{	cout << " *** Bad character '" << nc << "' (integer code " << (int)nc << ") within free format numbers" << endl;
			return false;
		}
		
		// add to number string
		numstr[numOffset++] = nc;
				
		// go on
		offset++;
	}
	
	// decode last number now
	if(numOffset)
	{	numstr[numOffset]=0;
		sscanf(numstr,"%lf",&dval);
		values.push_back(dval*scaling);
	}
	
	return true;
}
