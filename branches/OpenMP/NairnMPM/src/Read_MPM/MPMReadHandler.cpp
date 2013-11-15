/********************************************************************************
    MPMReadHandler.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 23 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.        
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/ArchiveData.hpp"
#include "Read_MPM/MPMReadHandler.hpp"
#include "Custom_Tasks/CalcJKTask.hpp"
#include "Custom_Tasks/ReverseLoad.hpp"
#include "Custom_Tasks/VTKArchive.hpp"
#include "Custom_Tasks/HistoryArchive.hpp"
#include "Custom_Tasks/AdjustTimeStepTask.hpp"
#include "Custom_Tasks/CarnotCycle.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Global_Quantities/GlobalQuantity.hpp"
#include "MPM_Classes/MatPoint2D.hpp"
#include "MPM_Classes/MatPointAS.hpp"
#include "MPM_Classes/MatPoint3D.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "Boundary_Conditions/MatPtHeatFluxBC.hpp"
#include "Read_XML/NodesController.hpp"
#include "Read_XML/ElementsController.hpp"
#include "Read_XML/ParseController.hpp"
#include "Read_XML/MaterialController.hpp"
#include "Read_MPM/CrackController.hpp"
#include "Read_MPM/MpsController.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Materials/MaterialBase.hpp"
#include "Materials/RigidMaterial.hpp"
#include "Materials/PressureLaw.hpp"
#include "Custom_Tasks/PropagateTask.hpp"
#include "Read_XML/ShapeController.hpp"

// Element types
#include "Elements/FourNodeIsoparam.hpp"

/********************************************************************************
	MPMReadHandler: Constructors and Destructor
********************************************************************************/

MPMReadHandler::MPMReadHandler()
{
	crackCtrl=new CrackController();
    currentTask=NULL;
}

MPMReadHandler::~MPMReadHandler()
{
	// set material propagate
	MaterialBase *obj=(MaterialBase *)matCtrl->firstObject;
	while(obj!=NULL)
		obj=obj->SetFinalPropagate();
		
	// locate first crack
	firstCrack=(CrackHeader *)crackCtrl->firstObject;
	delete crackCtrl;
}

/********************************************************************************
	MPMReadHandler: Constructors and Destructor
********************************************************************************/

// Custom MPM element start
bool MPMReadHandler::myStartElement(char *xName,const Attributes& attrs)
{
    char *aName,*value;
    char quantityName[100];
    unsigned int i,numAttr;
    
    //-------------------------------------------------------
    // <MPMHeader> section
    if(strcmp(xName,"MPMHeader")==0)
    	block=MPMHEADER;
        
    else if(strcmp(xName,"MPMMethod")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=INT_NUM;
        inputPtr=(char *)&fmobj->mpmApproach;
    }
    
    else if(strcmp(xName,"MatlPtsPerElement")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=INT_NUM;
        inputPtr=(char *)&fmobj->ptsPerElement;
    }
    
    else if(strcmp(xName,"TimeFactor")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)fmobj->GetCFLPtr();
    }
    
    else if(strcmp(xName,"MaxTime")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&fmobj->maxtime;
        gScaling=ReadUnits(attrs,TIME_UNITS);
    }
    
    else if(strcmp(xName,"TimeStep")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&timestep;
        gScaling=ReadUnits(attrs,TIME_UNITS);
    }
    
    else if(strcmp(xName,"ArchiveTime")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&archiver->archTime;
        gScaling=ReadUnits(attrs,TIME_UNITS);
    }
    
    else if(strcmp(xName,"FirstArchiveTime")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&archiver->firstArchTime;
        gScaling=ReadUnits(attrs,TIME_UNITS);
    }
    
    else if(strcmp(xName,"GlobalArchiveTime")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&archiver->globalTime;
        gScaling=ReadUnits(attrs,TIME_UNITS);
    }
    
    else if(strcmp(xName,"GlobalArchive")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	numAttr=attrs.getLength();
		int setWhichMat=0;
		quantityName[0]=0;
		char whichMat[200];
		whichMat[0]=0;
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"type")==0)
				strcpy(quantityName,value);
            else if(strcmp(aName,"material")==0 || strcmp(aName,"mat")==0)
                sscanf(value,"%d",&setWhichMat);
			else if(strcmp(aName,"matname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(whichMat,value);
			}
            delete [] aName;
            delete [] value;
        }
		// if gave a matname, it takes precedence over mat number
		if(strlen(whichMat)>0)
			setWhichMat = matCtrl->GetIDFromNewName(whichMat);
		new GlobalQuantity(quantityName,setWhichMat);
    }
	
    else if(strcmp(xName,"Damping")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&bodyFrc.damping;
		bodyFrc.useDamping=TRUE;
		
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"function")==0)
			{	bodyFrc.SetGridDampingFunction(value);
			}
            delete [] aName;
            delete [] value;
        }
    }
    
    else if(strcmp(xName,"FeedbackDamping")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
		bodyFrc.useFeedback=TRUE;
        inputPtr=(char *)&bodyFrc.dampingCoefficient;
		
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"target")==0)
			{	bodyFrc.SetTargetFunction(value);
			}
            else if(strcmp(aName,"max")==0)
            {   double maxAlpha;
                sscanf(value,"%lf",&maxAlpha);
                bodyFrc.SetMaxAlpha(maxAlpha);
            }
            delete [] aName;
            delete [] value;
        }
    }
    
    else if(strcmp(xName,"Diffusion")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		DiffusionTask::active=TRUE;
		DiffusionTask::reference=ReadNumericAttribute("reference",attrs,(double)0.0);
		if(DiffusionTask::reference<0.) DiffusionTask::reference=0.;
		if(DiffusionTask::reference>1.) DiffusionTask::reference=1.;
    }
	
	else if(strcmp(xName,"DefGradTerms")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		input=INT_NUM;
        inputPtr=(char *)&MaterialBase::incrementalDefGradTerms;
	}

    else if(strcmp(xName,"GIMP")==0)
    {   // no attribute or empty implies uGIMP (backward compatibility) or look for key words
		ValidateCommand(xName,MPMHEADER,ANY_DIM);
		ElementBase::useGimp = UNIFORM_GIMP;
		ElementBase::analysisGimp = UNIFORM_GIMP;
		maxShapeNodes = fmobj->np==THREED_MPM ? 28 : 10 ;
        numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"type")==0)
			{	value=XMLString::transcode(attrs.getValue(i));
                if(strcmp(value,"Dirac")==0 || strcmp(value,"Classic")==0 || strcmp(value,"0")==0)
                {   ElementBase::useGimp = POINT_GIMP;
                    ElementBase::analysisGimp = POINT_GIMP;
					maxShapeNodes = fmobj->np==THREED_MPM ? 9 : 5 ;
                }
                else if(strcmp(value,"uGIMP")==0 || strcmp(value,"GIMP")==0 || strcmp(value,"1")==0)
                {   ElementBase::useGimp = UNIFORM_GIMP;
                    ElementBase::analysisGimp = UNIFORM_GIMP;
					maxShapeNodes = fmobj->np==THREED_MPM ? 28 : 10 ;
                }
                else if(strcmp(value,"lCPDI")==0 || strcmp(value,"CPDI")==0 || strcmp(value,"2")==0)
                {   ElementBase::useGimp = LINEAR_CPDI;
                    ElementBase::analysisGimp = LINEAR_CPDI;
					maxShapeNodes = fmobj->np==THREED_MPM ? 40 : 17 ;		// 3D could need 65
                }
                else if(strcmp(value,"qCPDI")==0 || strcmp(value,"3")==0)
                {   ElementBase::useGimp = QUADRATIC_CPDI;
                    ElementBase::analysisGimp = QUADRATIC_CPDI;
					maxShapeNodes = fmobj->np==THREED_MPM ? 40 : 17 ;		// 3D not allowed for qCPDI
                }
                else
                    throw SAXException("GIMP type must be Classic, uGIMP, lCPDI, or qCPDI.");
				delete [] value;
			}
			delete [] aName;
        }
    }

    else if(strcmp(xName,"MultiMaterialMode")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		fmobj->multiMaterialMode=true;
		block=MULTIMATERIAL;
        numAttr=attrs.getLength();
		double scanInput;
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
			sscanf(value,"%lf",&scanInput);
            if(strcmp(aName,"Vmin")==0)
			{	contact.materialContactVmin=scanInput;
			}
            else if(strcmp(aName,"Dcheck")==0)
            {   if(scanInput<0.5 && scanInput>-0.5)
                    contact.displacementCheck=FALSE;                // 0
				else if(scanInput<1.5)
					contact.displacementCheck=TRUE;                 // 1  
				else
                    throw SAXException("Dcheck attribute on MultiMaterialMode must be 0 or 1.");
			}
            else if(strcmp(aName,"Normals")==0)
 			{	if(scanInput<0.5 && scanInput>-0.5)
					contact.materialNormalMethod=MAXIMUM_VOLUME_GRADIENT;           // 0 - MAXG
				else if(scanInput<1.5)
					contact.materialNormalMethod=MAXIMUM_VOLUME;                    // 1 - MAXV
				else if(scanInput<2.5)
					contact.materialNormalMethod=AVERAGE_MAT_VOLUME_GRADIENTS;      // 2 - AVGG
				else if(scanInput<3.5)
					contact.materialNormalMethod=EACH_MATERIALS_MASS_GRADIENT;		// 3 - OWNG
                else
                    throw SAXException("Normals attribute on MultiMaterialMode must be 0 to 3.");
			}
            else if(strcmp(aName,"RigidBias")==0)
			{	contact.rigidGradientBias=scanInput;
			}
			delete [] aName;
            delete [] value;
        }
    }
	
    else if(strcmp(xName,"ArchiveRoot")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
        input=TEXT_BLOCK;
		inputID=ARCHIVEROOT_NAME;
        numAttr=attrs.getLength();
		int scanInt;
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"unique")==0)
			{	value=XMLString::transcode(attrs.getValue(i));
				sscanf(value,"%d",&scanInt);
				if(scanInt==1) inputID=UNIQUE_ARCHIVEROOT_NAME;
				delete [] value;
			}
			delete [] aName;
        }
    }
    
    else if(strcmp(xName,"MPMArchiveOrder")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
        input=MPMORDER_BLOCK;
    }
    
    else if(strcmp(xName,"CrackArchiveOrder")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
        input=CRACKORDER_BLOCK;
    }

    else if(strcmp(xName,"StressFreeTemp")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&thermal.reference;
    }

    else if(strcmp(xName,"LeaveLimit")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=INT_NUM;
        inputPtr=(char *)&fmobj->warnParticleLeftGrid;
    }

    // Cracks in MPM Header
    else if(strcmp(xName,"Cracks")==0)
	{	ValidateCommand(xName,MPMHEADER,MUST_BE_2D);
        block=CRACKHEADER;
    }
    
    else if(strcmp(xName,"Friction")==0)
    {	if(block!=CRACKHEADER && block!=MULTIMATERIAL && block!=MATERIAL)
			ThrowCompoundErrorMessage(xName," command found at invalid location.","");
    	input=DOUBLE_NUM;
		double *theDn,*theDnc,*theDt;
		char othername[200];
		othername[0]=0;
		if(block==CRACKHEADER)
        {   // Commmand in the <Cracks> element of the <MPMHeader>
			inputPtr=(char *)&contact.friction;
			theDn=&contact.Dn;
			theDnc=&contact.Dnc;
			theDt=&contact.Dt;
		}
		else if(block==MULTIMATERIAL)
        {   // Commmand in the <MultiMaterialMode> element of the <MPMHeader>
			inputPtr=(char *)&contact.materialFriction;
			theDn=&contact.materialDn;
			theDnc=&contact.materialDnc;
			theDt=&contact.materialDt;
		}
		else
        {   // Command in a material to set custom material-to-material contact
			matCtrl->friction=0.;
			inputPtr=(char *)&(matCtrl->friction);
			matCtrl->Dn=-1.;
			theDn=&(matCtrl->Dn);
			matCtrl->Dnc=-101.;
			theDnc=&(matCtrl->Dnc);
			matCtrl->Dt=-1.;
			theDt=&(matCtrl->Dt);
			matCtrl->otherMatID=-1;
		}
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"Dn")==0 || strcmp(aName,"Dnt")==0)
			{	sscanf(value,"%lf",theDn);
			}
            else if(strcmp(aName,"Dnc")==0)
			{	sscanf(value,"%lf",theDnc);
			}
			else if(strcmp(aName,"Dt")==0)
			{	sscanf(value,"%lf",theDt);
			}
			else if(strcmp(aName,"mat")==0)
            {   // only relevant in material definition
				sscanf(value,"%d",&(matCtrl->otherMatID));
			}
			else if(strcmp(aName,"matname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(othername,value);
			}
			delete [] value;
            delete [] aName;
        }
		// if gave a mat name, it takes precedence over mat number (only relevant in material definition
		if(strlen(othername)>0)
			matCtrl->otherMatID = matCtrl->GetIDFromNewName(othername);
        
        // the inputPtr above will read in the friction coefficient (or numerical flag for other types of contact
        
        // the end of of friction element will transfer the material controller variable to the material
    }
    
    else if(strcmp(xName,"Propagate")==0 || strcmp(xName,"AltPropagate")==0)
    {	if(block!=CRACKHEADER && block!=MATERIAL)
			ThrowCompoundErrorMessage(xName," command found at invalid location.","");
		ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
    	numAttr=attrs.getLength();
		
		// get which to set
		int setIndex = strcmp(xName,"Propagate")==0 ? 0 : 1 ;
		
		// read attributes
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"criterion")==0)
            {	value=XMLString::transcode(attrs.getValue(i));
                if(block==CRACKHEADER)
                    sscanf(value,"%d",&fmobj->propagate[setIndex]);
                else
				{	int matCriterion;
                    sscanf(value,"%d",&matCriterion);
					matCtrl->SetCriterion(matCriterion,setIndex);
				}
				delete [] value;
            }
            else if(strcmp(aName,"direction")==0)
            {	value=XMLString::transcode(attrs.getValue(i));
                if(block==CRACKHEADER)
                    sscanf(value,"%d",&fmobj->propagateDirection[setIndex]);
                else
				{	int matDirection;
                    sscanf(value,"%d",&matDirection);
					matCtrl->SetDirection(matDirection,setIndex);
				}
				delete [] value;
            }
            else if(strcmp(aName,"traction")==0)
            {	value=XMLString::transcode(attrs.getValue(i));
                if(block==CRACKHEADER)
                    sscanf(value,"%d",&fmobj->propagateMat[setIndex]);
                else
				{	int tractionMat;
                    sscanf(value,"%d",&tractionMat);
					matCtrl->SetTractionMat(tractionMat,setIndex);
				}
				delete [] value;
            }
            delete [] aName;
        }
    }
    
    else if(strcmp(xName,"PropagateLength")==0)
	{	ValidateCommand(xName,CRACKHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&PropagateTask::cellsPerPropagationStep;
    }
    
    else if(strcmp(xName,"MovePlane")==0)
	{	ValidateCommand(xName,CRACKHEADER,MUST_BE_2D);
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"type")==0)
			{	if(strcmp(value,"cm")==0 || strcmp(value,"CM")==0)
					contact.SetMoveOnlySurfaces(FALSE);
				else
					contact.SetMoveOnlySurfaces(TRUE);				
            }
			else if(strcmp(aName,"prevent")==0)
			{	if(strcmp(value,"no")==0 || strcmp(value,"No")==0 || strcmp(value,"NO")==0)
					contact.SetPreventPlaneCrosses(FALSE);
				else
					contact.SetPreventPlaneCrosses(TRUE);				
            }
			delete [] value;
            delete [] aName;
        }
    }

    else if(strcmp(xName,"ContactPosition")==0)
    {	if(block!=CRACKHEADER && block!=MULTIMATERIAL)
			ThrowCompoundErrorMessage(xName," command found at invalid location.","");
		mpmgrid.SetContactByDisplacements(FALSE);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&mpmgrid.positionCutoff;
    }
	
    else if(strcmp(xName,"JContour")==0)
	{	ValidateCommand(xName,CRACKHEADER,MUST_BE_2D);
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"type")==0)
                sscanf(value,"%d",&JContourType);
            else if(strcmp(aName,"size")==0)
                sscanf(value,"%d",&JGridSize);
            else if(strcmp(aName,"terms")==0)
                sscanf(value,"%d",&JTerms);
            delete [] aName;
            delete [] value;
        }
    }

    //-------------------------------------------------------
    // <Mesh> section
	
	// Start Mesh section
    else if(strcmp(xName,"Mesh")==0)
	{	ValidateCommand(xName,NO_BLOCK,ANY_DIM);
		block=MESHBLOCK;
		archiver->SetArchiveMesh(FALSE);
		value=ReadTagValue("output",attrs);
		if(value!=NULL)
		{	if(strcmp(value,"file")==0) archiver->SetArchiveMesh(TRUE);
			delete [] value;
		}
	}
	
    else if(strcmp(xName,"elem")==0)
	{	ValidateCommand(xName,ELEMENTLIST,ANY_DIM);
    	input=NODE_BLOCK;
		value=ReadTagValue("type",attrs);
		if(value==NULL)
            throw SAXException("<elem> does not specify element type.");
		if(!theElems->SetElemIDStr(value))
            throw SAXException("Invalid or incompatible element type.");
        delete [] value;
    }
    
    //-------------------------------------------------------
    // Material points block
    else if(strcmp(xName,"MaterialPoints")==0)
	{	ValidateCommand(xName,NO_BLOCK,ANY_DIM);
    	if(nmpms!=0)
            throw SAXException("<MaterialPoints> section not allowed after material points allocated.");
		if(nelems==0)
            throw SAXException("<MaterialPoints> must come after <Mesh> creation.");
        block=POINTSBLOCK;
		mpCtrl=new MpsController();
    }
	
	else if(strcmp(xName,"MatlPts")==0)
	{   throw SAXException("<MatlPts> should be replaced by <PointList>.");
	}

    else if(strcmp(xName,"PointList")==0)
	{	ValidateCommand(xName,POINTSBLOCK,ANY_DIM);
    	block=MATLPTS;
    }
	
	// <GRIDBLOCK is to distiguish from pt used in generators called below (must be below)
    else if((strcmp(xName,"pt")==0 || strcmp(xName,"vel")==0) && (block<GRIDBLOCK))
	{	ValidateCommand(xName,NO_BLOCK,ANY_DIM);
		double aScaling;
		if(strcmp(xName,"pt")==0)
		{	aScaling=ReadUnits(attrs,LENGTH_UNITS);
			// <Intensity> may have <vel>, but cannot have <pt>
			if(block==INTENSITYBLOCK)
				throw SAXException("<pt> found at invalid location.");
		}
		else
		{	aScaling=ReadUnits(attrs,VELOCITY_UNITS);
			// <NodeList> and <CrackList> may have <pt>, but may not have <vel>
			if(block==NODELIST || block==CRACKLIST)
				throw SAXException("<vel> found at invalid location.");
		}
		Vector xp;
		ZeroVector(&xp);
        int tipMatnum=-1;
		int matid=0;
    	numAttr=attrs.getLength();
		char matname[200];
		matname[0]=0;
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"x")==0)
            	sscanf(value,"%lf",&xp.x);
            else if(strcmp(aName,"y")==0)
            	sscanf(value,"%lf",&xp.y);
            else if(strcmp(aName,"z")==0)
            	sscanf(value,"%lf",&xp.z);
            else if(strcmp(aName,"tip")==0)
            	sscanf(value,"%d",&tipMatnum);		// <CrackList> only
            else if(strcmp(aName,"mat")==0)
            	sscanf(value,"%d",&matid);			// <CrackList> only, 1-based ID or following name use
			else if(strcmp(aName,"matname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(matname,value);
			}
            delete [] aName;
            delete [] value;
        }
		// if gave a matname, it takes precedence over mat number
		if(strlen(matname)>0)
			matid = matCtrl->GetIDFromNewName(matname);
		ScaleVector(&xp,aScaling);
        
        if(block==NODELIST)
			theNodes->AddNode(xp.x,xp.y,xp.z);
        else if(block==MATLPTS)
		{	if(!mpCtrl->SetPtOrVel(xName,&xp))
				throw SAXException("A material point in an <mp> command is not within the grid.");
        }
        else if(block==CRACKLIST)
		{	if(!crackCtrl->AddSegment(new CrackSegment(xp.x,xp.y,tipMatnum,matid)))
                throw SAXException("Crack not in the mesh or out of memory adding a crack segment.");
        }
		else if(block==INTENSITYBLOCK)
			SetLevelVelocity(xp.x,xp.y,xp.z);
        else
            throw SAXException("<pt> or <vel> found at invalid location.");
    }
    
    else if(strcmp(xName,"mp")==0)
	{	ValidateCommand(xName,MATLPTS,ANY_DIM);
		int elemNum=1;
    	int matl=0;
        double angle=0.,pConcentration=0.,dval;
		double thick=mpmgrid.GetDefaultThickness();
		double pTempInitial=thermal.reference;
    	numAttr=attrs.getLength();
		char mpmat[200];
		mpmat[0]=0;
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            sscanf(value,"%lf",&dval);
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"elem")==0)
			{	elemNum=(int)(dval+0.5);
				if(elemNum<1 || elemNum>nelems) elemNum=1;
			}
            else if(strcmp(aName,"matl")==0 || strcmp(aName,"mat")==0)
                matl=(int)(dval+0.5);
			else if(strcmp(aName,"matname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(mpmat,value);
			}
            else if(strcmp(aName,"thick")==0)
                thick=dval;
            else if(strcmp(aName,"conc")==0)
			{	pConcentration=dval;
				if(pConcentration<0. || pConcentration>1.)
					throw SAXException("Material point concentration potential must be from 0 to 1");
			}
            else if(strcmp(aName,"wtconc")==0)
			{	pConcentration=-dval;
				if(pConcentration>0.)
					throw SAXException("Material point weight fraction concentration must be >= 0");
			}
            else if(strcmp(aName,"temp")==0)
                pTempInitial=dval;
            else
                angle=dval;
            delete [] aName;
            delete [] value;
        }
		// if gave a matname, it takes precedence over mat number
		if(strlen(mpmat)>0)
			matl = matCtrl->GetIDFromNewName(mpmat);
        
        // create object for 3D or 2D analysis
		if(fmobj->IsThreeD())
		{	mpCtrl->AddMaterialPoint(new MatPoint3D(elemNum,matl,angle),pConcentration,pTempInitial);
		}
		else if(fmobj->IsAxisymmetric())
		{	// thickness set to x position in pt command by SetPtOrVel() when it calls SetOrigin()
			mpCtrl->AddMaterialPoint(new MatPointAS(elemNum,matl,angle,1.),pConcentration,pTempInitial);
		}
		else
		{	mpCtrl->AddMaterialPoint(new MatPoint2D(elemNum,matl,angle,thick),pConcentration,pTempInitial);
		}
    }
 
	// particle mass (usually when reading points from a file)
    else if(strcmp(xName,"mass")==0)
	{	ValidateCommand(xName,MATLPTS,ANY_DIM);
		double aScaling=ReadUnits(attrs,MASS_UNITS);
		double pmass=-1.;
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"m")==0)
			{   value=XMLString::transcode(attrs.getValue(i));
            	sscanf(value,"%lf",&pmass);
				delete [] value;
			}
            delete [] aName;
        }
		pmass*=aScaling;
		if(pmass>0.) mpCtrl->SetPtMass(pmass);
    }

    //-------------------------------------------------------
    // Crack List
    else if(strcmp(xName,"CrackList")==0)
	{	ValidateCommand(xName,NO_BLOCK,MUST_BE_2D);
    	block=CRACKLIST;
		CrackHeader *newCrack=new CrackHeader();
		crackCtrl->AddCrack(newCrack);
		newCrack->SetContact(contact.friction,contact.Dn,contact.Dnc,contact.Dt);
		double gridThickness=mpmgrid.GetThickness();
		if(gridThickness>0.) newCrack->SetThickness(gridThickness);
		
		// read crack attributes
		double dval;
		numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
			if(strcmp(aName,"type")==0)
			{	if(strcmp(value,"fixed")==0)
					newCrack->SetFixedCrack(TRUE);
				else
					newCrack->SetFixedCrack(FALSE);
			}
			else if(strcmp(aName,"friction")==0)
            {	sscanf(value,"%lf",&dval);
				newCrack->SetFriction(dval);
			}
            else if(strcmp(aName,"Dn")==0 || strcmp(aName,"Dnt")==0)
            {	sscanf(value,"%lf",&dval);
				newCrack->SetDn(dval);
			}
            else if(strcmp(aName,"Dnc")==0)
            {	sscanf(value,"%lf",&dval);
				newCrack->SetDnc(dval);
			}
            else if(strcmp(aName,"Dt")==0)
            {	sscanf(value,"%lf",&dval);
				newCrack->SetDt(dval);
			}
            delete [] aName;
            delete [] value;
        }
    }

	// crack thickness (not needed if grid thickness correctly set in structure grid command)
    else if(strcmp(xName,"Thickness")==0 && block==CRACKLIST)
	{	ValidateCommand(xName,CRACKLIST,MUST_BE_2D);
    	input=DOUBLE_NUM;
		CrackHeader *newCrack=(CrackHeader *)crackCtrl->currentObject();
        inputPtr=(char *)newCrack->GetThicknessPtr();
        gScaling=ReadUnits(attrs,LENGTH_UNITS);
    }
    //-----------------------------------------------------------
	// <GridBCs> section
	
    // Prepares to read grid-based BCs
    else if(strcmp(xName,"GridBCs")==0)
	{	ValidateCommand(xName,NO_BLOCK,ANY_DIM);
	    block=GRIDBCHEADER;
		velocityBCs=new ParseController();
		concBCs=new ParseController();
		tempBCs=new ParseController();
    }
 
	// a velocity BC
    else if(strcmp(xName,"fix")==0)
	{	ValidateCommand(xName,FIXEDNODES,ANY_DIM);
    	int node=0;
		int dof=0,style=CONSTANT_VALUE,velID=0;;
        double ftime=0.,angle=0.,angle2=0.;
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"node")==0)
            	sscanf(value,"%d",&node);
            else if(strcmp(aName,"dof")==0)
            	sscanf(value,"%d",&dof);
            else if(strcmp(aName,"style")==0)
            	sscanf(value,"%d",&style);
            else if(strcmp(aName,"time")==0)
                sscanf(value,"%lf",&ftime);
            else if(strcmp(aName,"angle")==0)
                sscanf(value,"%lf",&angle);
            else if(strcmp(aName,"angle2")==0)
                sscanf(value,"%lf",&angle2);
            else if(strcmp(aName,"id")==0)
                sscanf(value,"%d",&velID);
            delete [] aName;
            delete [] value;
        }
		if(fmobj->IsThreeD())
		{	if(dof!=1 && dof!=2 && dof!=3 && dof!=12 && dof!=13 && dof!=23 && dof!=123)
				throw SAXException("'dir' in fix element must be 1, 2, 3, 12, 13, 23, or 123 for 3D analyses.");
		}
        else if(dof!=1 && dof!=2 && dof!=12)
            throw SAXException("'dir' in fix element must be 1, 2, or 12 for 2D analyses.");
        if(velID>0)
            throw SAXException("'id' for velocity boundary conditions must be <= 0.");
        
        // create object and get input
        // note that dof is input style (1,2,3,12,13,23, or 123)
        NodalVelBC *newVelBC=new NodalVelBC(node,dof,style,(double)0.,ftime,angle,angle2);
        newVelBC->SetID(velID);
		velocityBCs->AddObject(newVelBC);
		
		if(style!=FUNCTION_VALUE)
		{	input=DOUBLE_NUM;
			inputPtr=(char *)&newVelBC->value;
		}
		else
		{	input=FUNCTION_BLOCK;
			inputPtr=(char *)newVelBC;
		}
    }

    //-----------------------------------------------------------
    // <ParticleBCs> section
	
	// Prepare to read particle BCs
    else if(strcmp(xName,"ParticleBCs")==0)
	{	ValidateCommand(xName,NO_BLOCK,ANY_DIM);
		block=PARTICLEBCHEADER;
		mpLoadCtrl=new ParseController();
		mpTractionCtrl=new ParseController();
		mpConcFluxCtrl=new ParseController();
		mpHeatFluxCtrl=new ParseController();
    }

    // Loads on material points
    else if(strcmp(xName,"LoadBCs")==0)
	{	ValidateCommand(xName,PARTICLEBCHEADER,ANY_DIM);
	    block=LOADEDPOINTS;
    }
	
	// A particle load
    else if(strcmp(xName,"load")==0)
	{	ValidateCommand(xName,LOADEDPOINTS,ANY_DIM);
    	int ptNum=0;
		int dof=0,style=1;
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   value=XMLString::transcode(attrs.getValue(i));
            aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"pt")==0)
				sscanf(value,"%d",&ptNum);
            else if(strcmp(aName,"style")==0)
				sscanf(value,"%d",&style);
            else if(strcmp(aName,"dof")==0)
				sscanf(value,"%d",&dof);
            delete [] aName;
            delete [] value;
        }
		if(fmobj->IsThreeD())
		{	if(dof<1 || dof>3)
				throw SAXException("'dir' in <load> must be 1, 2, or 3 for 3D analyses.");
		}
        else if(dof>2 || dof<1)
            throw SAXException("'dir' in <load> must be 1 or 2 for 2D analyses.");
        
        // create object and prepare for data
		mpLoadCtrl->AddObject(new MatPtLoadBC(ptNum,dof,style));
		if(style!=FUNCTION_VALUE)
			input=BC_BLOCK;
		else
		{	input=FUNCTION_BLOCK;
			inputPtr=(char *)mpLoadCtrl->lastObject;
		}
    }

    //-------------------------------------------------------
    // Mesh, material point, BC, and more generator commands
    else if(GenerateInput(xName,attrs))
    {
    }

    //-------------------------------------------------------
    // BMP file generator commands
    else if(BMPFileInput(xName,attrs))
    {
    }

    //-------------------------------------------------------
    // <Thermal> section
	
	// thermal ramp
    else if(strcmp(xName,"Isothermal")==0)
	{	ValidateCommand(xName,THERMAL,ANY_DIM);
		thermal.Activate();
    	input=DOUBLE_NUM;
        inputPtr=(char *)&thermal.isoDeltaT;
		numAttr=attrs.getLength();
		for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
			if(strcmp(aName,"time")==0)
				sscanf(value,"%lf",&thermal.isoRampTime);
			else if(strcmp(aName,"start")==0)
				sscanf(value,"%lf",&thermal.rampStart);
			delete [] aName;
			delete [] value;
		}
    }
    
	// Turn on conduction analysis
    else if(strcmp(xName,"Conduction")==0)
	{	ValidateCommand(xName,THERMAL,ANY_DIM);
		ConductionTask::active=TRUE;
	}

	// Turn on crack tip heating
    else if(strcmp(xName,"CrackTipHeating")==0)
	{	ValidateCommand(xName,THERMAL,MUST_BE_2D);
		ConductionTask::crackTipHeating=TRUE;
	}
	
	// Turn on crack tip heating
    else if(strcmp(xName,"CrackContactHeating")==0)
	{	ValidateCommand(xName,THERMAL,MUST_BE_2D);
		ConductionTask::crackContactHeating=TRUE;
	}
	
	// Turn on crack tip heating
    else if(strcmp(xName,"ContactHeating")==0)
	{	ValidateCommand(xName,THERMAL,ANY_DIM);
		ConductionTask::matContactHeating=TRUE;
	}
	
	// Turn on adiabatic mode
    else if(strcmp(xName,"EnergyCoupling")==0)
	{	ValidateCommand(xName,THERMAL,ANY_DIM);
		ConductionTask::energyCoupling=TRUE;
	}
	
	// Turn off artificial viscosity heating (only on when adibatic)
    else if(strcmp(xName,"NoAVHeating")==0)
	{	ValidateCommand(xName,THERMAL,ANY_DIM);
		ConductionTask::AVHeating=FALSE;
	}
    //-------------------------------------------------------
    // <Gravity> section
	
	// begin Gravity section
    else if(strcmp(xName,"Gravity")==0)
	{	ValidateCommand(xName,NO_BLOCK,ANY_DIM);
    	block=GRAVITY;
		bodyFrc.Activate();
    }
    
	// Gravity in X direction
    else if(strcmp(xName,"BodyXForce")==0)
	{	ValidateCommand(xName,GRAVITY,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&bodyFrc.gforce.x;
    }
    
	// Gravity in Y direction
    else if(strcmp(xName,"BodyYForce")==0)
	{	ValidateCommand(xName,GRAVITY,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&bodyFrc.gforce.y;
    }

	// Gravity in Y direction
    else if(strcmp(xName,"BodyZForce")==0)
	{	ValidateCommand(xName,GRAVITY,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&bodyFrc.gforce.z;
    }

    //-------------------------------------------------------
    // <CustomTasks> section
	
	// begin CustomTasks section
    else if(strcmp(xName,"CustomTasks")==0)
    	block=CUSTOMTASKS;

	// Schedule a custom task
    else if(strcmp(xName,"Schedule")==0)
	{	CustomTask *nextTask=NULL;
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"name")==0)
            {   value=XMLString::transcode(attrs.getValue(i));
                if(strcmp(value,"ReverseLoad")==0)
                {   nextTask=(CustomTask *)(new ReverseLoad());
                    if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
                }
				else if(strcmp(value,"VTKArchive")==0)
                {   nextTask=(CustomTask *)(new VTKArchive());
                    if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
                }
				else if(strcmp(value,"HistoryArchive")==0)
                {   nextTask=(CustomTask *)(new HistoryArchive());
                    if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
                }
				else if(strcmp(value,"AdjustTimeStep")==0)
                {   nextTask=(CustomTask *)(new AdjustTimeStepTask());
                    if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
                }
				else if(strcmp(value,"CarnotCycle")==0)
                {   nextTask=(CustomTask *)(new CarnotCycle());
                    if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
                }
                else
                    throw SAXException("Unknown custom task requested for scheduling.");
                delete [] value;
            }
            delete [] aName;
            if(nextTask!=NULL) break;
        }
        if(nextTask!=NULL)
        {   block=TASKPARAMETERS;
            if(currentTask==NULL)
                theTasks=nextTask;
            else
                ((CustomTask *)currentTask)->nextTask=nextTask;
            currentTask=(char *)nextTask;
        }
    }
    
	// Set Custom Task parameter
    else if(strcmp(xName,"Parameter")==0)
	{	inputPtr=NULL;
    	numAttr=attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"name")==0)
            {   value=XMLString::transcode(attrs.getValue(i));
                inputPtr=((CustomTask *)currentTask)->InputParam(value,input);
                delete [] value;
            }
            delete [] aName;
            if(inputPtr!=NULL) break;
        }
        if(inputPtr==NULL)
            throw SAXException("Unrecognized task parameter was found.");
    }

	//-------------------------------------------------------
	// Unknown element
	else
		return FALSE;
    
    return TRUE;
}

// End an element
void MPMReadHandler::myEndElement(char *xName)
{
    if(strcmp(xName,"CustomTasks")==0 || strcmp(xName,"Gravity")==0)
    {	block=NO_BLOCK;
    }
    
    else if(strcmp(xName,"CrackList")==0)
    {   crackCtrl->FinishCrack();
        block = NO_BLOCK;
    }
	
	else if(strcmp(xName,"MPMHeader")==0)
	{	block=NO_BLOCK;
	
		// create new folder now, to be sure it works
		if(archiver->GetArchiveRoot()==NULL)
			throw SAXException("No archive root file name was defined.");
		else if(!archiver->MakeArchiveFolder())
			throw SAXException("The system('mkdir') command failed");
	}
	
	else if(strcmp(xName,"GridBCs")==0)
	{	firstVelocityBC=(NodalVelBC *)velocityBCs->firstObject;
		lastVelocityBC=(NodalVelBC *)velocityBCs->lastObject;
		firstConcBC=(NodalConcBC *)concBCs->firstObject;
		lastConcBC=(NodalConcBC *)concBCs->lastObject;
		firstTempBC=(NodalTempBC *)tempBCs->firstObject;
		lastTempBC=(NodalTempBC *)tempBCs->lastObject;
		delete velocityBCs;
		delete concBCs;
		delete tempBCs;
	}
	
	else if(strcmp(xName,"ParticleBCs")==0)
	{	firstLoadedPt=(MatPtLoadBC *)mpLoadCtrl->firstObject;
		firstTractionPt=(MatPtTractionBC *)mpTractionCtrl->firstObject;
		firstFluxPt=(MatPtFluxBC *)mpConcFluxCtrl->firstObject;
		firstHeatFluxPt=(MatPtHeatFluxBC *)mpHeatFluxCtrl->firstObject;
		delete mpLoadCtrl;
		delete mpTractionCtrl;
		delete mpConcFluxCtrl;
		delete mpHeatFluxCtrl;
	}
	
    else if(strcmp(xName,"MaterialPoints")==0)
    {	block=NO_BLOCK;
		if(!mpCtrl->SetMPArray())
			throw SAXException("Memory error or no material points were defined in <MaterialPoints> block");
		delete mpCtrl;
    }
    
	else if(strcmp(xName,"DisplacementBCs")==0)
	{   block=GRIDBCHEADER;
	}
	
    else if(strcmp(xName,"PointList")==0)
    {	block=POINTSBLOCK;
    }
    
    else if(strcmp(xName,"Cracks")==0)
    {	block=MPMHEADER;
    }
    
    else if(strcmp(xName,"MultiMaterialMode")==0)
    {	block=MPMHEADER;
    }
    
    else if(strcmp(xName,"Schedule")==0)
    {	block=CUSTOMTASKS;
    }
    
    else if(strcmp(xName,"Friction")==0)
    {	if(block==MATERIAL)
			matCtrl->SetMaterialFriction();
    }
    
    else if(EndGenerator(xName))
    {   
    }
    
    else if(EndBMPInput(xName,POINTSBLOCK))
    {   
    }
	
	else if(strcmp(xName,"JANFEAInput")==0)
	{	// all done
		CreateSymmetryBCs();
	}
}

// Decode block of characters if input!=NO_INPUT
void MPMReadHandler::myCharacters(char *xData,const unsigned int length)
{
	MatPtLoadBC *newBC;
    
    switch(input)
	{	case TEXT_BLOCK:
			switch(inputID)
			{	case ARCHIVEROOT_NAME:
					if(!archiver->SetArchiveRoot(xData,false))
						throw SAXException("Cannot set two <ArchiveRoot> paths.");
					break;
				case UNIQUE_ARCHIVEROOT_NAME:
					if(!archiver->SetArchiveRoot(xData,true))
						throw SAXException("Cannot set two <ArchiveRoot> paths.");
					break;
				default:
					break;
			}
			break;
			
        case NODE_BLOCK:
			if(!theElems->CreateElement(xData))
				throw SAXException("Unknown element type was found in <elem> command.");
            break;
        
        case BC_BLOCK:
			newBC=(MatPtLoadBC *)mpLoadCtrl->lastObject;
            sscanf(xData,"%lf,%lf",&newBC->value,&newBC->ftime);
            break;
		
		case FUNCTION_BLOCK:
			((BoundaryCondition *)inputPtr)->SetFunction(xData);
			break;
        
		case SETTING_FUNCTION_BLOCK:
		case SETTING_FUNCTION2_BLOCK:
		case SETTING_FUNCTION3_BLOCK:
		case VALUE_FUNCTION_BLOCK:
			((RigidMaterial *)inputPtr)->SetSettingFunction(xData,input);
			break;
		
		case STRESS_FUNCTION_BLOCK:
			((PressureLaw *)inputPtr)->SetStressFunction(xData);
			break;
        
        case MPMORDER_BLOCK:
            // mpm archive order only
			archiver->SetMPMOrder(xData);
            break;
        
        case CRACKORDER_BLOCK:
            // crack archive order only
			archiver->SetCrackOrder(xData);
            break;
		
		case POLYHEDRON_BLOCK:
			// must be in active body controller
			theShape->SetProperty(xData,this);
			break;
        
        case HARDENING_LAW_SELECTION:
            ((MaterialBase *)inputPtr)->SetHardeningLaw(xData);
            break;
            
        default:
            break;
    }
}
