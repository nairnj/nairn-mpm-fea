/********************************************************************************
    MPMReadHandler.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 23 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.        
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/ArchiveData.hpp"
#include "Read_MPM/MPMReadHandler.hpp"
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
#include "Materials/TaitLiquid.hpp"
#include "Materials/PressureLaw.hpp"
#include "Custom_Tasks/PropagateTask.hpp"
#include "Read_XML/ShapeController.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "System/UnitsController.hpp"
#include "Materials/ContactLaw.hpp"
#include "Elements/FourNodeIsoparam.hpp"
#include "Boundary_Conditions/InitialCondition.hpp"
#include "Custom_Tasks/PhaseFieldDiffusion.hpp"

int cracksDim = MUST_BE_2D;

#pragma mark MPMReadHandler: Constructors and Destructor

// throws std::bad_alloc
MPMReadHandler::MPMReadHandler()
{
	crackCtrl=new CrackController();
    currentTask=NULL;
	currentContact=NULL;
}

MPMReadHandler::~MPMReadHandler()
{
	// set material propagate
	MaterialBase *obj=(MaterialBase *)matCtrl->firstObject;
	while(obj!=NULL)
		obj=obj->SetFinalPropagate();
		
	// locate first crack, count them, and populate an array
	firstCrack=(CrackHeader *)crackCtrl->firstObject;
    crackCtrl->SetCracksArray();
	delete crackCtrl;
}

#pragma mark MPMReadHandler: Methods

// Return NairnFEA class object
CommonAnalysis *MPMReadHandler::GetCommonAnalysis(void) { return fmobj; }

// Custom MPM element start
// throws std::bad_alloc, SAXException()
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
	{	// XML expects total points per element, the end element code converts to linear
		if(block==MPMHEADER)
		{	input=INT_NUM;
        	inputPtr=(char *)&fmobj->ptsPerElement;
		}
		else if(block==BODYPART)
		{	if(fmobj->customPtsPerElement<0)
				ThrowSAXException("MatlPtsPerElement cannot be used in Hole blocks");
			input=INT_NUM;
			inputPtr=(char *)&fmobj->customPtsPerElement;
		}
        else if(block==BMPBLOCK)
        {   input=INT_NUM;
            inputPtr=(char *)&bmpCustomPtsPerElement;      // but supplied per element
        }
		else
			ThrowCompoundErrorMessage(xName,"command found at invalid location","");
    }
	
	else if(strcmp(xName,"Timing")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
        numAttr=(int)attrs.getLength();
		gScaling=ReadUnits(attrs,SEC_UNITS);
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"CFL")==0)
			{	sscanf(value,"%lf",fmobj->GetCFLPtr());
            }
			else if(strcmp(aName,"TransCFL")==0)
			{	sscanf(value,"%lf",fmobj->GetTransCFLPtr());
			}
            else if(strcmp(aName,"step")==0)
			{	sscanf(value,"%lf",&timestep);
				timestep *= gScaling;					// convert to sec
            }
            else if(strcmp(aName,"max")==0)
			{	sscanf(value,"%lf",&fmobj->maxtime);
				fmobj->maxtime *= gScaling;			// convert to sec
            }
            delete [] aName;
            delete [] value;
		}
    }
    
    else if(strcmp(xName,"TimeFactor")==0)
	{	// Deprecated - use Timing instead
		ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)fmobj->GetCFLPtr();
    }
    
	else if(strcmp(xName,"TransTimeFactor")==0)
	{	// Deprecated - use Timing instead
		ValidateCommand(xName,MPMHEADER,ANY_DIM);
		input=DOUBLE_NUM;
		inputPtr=(char *)fmobj->GetTransCFLPtr();
	}
	
    else if(strcmp(xName,"MaxTime")==0)
	{	// Deprecated - use Timing instead
		ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&fmobj->maxtime;
        gScaling=ReadUnits(attrs,SEC_UNITS);
    }
    
    else if(strcmp(xName,"TimeStep")==0)
	{	// Deprecated - use Timing instead
		ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&timestep;
        gScaling=ReadUnits(attrs,SEC_UNITS);
    }
    
#ifdef RESTART_OPTION
    else if(strcmp(xName,"RestartScaling")==0)
    {   // If 1.e-6 < fabs(val) < 1, restart if acceleration too high with new time step
        // If >0 warn once, but then silent. If <0 print each time restared
        ValidateCommand(xName,MPMHEADER,ANY_DIM);
        input=DOUBLE_NUM;
        inputPtr=(char *)&fmobj->restartScaling;
        numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"CFL")==0)
            {   value=XMLString::transcode(attrs.getValue(i));
                sscanf(value,"%lf",&fmobj->restartCFL);
                delete [] value;
            }
            delete [] aName;
        }
    }
#else
    else if(strcmp(xName,"RestartScaling")==0)
    {    throw SAXException("<RestartScaling> requires code compiled with RESTART_OPTION");
    }
#endif
    
	else if(strcmp(xName,"ArchiveGroup")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
        numAttr=(int)attrs.getLength();
		gScaling=ReadUnits(attrs,SEC_UNITS);
		double *atptr = archiver->GetArchTimePtr();			// creates the next group
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"time")==0)
			{	sscanf(value,"%lf",atptr);
				*atptr *= gScaling;					// convert to sec
            }
            else if(strcmp(aName,"start")==0)
			{	double *ftptr = archiver->GetFirstArchTimePtr();
				sscanf(value,"%lf",ftptr);
				*ftptr *= gScaling;					// convert to sec
            }
            else if(strcmp(aName,"maxProps")==0)
            {	int maxProps;
				sscanf(value,"%d",&maxProps);
				archiver->SetMaxiumPropagations(maxProps);
            }
            delete [] aName;
            delete [] value;
		}
    }
	
    else if(strcmp(xName,"ArchiveTime")==0)
	{	// Deprecated - use ArchiveGroup instead
		ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
		inputPtr=(char *)archiver->GetArchTimePtr();
        gScaling=ReadUnits(attrs,SEC_UNITS);
        numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"maxProps")==0)
            {	int maxProps;
				sscanf(value,"%d",&maxProps);
				archiver->SetMaxiumPropagations(maxProps);
            }
            delete [] aName;
            delete [] value;
		}
    }
    
    else if(strcmp(xName,"FirstArchiveTime")==0)
	{	// Deprecated - use ArchiveGroup instead
		ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)archiver->GetFirstArchTimePtr();
        gScaling=ReadUnits(attrs,SEC_UNITS);
    }
    
    else if(strcmp(xName,"GlobalArchiveTime")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        inputPtr=(char *)archiver->GetGlobalTimePtr();
		gScaling=ReadUnits(attrs,SEC_UNITS);
    }
    
    else if(strcmp(xName,"GlobalArchive")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	numAttr=(int)attrs.getLength();
		int setWhichMat=0;
		quantityName[0]=0;
		char whichMat[200];
		whichMat[0]=0;
		Vector *ptLoc = NULL;
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
			else if(strcmp(aName,"pt")==0)
			{	if(ptLoc!=NULL)
					ThrowCompoundErrorMessage(xName,"has more than one pts attribute","");
				vector<double> pts;
				if(!CommonReadHandler::GetFreeFormatNumbers(value,pts,1.0))
					ThrowCompoundErrorMessage(xName,"pts attribute has invalid number format","");
				int numpts = (int)pts.size();
				ptLoc = new Vector;
				if(numpts==2)
					*ptLoc = MakeVector(pts[0],pts[1],0.);
				else if(numpts==3)
					*ptLoc = MakeVector(pts[0],pts[1],pts[2]);
				else
					ThrowCompoundErrorMessage(xName,"pts attribute requires 2 or three numbers","");
			}
            delete [] aName;
            delete [] value;
        }
		// if gave a matname, it takes precedence over mat number
		if(ptLoc!=NULL)
		{	new GlobalQuantity(quantityName,ptLoc);
		}
		else
		{	if(strlen(whichMat)>0)
				setWhichMat = matCtrl->GetIDFromNewName(whichMat);
			new GlobalQuantity(quantityName,setWhichMat);
		}
    }
	
	// Damping in MPMHEADER only, PDamping in MPMHEADER or in MATERIAL
    else if(strcmp(xName,"Damping")==0 || strcmp(xName,"PDamping")==0)
	{	input=DOUBLE_NUM;
        bool gridDamp = true,matDamp = false;
        if(strcmp(xName,"Damping")==0)
		{	if(block!=MPMHEADER)
				ThrowCompoundErrorMessage(xName,"command found at invalid location","");
            inputPtr=(char *)&bodyFrc.damping;
            bodyFrc.useDamping = true;
        }
        else if(block==MATERIAL)
		{   // material properties
			inputPtr=(char *)&(matCtrl->matPdamping);
			matCtrl->matPdamping = -1.e12;
			matDamp = true;
		}
		else
        {	// particle damping in the header
			if(block!=MPMHEADER)
				ThrowCompoundErrorMessage(xName,"command found at invalid location","");
			inputPtr=(char *)&bodyFrc.pdamping;
            bodyFrc.usePDamping = true;
            gridDamp = false;
        }
		
        numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"function")==0)
			{	if(matDamp)
					ThrowSAXException("Material damping function of time is not allowed");
				else
					bodyFrc.SetGridDampingFunction(value,gridDamp);
			}
            else if(strcmp(aName,"PIC")==0)
				ThrowSAXException("PIC attribute on damping remove; use PeriodicXPIC task intead");
            delete [] aName;
            delete [] value;
        }
    }
    
    else if(strcmp(xName,"FeedbackDamping")==0 || strcmp(xName,"PFeedbackDamping")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=DOUBLE_NUM;
        bool gridDamp = true;
        if(strcmp(xName,"FeedbackDamping")==0)
		{	// on the grid
			bodyFrc.useFeedback = true;
            inputPtr=(char *)&bodyFrc.dampingCoefficient;
            bodyFrc.useDamping = true;
        }
        else
		{	// for particles
			bodyFrc.usePFeedback = true;
            inputPtr=(char *)&bodyFrc.pdampingCoefficient;
            gridDamp = false;
            bodyFrc.usePDamping = true;
        }
		
        numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"target")==0)
			{	bodyFrc.SetTargetFunction(value,gridDamp);
			}
            else if(strcmp(aName,"max")==0)
            {   double maxAlpha;
                sscanf(value,"%lf",&maxAlpha);
                bodyFrc.SetMaxAlpha(maxAlpha,gridDamp);
            }
            delete [] aName;
            delete [] value;
        }
    }

	// XPIC option - get order (1=PIC, 2 is XPIC, <1 invalid)
    else if(strcmp(xName,"XPIC")==0)
	{	throw SAXException("<XPIC> command no longer allowed; use PeriodicXPIC custom task.");
	}

    else if(strcmp(xName,"Diffusion")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		
		// read attributes to get style and reference
		numAttr=(int)attrs.getLength();
		double refCon=0.;
		int diffStyle=MOISTURE_DIFFUSION;
        int noLimit=0;
		for(i=0;i<numAttr;i++)
		{   aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
			if(strcmp(aName,"reference")==0)
				sscanf(value,"%lf",&refCon);
			else if(strcmp(aName,"style")==0)
				sscanf(value,"%d",&diffStyle);
            else if(strcmp(aName,"nolimit")==0)
                sscanf(value,"%d",&noLimit);
			delete [] aName;
			delete [] value;
		}
		
		if(diffStyle==MOISTURE_DIFFUSION)
		{	if(fmobj->HasDiffusion())
				throw SAXException("Only one 'Diffusion' task allowed in an analysis.");
#ifdef POROELASTICITY
			else if(fmobj->HasPoroelasticity())
				throw SAXException("Cannot use both 'Diffusion' and 'Poroelasticity' in the same analysis.");
#endif
			// create task
			diffusion = new DiffusionTask(fmax(refCon,0.),(double)noLimit,MOISTURE_DIFFUSION);
		}
        else if(diffStyle==POROELASTICITY_DIFFUSION)
        {
#ifdef POROELASTICITY
            if(fmobj->HasPoroelasticity())
                throw SAXException("Only one 'Poroelasticity' command allowed in an analysis.");
            else if(fmobj->HasDiffusion())
                throw SAXException("Cannot use both 'Poroelasticity' and 'Diffusion' in the same analysis.");
            
            // Reference pressure: legacy units are MPa - convert to Pa and must be positive
            double refPressure = UnitsController::Scaling(1.e6)*refCon;
            
            // Viscosity: Legacy units are cP - 1 cP = 0.001 Pa-sec
            // Default: .001 Pa-sec in Legacy and 1 in CU
            double ppViscosity = UnitsController::Scaling(1.e-3)*ReadNumericAttribute("viscosity",attrs,(double)1.0);
            
            // create task
            diffusion = new DiffusionTask(refPressure,ppViscosity,POROELASTICITY_DIFFUSION);
#else
            throw SAXException("Unrecognized 'Diffusion' style.");
#endif
        }
		else if(diffStyle<NUM_DUFFUSION_OPTIONS)
		{	// create task for fracture, battery, or conduction modeling
            PhaseFieldDiffusion *phaseTask = new PhaseFieldDiffusion(diffStyle);
			
			// tack on to other diffustion tasks list
			if(otherDiffusion==NULL)
				otherDiffusion = phaseTask;
			else
			{	// insert at end of the list
				TransportTask *currentTask = (TransportTask *)otherDiffusion;
				while(currentTask!=NULL)
				{	if(currentTask->nextTask==NULL)
					{	currentTask->nextTask = phaseTask;
						break;
					}
					currentTask = currentTask->GetNextTransportTask();
				}
			}
		}
		else
		{	throw SAXException("An unrecognized 'Diffusion' style was requested.");
		}
    }
	
#ifdef POROELASTICITY
	else if(strcmp(xName,"Poroelasticity")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		if(fmobj->HasPoroelasticity())
			throw SAXException("Only one 'Poroelasticity' command allowed in an analysis.");
		else if(fmobj->HasDiffusion())
			throw SAXException("Cannot use both 'Poroelasticity' and 'Diffusion' in the same analysis.");
		
		// Reference pressure: legacy units are MPa - convert to Pa and must be positive
		double refPressure = UnitsController::Scaling(1.e6)*ReadNumericAttribute("reference",attrs,(double)0.0);
		
		// Viscosity: Legacy units are cP - 1 cP = 0.001 Pa-sec
		// Default: .001 Pa-sec in Legacy and 1 in CU
		double ppViscosity = UnitsController::Scaling(1.e-3)*ReadNumericAttribute("viscosity",attrs,(double)1.0);
		
		// create task
		diffusion = new DiffusionTask(refPressure,ppViscosity,POROELASTICITY_DIFFUSION);
	}
#else
	else if(strcmp(xName,"Poroelasticity")==0)
	{	throw SAXException("<Poroelasticity> command requires POROELASTICITY defined in MPMPrefix.hpp.");
	}
#endif
	
	else if(strcmp(xName,"DefGradTerms")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		input=INT_NUM;
        inputPtr=(char *)&MaterialBase::incrementalDefGradTerms;
	}

	else if(strcmp(xName,"ExtrapolateRigid")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		MaterialBase::extrapolateRigidBCs = true;
	}

 	else if(strcmp(xName,"SkipPostExtrapolation")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		fmobj->skipPostExtrapolation = true;
	}
	
	else if(strcmp(xName,"TrackParticleSpin")==0 || strcmp(xName,"TrackGradV")==0)
	{
		throw SAXException("<TrackParticleSpin> requires OSParticulas compiled with ADD_PARTICLE_SPIN");
	}

	else if(strcmp(xName,"TransportOnly")==0)
	{	// future may want to programmatically active this mode
		throw SAXException("<TransportOnly> setting requires OSParticulas compiled with TRANSPORT_ONLY");
	}
	
	else if(strcmp(xName,"NeedsMechanics")==0)
	{	// future may want to programmatically active this mode
	}

	else if(strcmp(xName,"ExactTractions")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		fmobj->exactTractions = true;
	}

	else if(strcmp(xName,"GIMP")==0)
    {   // no attribute or empty implies uGIMP (backward compatibility) or look for key words
		ValidateCommand(xName,MPMHEADER,ANY_DIM);
		
		// default to uGIMP
		ElementBase::useGimp = UNIFORM_GIMP;
		if(fmobj->np==THREED_MPM)
		{	ElementBase::gridNiNodes = 8;
			maxShapeNodes = 28;
		}
		else
		{	ElementBase::gridNiNodes = 4;
			maxShapeNodes = 10;
		}
		
		// look at attributes
        numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            if(strcmp(aName,"type")==0)
			{	value=XMLString::transcode(attrs.getValue(i));
                if(strcmp(value,"Dirac")==0 || strcmp(value,"Classic")==0 || strcmp(value,"0")==0)
                {   ElementBase::useGimp = POINT_GIMP;
					maxShapeNodes = fmobj->np==THREED_MPM ? 9 : 5 ;
					ElementBase::gridNiNodes = fmobj->np==THREED_MPM ? 8 : 4 ;
                }
                else if(strcmp(value,"uGIMP")==0 || strcmp(value,"GIMP")==0 || strcmp(value,"1")==0)
                {   ElementBase::useGimp = UNIFORM_GIMP;
					maxShapeNodes = fmobj->np==THREED_MPM ? 28 : 10 ;
					ElementBase::gridNiNodes = fmobj->np==THREED_MPM ? 8 : 4 ;
                }
                else if(strcmp(value,"lCPDI")==0 || strcmp(value,"CPDI")==0 || strcmp(value,"2")==0)
                {   ElementBase::useGimp = LINEAR_CPDI;
					maxShapeNodes = fmobj->np==THREED_MPM ? 50 : 17 ;
					ElementBase::gridNiNodes = fmobj->np==THREED_MPM ? 8 : 4 ;
                }
                else if(strcmp(value,"qCPDI")==0 || strcmp(value,"3")==0)
                {   ElementBase::useGimp = QUADRATIC_CPDI;
					maxShapeNodes = fmobj->np==THREED_MPM ? 40 : 37 ;
					ElementBase::gridNiNodes = fmobj->np==THREED_MPM ? 8 : 4 ;
                }
				else if(strcmp(value,"Finite")==0 || strcmp(value,"4")==0)
				{   ElementBase::useGimp = FINITE_GIMP;
					maxShapeNodes = fmobj->np==THREED_MPM ? 40 : 20 ;		// 3D not allowed yet
					ElementBase::gridNiNodes = fmobj->np==THREED_MPM ? 8 : 4 ;
				}
				else if(strcmp(value,"B2GIMP")==0 || strcmp(value,"5")==0)
				{   ElementBase::useGimp = BSPLINE_GIMP;
					maxShapeNodes = fmobj->np==THREED_MPM ? 65 : 17 ;
					ElementBase::gridNiNodes = fmobj->np==THREED_MPM ? 27 : 9 ;
				}
				else if(strcmp(value,"B2SPLINE")==0 || strcmp(value,"6")==0)
				{   ElementBase::useGimp = BSPLINE;
					maxShapeNodes = fmobj->np==THREED_MPM ? 28 : 10 ;
					ElementBase::gridNiNodes = fmobj->np==THREED_MPM ? 27 : 9 ;
				}
				else if(strcmp(value,"B2CPDI")==0 || strcmp(value,"7")==0)
				{   ElementBase::useGimp = BSPLINE_CPDI;
					maxShapeNodes = fmobj->np==THREED_MPM ? 100 : 37 ;
					ElementBase::gridNiNodes = fmobj->np==THREED_MPM ? 27 : 9 ;
				}
                else
                    throw SAXException("GIMP type must be Classic, uGIMP, lCPDI, qCPDI, B2SPLINE, B2GIMP, B2CPDI, or Finite");
				delete [] value;
			}
			delete [] aName;
        }
    }

	else if(strcmp(xName,"CPDIrcrit")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		input=DOUBLE_NUM;
		inputPtr=(char *)&ElementBase::rcrit;
	}

    else if(strcmp(xName,"CPDIrcritdynamic")==0)
    {
        throw SAXException("<CPDIrcritdynamic> requires OSParticulas compiled with RCRIT_VERIFICATION");
    }

	else if(strcmp(xName,"MultiMaterialMode")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		fmobj->multiMaterialMode = true;
		block=MULTIMATERIAL;
        numAttr=(int)attrs.getLength();
		double scanInput,polarAngle=0.,azimuthAngle=0.;
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
			sscanf(value,"%lf",&scanInput);
            if(strcmp(aName,"Normals")==0)
 			{	if(scanInput<0.5 && scanInput>-0.5)
					mpmgrid.materialNormalMethod=MAXIMUM_VOLUME_GRADIENT;           // 0 - MAXG
				else if(scanInput<1.5)
					mpmgrid.materialNormalMethod=MAXIMUM_VOLUME;                    // 1 - MAXV
				else if(scanInput<2.5)
					mpmgrid.materialNormalMethod=AVERAGE_MAT_VOLUME_GRADIENTS;      // 2 - AVGG
				else if(scanInput<3.5)
					mpmgrid.materialNormalMethod=EACH_MATERIALS_MASS_GRADIENT;		// 3 - OWNG
				else if(scanInput<4.5)
					mpmgrid.materialNormalMethod=SPECIFIED_NORMAL;					// 4 - SN
				else if(scanInput<5.5)
					mpmgrid.materialNormalMethod=LINEAR_REGRESSION;					// 5 - LINR
				else if(scanInput<6.5)
					mpmgrid.materialNormalMethod=LOGISTIC_REGRESSION;				// 6 - LOGR
                else
                    throw SAXException("Normals attribute on MultiMaterialMode must be 0 to 6.");
			}
            else if(strcmp(aName,"RigidBias")==0)
			{	mpmgrid.rigidGradientBias=scanInput;
			}
            else if(strcmp(aName,"Polar")==0)
			{	polarAngle=scanInput;
			}
            else if(strcmp(aName,"Azimuth")==0)
			{	azimuthAngle=scanInput;
			}
			else if(strcmp(aName,"Lumping")==0)
			{	mpmgrid.lumpingMethod=(int)(scanInput+0.5);
			}
			// Note old Vmin and Dcheck are now ignored
			delete [] aName;
            delete [] value;
        }
		
		// set normal if specified
		if(mpmgrid.materialNormalMethod==SPECIFIED_NORMAL)
			mpmgrid.SetContactNormal(polarAngle,azimuthAngle);
    }

    else if(strcmp(xName,"ArchiveRoot")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
        input=TEXT_BLOCK;
		inputID=ARCHIVEROOT_NAME;
        numAttr=(int)attrs.getLength();
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
    
	else if(strcmp(xName,"MPMArchive")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		int archByte = 0;				// 0 fails, <0 crack archive, >0 mpm archive
		char setting = 'Y';
		int historyNum = 0;				// 1 to 4 on history#
		numAttr=(int)attrs.getLength();
		for(i=0;i<numAttr;i++)
		{   aName=XMLString::transcode(attrs.getLocalName(i));
			if(strcmp(aName,"result")==0)
			{	value=XMLString::transcode(attrs.getValue(i));
				if(CIstrcmp(value,"velocity")==0)
					archByte = ARCH_Velocity;
				else if(CIstrcmp(value,"stress")==0)
					archByte = ARCH_Stress;
				else if(CIstrcmp(value,"strain")==0)
					archByte = ARCH_Strain;
				else if(CIstrcmp(value,"rotstrain")==0)
					archByte = ARCH_RotStrain;
				else if(CIstrcmp(value,"plasticstrain")==0)
					archByte = ARCH_PlasticStrain;
				else if(CIstrcmp(value,"strainenergy")==0)
					archByte = ARCH_StrainEnergy;
				else if(CIstrcmp(value,"workenergy")==0)
					archByte = ARCH_WorkEnergy;
				else if(CIstrcmp(value,"plasticenergy")==0)
					archByte = ARCH_PlasticEnergy;
				else if(CIstrcmp(value,"heatenergy")==0)
					archByte = ARCH_HeatEnergy;
				else if(CIstrcmp(value,"lp")==0)
					archByte = ARCH_SpinMomentum;
				else if(CIstrcmp(value,"wp")==0)
					archByte = ARCH_SpinVelocity;
				else if(CIstrcmp(value,"temperature")==0)
					archByte = ARCH_DeltaTemp;
				else if(CIstrcmp(value,"concentration")==0 || CIstrcmp(value,"porepressure")==0)
					archByte = ARCH_Concentration;
				else if( CIstrcmp(value,"history1")==0 || CIstrcmp(value,"history2")==0
							|| CIstrcmp(value,"history3")==0 || CIstrcmp(value,"history4")==0 )
				{	archByte = ARCH_History;
					historyNum = value[7]&0x0F;				// bit 1,2,3,4 (e.g., 0x31 to 0x01)
				}
				else if( CIstrcmp(value,"history5")==0 || CIstrcmp(value,"history6")==0 || CIstrcmp(value,"history7")==0
						|| CIstrcmp(value,"history8")==0 || CIstrcmp(value,"history9")==0 )
				{	archByte = ARCH_History59;
					historyNum = value[7]&0x0F - 4;			// bit 1,2,3,4,5 (e.g., 0x35 to 0x05)
				}
				else if( CIstrcmp(value,"history10")==0 || CIstrcmp(value,"history11")==0 || CIstrcmp(value,"history12")==0
						|| CIstrcmp(value,"history13")==0 || CIstrcmp(value,"history14")==0 )
				{	archByte = ARCH_History1014;
					historyNum = value[8]&0x0F + 1;				// bit 1,2,3,4,5 (e.g., 0x31 to 0x01)
				}
				else if( CIstrcmp(value,"history15")==0 || CIstrcmp(value,"history16")==0 || CIstrcmp(value,"history17")==0
						|| CIstrcmp(value,"history18")==0 || CIstrcmp(value,"history19")==0 )
				{	archByte = ARCH_History1519;
					historyNum = value[8]&0x0F - 4;			// bit 1,2,3,4,5 (e.g., 0x35 to 0x05)
				}
				else if(CIstrcmp(value,"elementcrossings")==0)
					archByte = ARCH_ElementCrossings;
				else if(CIstrcmp(value,"damagenormal")==0)
					archByte = ARCH_DamageNormal;
				else if(CIstrcmp(value,"size")==0)
					archByte = ARCH_Size;

				else if(CIstrcmp(value,"jintegral")==0)
					archByte = -ARCH_JIntegral;
				else if(CIstrcmp(value,"stressintensity")==0)
					archByte = -ARCH_StressIntensity;
				else if(CIstrcmp(value,"energybalance")==0 || CIstrcmp(value,"czmdisp")==0)
					archByte = -ARCH_CZMDISP;
                else if( CIstrcmp(value,"traction1")==0 || CIstrcmp(value,"traction2")==0 || CIstrcmp(value,"traction3")==0
                        || CIstrcmp(value,"traction4")==0 || CIstrcmp(value,"traction5")==0 )
                {   archByte = ARCH_Traction15;
                    historyNum = value[8]&0x0F;            // bit 1,2,3,4,5 (e.g., 0x31 to 0x01, etc.)
                }
                else if( CIstrcmp(value,"traction6")==0 || CIstrcmp(value,"traction7")==0 || CIstrcmp(value,"traction8")==0
                        || CIstrcmp(value,"traction9")==0 || CIstrcmp(value,"traction10")==0 )
                {   archByte = ARCH_Traction610;
                    historyNum = value[8]&0x0F-5;            // bit 1,2,3,4,5 (e.g., 0x36 to 0x06 to 0x01, etc.)
                    if(historyNum<0) historyNum = 5;
                }
				delete [] value;
			}
			else if(strcmp(aName,"setting")==0)
			{	value=XMLString::transcode(attrs.getValue(i));
				if(strlen(value)>0) setting = value[0];
				delete [] value;
			}
			delete [] aName;
		}
		
		// set a byte
		if(archByte==0)
			throw SAXException("<MPMArchive> command has missing on unrecognized 'result' attribute.");
		else if(archByte == ARCH_History || archByte == ARCH_History59 || archByte == ARCH_History1014 || archByte == ARCH_History1519)
		{	// get current byte and switch 'Y' and 'N' to 1 bit or no bits
			char history = archiver->GetMPMOrderByte(archByte);
			if(history=='Y')
				history = archByte == ARCH_History ? 0x31 : 0x21 ;
			else if(history=='N')
				history = archByte == ARCH_History ? 0x30 : 0x20 ;
			
			// convert bit number from above to Hex value (1, 2, 4, 8, or 16 for 1 to 5)
			char historyBit = 1 << (historyNum-1);
			
			// Add or remove this new bit (remove only if was set)
			if(setting=='Y')
				history |= historyBit;
			else if(history&historyBit)
				history -= historyBit;
			
			// convert new setting of no bits to 'N', 1 bit left as is instead of 'Y'
			if(archByte == ARCH_History)
			{	if(history == 0x30) history = 'N';
			}
			else if(history == 0x20)
				history = 'N';
			
			// replace previous character with new one
			archiver->SetMPMOrderByte(archByte,history);
		}
        else if(archByte == ARCH_Traction15 || archByte == ARCH_Traction610)
        {    // get current byte and switch 'Y' and 'N' to 1 bit or no bits
            char history = archiver->GetCrackOrderByte(archByte);
            if(history=='Y')
                history = 0x21 ;
            else if(history=='N')
                history = 0x20 ;
            
            // convert bit number from above to Hex value (1, 2, 4, 8, or 16 for 1 to 5)
            char historyBit = 1 << (historyNum-1);
            
            // Add or remove this new bit (remove only if was set)
            if(setting=='Y')
                history |= historyBit;
            else if(history&historyBit)
                history -= historyBit;
            
            // convert new setting of no bits to 'N', 1 bit left as is instead of 'Y'
            if(history == 0x20) history = 'N';
            
            // replace previous character with new one
            archiver->SetCrackOrderByte(archByte,history);
        }
		else if(archByte>0)
			archiver->SetMPMOrderByte(archByte,setting);
		else
			archiver->SetCrackOrderByte(-archByte,setting);
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
	{	// <0 to allow delete that many, >0 to push back then many, 0 (default) set to 1% of particles
		ValidateCommand(xName,MPMHEADER,ANY_DIM);
    	input=INT_NUM;
        inputPtr=(char *)&fmobj->warnParticleLeftGrid;
    }

    else if(strcmp(xName,"DeleteLimit")==0)
    {	// >1 deletes and continues until that many nan particles, <=1 aborts on first nan partiles
		ValidateCommand(xName,MPMHEADER,ANY_DIM);
        input=INT_NUM;
        inputPtr=(char *)&fmobj->warnParticleDeleted;
    }
	
	// Custom patching in MPM Header
	else if(strcmp(xName,"PatchGrid")==0)
	{	ValidateCommand(xName,MPMHEADER,ANY_DIM);
		block=PATCHGRID;
		int xnum=1,ynum=1,znum=1;
		numAttr=(int)attrs.getLength();
		for(i=0;i<numAttr;i++)
		{   aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
			if(strcmp(aName,"x")==0)
				sscanf(value,"%d",&xnum);
			else if(strcmp(aName,"y")==0)
				sscanf(value,"%d",&ynum);
			else if(strcmp(aName,"z")==0)
				sscanf(value,"%d",&znum);
			delete [] value;
			delete [] aName;
		}
		mpmgrid.SetCustomPatching(xnum,ynum,znum);
	}
    
	else if(strcmp(xName,"Xpatches")==0)
	{	ValidateCommand(xName,PATCHGRID,ANY_DIM);
		mpmgrid.SetSizeDirection(1);
		input = PATCH_SIZES;
	}
	
	else if(strcmp(xName,"Ypatches")==0)
	{	ValidateCommand(xName,PATCHGRID,ANY_DIM);
		mpmgrid.SetSizeDirection(2);
		input = PATCH_SIZES;
	}
	
	else if(strcmp(xName,"Zpatches")==0)
	{	ValidateCommand(xName,PATCHGRID,ANY_DIM);
		mpmgrid.SetSizeDirection(3);
		input = PATCH_SIZES;
	}
	
    // Cracks in MPM Header
    else if(strcmp(xName,"Cracks")==0)
	{	ValidateCommand(xName,MPMHEADER,cracksDim);
        block=CRACKHEADER;
    }
    
    else if(strcmp(xName,"Friction")==0)
    {	if(block!=CRACKHEADER && block!=MULTIMATERIAL && block!=MATERIAL)
			ThrowCompoundErrorMessage(xName," command found at invalid location.","");
		if(currentContact!=NULL)
			ThrowCompoundErrorMessage(xName," command nested in another <Friction> command.","");
		
		char tempName[20];
		strcpy(tempName,"Input Only");
		currentContact = new ContactLaw(tempName,CONTACTLAW);
    	input=DOUBLE_NUM;
		inputPtr=(char *)&(currentContact->contactProps.friction);
		char othername[200],lawname[200];
		othername[0]=0;
		lawname[0]=0;
		
    	numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"Dn")==0 || strcmp(aName,"Dnt")==0)
			{	sscanf(value,"%lf",&(currentContact->contactProps.Dn));
				currentContact->contactProps.Dn *= UnitsController::Scaling(1.e6);
			}
            else if(strcmp(aName,"Dnc")==0)
			{	sscanf(value,"%lf",&(currentContact->contactProps.Dnc));
				currentContact->contactProps.Dnc *= UnitsController::Scaling(1.e6);
			}
			else if(strcmp(aName,"Dt")==0)
			{	sscanf(value,"%lf",&(currentContact->contactProps.Dt));
				currentContact->contactProps.Dt *= UnitsController::Scaling(1.e6);
			}
			else if(strcmp(aName,"mat")==0)
            {   // only relevant in material definition
				sscanf(value,"%d",&(currentContact->contactProps.otherMatID));
			}
			else if(strcmp(aName,"law")==0)
			{	sscanf(value,"%d",&(currentContact->contactProps.contactLawID));
			}
			else if(strcmp(aName,"matname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(othername,value);
			}
			else if(strcmp(aName,"lawname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(lawname,value);
			}
			delete [] value;
            delete [] aName;
        }
		
		// if gave a mat name, it takes precedence over law number
		if(strlen(lawname)>0)
			currentContact->contactProps.contactLawID = matCtrl->GetIDFromNewName(lawname);
		
		// if gave a mat name, it takes precedence over mat number (only relevant in material definition)
		if(strlen(othername)>0)
			currentContact->contactProps.otherMatID = matCtrl->GetIDFromNewName(othername);
        
        // the inputPtr above will read in the friction coefficient (or numerical flag for other types of contact
		// only used in old style when law name is not used
        
        // the end of of friction element will transfer properties to appropriate contact properties
    }
    
    else if(strcmp(xName,"Propagate")==0 || strcmp(xName,"AltPropagate")==0)
    {	if(block!=CRACKHEADER && block!=MATERIAL)
			ThrowCompoundErrorMessage(xName," command found at invalid location.","");
		ValidateCommand(xName,NO_BLOCK,cracksDim);
        if(fmobj->IsThreeD())
            throw SAXException("3D cracks not yet programming for crack propagation.");
    	numAttr=(int)attrs.getLength();
		
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
	{	ValidateCommand(xName,CRACKHEADER,cracksDim);
    	input=DOUBLE_NUM;
        inputPtr=(char *)&PropagateTask::cellsPerPropagationStep;
    }
    
    else if(strcmp(xName,"MovePlane")==0)
	{	ValidateCommand(xName,CRACKHEADER,cracksDim);
    	numAttr=(int)attrs.getLength();
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
	{	if(block==CRACKHEADER)
		{	contact.crackContactByDisplacements = false;
			input=DOUBLE_NUM;
			inputPtr=(char *)&contact.crackPositionCutoff;
		}
		else if(block==MULTIMATERIAL)
		{	mpmgrid.contactByDisplacements = false;
			input=DOUBLE_NUM;
			inputPtr=(char *)&mpmgrid.positionCutoff;
		}
		else
			ThrowCompoundErrorMessage(xName," command found at invalid location.","");
    }
	
	else if(strcmp(xName,"CrackParticleSize")==0)
	{
		throw SAXException("<CrackParticleSize> feature requires OSParticulas compiled with CRACK_GIMP");
	}
	
    else if(strcmp(xName,"JContour")==0)
	{	ValidateCommand(xName,CRACKHEADER,cracksDim);
    	numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
        {   aName=XMLString::transcode(attrs.getLocalName(i));
            value=XMLString::transcode(attrs.getValue(i));
            if(strcmp(aName,"type")==0)
                sscanf(value,"%d",&JContourType);
            else if(strcmp(aName,"size")==0)
                sscanf(value,"%d",&JGridSize);
            else if(strcmp(aName,"terms")==0)
                sscanf(value,"%d",&JTerms);
            else if(strcmp(aName,"gridenergy")==0)
                sscanf(value,"%d",&JGridEnergy);
            delete [] aName;
            delete [] value;
        }
    }

	else if(strcmp(xName,"ShiftCracks")==0 )
	{	throw SAXException("<ShiftCracks> feature requires OSParticulas compiled with SHIFT_CRACKS");
	}
	
    //-------------------------------------------------------
    // <Mesh> section
	
	// Start Mesh section - may be at root level or in a 3D crack list
	else if(strcmp(xName,"Mesh")==0)
	{	ValidateCommand(xName,NO_BLOCK,ANY_DIM);
		block = MESHBLOCK;
		archiver->SetArchiveMesh(false);
		value=ReadTagValue("output",attrs);
		if(value!=NULL)
		{	if(strcmp(value,"file")==0) archiver->SetArchiveMesh(true);
			delete [] value;
		}
	}
	
    else if(strcmp(xName,"elem")==0)
	{
		ValidateCommand(xName,ELEMENTLIST,ANY_DIM);
    	input=NODE_BLOCK;
		numAttr=(int)attrs.getLength();
		double segLength = -1.;
		int elemType = -1;
		int czmID=0;
		char czmName[200];
		czmName[0]=0;
		for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
			if(strcmp(aName,"type")==0)
				sscanf(value,"%d",&elemType);
			else if(strcmp(aName,"length")==0)
				sscanf(value,"%lf",&segLength);
			else if(strcmp(aName,"czm")==0)
				sscanf(value,"%d",&czmID);
			else if(strcmp(aName,"czmname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(czmName,value);
			}
			delete [] aName;
			delete [] value;
		}
		// if gave a matname, it takes precedence over mat number
		if(strlen(czmName)>0)
			czmID = matCtrl->GetIDFromNewName(czmName);
		if(elemType<0)
            throw SAXException("<elem> does not specify an element type.");
		if(!theElems->SetElemIDCodes(elemType,block))
			throw SAXException("Invalid or incompatible element type.");
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
	
	// <pt> may be in:
	//  a. <NodeList> for node in a mesh (which might be 3D crack mesh too)
	//  b. MATLPTS to define a material point (XML only and subordinate to <mp> defining only point location)
	//  c. 2D <CrackList> for point in 2D crack
	// <pt> also be in shape commands, but those are in generators (and <GRIDBLOCK skips them here)
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
    	numAttr=(int)attrs.getLength();
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
            	sscanf(value,"%d",&tipMatnum);		// 2D <CrackList> only
            else if(strcmp(aName,"mat")==0)
            	sscanf(value,"%d",&matid);			// 2D <CrackList> only, 1-based ID or following name use
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
		{	// for mesh only
			theNodes->AddNode(xp.x,xp.y,xp.z);
		}
        else if(block==MATLPTS)
		{	if(!mpCtrl->SetPtOrVel(xName,&xp))
				throw SAXException("A material point in an <mp> command is not within the grid.");
        }
        else if(block==CRACKLIST)
		{
            if(!crackCtrl->AddSegment(new CrackSegment(&xp,tipMatnum,matid),false))
                throw SAXException("Crack point not in the mesh or out of memory adding a crack segment.");
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
		double angle=0.,dval;
		double pConcInitial = diffusion!=NULL ? diffusion->reference : 0.;
		double thick=mpmgrid.GetDefaultThickness();
		double pTempInitial=thermal.reference;
    	numAttr=(int)attrs.getLength();
		char mpmat[200];
		mpmat[0]=0;
		Vector lp = MakeVector(0.5, 0.5, 0.5);
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
			{	// only used if diffusion is active
				if(fmobj->HasDiffusion())
				{	pConcInitial=dval;
					if(pConcInitial<0. || pConcInitial>1.)
						throw SAXException("Material point concentration potential must be from 0 to 1");
				}
			}
			else if(strcmp(aName,"wtconc")==0)
			{	// only used if diffusion is active
				if(fmobj->HasDiffusion())
				{	pConcInitial=-dval;
					if(pConcInitial>0.)
						throw SAXException("Material point weight fraction concentration must be >= 0");
				}
			}
#ifdef POROELASTICITY
			else if(strcmp(aName,"pp")==0)
			{	// only used if poroelasticity is active
				if(fmobj->HasPoroelasticity())
				{	pConcInitial=dval;
					pConcInitial *= UnitsController::Scaling(1.e6);
				}
			}
#endif
            else if(strcmp(aName,"temp")==0)
                pTempInitial=dval;
            else if(strcmp(aName,"angle")==0)
                angle=dval;
            else if(strcmp(aName,"lpx")==0)
                lp.x = dval;
			else if(strcmp(aName,"lpy")==0)
                lp.y = dval;
            else if(strcmp(aName,"lpz")==0)
                lp.z = dval;
            delete [] aName;
            delete [] value;
        }
		// if gave a matname, it takes precedence over mat number
		if(strlen(mpmat)>0)
			matl = matCtrl->GetIDFromNewName(mpmat);
        
        // create object for 3D or 2D analysis
		if(fmobj->IsThreeD())
		{	MatPoint3D *newMpt=new MatPoint3D(elemNum,matl,angle);
			mpCtrl->AddMaterialPoint(newMpt,pConcInitial,pTempInitial);
			newMpt->SetDimensionlessSize(&lp);
		}
		else if(fmobj->IsAxisymmetric())
		{	// thickness set to x position in pt command by SetPtOrVel() when it calls SetOrigin()
			MatPointAS *newMpt=new MatPointAS(elemNum,matl,angle,1.);
			mpCtrl->AddMaterialPoint(newMpt,pConcInitial,pTempInitial);
			newMpt->SetDimensionlessSize(&lp);
		}
		else
		{	MatPoint2D *newMpt=new MatPoint2D(elemNum,matl,angle,thick);
			mpCtrl->AddMaterialPoint(newMpt,pConcInitial,pTempInitial);
			newMpt->SetDimensionlessSize(&lp);
		}
    }
 
	// particle mass (usually when reading points from a file)
    else if(strcmp(xName,"mass")==0)
	{	ValidateCommand(xName,MATLPTS,ANY_DIM);
		double aScaling=ReadUnits(attrs,MASS_UNITS);
		double pmass=-1.;
    	numAttr=(int)attrs.getLength();
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
	{	ValidateCommand(xName,NO_BLOCK,cracksDim);
    	block=CRACKLIST;
		CrackHeader *newCrack = new CrackHeader();
		crackCtrl->AddCrack(newCrack,fmobj->IsThreeD());
		
		// only needed for 2D
		double gridThickness=mpmgrid.GetThickness();
		if(gridThickness>0.) newCrack->SetThickness(gridThickness);
		
		// to hold custom contact law
		char tempName[20];
		strcpy(tempName,"Input Only");
		currentContact = new ContactLaw(tempName,CONTACTLAW);
		char lawname[200],Tpropname[200];
		lawname[0]=0;
		Tpropname[0]=0;
		bool hasCustomContact = false;
		int tractionPropID=-1;

		// read crack attributes
		double dval;
		numAttr=(int)attrs.getLength();
        for(i=0;i<numAttr;i++)
		{	aName=XMLString::transcode(attrs.getLocalName(i));
			value=XMLString::transcode(attrs.getValue(i));
			if(strcmp(aName,"type")==0)
			{	if(strcmp(value,"fixed")==0)
					newCrack->SetFixedCrack(TRUE);
				else
					newCrack->SetFixedCrack(FALSE);
			}
			else if(strcmp(aName,"law")==0)
			{	sscanf(value,"%d",&(currentContact->contactProps.contactLawID));
			}
			else if(strcmp(aName,"lawname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(lawname,value);
			}
			else if(strcmp(aName,"Tprop")==0)
			{	sscanf(value,"%d",&tractionPropID);
			}
			else if(strcmp(aName,"Tpropname")==0)
			{	if(strlen(value)>199) value[200]=0;
				strcpy(Tpropname,value);
			}
			else if(strcmp(aName,"friction")==0)
            {	// deprecated - use law or lawname instead
				sscanf(value,"%lf",&dval);
				currentContact->contactProps.friction = dval;
				hasCustomContact = true;
			}
            else if(strcmp(aName,"Dn")==0 || strcmp(aName,"Dnt")==0)
            {	// deprecated - use law or lawname instead
				sscanf(value,"%lf",&dval);
				currentContact->contactProps.Dn = dval*UnitsController::Scaling(1.e6);
				currentContact->contactProps.friction = 11.;
				hasCustomContact = true;
			}
            else if(strcmp(aName,"Dnc")==0)
            {	// deprecated - use law or lawname instead
				sscanf(value,"%lf",&dval);
				currentContact->contactProps.Dnc = dval*UnitsController::Scaling(1.e6);
				currentContact->contactProps.friction = 11.;
				hasCustomContact = true;
			}
            else if(strcmp(aName,"Dt")==0)
            {	// deprecated - use law or lawname instead
				sscanf(value,"%lf",&dval);
				currentContact->contactProps.Dt = dval*UnitsController::Scaling(1.e6);
				currentContact->contactProps.friction = 11.;
				hasCustomContact = true;
			}
            delete [] aName;
            delete [] value;
        }
		
		// if gave a mat name, it takes precedence over law number
		if(lawname[0]!=0)
			currentContact->contactProps.contactLawID = matCtrl->GetIDFromName(lawname);
		
		// If use old method to create contact or interface properties, convert to contact law now
		int finalLawID = currentContact->contactProps.contactLawID;
		if(finalLawID<0 && hasCustomContact)
		{	char ccName[100];
			size_t ccSize=100;
			snprintf(ccName,ccSize,"Crack #%d",newCrack->GetNumber());
			finalLawID = ContactLaw::ConvertOldStyleToContactLaw(matCtrl,&(currentContact->contactProps),0,ccName);
		}
		
		// If has a law, put on the crack
		if(finalLawID>=0)
			newCrack->SetContactLawID(finalLawID);
		
		// if gave traction mat name, it takes precedence over traction number
		if(Tpropname[0]!=0)
			tractionPropID = matCtrl->GetIDFromName(Tpropname);
		if(tractionPropID>0)
			newCrack->SetTractionPropID(tractionPropID);
		
		// done with currentContact
		delete currentContact;
		currentContact = NULL;
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
    	numAttr=(int)attrs.getLength();
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
        NodalVelBC *newVelBC = new NodalVelBC(node,dof,style,(double)0.,ftime,angle,angle2);
        newVelBC->SetID(velID);
		velocityBCs->AddObject(newVelBC);
		
		if(style!=FUNCTION_VALUE)
		{	input=DOUBLE_NUM;
			if(style==LINEAR_VALUE) gScaling = UnitsController::Scaling(1.e3);		// convert Legacy 1/ms to 1/s
			inputPtr = (char *)newVelBC->GetBCValuePtr();
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
        damageICCtrl=new ParseController();
    }

    // Loads on material points (deprecated method)
    else if(strcmp(xName,"LoadBCs")==0)
	{	ValidateCommand(xName,PARTICLEBCHEADER,ANY_DIM);
	    block=LOADEDPOINTS;
    }
	
	// A particle load
    else if(strcmp(xName,"load")==0)
	{	ValidateCommand(xName,LOADEDPOINTS,ANY_DIM);
    	int ptNum=0;
		int dof=0,style=1;
    	numAttr=(int)attrs.getLength();
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
	{	throw SAXException("The <Isothermal> command has been replaced by the PropertyRamp custom task.");
    }
    
	// Turn on conduction analysis
    else if(strcmp(xName,"Conduction")==0)
	{	ValidateCommand(xName,THERMAL,ANY_DIM);
		ConductionTask::active=true;
	}

	// Turn on crack tip heating
    else if(strcmp(xName,"CrackTipHeating")==0)
	{	ValidateCommand(xName,THERMAL,cracksDim);
		ConductionTask::crackTipHeating = true;
	}
	
	// Turn on crack contact heating by friction
    else if(strcmp(xName,"CrackContactHeating")==0)
	{	ValidateCommand(xName,THERMAL,cracksDim);
		ConductionTask::crackContactHeating = true;
	}
	
	// Turn on material contact heating by friction
    else if(strcmp(xName,"ContactHeating")==0)
	{	ValidateCommand(xName,THERMAL,ANY_DIM);
		ConductionTask::matContactHeating=true;
	}
	
	// Turn on adiabatic mode
    else if(strcmp(xName,"EnergyCoupling")==0)
	{	ValidateCommand(xName,THERMAL,ANY_DIM);
		ConductionTask::adiabatic = true;
	}
	
	else if(strcmp(xName,"ContactHeatFlow")==0)
	{	throw SAXException("<ContactHeatFlow> command requires OSParticulas.");
	}
	
    //-------------------------------------------------------
    // <Gravity> section
	
	// begin Gravity section - starts block
	// Use Body(XYZ)Force to set gravity or GridBody(XYZ)Force to set using functions
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
	
	// Grid Body Force in X direction
    else if(strcmp(xName,"GridBodyXForce")==0)
	{	ValidateCommand(xName,GRAVITY,ANY_DIM);
    	input=GRID_X_BODY_FORCE_FUNCTION_BLOCK;
        inputPtr=(char *)&bodyFrc;		// although not needed
    }
    
	// Grid Body Force in Y direction
    else if(strcmp(xName,"GridBodyYForce")==0)
	{	ValidateCommand(xName,GRAVITY,ANY_DIM);
    	input=GRID_Y_BODY_FORCE_FUNCTION_BLOCK;
        inputPtr=(char *)&bodyFrc;		// although not needed
    }
    
	// Grid Body Force in Z direction
    else if(strcmp(xName,"GridBodyZForce")==0)
	{	ValidateCommand(xName,GRAVITY,ANY_DIM);
    	input=GRID_Z_BODY_FORCE_FUNCTION_BLOCK;
        inputPtr=(char *)&bodyFrc;		// although not needed
    }
    
    //-------------------------------------------------------
    // <CustomTasks> section
	
	// begin CustomTasks section
    else if(strcmp(xName,"CustomTasks")==0)
    	block=CUSTOMTASKS;

	// Schedule a custom task
    else if(strcmp(xName,"Schedule")==0)
	{	ScheduleCustomTask(attrs);
    }
    
	// Set Custom Task parameter
    else if(strcmp(xName,"Parameter")==0)
	{	SetCustomTasksParameter(attrs);
    }

	//-------------------------------------------------------
	// Unknown element
	else
		return FALSE;
    
    return TRUE;
}

// End an element
// throws SAXException()
void MPMReadHandler::myEndElement(char *xName)
{
    if(strcmp(xName,"CustomTasks")==0 || strcmp(xName,"Gravity")==0)
    {	block=NO_BLOCK;
    }
    
    else if(strcmp(xName,"CrackList")==0)
    {   if(!crackCtrl->FinishCrack())
            throw SAXException("All cracks must have at least one segement or surface element");
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
	{	firstVelocityBC = (NodalVelBC *)velocityBCs->firstObject;
		lastVelocityBC = (NodalVelBC *)velocityBCs->lastObject;
		firstConcBC = (NodalConcBC *)concBCs->firstObject;
		lastConcBC = (NodalConcBC *)concBCs->lastObject;
		firstTempBC = (NodalTempBC *)tempBCs->firstObject;
		lastTempBC = (NodalTempBC *)tempBCs->lastObject;
		delete velocityBCs;
		delete concBCs;
		delete tempBCs;
		block=NO_BLOCK;
	}
	
	else if(strcmp(xName,"ParticleBCs")==0)
	{	firstLoadedPt=(MatPtLoadBC *)mpLoadCtrl->firstObject;
		firstTractionPt=(MatPtTractionBC *)mpTractionCtrl->firstObject;
		firstFluxPt=(MatPtFluxBC *)mpConcFluxCtrl->firstObject;
		firstHeatFluxPt=(MatPtHeatFluxBC *)mpHeatFluxCtrl->firstObject;
        firstDamagedPt=(InitialCondition *)damageICCtrl->firstObject;
		delete mpLoadCtrl;
		delete mpTractionCtrl;
		delete mpConcFluxCtrl;
		delete mpHeatFluxCtrl;
		block=NO_BLOCK;
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
    
	else if(strcmp(xName,"PatchGrid")==0)
	{	block=MPMHEADER;
	}
	
	else if(strcmp(xName,"Cracks")==0)
    {	// install frictionless if not provded
        // but only if Cracks element is in the MPMHeader
		if(contact.crackContactLawID<0)
		{	contact.crackContactLawID = ContactLaw::ConvertOldStyleToContactLaw(matCtrl,NULL,0.,"Cracks default");
		}
    	block=MPMHEADER;
    }
    
    else if(strcmp(xName,"MultiMaterialMode")==0)
    {	// install frictionless if not provded
        // but only if MultiMaterialMode element is in the MPMHeader
        if(mpmgrid.materialContactLawID<0)
		{	mpmgrid.materialContactLawID = ContactLaw::ConvertOldStyleToContactLaw(matCtrl,NULL,0.,"MultiMaterialMode default");
		}
		block=MPMHEADER;
    }
    
    else if(strcmp(xName,"Schedule")==0)
	{	((CustomTask *)currentTask)->Finalize();
		block=CUSTOMTASKS;
	}
    
    else if(strcmp(xName,"Friction")==0)
    {	// Transfer contact law to appropriate location
		
		// convert to the final contact law
		int finalLawID = currentContact->contactProps.contactLawID;
		if(finalLawID<0)
			finalLawID = ContactLaw::ConvertOldStyleToContactLaw(matCtrl,&(currentContact->contactProps),0.,"deprecated Friction");
		
		// assign law to appropriate place
		if(block==CRACKHEADER)
		{	contact.crackContactLawID = finalLawID;
		}
		else if(block==MULTIMATERIAL)
		{	mpmgrid.materialContactLawID = finalLawID;
		}
		else
		{	// custom material/material contact
			matCtrl->SetMaterialFriction(finalLawID,currentContact->contactProps.otherMatID);
		}
		
		// done with currentContact
		delete currentContact;
		currentContact = NULL;
    }
    
    else if(strcmp(xName,"PDamping")==0)
    {	if(block==MATERIAL)
			matCtrl->SetMaterialDamping();
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
	
	else if(strcmp(xName,"MatlPtsPerElement")==0)
	{	// This guarantees these entries are always valid
        int totalPoints;
        if(block==BMPBLOCK)
            totalPoints = bmpCustomPtsPerElement;
        else if(block==MPMHEADER)
            totalPoints = fmobj->ptsPerElement;
        else
            totalPoints = abs(fmobj->customPtsPerElement);
		if(totalPoints<0)
			throw SAXException("Material points per element must be positive.");
		int oneSide;
		if(fmobj->IsThreeD())
		{	oneSide = pow((double)totalPoints,1./3.)+0.1;
			if(oneSide*oneSide*oneSide!=totalPoints)
				throw SAXException("Material points per element in 3D must be a cubed number.");
		}
		else
		{	oneSide = sqrt((double)totalPoints)+0.1;
			if(oneSide*oneSide!=totalPoints)
				throw SAXException("Material points per element in 2D must be a squared number.");
		}
		// align input number with number per side
		if(block==MPMHEADER)
			fmobj->ptsPerSide = oneSide;
        if(block==BMPBLOCK)
            bmpCustomPtsPerSide = oneSide;
		else
        {   // if provided negative number, set per side to negative (keep total points positive)
            fmobj->customPtsPerSide = fmobj->customPtsPerElement>0 ? oneSide : -oneSide;
            fmobj->customPtsPerElement = totalPoints;
        }
	}
}

// Decode block of characters if input!=NO_INPUT
// throws SAXException()
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
			double bcvalue,bcftime;
            sscanf(xData,"%lf,%lf",&bcvalue,&bcftime);
			newBC->SetBCValue(bcvalue);
			newBC->SetBCFirstTime(bcftime);
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
		
		case PRESSURE_FUNCTION_BLOCK:
			((TaitLiquid *)inputPtr)->SetPressureFunction(xData);
			break;
		
		case GRID_X_BODY_FORCE_FUNCTION_BLOCK:
		case GRID_Y_BODY_FORCE_FUNCTION_BLOCK:
		case GRID_Z_BODY_FORCE_FUNCTION_BLOCK:
			bodyFrc.SetGridBodyForceFunction(xData,input);
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
		
		case PATCH_SIZES:
			mpmgrid.SetPatchSizes(xData);
			break;
        
        case HARDENING_LAW_SELECTION:
            ((MaterialBase *)inputPtr)->SetHardeningLaw(xData);
            break;

        case INITIATION_LAW_SELECTION:
            ((MaterialBase *)inputPtr)->SetInitiationLaw(xData);
            break;
			
        case SOFTI_LAW_SELECTION:
		case SOFTII_LAW_SELECTION:
            ((MaterialBase *)inputPtr)->SetSofteningLaw(xData,input);
            break;
			
		case TEXT_PARAMETER:
			// must be in active custom tasks
			((CustomTask *)currentTask)->SetTextParameter(xData,inputPtr);
			break;
            
        default:
            break;
    }
}
