/********************************************************************************
	CustomThermalRamp.cpp
	nairn-mpm-fea

	Created by John Nairn on 9/14/16.
	Copyright (c) 2016 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/CustomThermalRamp.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "System/UnitsController.hpp"
#include "Read_XML/mathexpr.hpp"
#include "System/ArchiveData.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Exceptions/CommonException.hpp"
#include "Read_XML/BMPLevel.hpp"

extern double timestep;

double CustomThermalRamp::varTime=0.;
double CustomThermalRamp::varXValue=0.;
double CustomThermalRamp::varYValue=0.;
double CustomThermalRamp::varZValue=0.;
PRVar rampVarArray[4] = { NULL, NULL, NULL, NULL };

#pragma mark Constructors and Destructors

// Constructors
CustomThermalRamp::CustomThermalRamp()
{
	isoRampTime = -1.;
	rampStart = 0.;
	isoDeltaT = 0.;
	sigmoidal = 0;
	scaleFxn = NULL;
	
	// for imaged ramp
	bmpFile = NULL;
	orig = MakeVector(0.,0.,-1.e-9);	// <-1.-8 means zlevel not supplied
	width = -1.e-9;			// <-1.-8 means not supplied
	height = -1.e-9;
	flipped = false;
	deltaTmin = 0.;
	deltaTmax = 0.;
	firstLevel = NULL;
	currentLevel = NULL;
	firstMatPt = NULL;
	currentMatPt = NULL;
}

// Destructor (and it is virtual)
CustomThermalRamp::~CustomThermalRamp()
{	if(scaleFxn!=NULL)
	{	delete rampVarArray[0];
		delete rampVarArray[1];
		delete rampVarArray[2];
		delete rampVarArray[3];
		delete scaleFxn;
	}
	if(bmpFile!=NULL) delete [] bmpFile;
	
	// delete defined levels
	while(firstLevel!=NULL)
	{	currentLevel=(BMPLevel *)firstLevel->GetNextObject();
		delete firstLevel;
		firstLevel=currentLevel;
	}
	currentLevel=NULL;
	
	// delete mapped material points
	while(firstMatPt!=NULL)
	{	currentMatPt=(BMPLevel *)firstMatPt->GetNextObject();
		delete firstMatPt;
		firstMatPt=currentMatPt;
	}
	currentMatPt=NULL;
}

// Return name of this task
const char *CustomThermalRamp::TaskName(void) { return "Ramp particle temperatures"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *CustomThermalRamp::InputParam(char *pName,int &input,double &gScaling)
{
	if(strcmp(pName,"time")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&isoRampTime,gScaling,1.e-3);
	}
	
	else if(strcmp(pName,"DeltaT")==0)
	{	input=DOUBLE_NUM;
		return (char *)&isoDeltaT;
	}
	
	else if(strcmp(pName,"start")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&rampStart,gScaling,1.e-3);
	}
	
	else if(strcmp(pName,"sigmoidal")==0)
	{	input=INT_NUM;
		return (char *)&sigmoidal;
	}
	
	else if(strcmp(pName,"scale")==0)
	{	input=TEXT_PARAMETER;
		return (char *)&scaleFxn;
	}
	
	else if(strcmp(pName,"file")==0)
	{	input=TEXT_PARAMETER;
		return (char *)&bmpFile;
	}
	
	else if(strcmp(pName,"xorigin")==0)
	{	input=DOUBLE_NUM;
		return (char *)&orig.x;
	}
	
	else if(strcmp(pName,"yorigin")==0)
	{	input=DOUBLE_NUM;
		return (char *)&orig.y;
	}
	
	else if(strcmp(pName,"zlevel")==0)
	{	input=DOUBLE_NUM;
		return (char *)&orig.z;
	}
	
	else if(strcmp(pName,"width")==0)
	{	input=DOUBLE_NUM;
		return (char *)&width;
	}
	
	else if(strcmp(pName,"height")==0)
	{	input=DOUBLE_NUM;
		return (char *)&height;
	}
	
	else if(strcmp(pName,"DeltaTmin")==0)
	{	input=DOUBLE_NUM;
		return (char *)&deltaTmin;
	}
	
	else if(strcmp(pName,"DeltaTmax")==0)
	{	input=DOUBLE_NUM;
		return (char *)&deltaTmax;
	}
	
	else if(strcmp(pName,"flipped")==0)
	{	input=INT_NUM;
		return (char *)&flipped;
	}
	
	else
	{	// look for "map # #"
		char* check4function = strstr(pName,"map ");
		if(check4function==&pName[0])
		{	// read two numbers
			vector<double> range;
			if(!CommonReadHandler::GetFreeFormatNumbers(&pName[4],range,1.0)) return NULL;
			if(range.size()!=2)
			{	cout << "*** map parameter must supply exactly two numbers" << endl;
				return NULL;
			}
			
			// create level for provided deltT
			int imin = int(range[0]+.5);
			int imax = int(range[1]+.5);
			if(imin<0 || imax>255 || imax<imin)
			{	cout << "*** map embedded numbers mut be 0 to 255 and first < second" << endl;
				return NULL;
			}
			BMPLevel *newLevel=new BMPLevel(0,imin,imax);
			if(currentLevel==NULL)
				firstLevel=newLevel;
			else
				currentLevel->SetNextObject(newLevel);
			currentLevel=newLevel;
			input=DOUBLE_NUM;
			return (char *)&(currentLevel->temperature);
		}
	}
	
	// check remaining commands
	return CustomTask::InputParam(pName,input,gScaling);
}

// Actually sets load function
// throws std::bad_alloc, SAXException()
void CustomThermalRamp::SetTextParameter(char *fxn,char *ptr)
{
	if(ptr == (char *)&bmpFile)
	{	// bmp file name
		if(bmpFile!=NULL)
			ThrowSAXException("Duplicate bit mapped file name for thermal ramp");
		if(fxn==NULL)
			ThrowSAXException("Bit mapped file name is missing");
		if(strlen(fxn)==0)
			ThrowSAXException("Bit mapped file name is missing");
		
		// save is
		bmpFile = new char[strlen(fxn)+1];
		strcpy(bmpFile,fxn);
	}
	
	else
	{	// scale function
		if(scaleFxn!=NULL)
			ThrowSAXException("Duplicate ramp scale function was entered");
		if(fxn==NULL)
			ThrowSAXException("Ramp scale function is missing");
		if(strlen(fxn)==0)
			ThrowSAXException("Ramp scale function is missing");
		
		// create variable
		if(rampVarArray[0]==NULL)
		{	rampVarArray[0]=new RVar("t",&varTime);
			rampVarArray[1]=new RVar("x",&varXValue);
			rampVarArray[2]=new RVar("y",&varYValue);
			rampVarArray[3]=new RVar("z",&varZValue);
		}
		
		// create function
		scaleFxn = new ROperation(fxn,4,rampVarArray);
		if(scaleFxn->HasError())
			ThrowSAXException("Ramp scale function is not valid");
	}
}

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *CustomThermalRamp::Initialize(void)
{
	// Ramp should be at least one time step
	int nsteps = (int)(isoRampTime/timestep);
	if(nsteps<1)
	{	nsteps = 1;
		isoRampTime = timestep;
	}
	
	// start after zero
	if(rampStart<0.) rampStart = 0.;
	
	// official end time
	endTime = rampStart+isoRampTime+timestep;
	
	// current Delta T
	currentDeltaT = 0.;
	
	cout << "Ramp particle temperatures." << endl;
	
	// Two options
	// 1. Ramp to isoDeltaT, optionally scalled by scaleFxn (when bmpFile is NULL)
	// 2. Ramp to bit mapped file settings
	
	char hline[200];
	
	if(bmpFile==NULL)
	{	sprintf(hline,"   Final isothermal temperature difference: %g C",isoDeltaT);
		cout << hline << endl;
		sprintf(hline,"   Ramped between %g and %g %s",rampStart*UnitsController::Scaling(1.e3),
					(rampStart+isoRampTime)*UnitsController::Scaling(1.e3),UnitsController::Label(ALTTIME_UNITS));
		cout << hline << endl;
		cout << "      (which covers " << nsteps << " time steps)" << endl;
		if(sigmoidal)
			cout << "      (use sigmoidal ramp)" << endl;
		if(scaleFxn!=NULL)
		{	char *expr = scaleFxn->Expr('#');
			cout << "   Scaling function =  " << expr << endl;
			delete [] expr;
		}
	}
	else
	{	// Read the File
		XYInfoHeader info;
		unsigned char **rows;
		char *bmpFullPath = archiver->ExpandOutputPath(bmpFile);
		CommonReadHandler::ReadBMPFile(bmpFullPath,info,&rows);
		delete [] bmpFullPath;
		
		// validate
		if(firstLevel==NULL && DbleEqual(deltaTmin,0.) && DbleEqual(deltaTmax,0.))
			throw CommonException("ThermalRamp for bit mapped file does not define any ramp values.","CustomThermalRamp::Initialize");
		deltaTscale = (deltaTmax-deltaTmin)/255.;
		
		// get final width an height
		Vector pw;
		const char *msg = CommonReadHandler::DecodeBMPWidthAndHeight(info,width,height,orig.z,pw,fmobj->IsThreeD());
		if(msg != NULL)
			throw CommonException("ThermalRamp from bit mapped file specify width and/or height as size or pixels per mm.","CustomThermalRamp::Initialize");
		
		// print settings
		cout << "   Ramp for bit mapped file: " << bmpFile << endl;
		if(firstLevel!=NULL)
		{	BMPLevel *nextLevel = firstLevel;
			while(nextLevel!=NULL)
			{	nextLevel->OutputLevel();
				nextLevel = (BMPLevel *)nextLevel->GetNextObject();
			}
		}
		else
		{	cout << "      0 to 255 linearly mapped to dT = " << deltaTmin << " to " << deltaTmax << endl;
		}
		
		// loop over nonrigid material points and find those needed thermal ramp
		DomainMap map;
		for(int p=0;p<nmpmsNR;p++)
		{	// get box
			double r1x,r2y,r3z;
			Vector del;
			mpm[p]->GetUndeformedSemiSides(&r1x,&r2y,&r3z);
			if(fmobj->IsThreeD())
				del = MakeVector(r1x,r2y,r3z);
			else
				del = MakeVector(r1x,r2y,-1.);
			
			// if point in the view area, then check it
			if(CommonReadHandler::MapDomainToImage(info,mpm[p]->pos,orig,del,pw,width,height,map))
			{	// find maximum level or bit value
				double particleDeltaT = 0.;
				if(firstLevel!=NULL)
				{	BMPLevel *nextLevel = CommonReadHandler::FindBMPLevel(firstLevel,map,rows);
					if(nextLevel!=NULL) particleDeltaT = nextLevel->temperature;
				}
				else
				{	double totalIntensity = CommonReadHandler::FindAverageValue(map,rows);
					if(totalIntensity>0.)
						particleDeltaT = deltaTmin+totalIntensity*deltaTscale;
				}
				
				if(!DbleEqual(particleDeltaT,0.))
				{	// holder for this material point pointer
					BMPLevel *newMatPt=new BMPLevel();
					newMatPt->temperature = particleDeltaT;
					newMatPt->SetContextInfo((void *)mpm[p]);
					if(currentMatPt==NULL)
						firstMatPt=newMatPt;
					else
						currentMatPt->SetNextObject(newMatPt);
					currentMatPt=newMatPt;
				}
			}
		}
		
		// any particle found
		if(firstMatPt==NULL)
			throw CommonException("ThermalRamp from bit mapped file does not overlap any material points.","CustomThermalRamp::Initialize");
	}
	
	return nextTask;
}

// called when MPM step is getting ready to do custom tasks
// never uses extrapolations so no need to set needExtraps
CustomTask *CustomThermalRamp::PrepareForStep(bool &needExtraps)
{
	double effTime = mtime+timestep;
	if(effTime<rampStart || effTime>endTime)
		doRamp = false;
	else
		doRamp = true;
	
	return nextTask;
}

// Adjust time step now
CustomTask *CustomThermalRamp::StepCalculation(void)
{
	// exit when not needed
	if(!doRamp) return nextTask;
	
	// get new temperature
	double rampFraction = (mtime+timestep-rampStart)/isoRampTime;
	if(rampFraction>1.)
	{	rampFraction = 1.;
	}
	else if(sigmoidal)
	{	rampFraction = 1./(1+exp(-12.*(rampFraction-0.5)));
	}
	
	// Adjust particle temperatures
	double newDeltaT,deltaT;
	if(firstMatPt == NULL)
	{	// to current particles
		newDeltaT = rampFraction*isoDeltaT;
		
		// get temperature change
		deltaT = newDeltaT - currentDeltaT;
		currentDeltaT = newDeltaT;
		
		// loop over nonrigid material points and update temperature
		for(int p=0;p<nmpmsNR;p++)
		{	if(scaleFxn==NULL)
				mpm[p]->pTemperature += deltaT;
			else
			{	varTime = mtime*UnitsController::Scaling(1000.);
				varXValue = mpm[p]->pos.x;
				varYValue = mpm[p]->pos.y;
				varZValue = mpm[p]->pos.y;
				mpm[p]->pTemperature += deltaT*scaleFxn->Val();
			}
		}
	}
	
	else
	{	// to list of material points
		BMPLevel *nextMatPt = firstMatPt;
		while(nextMatPt!=NULL)
		{	newDeltaT = rampFraction*nextMatPt->temperature;
			
			// get temperature change on particle
			deltaT = newDeltaT - nextMatPt->currentDeltaT;
			nextMatPt->currentDeltaT = newDeltaT;
			MPMBase *mptr = (MPMBase *)nextMatPt->GetContextInfo();
			mptr->pTemperature += deltaT;
			
			// next one
			nextMatPt = (BMPLevel *)nextMatPt->GetNextObject();
		}
	}
	
	return nextTask;
}

// Called when custom tasks are all done on a step
CustomTask *CustomThermalRamp::FinishForStep(bool &removeMe)
{	double effTime = mtime+timestep;
	if(effTime>endTime) removeMe = true;
	return nextTask;
}


