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
#include "System/ArchiveData.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Exceptions/CommonException.hpp"
#include "Read_XML/BMPLevel.hpp"
#include "Materials/MaterialBase.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Read_XML/Expression.hpp"

extern double timestep;

#pragma mark Constructors and Destructors

// Constructors
CustomThermalRamp::CustomThermalRamp()
{
	isoRampTime = -1.;
	property = RAMP_TEMP;
	rampStart = 0.;
	isoDeltaT = 0.;
	sigmoidal = 0;
	scaleFxn = NULL;
	
	// for imaged ramp
	bmpFile = NULL;
	orig = MakeVector(0.,0.,-1.e-9);	// <-1.-8 means zlevel not supplied
	width = -1.e9;			// <-1.-8 means not supplied
	height = -1.e9;
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
{
	if(scaleFxn!=NULL) delete scaleFxn;
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

// Return name of this task (not user)
const char *CustomThermalRamp::TaskName(void)
{	return "Ramp particle temperature, concentration, or other property";
}

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *CustomThermalRamp::InputParam(char *pName,int &input,double &gScaling)
{
	if(strcmp(pName,"time")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&isoRampTime,gScaling,1.e-3);
	}
	
	else if(strcmp(pName,"property")==0)
	{	input=INT_NUM;
		return (char *)&property;
	}
	
	else if(strcmp(pName,"DeltaT")==0 || strcmp(pName,"Delta")==0)
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
	
	else if(strcmp(pName,"DeltaTmin")==0 || strcmp(pName,"DeltaMin")==0)
	{	input=DOUBLE_NUM;
		return (char *)&deltaTmin;
	}
	
	else if(strcmp(pName,"DeltaTmax")==0 || strcmp(pName,"DeltaMax")==0)
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
		{	ThrowSAXException("Duplicate bit mapped file name for ramp");
			return;
		}
		if(fxn==NULL)
		{	ThrowSAXException("Bit mapped file name is missing");
			return;
		}
		if(strlen(fxn)==0)
		{	ThrowSAXException("Bit mapped file name is missing");
			return;
		}
		
		// save is
		bmpFile = new char[strlen(fxn)+1];
		strcpy(bmpFile,fxn);
	}
	
	else if(ptr == (char *)&scaleFxn)
	{	// scale function
		if(scaleFxn!=NULL)
		{	ThrowSAXException("Duplicate ramp scale function was entered");
			return;
		}
		if(fxn==NULL)
		{	ThrowSAXException("Ramp scale function is missing");
			return;
		}
		if(strlen(fxn)==0)
		{	ThrowSAXException("Ramp scale function is missing");
			return;
		}

		scaleFxn =  Expression::CreateExpression(fxn,"Ramp scale function is not valid");
	}
	
	else
		CustomTask::SetTextParameter(fxn,ptr);
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
        endTime = 1.9*timestep;
	}
    else
        endTime = isoRampTime+timestep;
	
	// start after zero
	if(rampStart<0.) rampStart = 0.;
	endTime += rampStart;
	
	// current Delta T
	currentDeltaT = 0.;
	
	// convert to ramp type
    char rampValue[100],rampUnits[20];
    switch(property)
    {   case RAMP_TEMP:
            strcpy(rampValue,"isothermal temperature");
            cout << "Ramp particle temperatures:" << endl;
            strcpy(rampUnits,"K");
            break;
        
        case RAMP_FLUID:
            if(fmobj->HasDiffusion())
            {   property = RAMP_CONC;
                strcpy(rampValue,"concentration");
                cout << "Ramp particle concentration:" << endl;
                rampUnits[0] = 0;
            }
#ifdef POROELASTICITY
            else if(fmobj->HasPoroelasticity())
            {   property = RAMP_PP;
                strcpy(rampValue,"pore pressure");
                cout << "Ramp particle pore pressure:" << endl;
                strcpy(rampUnits,UnitsController::Label(PRESSURE_UNITS));
            }
#endif
            else
            {    throw CommonException(
                            "Cannot ramp particle fluid unless doing diffusion or poroelasticity calculations.",
                            "CustomThermalRamp::Initialize");
            }
            break;
            
        case RAMP_OOPSE:
            if(fmobj->np==PLANE_STRESS_MPM)
            {   property = RAMP_SZZ;
                strcpy(rampValue,"sigma(zz)");
                cout << "Ramp particle sigma(zz):" << endl;
                strcpy(rampUnits,UnitsController::Label(PRESSURE_UNITS));
            }
            else if(fmobj->np==PLANE_STRAIN_MPM)
            {   property = RAMP_EZZ;
                strcpy(rampValue,"strain(zz)");
                cout << "Ramp particle strain(zz):" << endl;
                strcpy(rampUnits,"");
            }
            else
            {    throw CommonException(
                                      "Cannot ramp out-of-plane stress or strain unless doing plane stress or strain calculations.",
                                      "CustomThermalRamp::Initialize");
            }
            break;

        case RAMP_ADAM_STRENGTH:
            property = RAMP_STRENGTH;
            // convert to first step
            nsteps = 1;
            isoRampTime = timestep;
            endTime = rampStart+1.9*timestep;
            strcpy(rampValue,"relative initiation stress");
            cout << "Relative particle initiation stress:" << endl;
            rampUnits[0] = 0;
            break;
            
        case RAMP_ADAM_TOUGHNESS:
            property = RAMP_TOUGHNESS;
            // convert to first step
            nsteps = 1;
            isoRampTime = timestep;
            endTime = rampStart+1.9*timestep;
            strcpy(rampValue,"scaled damge toughness");
            cout << "Relative particle damage toughness:" << endl;
            rampUnits[0] = 0;
            break;
        
        default:
            throw CommonException("Invalid property select for custom thermal ramp task.","CustomThermalRamp::Initialize");
            break;
    }
	
	// Two options
	// 1. Ramp to isoDeltaT, optionally scaled by scaleFxn (when bmpFile is NULL)
	// 2. Ramp to bit mapped file settings
	
	char hline[200];
	
	if(bmpFile==NULL)
	{	sprintf(hline,"   Final %s difference: %g %s",rampValue,isoDeltaT,rampUnits);
		cout << hline << endl;
		sprintf(hline,"   Ramped between %g and %g %s",rampStart*UnitsController::Scaling(1.e3),
					(rampStart+isoRampTime)*UnitsController::Scaling(1.e3),UnitsController::Label(ALTTIME_UNITS));
		cout << hline << endl;
		cout << "      (which covers " << nsteps << " time steps)" << endl;
		if(sigmoidal && nsteps>1)
			cout << "      (use sigmoidal ramp)" << endl;
		if(scaleFxn!=NULL)
			cout << "   Scaling function =  " << scaleFxn->GetString() << endl;
	}
	else
	{	// Read the File
		XYInfoHeader info;
		unsigned char **rows;
		char *bmpFullPath = archiver->ExpandOutputPath(bmpFile);
		rows = (unsigned char **)CommonReadHandler::ReadXYFile(bmpFullPath,info,BYTE_DATA,false);
		delete [] bmpFullPath;
		
		// validate
		if(firstLevel==NULL && DbleEqual(deltaTmin,0.) && DbleEqual(deltaTmax,0.))
			throw CommonException("PropertyRamp using bit mapped file does not define any ramp values.","CustomThermalRamp::Initialize");
		deltaTscale = (deltaTmax-deltaTmin)/255.;
		
		// get final width an height
		Vector pw;
		const char *msg = CommonReadHandler::DecodeBMPWidthAndHeight(info,width,height,orig.z,pw,fmobj->IsThreeD());
		if(msg != NULL)
			throw CommonException("PropertyRamp using bit mapped file must specify width and/or height as size or pixels per mm.","CustomThermalRamp::Initialize");
        if(flipped!=0)
            pw.z = info.topDown==0 ? -1. : 1. ;
        else
            pw.z = info.topDown==0 ? 1. : -1. ;

		// print settings
		cout << "   Ramp using bit mapped file: " << bmpFile << endl;
		if(firstLevel!=NULL)
		{	BMPLevel *nextLevel = firstLevel;
			while(nextLevel!=NULL)
			{	nextLevel->OutputLevel();
				nextLevel = (BMPLevel *)nextLevel->GetNextObject();
			}
		}
		else
		{	cout << "      0 to 255 linearly mapped to " << rampValue << " = " << deltaTmin << " to " << deltaTmax << endl;
		}
		
		// loop over nonrigid material points and find those needing the ramp
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
			throw CommonException("PropertyRamp using bit mapped file does not overlap any material points.","CustomThermalRamp::Initialize");
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
	{	// This section is calculated ramp applied to all particles
		
		// to current particles
		newDeltaT = rampFraction*isoDeltaT;
		
		// get temperature change
		deltaT = newDeltaT - currentDeltaT;
		currentDeltaT = newDeltaT;
		
		// loop over nonrigid material points and update temperature
		for(int p=0;p<nmpmsNR;p++)
        {   double dTp = deltaT;
			if(scaleFxn!=NULL)
				dTp = deltaT*scaleFxn->XYZTValue(&mpm[p]->pos,mtime*UnitsController::Scaling(1000.));
			
			if(property==RAMP_TEMP)
			{	mpm[p]->pTemperature += dTp;
				// heat and entropy (per unit mass)
				const MaterialBase *matRef=theMaterials[mpm[p]->MatID()];
				double dqp = matRef->GetHeatCapacity(mpm[p])*dTp;
				mpm[p]->AddHeatEnergy(dqp);
				mpm[p]->AddEntropy(dqp,mpm[p]->pPreviousTemperature);
			}
			else if(property==RAMP_SZZ)
			{	// for generalized plane stress convert to specific stress units
				dTp *= UnitsController::Scaling(1.e6)/mpm[p]->GetRho();
			
				// set particle sigma(zz) and store increment for res.doopse
				mpm[p]->StoreThicknessStressIncrement(dTp);
			}
			else if(property==RAMP_EZZ)
			{	// for generalized plane strain convert to absolute strain
				dTp *= 0.01;
				
				// set particle sigma(zz) and store increment for res.doopse
				mpm[p]->StoreThicknessStrainIncrement(dTp);
			}
			else if(property==RAMP_PP)
			{	// convert to MPa to Pa in Legacy
				dTp*=UnitsController::Scaling(1.e6);
				mpm[p]->pConcentration += dTp;
			}
			else if(property==RAMP_CONC)
			{	// add concentration, contrained to 0 to 1
				mpm[p]->pConcentration += dTp;
				if(mpm[p]->pConcentration<0.)
					mpm[p]->pConcentration = 0.;
				else if(mpm[p]->pConcentration>1.)
					mpm[p]->pConcentration = 1.;
			}
            else if(property==RAMP_STRENGTH)
            {   mpm[p]->SetRelativeStrength(newDeltaT);
            }
            else if(property==RAMP_TOUGHNESS)
            {   mpm[p]->SetRelativeToughness(newDeltaT);
            }
		}
	}
	
	else
	{	// This section maps image to particle change
		
		// to list of material points
		BMPLevel *nextMatPt = firstMatPt;
		while(nextMatPt!=NULL)
		{	newDeltaT = rampFraction*nextMatPt->temperature;
 			
			// get temperature change on particle
			deltaT = newDeltaT - nextMatPt->currentDeltaT;
			nextMatPt->currentDeltaT = newDeltaT;
			MPMBase *mptr = (MPMBase *)nextMatPt->GetContextInfo();
			
			if(property==RAMP_TEMP)
			{	mptr->pTemperature += deltaT;
				// heat and entropy (per unit mass)
				const MaterialBase *matRef=theMaterials[mptr->MatID()];
				double dqp = matRef->GetHeatCapacity(mptr);
				mptr->AddHeatEnergy(dqp);
				mptr->AddEntropy(dqp,mptr->pPreviousTemperature);
			}
			else if(property==RAMP_SZZ)
			{	// for generalized plane stress convert to specific stress units
				deltaT *= UnitsController::Scaling(1.e6)/mptr->GetRho();
				
				// set particle sigma(zz) and store increment for res.doopse
				mptr->StoreThicknessStressIncrement(deltaT);
			}
			else if(property==RAMP_EZZ)
			{	// for generalized plane strain convert to absolute strain
				deltaT *= 0.01;
				
				// set particle sigma(zz) and store increment for res.doopse
				mptr->StoreThicknessStrainIncrement(deltaT);
			}
			else if(property==RAMP_PP)
			{	// Convert MPa to Pa in Legacy
				deltaT *= UnitsController::Scaling(1.e6);
				mptr->pConcentration += deltaT;
			}
			else if(property==RAMP_CONC)
			{	// add temperature contrained to 0 to 1
				mptr->pConcentration += deltaT;
				if(mptr->pConcentration<0.)
					mptr->pConcentration = 0.;
				else if(mptr->pConcentration>1.)
					mptr->pConcentration = 1.;
			}
            else if(property==RAMP_STRENGTH)
            {   mptr->SetRelativeStrength(newDeltaT);
            }
            else if(property==RAMP_TOUGHNESS)
            {   mptr->SetRelativeToughness(newDeltaT);
            }

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


