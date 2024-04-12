/********************************************************************************
    LoadControl.cpp
    nairn-mpm-fea
    
    Created by Chad Hammerquist April 2016
    Copyright (c) 2016 John A. Nairn, All rights reserved.
    Revised 2023 by John A. Nairn
	
	PID controller to update velocity
********************************************************************************/

#include "stdafx.h"
#include <fstream>

#include "Custom_Tasks/LoadControl.hpp"
#include "Custom_Tasks/PropagateTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/RigidMaterial.hpp"
#include "Exceptions/CommonException.hpp"
#include "Global_Quantities/GlobalQuantity.hpp"
#include "System/ArchiveData.hpp"
#include "System/UnitsController.hpp"
#include "Read_XML/Expression.hpp"

#define DYNAMIC_ADJUST

#pragma mark INITIALIZE

// Constructors
LoadControl::LoadControl() : CustomTask()
{
	// input required
	material = -1;
    matName = NULL;
    Load_Function=NULL;

    // defdaults
	direction = 1;
    velocity = 1;
    minVelocity = 0.;
    smoothF = 0.9;
	smoothA = 1.;     // default is fixed value
    At=0.0;
	smoothErr = 0.0;
    Kp= 0.1;
    Ki = 0.;
    Kd = 0.;
    archive = NULL;

	// internal variables
    prevUpdateTime = 0.;
    startTime = 0.;
	numStartup=0;
    sumX=0.;
    sumY=0.;
    sumX2=0.;
    sumXY=0.;
	error=0.0;
    intError=0.0;
}

// Return name of this task
const char *LoadControl::TaskName(void) { return "Load Control Task"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
// throws CommonException()
char *LoadControl::InputParam(char *pName,int &input,double &gScaling)
{
   	// which material for global quantity
    if(strcmp(pName,"material")==0)
    {	input=INT_NUM;
        return (char *)&material;
    }
    else if(strcmp(pName,"matname") == 0)
    {   input = TEXT_PARAMETER;
        return (char *)&matName;
    }

    // direction
    else if(strcmp(pName,"direction")==0)
    {	input=INT_NUM;
        return (char *)&direction;
    }
    
    // starting velocity
    else if(strcmp(pName,"velocity")==0)
    {	input=DOUBLE_NUM;
        return (char *)&velocity;
    }
    
    // starting velocity
    else if(strcmp(pName,"minVelocity")==0)
    {   input=DOUBLE_NUM;
        return (char *)&minVelocity;
    }
    
   // starting velocity
    else if(strcmp(pName,"startupTime")==0)
    {   input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&startTime,gScaling,0.001);
    }
    
    // get load function
    else if(strcmp(pName,"Load")==0)
    {   input=TEXT_PARAMETER;
        return (char *)&Load_Function;
    }
    
    // Gain Factor 1
    else if(strcmp(pName,"Kp")==0)
    {	input=DOUBLE_NUM;
		return (char *)&Kp;
    }
	
	// Gain Factor 2
    else if(strcmp(pName,"Ki")==0)
    {	input=DOUBLE_NUM;
		return (char *)&Ki;
    }
	
    // Gain Factor 3
    else if(strcmp(pName,"Kd")==0)
    {	input=DOUBLE_NUM;
        return (char *)&Kd;
    }
    
    // archive PID calcs to file
    else if(strcmp(pName,"Archive")==0)
    {   input=TEXT_PARAMETER;
        return (char *)&archive;
    }
    
    // Smoothing factor
    else if(strcmp(pName,"smoothF")==0)
    {   input=DOUBLE_NUM;
        return (char *)&smoothF;
    }

    // Smoothing factor
    else if(strcmp(pName,"smoothA")==0)
    {   input=DOUBLE_NUM;
        return (char *)&smoothA;
    }

    //  Initial slope offset
    else if(strcmp(pName,"At")==0)
    {   input=DOUBLE_NUM;
        return (char *)&At;
    }
    
    // Smoothing factor
    else if(strcmp(pName,"smoothErr")==0)
    {    input=DOUBLE_NUM;
       return (char *)&smoothErr;
    }
   
	// not recognized - see if parent can handle
	return CustomTask::InputParam(pName,input,gScaling);
}

// Actually sets t4xt-based parameters
// throws std::bad_alloc, SAXException()
void LoadControl::SetTextParameter(char *fxn,char *ptr)
{
    // Load_Function
    if(ptr == (char *)&Load_Function)
    {    if(Load_Function!=NULL)
            ThrowSAXException("Duplicate load control function was supplied");
        if(fxn==NULL)
            ThrowSAXException("Load control function is missing");
        if(strlen(fxn)==0)
            ThrowSAXException("Load control function is empty");

        Load_Function = Expression::CreateExpression(fxn,"Load control function is not valid");
    }
    
    else if(ptr == (char *)&archive)
    {    if(archive!=NULL)
            ThrowSAXException("Duplicate archive file name was supplied");
        if(fxn==NULL)
            ThrowSAXException("Archive file name is missing");
        if(strlen(fxn)==0)
            ThrowSAXException("Archive file name is empty");
        
        // create file and store name
        archive = archiver->CreateFileInArchiveFolder(fxn);
        
        // save base name
        archiveBase = new char[strlen(fxn)+1];
        strcpy(archiveBase,fxn);
    }
        
    else if(ptr == (char *)&matName)
    {    // material name needed
        if(matName!=NULL)
            ThrowSAXException("Duplicate material name supplied to task");
        if(fxn==NULL)
            ThrowSAXException("Material name is missing");
        if(strlen(fxn)==0)
            ThrowSAXException("Material name is missing");
        
        // save is
        matName = new char[strlen(fxn)+1];
        strcpy(matName,fxn);
    }
    else
        CustomTask::SetTextParameter(fxn,ptr);
}

#pragma mark GENERIC TASK METHODS

// at beginning of analysis
// throws CommonException()
CustomTask *LoadControl::Initialize(void)
{
	int just_some_int = 0; // need this for a function input
    int whichMat;
	char load_we_want[11],disp_we_want[11];    
	
    // Check materials
    if(material <= 0 || material > nmat)
    {    if(matName == NULL)
            throw CommonException("No material defined", "LoadControl::Initialize");
        
        for(int i=0;i<nmat;i++)
        {    if(strcmp(matName,theMaterials[i]->name)==0)
            {    if(material>0)
                    throw CommonException("More than one material has the requested name", "LoadControl::Initialize");
                material = i+1;
            }
        }
        
        if(material <= 0 || material > nmat)
            throw CommonException("No material matching the specified name", "LoadControl::Initialize");
    }

    // convert direction
	RigidMaterial *ControlledMaterial = (RigidMaterial *)theMaterials[material-1];
	if(ControlledMaterial->IsRigidContact() == 0)
    {   if(direction==1)
        {   strcpy(load_we_want,"reactionx");
            globalLoadQuantity = GlobalQuantity::DecodeGlobalQuantity(load_we_want,&just_some_int,&whichMat);
		}
		else if(direction == 2)
        {   strcpy(load_we_want,"reactiony");
            globalLoadQuantity = GlobalQuantity::DecodeGlobalQuantity(load_we_want,&just_some_int,&whichMat);
        }
		else if(direction == 3)
        {   strcpy(load_we_want,"reactionz");
            globalLoadQuantity = GlobalQuantity::DecodeGlobalQuantity(load_we_want,&just_some_int,&whichMat);
		}
		else
        {   throw CommonException("Direction not valid","LoadControl::Initialize");
		}
	}
	else
    {   if(direction==1)
        {   strcpy(load_we_want,"contactx");
            globalLoadQuantity = GlobalQuantity::DecodeGlobalQuantity(load_we_want,&just_some_int,&whichMat);
		}
		else if(direction == 2)
        {   strcpy(load_we_want,"contacty");
            globalLoadQuantity = GlobalQuantity::DecodeGlobalQuantity(load_we_want,&just_some_int,&whichMat);
        }
		else if(direction == 3)
        {   strcpy(load_we_want,"contactz");
            globalLoadQuantity = GlobalQuantity::DecodeGlobalQuantity(load_we_want,&just_some_int,&whichMat);
		}
		else
        {   throw CommonException("Direction not valid","LoadControl::Initialize");
		}
	}
    ControlledMaterial->SetControlVelocity(velocity,direction);

    // get initial velocity
    if(DbleEqual(minVelocity,0.))
        minVelocity = fabs(0.05*velocity);
    else
        minVelocity = fabs(minVelocity);
    
    // startup time
    if(DbleEqual(startTime,0.))
        throw CommonException("The start up time cannot be zero","LoadControl::Initialize");
    if(startTime<0.)
    {   // use provided At, but cannot be zero
        if(DbleEqual(At,0.))
            throw CommonException("Non-zero At required when startupTime < 0","LoadControl::Initialize");
        startTime = -startTime;
        numStartup = -1;
    }
		
    // check F smoothing [0,1)
    if(smoothF<0.0||smoothF>=1.0)
    {   throw CommonException("Parmeter SmoothF is not in interval [0,1)","LoadControl::Initialize");
    }
    
	// Check At smoothing [0,1]
	if(smoothA<0.0||smoothA>1.0)
    {   throw CommonException("Parmeter SmoothA is not in interval [0,1]","LoadControl::Initialize");
	}
    
    // check error smoothing [0,1)
	if(smoothErr<0.0||smoothErr>=1.0)
    {   throw CommonException("Parmeter SmoothErr is not in interval [0,1)","LoadControl::Initialize");
	}
    
	// Find index of stored load.
	int LoadIndex=0;
	GlobalQuantity *nextGlobal=firstGlobal;
	while(nextGlobal!=NULL)
	{	if(nextGlobal->IsSameQuantity(globalLoadQuantity,just_some_int,material)) break;
		LoadIndex++;
		nextGlobal=nextGlobal->GetNextGlobal();
	}
	
	// error if not found
	if(nextGlobal==NULL)
    {   throw CommonException("LoadControl task's global quantity (Contact/Reaction) is not available in the current global archives","LoadControl::Initialize");
    }
	// store in quantity
    globalLoadQuantity = LoadIndex;
	
	// get displacement
	if(direction==1)
    {   strcpy(disp_we_want,"dispx");
        globalDistanceQuantity = GlobalQuantity::DecodeGlobalQuantity(disp_we_want,&just_some_int,&whichMat);
	}
	else if(direction == 2)
    {   strcpy(disp_we_want,"dispy");
        globalDistanceQuantity = GlobalQuantity::DecodeGlobalQuantity(disp_we_want,&just_some_int,&whichMat);
    }
	else if(direction == 3)
    {   strcpy(disp_we_want,"dispz");
        globalDistanceQuantity = GlobalQuantity::DecodeGlobalQuantity(disp_we_want,&just_some_int,&whichMat);
	}
		
	int DispIndex=0;
	nextGlobal=firstGlobal;
	while(nextGlobal!=NULL)
	{	if(nextGlobal->IsSameQuantity(globalDistanceQuantity,just_some_int,material)) break;
		DispIndex++;
		nextGlobal=nextGlobal->GetNextGlobal();
	}
	if(nextGlobal==NULL)
    {   throw CommonException("LoadControl task's global quantity for displacement is not available in the current global archives","LoadControl::Initialize");
    }
	// store in quantity
    globalDistanceQuantity= DispIndex;
	
	// Print the task parameters
	cout << "Control load in direction " << direction << " using " << load_we_want << endl;
    cout << "  Material #" << material << " (" << theMaterials[material-1]->name << ")" << endl;
	cout << "  Target load function =  " << Load_Function->GetString() << endl;
	cout << "  Starting Velocity = " << velocity << " " << UnitsController::Label(CUVELOCITY_UNITS);
    cout << " for " << startTime*UnitsController::Scaling(1000.)
                    << " " << UnitsController::Label(ALTTIME_UNITS) << endl;
    if(numStartup<0) cout << "  Start up At = " << At << endl;
#ifdef DYNAMIC_ADJUST
    cout << "  Minimum Velocity = " << minVelocity << " " << UnitsController::Label(CUVELOCITY_UNITS) << endl;
#endif
 	cout << "  PID parameters: Kp = " << Kp << "  Ki = " << Ki << "  Kd = "<< Kd<<endl;
    cout << "  Smoothing: F = " << smoothF << " A = " << smoothA << " Err = " << smoothErr << endl;
 	
	// archive control variables
	if(archive!=NULL)
	{	cout << "  Control Archive: " << archiver->GetArchiveRoot() << archiveBase << endl;
		
		// write header of the file
		ofstream carch;
		try
		{	carch.open(archive,ios::out | ios::app);
			if(!carch.is_open())
				archiver->FileError("File error opening PID control archive file",archive,"LoadControl::Initialize");
			
			carch << "# Control load in direction " << direction << " using " << load_we_want << endl;
            carch << "#   Material #" << material << " (" << theMaterials[material-1]->name << ")" << endl;
			carch << "#   Target load function =  " << Load_Function->GetString() << endl;
			carch << "#   Initial Velocity = " << velocity << " " << UnitsController::Label(CUVELOCITY_UNITS);
            carch << " for " << startTime*UnitsController::Scaling(1000.)
                            << " " << UnitsController::Label(ALTTIME_UNITS) << endl;
            if(numStartup<0) carch << "#   Start up At = " << At << endl;
#ifdef DYNAMIC_ADJUST
            carch << "#   Minimum Velocity = " << minVelocity << " " << UnitsController::Label(CUVELOCITY_UNITS) << endl;
#endif
			carch << "#   PID parameters: Kp = " << Kp << "  Ki = " << Ki << "  Kd = " << Kd<< endl;
            carch << "#   Smoothing: F = " << smoothF << " A = " << smoothA << " Err = " << smoothErr << endl;
            carch << "#setColor\tblack\tred\tblue\tgreen\tbrown\tpurple" << endl;
            carch << "#setName\tAt\terror\tInt_error\tD_error\tvelocity\tforce" << endl;
 
			if(carch.bad())
				archiver->FileError("File error writing PID control archive file",archive,"LoadControl::Initialize");
			carch.close();
			if(carch.bad())
				archiver->FileError("File error closing PID control archive file",archive,"LoadControl::Initialize");
		}
		
		catch(CommonException& err)
		{   // divert to standard output and try to continue
			cout << "# PID control file error - check disk for amount of free space" << endl;
			cout << "# " << err.Message() << endl;
			if(carch.is_open()) carch.close();
		}
	}

    return nextTask;
}

// Update Velocity
CustomTask *LoadControl::StepCalculation(void)
{
	// Only control last step ended with a global archive update
    if(fmobj->mstep!=archiver->GetLastArchivedStep()+1) return nextTask;
    
    // skip the first one with zero data
    if(fmobj->mstep<=2) return nextTask;
    
    // time since last update
    double dTime = mtime-prevUpdateTime;
    prevUpdateTime = mtime;
	
    // Get Load and update distance
	// Get the load we want
	double Load_needed = Load_Function->TValue(mtime*UnitsController::Scaling(1000.));          // get desired load
	double currentLoad = archiver->GetLastArchived(globalLoadQuantity);            // get current load
    prevLoad = load;
    load = smoothF*prevLoad + (1.-smoothF)*currentLoad;
	double distance = archiver->GetLastArchived(globalDistanceQuantity);    // get current displacement

	// don't do anything for first 20 iterations to get an estimate for parameters
	if(mtime<startTime)
    {   // keep going until get some load
		if(fabs(load)>1e-7 && numStartup>=0)
        {   numStartup++;
            sumX += distance;
            sumX2 += distance*distance;
            sumY += load;
            sumXY += load*distance;
		}
        
		// Put starting velocity in Rigid Material (Function has to be NULL)
		//RigidMaterial *ControlledMaterial = (RigidMaterial *)theMaterials[material-1];
		//ControlledMaterial->SetControlVelocity(velocity,direction);
		return nextTask;
	}
    else if(numStartup>=0)
    {   // note that zero is an error here
        if(numStartup==0)
        {   throw CommonException("Load control start up time did not include and global archives",
                                  "LoadControl::StepCalculation");
        }
        double n = (double)numStartup;
        At = (n*sumXY - sumX*sumY)/(n*sumX2 - sumX*sumX);
        //double b = (sumY - At*sumX)/n;
        // cout << "# First n = " << numStartup << ", At = " << At << ", b = " << b << endl;
        numStartup=-1;
    }
    else
    {   // adjust slope and intercept to fit F = At*(distance-Bt)
        // .... but just need dF/d(distance) = At
#ifdef DYNAMIC_ADJUST
        // adjust only if velocity not too small
        if(fabs(velocity)>minVelocity)
        {   double localAt = (load-prevLoad)/(velocity*dTime);
            At = smoothA*At + (1.0-smoothA)*localAt;
        }
#else
        // At from load and displacement and then smooth
        At = smoothA*At + (1.0-smoothA)*(load/distance);
#endif
    }
	
	// Calculate error (mm or L) with exponential smoothing
	double raw_error = (Load_needed-load)/At;
    prevError = error;
	error = (1.0-smoothErr)*raw_error + smoothErr*prevError;
    
    // Integrated error divided by dTime (mm-s or L-t)
    intError += error;
	
    // Derivative (time units are mm/sec or L/t) * dTime
    double dedt = error-prevError;

	// calculate scaled error for PI controller
    velocity = (Kp*error + Ki*intError + Kd*dedt)/dTime;		// control output (velocity (mm/s)

	// archive control variables
	if(archive!=NULL)
	{	// write to file
		ofstream carch;
		try
		{	carch.open(archive,ios::out | ios::app);
			if(!carch.is_open())
				archiver->FileError("File error opening PID control archive file",archive,"LoadControl::Initialize");
			
            // for Legacy, time in ms, At in F/distance
            // error in mm, integrated error in mm-msec, derivative in mm/ms
            // velocity in mm/sec
            double tscale = UnitsController::Scaling(1000.);
			carch << mtime*tscale;
            carch << "\t" << At << "\t" << error;
            carch << "\t" << intError << "\t" << dedt << "\t" << velocity << "\t" << load << endl;
			if(carch.bad())
				archiver->FileError("File error writing PID control archive file",archive,"LoadControl::Initialize");
			
			carch.close();
			if(carch.bad())
				archiver->FileError("File error closing PID control archive file",archive,"LoadControl::Initialize");
		}
		
		catch(CommonException& err)
		{   // divert to standard output and try to continue
			cout << "# PID control file error - check disk for amount of free space" << endl;
			cout << "# " << err.Message() << endl;
			if(carch.is_open()) carch.close();
		}
	}
	
	// update velocity in rigid material
	RigidMaterial *ControlledMaterial = (RigidMaterial *)theMaterials[material-1];
	ControlledMaterial->SetControlVelocity(velocity,direction);
	
	// done
	return nextTask;
}

