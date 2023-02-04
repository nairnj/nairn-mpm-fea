/********************************************************************************
	DeleteDamaged.cpp

	Created by Chad Hammerquist October 2017

	This task deletes damaged particles
	Parameters
*********************************************************************************/

#include "stdafx.h"
#include "MPM_Classes/MPMBase.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "NairnMPM_Class/Reservoir.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/DeleteDamaged.hpp"
#include "System/UnitsController.hpp"
#include "Materials/MaterialBase.hpp"
#include "Materials/Elastic.hpp"

extern double damageState;

// Constructors
DeleteDamaged::DeleteDamaged() : CustomTask()
{
	//  Initial values
	material = -1;
	matName = NULL;
	
	// direction and optional cod requirement
	damage_direction = TOTAL_COD;
	minRelativeCod = -1.;
    
    // track the number
    numDeleted = 0;
    
    // no deleted by time
    deleteTime = -1.;
}

// Return name of this task
const char *DeleteDamaged::TaskName(void) { return "Delete damaged particles"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *DeleteDamaged::InputParam(char *pName, int &input, double &gScaling)
{
 	if(strcmp(pName,"material") == 0)
	{
		input = INT_NUM;
		return (char *)&material;
	}
	else if(strcmp(pName,"matname") == 0)
	{
		input = TEXT_PARAMETER;
		return (char *)&matName;
	}
	else if(strcmp(pName,"direction") == 0)
	{	// must be 1,2,3,4,5 (see constants), anything else changed to TOTAL_COD=4
		input = INT_NUM;
		return (char *)&damage_direction;
	}
	else if(strcmp(pName,"minCOD") == 0)
	{	// inactive unless entered > 0
		input = DOUBLE_NUM;
		return (char *)&minRelativeCod;
	}
    else if(strcmp(pName,"deleteTime") == 0)
    {   // delete after this time (inactive unless entered > 0)
        input = DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&deleteTime,gScaling,1.e-3);
    }
	
	return CustomTask::InputParam(pName,input,gScaling);

}

// Get material by name instead of number
// throws std::bad_alloc, SAXException()
void DeleteDamaged::SetTextParameter(char *tdata,char *ptr)
{
	if(ptr == (char *)&matName)
	{	// material name needed
		if(matName!=NULL)
			ThrowSAXException("Duplicate material name supplied to task");
		if(tdata==NULL)
			ThrowSAXException("Material name is missing");
		if(strlen(tdata)==0)
			ThrowSAXException("Material name is missing");
		
		// save is
		matName = new char[strlen(tdata)+1];
		strcpy(matName,tdata);
	}
	
	else
		CustomTask::SetTextParameter(tdata,ptr);
}

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *DeleteDamaged::Initialize(void)
{
	// requires reservoir
	if(mpmReservoir==NULL)
		throw CommonException("DeleteDamaged custom task requires an available reservoir (i.e., a structured grid)",
							  "DeleteDamaged::Initialize");
    
	// Check materials
	if(material <= 0 || material > nmat)
	{	if(matName == NULL)
			throw CommonException("No material defined", "DeleteDamaged::Initialize");
		
		for(int i=0;i<nmat;i++)
		{	if(strcmp(matName,theMaterials[i]->name)==0)
			{	if(material>0)
					throw CommonException("More than one material has the requested name", "DeleteDamaged::Initialize");
				material = i+1;
			}
		}
		
		if(material <= 0 || material > nmat)
			throw CommonException("No material matching the specified name", "DeleteDamaged::Initialize");
	}
	
	// validate softening material (but may by NULL here)
	
	// output
	cout << "Delete particles after decohesion or specified time" << endl;
	cout << "   Removes failed particles of material: " << theMaterials[material-1]->name << " (" << material << ")" << endl;
    if(minRelativeCod>0.)
    {   cout << "   ... but not removed until ";
        if (damage_direction == TENSILE_COD)
        {	cout << "normal direction COD";
        }
        else if (damage_direction == SHEARXY_COD)
        {	cout << "xy shear slippage relative to crack plane";
        }
        else if (damage_direction == SHEARXZ_COD)
        {	cout << "xz shear slippage relative to crack plane";
        }
        else if (damage_direction == SHEAR_COD)
        {	cout << "tangential direction COD";
        }
        else
        {	// convert all other to any direction
            damage_direction = TOTAL_COD;
            cout << "magnitude of COD";
        }
        const char *label = UnitsController::Label(CULENGTH_UNITS);
		cout << " is at least " << minRelativeCod << " " << label << endl;
	}
    if(deleteTime>0.)
    {   cout << "   Removes particles at time " << deleteTime*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS) << endl;
    }

	cout << endl;

	return nextTask;
}


// do custom calculation
CustomTask *DeleteDamaged::StepCalculation(void)
{
#pragma omp parallel for
	for(int p = 0; p<nmpms; p++)
    {   MPMBase *point = mpm[p];
		if(point->InReservoir()) continue;
		
		// only look a particles of the selected material (convert to 1 based)
		if((point->MatID() + 1) != material) continue;
        
        // is it time to delete?
        if((deleteTime>0.) && (mtime>=deleteTime))
        {
#pragma omp critical (delparticle)
			{
				mpmReservoir->DeleteParticle(point);
				numDeleted++;
			}
            continue;
        }
		
		// get softening data, but skip if not softening material
        // This code assume all softening materials use the same soft history
        //    variable to make when failre
        int failureSurface = theMaterials[material-1]->GetTractionFailureSurface();
        if(failureSurface<0) continue;
		double *soft = theMaterials[material-1]->GetSoftHistoryPtr(point);
		
        // is it failed
        if(soft[SOFT_DAMAGE_STATE]>damageState)
		{	// do we care about cod too?
            double magCod = minRelativeCod+1.;
			if(minRelativeCod>0.)
			{	// Optional delay of deletion unless cod exceeds input value
				const MaterialBase *matref = theMaterials[point->MatID()];
				Vector cod = matref->GetCrackingCOD(point,!(fmobj->IsThreeD()));
				
				if (damage_direction == TENSILE_COD)
					magCod = cod.x;
				else if (damage_direction == SHEARXY_COD)
					magCod = fabs(cod.y);
				else if (damage_direction == SHEARXZ_COD)
					magCod = fabs(cod.z);
				else if (damage_direction == SHEAR_COD)
					magCod = sqrt(cod.y*cod.y+cod.z+cod.z);
				else
					magCod = sqrt(DotVectors(&cod,&cod));
			}
			
			// delete the particle
			if(magCod>minRelativeCod)
			{
#pragma omp critical (delparticle)
				{
					mpmReservoir->DeleteParticle(point);
					numDeleted++;
				}
 			}
		}
	}

	// next
	return nextTask;
}

// called after the analysis is done. A task can output a report to the main.mpm
// file. It should end with an empty line
bool DeleteDamaged::HasReport(void) { return true; }
CustomTask *DeleteDamaged::Report(void)
{
    cout << TaskName() << endl;
    cout << "    Number of particles deleted: " << numDeleted << endl;
    cout << endl;
    return nextTask;
}
