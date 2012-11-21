/*********************************************************************
    OutputSetup.cpp
    Nairn Research Group MPM Code
    
    Created by jnairn on Wed Jan 23 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
*********************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "System/ArchiveData.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Exceptions/CommonException.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "Cracks/CrackHeader.hpp"

/*********************************************************************
    Print setup information
*********************************************************************/

// print title of results file
void NairnMPM::PrintAnalysisTitle(void)
{
	cout << "MPM ANALYSIS BY ";
	CoutCodeVersion();
}

// print analysis type
void NairnMPM::PrintAnalysisType(void)
{
	// new section
    PrintSection("MPM ANALYSIS");
    
	// MPM method
    cout << "MPM Analysis Method: ";
    switch(mpmApproach)
    {	case USF_METHOD:
            cout << "USF";
            break;
        case USAVG_METHOD:
            cout << "USAVG";
            break;
        case SZS_METHOD:
            cout << "SZS";
            break;
        default:
            throw CommonException("Invalid MPM analysis method provided.","NairnMPM::PrintAnalysisType");
    }
	
    switch(ElementBase::useGimp)
    {   case UNIFORM_GIMP:
		case UNIFORM_GIMP_AS:
            cout << " / GIMP";
            break;
        case LINEAR_CPDI:
		case LINEAR_CPDI_AS:
            cout << " / Linear CPDI";
            break;
        case QUADRATIC_CPDI:
            cout << " / Quadratric CPDI";
            break;
    }
	cout << endl;
    
	// time step and max time
    cout << "Time step: min(" << 1000.*timestep << " ms, "
        << fmobj->GetCFLCondition() << " time for wave to cross one cell)\n";
    cout << "Maximum time: " << 1000.*maxtime << " ms\n\n";
}

// rest of MPM information
void NairnMPM::MyStartResultsOutput(void)
{
	BoundaryCondition *nextBC;
    
    //---------------------------------------------------
    // Gravity and Damping and Thermal Ramp
    PrintSection("GRAVITY, DAMPING, AND TRANSPORT");
	bodyFrc.Output();
	thermal.Output();
	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport=nextTransport->TransportOutput();
    cout << endl;

    //---------------------------------------------------
    // Cracks
    if(firstCrack!=NULL)
    {   PrintSection("CRACKS AND CRACK CONTACT");
		int i;
		for(i=0;i<=1;i++)
		{	if(i==0)
				cout << "Default propagation criterion: ";
			else
			{	if(propagate[i]==NO_PROPAGATION) break;
				cout << "Alternate propagation criterion: ";
			}
			switch(propagate[i])
			{   case NO_PROPAGATION:
					cout << "No propagation" << endl;
					break;
				case MAXHOOPSTRESS:
					cout << "Maximum hoop stress" << MaterialBase::PreferredDirection(propagateDirection[i]) << endl;
					break;
				case STEADYSTATEGROWTH:
					cout << "Constant crack speed" << MaterialBase::PreferredDirection(propagateDirection[i]) << endl;
					break;
				case CRITICALERR:
					cout << "Critical energy release rate" << MaterialBase::PreferredDirection(propagateDirection[i]) << endl;
					break;
				case TOTALENERGYBALANCE:
					cout << "Total energy balance" << MaterialBase::PreferredDirection(propagateDirection[i]) << endl;
					break;
				case STRAINENERGYDENSITY:
					cout << "Minimum strain energy density" << MaterialBase::PreferredDirection(propagateDirection[i]) << endl;
					break;
				case EMPIRICALCRITERION:
					cout << "Empirical criterion" << MaterialBase::PreferredDirection(propagateDirection[i]) << endl;
					break;
				case MAXCTODCRITERION:
					cout << "Maximum CTOD criterion " << MaterialBase::PreferredDirection(propagateDirection[i]) << endl;
					break;
				default:
					cout << "Unknown criterion" << endl;
					break;
			}
			if(propagate[i]!=NO_PROPAGATION && propagateMat[i]>0)
			{	cout << "   New crack surface has traction law material " << propagateMat[i] << endl;
				if(propagateMat[i]>nmat)
					throw CommonException("Propagation traction law material is not defined","NairnMPM::MyStartResultsOutput");
				if(!theMaterials[propagateMat[i]-1]->isTractionLaw())
					throw CommonException("Propagation traction law is not a traction law material","NairnMPM::MyStartResultsOutput");
			}
		}
		
		// default crack contact loaw
		int numberOfCracks=firstCrack->Count();
		contact.Output(numberOfCracks);
		
		// crack details
		cout << "Number of cracks = " << numberOfCracks << endl;
		CrackHeader *nextCrack=firstCrack;
		while(nextCrack!=NULL)
		{	nextCrack->Output();
			nextCrack=(CrackHeader *)nextCrack->GetNextObject();
		}
		cout << endl;
    }
    
    //---------------------------------------------------
    // Multimaterial Contact
	if(multiMaterialMode)
    {   PrintSection("MULTIMATERIAL CONTACT");
		contact.MaterialOutput();
		int i;
		for(i=0;i<nmat;i++)
			theMaterials[i]->ContactOutput(i+1);
		cout << endl;
	}
	
    //---------------------------------------------------
    // Fixed Displacements
    if(firstVelocityBC!=NULL)
    {   PrintSection("NODAL POINTS WITH FIXED DISPLACEMENTS");
        archiver->ArchiveVelocityBCs(firstVelocityBC);
    }
    
    //---------------------------------------------------
    // Loaded Material Points
    if(firstLoadedPt!=NULL)
    {   PrintSection("MATERIAL POINTS WITH EXTERNAL FORCES");
        cout << "Point   DOF ID     Load (N)     Arg (ms/ms^-1)  Function\n"
        << "----------------------------------------------------------\n";
        nextBC=(BoundaryCondition *)firstLoadedPt;
        while(nextBC!=NULL)
            nextBC=nextBC->PrintBC(cout);
        cout << endl;
    }
	
	//---------------------------------------------------
    // Traction Loaded Material Points
    if(firstTractionPt!=NULL)
    {   PrintSection("MATERIAL POINTS WITH TRACTIONS");
        cout << "Point   DOF Face ID   Stress (MPa)    Arg (ms/ms^-1)  Function\n"
        << "----------------------------------------------------------------\n";
        nextBC=(BoundaryCondition *)firstTractionPt;
        while(nextBC!=NULL)
            nextBC=nextBC->PrintBC(cout);
        cout << endl;
    }
	
   //---------------------------------------------------
    // Diffusion boundary conditions
	if(DiffusionTask::active)
	{   PrintSection("NODAL POINTS WITH FIXED CONCENTRATIONS");
		cout << " Node ID   Conc (/csat)   Arg (ms/ms^-1)  Function\n"
		<< "----------------------------------------------------\n";
		nextBC=firstConcBC;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
		
		//---------------------------------------------------
		// Concentration Flux Material Points
		PrintSection("MATERIAL POINTS WITH CONCENTRATION FLUX");
		cout << " Point DOF ID      Flux ( )     Arg (ms/ms^-1)\n"
		<< "-------------------------------------------------\n";
		nextBC=(BoundaryCondition *)firstFluxPt;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
	
	//---------------------------------------------------
    // Conduction boundary conditions
	if(ConductionTask::active)
	{   PrintSection("NODAL POINTS WITH FIXED TEMPERATURES");
		cout << " Node ID   Temp (-----)   Arg (ms/ms^-1)  Function\n"
		<< "----------------------------------------------------\n";
		nextBC=firstTempBC;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
}
