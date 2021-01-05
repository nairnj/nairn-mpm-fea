/*********************************************************************
    OutputSetup.cpp
    Nairn Research Group MPM Code
    
    Created by jnairn on Wed Jan 23 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
*********************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "System/ArchiveData.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Exceptions/CommonException.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Cracks/CrackHeader.hpp"
#include "System/UnitsController.hpp"
#include "Read_MPM/MPMReadHandler.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "Boundary_Conditions/MatPtHeatFluxBC.hpp"
#include "Boundary_Conditions/InitialCondition.hpp"

/*********************************************************************
    Print setup information
*********************************************************************/

// print title of results file
void NairnMPM::PrintAnalysisTitle(void)
{
	cout << "MPM ANALYSIS BY ";
	CoutCodeVersion();
}

// print analysis method details
// throws CommonException()
void NairnMPM::PrintAnalysisMethod(void)
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
        case USL_METHOD:
            if(mpmApproach==USL_METHOD) cout << "USL";
			if(fmobj->skipPostExtrapolation)
				cout << " (-2nd extrap)";
			else
				cout << " (+2nd extrap)";
            break;
        default:
            throw CommonException("Invalid MPM analysis method provided.","NairnMPM::PrintAnalysisMethod");
    }
	
    switch(ElementBase::useGimp)
	{	case POINT_GIMP:
			cout << " / Classic";
			break;
		case UNIFORM_GIMP:
		case UNIFORM_GIMP_AS:
            cout << " / GIMP";
            break;
        case LINEAR_CPDI:
		case LINEAR_CPDI_AS:
            cout << " / Linear CPDI";
			if(ElementBase::rcrit>0.) cout << " / rcrit = " << ElementBase::rcrit;
            break;
       case QUADRATIC_CPDI:
            cout << " / Quadratric CPDI";
			if(ElementBase::rcrit>0.) cout << " / rcrit = " << ElementBase::rcrit;
            break;
		case BSPLINE_GIMP:
		case BSPLINE_GIMP_AS:
			cout << " / B2GIMP";
			break;
		case BSPLINE_CPDI:
			cout << " / B2CPDI";
			if(ElementBase::rcrit>0.) cout << " / rcrit = " << ElementBase::rcrit;
			break;
		case BSPLINE:
			cout << " / B2SPLINE";
			break;
		default:
			cout << " / (unknown shape function method)";
    }
	cout << endl;

	// incremental F terms
	if(MaterialBase::incrementalDefGradTerms<0)
	{	MaterialBase::incrementalDefGradTerms = - MaterialBase::incrementalDefGradTerms;
		if(IsThreeD())  MaterialBase::incrementalDefGradTerms = 1;
	}
	cout << "Incremental F Terms: " << MaterialBase::incrementalDefGradTerms << endl;
	
	// time step and max time
    cout << "Time step: min(" << timestep*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS) << ", "
		<< fmobj->GetCFLCondition() << " time for wave to cross one cell)\n";
    cout << "Maximum time: " << maxtime*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS) << "\n\n";
}

// rest of MPM information
// throws CommonException()
void NairnMPM::CMStartResultsOutput(void)
{
    //---------------------------------------------------
    // Gravity and Damping and Thermal Ramp
    PrintSection("GRAVITY, DAMPING, AND THERMAL");
	bodyFrc.Output();
	thermal.Output();
    ConductionTask::ThermodynamicsOutput();
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
					throw CommonException("Propagation traction law material is not defined","NairnMPM::CMStartResultsOutput");
				if(theMaterials[propagateMat[i]-1]->MaterialStyle()!=TRACTION_MAT)
					throw CommonException("Propagation traction law is not a traction law material","NairnMPM::CMStartResultsOutput");
			}
		}
		
		// contact output and allocations
		contact.Output();
		
		// location for COD calculation - range is 0 to 3 segments from crack tip
		// future may want to read this as parameter
		// Only used to partition J or to implement COD-based criteria
		CrackHeader::SetCodLocation(2.);
		cout << "Crack COD found (when needed) " << 2 << " segments from crack tip" << endl;
		
		// crack details
		cout << "Number of cracks = " << numberOfCracks << endl;
		CrackHeader *nextCrack=firstCrack;
		while(nextCrack!=NULL)
		{	nextCrack->Output();
			nextCrack=(CrackHeader *)nextCrack->GetNextObject();
		}
		cout << endl;
    }
	else
	{	// turn off request for position extrapolation
		contact.crackContactByDisplacements = true;
	}
    
    //---------------------------------------------------
    // Multimaterial Contact
	if(multiMaterialMode)
    {   PrintSection("MULTIMATERIAL CONTACT");
		mpmgrid.MaterialOutput();
		int i;
		for(i=0;i<nmat;i++)
			theMaterials[i]->ContactOutput(i+1);
	}
	
    //---------------------------------------------------
    // Fixed Displacements
    if(firstVelocityBC!=NULL)
    {   PrintSection("NODAL POINTS WITH FIXED DISPLACEMENTS");
        archiver->ArchiveVelocityBCs(firstVelocityBC);
    }
    
}

// output boundary conditions (called in PreliminaryParticleCalcs())
void NairnMPM::OutputBCMassAndGrid(void)
{
	//---------------------------------------------------
	// Finish Results File
	BoundaryCondition *nextBC;
	
	//---------------------------------------------------
	// Loaded Material Points
	if(firstLoadedPt!=NULL)
	{   PrintSection("MATERIAL POINTS WITH EXTERNAL FORCES");
		cout << "Point   DOF ID     Load (" << UnitsController::Label(FEAFORCE_UNITS) << ")     Arg ("
				<< UnitsController::Label(BCARG_UNITS) << ")  Function\n"
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
		cout << "Point   DOF Face ID   Stress (" << UnitsController::Label(PRESSURE_UNITS) << ")    Arg ("
				<< UnitsController::Label(BCARG_UNITS) << ")  Function\n"
				<< "----------------------------------------------------------------\n";
		nextBC=(BoundaryCondition *)firstTractionPt;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
	
	//---------------------------------------------------
	// Diffusion boundary conditions
	if(fmobj->HasFluidTransport() && firstConcBC!=NULL)
	{	if(fmobj->HasDiffusion())
		{	PrintSection("NODAL POINTS WITH FIXED CONCENTRATIONS");
			cout << "  Node  ID    Conc (/csat)   Arg (" << UnitsController::Label(BCARG_UNITS) << ")  Function";
		}
#ifdef POROELASTICITY
		else
		{	PrintSection("NODAL POINTS WITH FIXED PORE PRESSURE");
			cout << "  Node  ID    Press. (" << UnitsController::Label(PRESSURE_UNITS) << ")    Arg ("
					<< UnitsController::Label(BCARG_UNITS) << ")  Function";
		}
#endif
		cout << "\n------------------------------------------------------\n";
		nextBC=firstConcBC;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
	
	//---------------------------------------------------
	// Concentration Flux Material Points
	if(fmobj->HasFluidTransport() && firstFluxPt!=NULL)
	{	if(fmobj->HasDiffusion())
		{	PrintSection("MATERIAL POINTS WITH CONCENTRATION FLUX");
			cout << " Point  DOF Face ID Flux (" << UnitsController::Label(BCCONCFLUX_UNITS) << ") Arg ("
					<< UnitsController::Label(BCARG_UNITS) << ")  Function";
		}
#ifdef POROELASTICITY
		else
		{	PrintSection("NODAL POINTS WITH PORE PRESSURE FLUX");
			cout << " Point  DOF Face ID Flux(v_F/(" << UnitsController::Label(CULENGTH_UNITS) << "^2-"
					<< UnitsController::Label(TIME_UNITS)<< ")) Arg ("
					<< UnitsController::Label(BCARG_UNITS) << ") Function";
		}
#endif
		cout << "\n---------------------------------------------------------------\n";
		nextBC=(BoundaryCondition *)firstFluxPt;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
	
	//---------------------------------------------------
	// Conduction boundary conditions
	if(ConductionTask::active && firstTempBC!=NULL)
	{   PrintSection("NODAL POINTS WITH FIXED TEMPERATURES");
		cout << " Node   ID   Temp (-----)   Arg (" << UnitsController::Label(BCARG_UNITS) << ")  Function\n"
				<< "------------------------------------------------------\n";
		nextBC=firstTempBC;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
	
	//---------------------------------------------------
	// Heat Flux Material Points
	if(ConductionTask::active && firstHeatFluxPt!=NULL)
	{	PrintSection("MATERIAL POINTS WITH HEAT FLUX");
		cout << " Point  DOF Face ID   Flux (" << UnitsController::Label(BCHEATFLUX_UNITS) << ")    Arg ("
				<< UnitsController::Label(BCARG_UNITS) << ")  Function\n"
				<< "---------------------------------------------------------------\n";
		nextBC=(BoundaryCondition *)firstHeatFluxPt;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
	
	//---------------------------------------------------
	// Damaged Material Points
	if(firstDamagedPt!=NULL)
	{   PrintSection("MATERIAL POINTS HAVE INITIAL DAMAGE");
		cout << "Plot history #1 to see initially damage particles" << endl;
		cout << "Those fully failed have history #1 = 3" << endl;
		cout << endl;
		
		// assign initial conditions and then delete
		InitialCondition *nextIC=firstDamagedPt,*prevIC;
		while(nextIC!=NULL)
		{   prevIC = nextIC;
			nextIC = prevIC->AssignInitialConditions(IsThreeD());
			delete prevIC;
		}
	}
	
	// Print particle information and other preliminary calc results
	PrintSection("FULL MASS MATRIX");
	
	char fline[200];
	sprintf(fline,"Number of Material Points: %d",nmpms);
	cout << fline << endl;
	
	// background grid info
	mpmgrid.Output(IsAxisymmetric());
	
	sprintf(fline,"Adjusted time step (%s): %.7e",UnitsController::Label(ALTTIME_UNITS),timestep*UnitsController::Scaling(1.e3));
	cout << fline << endl;
	
	// propagation time step and other settings when has cracks
	if(firstCrack!=NULL)
	{	if(propagate[0])
		{   sprintf(fline,"Propagation time step (%s): %.7e",UnitsController::Label(ALTTIME_UNITS),propTime*UnitsController::Scaling(1.e3));
			cout << fline << endl;
		}
	}
	
	// blank line
	cout << endl;
}

#pragma mark NairnMPM: Output accessors

// return name
const char *NairnMPM::CodeName(void) const
{
#ifdef OSPARTICULAS
	return "OSParticulas";
#else
	return "NairnMPM";
#endif
}

// MPM analysis type
void NairnMPM::GetAnalysisType(int np,char *fline)
{
	// analysis type
	switch(np)
	{	case PLANE_STRAIN_MPM:
			strcpy(fline,"2D Plane Strain MPM");
			break;
			
		case PLANE_STRESS_MPM:
			strcpy(fline,"2D Plane Stress MPM");
			break;
			
		case THREED_MPM:
			strcpy(fline,"3D MPM");
			break;
			
		case AXISYMMETRIC_MPM:
			strcpy(fline,"Axisymmetric MPM");
			break;
			
		default:
			throw CommonException("No FEA or MPM analysis type was provided.","CommonAnalysis::StartResultsOutput");
	}
	
	// Add augmentation
	strcat(fline,MPMAugmentation());
	strcat(fline," Analysis");
}

// string about MPM additions (when printing out analysis type
const char *NairnMPM::MPMAugmentation(void)
{
#ifdef IMPLICIT_STEP
	if (exactTractions) {
		if (implicit_step) {
			return " ET Implicit";
		}
		else {
			if (implicit_step)
				return " Implicit";
		}
	}
#endif
	
#ifdef ADD_PARTICLE_SPIN
	// With exact tracitons
	// options are none, +PS, +TGV +ET +PS+ET +TGV+ET
	if(plusParticleSpin)
	{	if(exactTractions)
			return " +PS(Bp) +ET";
		else
			return " +PS(Bp)";
	}
	else if(exactTractions)
		return " +ET";
#else
	if(exactTractions)
		return " +ET";
#endif
	
	// no features
	return "";
}

// title for nodes and elements section
const char *NairnMPM::NodesAndElementsTitle(void) const { return "NODES AND ELEMENTS (Background Grid)"; }

// Get file read handler
CommonReadHandler *NairnMPM::GetReadHandler(void)
{
	MPMReadHandler *handler=new MPMReadHandler();
	return handler;
}

// Explain usage of this program
void NairnMPM::Usage()
{
	CoutCodeVersion();
	cout << "    Expects Xerces 3.1.1 or newer" << endl;
	cout << "\nUsage:\n"
	"    NairnMPM [-options] <InputFile>\n\n"
	"This program reads the <InputFile> and does an MPM analysis.\n"
	"The results summary is directed to standard output. Errors are\n"
	"directed to standard error. The numerical results are saved\n"
	"to archived files as directed in <InputFile>.\n\n"
	"Options:\n"
	"    -a          Abort after setting up problem but before\n"
	"                   MPM steps begin. Initial conditions will\n"
	"                   be archived\n"
	"    -H          Show this help\n"
	"    -np #       Set number of processors for parallel code\n"
	"    -r          Reverse byte order in archive files\n"
	"                   (default is to not reverse the bytes)\n"
	"    -v          Validate input file if DTD is provided in !DOCTYPE\n"
	"                   (default is to skip validation)\n"
	"    -w          Output results to current working directory\n"
	"                   (default is output to folder of input file)\n"
	"    -?          Show this help\n\n"
	"See http://osupdocs.forestry.oregonstate.edu/index.php/Main_Page for documentation\n\n"
	<<  endl;
}

