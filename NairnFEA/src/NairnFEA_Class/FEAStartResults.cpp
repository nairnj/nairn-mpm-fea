/*********************************************************************
    BeginResults.cpp
    Nairn Research Group FEA Code
    
    Created by jnairn on Mon Feb 2 2003.
    Copyright (c) 2003, All rights reserved.    
*********************************************************************/

#include "stdafx.h"
#include "NairnFEA_Class/NairnFEA.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/NodalDispBC.hpp"
#include "Boundary_Conditions/NodalLoad.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Boundary_Conditions/EdgeBC.hpp"
#include "System/UnitsController.hpp"
#include "Read_XML/Expression.hpp"
#include "Read_FEA/FEAReadHandler.hpp"

/*********************************************************************
    Begin results information
*********************************************************************/

// print title of results file
void NairnFEA::PrintAnalysisTitle(void)
{
	cout << "FEA ANALYSIS BY ";
	CoutCodeVersion();
}

// finish start of FEA results file
// throws std::bad_alloc, CommonException()
void NairnFEA::CMStartResultsOutput(void)
{
    char hline[200];
    size_t hsize=200;
    int i;
    
    //---------------------------------------------------
    // Temperature
    PrintSection("THERMAL LOAD");
	if(temperatureExpr!=NULL)
	{	if(!Expression::CreateFunction(temperatureExpr,1))
			throw CommonException("The temperature expression is not a valid function","NairnFEA::CMStartResultsOutput");
		for(i=1;i<=nnodes;i++)
		{	nd[i]->gTemperature = Expression::FunctionValue(1,nd[i]->x,nd[i]->y,0.,0.,0.,0.)-stressFreeTemperature;
		}
		Expression::DeleteFunction(1);
	}
	else
	{	// unknown, but may have been set in explicit node commands
		temperatureExpr=new char[2];
		strcpy(temperatureExpr,"?");
	}
	
	snprintf(hline,hsize,"T0: %.2lf C\nT: %s C",stressFreeTemperature,temperatureExpr);
	cout << hline << endl;
	if(stressFreeTemperature!=0)
	{	snprintf(hline,hsize,"T-T0: %s - %.2lf C",temperatureExpr,stressFreeTemperature);
		cout << hline << endl;
	}
	cout << endl;
    
    //---------------------------------------------------
    // Fixed Displacements
	if(firstDispBC!=NULL)
	{	PrintSection("NODAL POINTS WITH FIXED DISPLACEMENTS");
		cout << " Node  DOF  Displacement (" << UnitsController::Label(CULENGTH_UNITS) << ")  Axis        Angle\n"
		<< "----------------------------------------------------\n";
		NodalDispBC *nextBC=firstDispBC;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
    
    //---------------------------------------------------
    // Loaded Nodes
	if(firstLoadBC!=NULL)
	{	PrintSection("NODAL POINTS WITH APPLIED LOADS");
		cout << " Node  DOF       Load (" << UnitsController::Label(FEAFORCE_UNITS) << ")\n"
			<< "------------------------------\n";
		NodalLoad *nextLoad=firstLoadBC;
		while(nextLoad!=NULL)
			nextLoad=nextLoad->PrintLoad();
		cout << endl;
	}

    //---------------------------------------------------
    // Stress element faces
	if(firstEdgeBC!=NULL)
	{	snprintf(hline,hsize,"FACES WITH APPLIED STRESS (%s)",UnitsController::Label(PRESSURE_UNITS));
		PrintSection(hline);
		cout << " Elem  Fc   Type      Nodal Stress     Nodal Stress     Nodal Stress\n"
			<< "----------------------------------------------------------------------\n";
		EdgeBC *nextEdge=firstEdgeBC;
		while(nextEdge!=NULL)
			nextEdge=nextEdge->PrintEdgeLoad();
		cout << endl;
	}
	
    //---------------------------------------------------
    // Periodic directions
	if(periodic.dof)
    {	PrintSection("PERIODIC DIRECTIONS");
		// x periodic, but y not (not allowed for axisymmetric)
		if(fmobj->periodic.dof==1)
		{	snprintf(hline,hsize,"x direction from %g to %g %s",periodic.xmin,periodic.xmax,UnitsController::Label(CULENGTH_UNITS));
			cout << hline << endl;
			if(periodic.fixDu)
			{	snprintf(hline,hsize,"   Displacement jump fixed at %g %s for exx = %g",periodic.du,
						UnitsController::Label(CULENGTH_UNITS),periodic.du/(periodic.xmax-periodic.xmin));
				cout << hline << endl;
			}
			if(periodic.fixDudy)
			{	snprintf(hline,hsize,"   Displacement jump slope fixed at %g",periodic.dudy);
				cout << hline << endl;
			}
		}
		
		// y periodic, but x not (only option for axisymmetrix - z periodic, but r not)
		else if(fmobj->periodic.dof==2)
		{	char pax = fmobj->IsAxisymmetric() ? 'z' : 'y' ;
            snprintf(hline,hsize,"%c direction from %g to %g %s",pax,periodic.ymin,periodic.ymax,UnitsController::Label(CULENGTH_UNITS));
			cout << hline << endl;
			if(periodic.fixDv)
			{	snprintf(hline,hsize,"   Displacement jump fixed at %g %s for e%c%c = %g",periodic.dv,
						UnitsController::Label(CULENGTH_UNITS),pax,pax,periodic.dv/(periodic.ymax-periodic.ymin));
				cout << hline << endl;
			}
			if(periodic.fixDvdx)
			{	snprintf(hline,hsize,"   Displacement jump slope fixed at %g",periodic.dvdx);
				cout << hline << endl;
			}
		}
		
		// x and y both periodic (not allowed for axisymmetric)
		else
		{	snprintf(hline,hsize,"x direction from %g to %g %s",periodic.xmin,periodic.xmax,UnitsController::Label(CULENGTH_UNITS));
			cout << hline << endl;
			if(periodic.fixDu)
			{	snprintf(hline,hsize,"   Displacement jump fixed at %g %s for exx = %g",periodic.du,
						UnitsController::Label(CULENGTH_UNITS),periodic.du/(periodic.xmax-periodic.xmin));
				cout << hline << endl;
			}
			if(periodic.fixDvdx)
			{	snprintf(hline,hsize,"   Displacement jump dv fixed at %g",periodic.dvdx);
				cout << hline << endl;
			}
            snprintf(hline,hsize,"y direction from %g to %g %s",periodic.ymin,periodic.ymax,UnitsController::Label(CULENGTH_UNITS));
			cout << hline << endl;
			if(periodic.fixDv)
			{	snprintf(hline,hsize,"   Displacement jump fixed at %g %s for eyy = %g",periodic.dv,
						UnitsController::Label(CULENGTH_UNITS),periodic.dv/(periodic.ymax-periodic.ymin));
				cout << hline << endl;
			}
			if(periodic.fixDudy)
			{	snprintf(hline,hsize,"   Displacement jump du fixed at %g",periodic.dudy);
				cout << hline << endl;
				
				// if both Dudy and Dvdx then global shear is specified
				if(periodic.fixDvdx)
				{	double b=periodic.dvdx/(periodic.xmax-periodic.xmin);
					double d=periodic.dudy/(periodic.ymax-periodic.ymin);
                    snprintf(hline,hsize,"    for gxy = %g and rotation = %g",b+d,(d-b)/2);
					cout << hline << endl;
				}
			}
		}
		cout << endl;
	}
}

#pragma mark NairnFEA: output accessors

// return name, caller should delete
const char *NairnFEA::CodeName(void) const { return "NairnFEA"; }

// FEA analysis type
void NairnFEA::GetAnalysisType(int np,char *fline)
{
	// analysis type
	switch(np)
	{	case PLANE_STRAIN:
			strcpy(fline,"2D Plane Strain Analysis");
			break;
			
		case PLANE_STRESS:
			strcpy(fline,"2D Plane Stress Analysis");
			break;
			
		case AXI_SYM:
			strcpy(fline,"Axisymmetric Analysis");
			break;
			
		default:
			throw CommonException("No FEA analysis type was provided.","CommonAnalysis::StartResultsOutput");
	}
}

// title for nodes and elements section
const char *NairnFEA::NodesAndElementsTitle(void) const { return "NODES AND ELEMENTS"; }

// Get file read handler
CommonReadHandler *NairnFEA::GetReadHandler(void)
{
	FEAReadHandler *handler=new FEAReadHandler();
	return handler;
}

// Explain usage of this program
void NairnFEA::Usage()
{
	CoutCodeVersion();
	cout << "    Expects Xerces 3.1.1 or newer" << endl;
	cout << "\nUsage:\n"
	"    NairnFEA [-options] <InputFile>\n\n"
	"This program reads the <InputFile> and does an FEA analysis.\n"
	"The results summary is directed to standard output.\n\n"
	"Options:\n"
	"    -a          Abort after setting up problem but before\n"
	"                   FEA Analysis\n"
	"    -H          Show this help\n"
	"    -np #       Set number of processors for parallel code\n"
	"    -v          Validate input file if DTD is provided in !DOCTYPE\n"
	"                   (default is to skip validation)\n"
	"    -w          Out results to current working directory\n"
	"                   (default isoutput to folder of input file)\n"
	"    -?          Show this help\n\n"
	"See http://osupdocs.forestry.oregonstate.edu/index.php/Main_Page for full documentation\n\n"
	<<  endl;
}


