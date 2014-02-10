/*********************************************************************
    BeginResults.cpp
    Nairn Research Group FEA Code
    
    Created by jnairn on Mon Feb 2 2003.
    Copyright (c) 2003, All rights reserved.    
*********************************************************************/

#include "NairnFEA_Class/NairnFEA.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/NodalDispBC.hpp"
#include "Boundary_Conditions/NodalLoad.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Boundary_Conditions/EdgeBC.hpp"
#include "Read_XML/mathexpr.hpp"

/*********************************************************************
    Begin results information
*********************************************************************/

// print title of results file
void NairnFEA::PrintAnalysisTitle(void)
{
	cout << "FEA ANALYSIS BY ";
	CoutCodeVersion();
}

// print analysis type
void NairnFEA::PrintAnalysisType(void) {}
    
// finish start of FEA results file
void NairnFEA::MyStartResultsOutput(void)
{
    char hline[200];
    int i;
    
    //---------------------------------------------------
    // Temperature
    PrintSection("THERMAL LOAD");
	if(temperatureExpr!=NULL)
	{	if(!CreateFunction(temperatureExpr))
			throw CommonException("The temperature expression is not a valid function","NairnFEA::MyStartResultsOutput");
		for(i=1;i<=nnodes;i++)
		{	nd[i]->gTemperature = FunctionValue(1,nd[i]->x,nd[i]->y,0.,0.,0.,0.)-stressFreeTemperature;
		}
		DeleteFunction();
	}
	else
	{	// unknown, but may have been set in explicit node commands
		temperatureExpr=new char[2];
		strcpy(temperatureExpr,"?");
	}
	
	sprintf(hline,"T0: %.2lf C\nT: %s C",stressFreeTemperature,temperatureExpr);
	cout << hline << endl;
	if(stressFreeTemperature!=0)
	{	sprintf(hline,"T-T0: %s - %.2lf C",temperatureExpr,stressFreeTemperature);
		cout << hline << endl;
	}
	cout << endl;
    
    //---------------------------------------------------
    // Fixed Displacements
	if(firstDispBC!=NULL)
	{	PrintSection("NODAL POINTS WITH FIXED DISPLACEMENTS");
		cout << " Node  DOF  Displacement (mm)  Axis        Angle\n"
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
		cout << " Node  DOF       Load (N)\n"
			<< "------------------------------\n";
		NodalLoad *nextLoad=firstLoadBC;
		while(nextLoad!=NULL)
			nextLoad=nextLoad->PrintLoad();
		cout << endl;
	}

    //---------------------------------------------------
    // Stress element faces
	if(firstEdgeBC!=NULL)
	{	PrintSection("FACES WITH APPLIED STRESS (MPa)");
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
		{	sprintf(hline,"x direction from %g to %g mm",periodic.xmin,periodic.xmax);
			cout << hline << endl;
			if(periodic.fixDu)
			{	sprintf(hline,"   Displacement jump fixed at %g mm for exx = %g",periodic.du,periodic.du/(periodic.xmax-periodic.xmin));
				cout << hline << endl;
			}
			if(periodic.fixDudy)
			{	sprintf(hline,"   Displacement jump slope fixed at %g",periodic.dudy);
				cout << hline << endl;
			}
		}
		
		// y periodic, but x not (only option for axisymmetrix - z periodic, but r not)
		else if(fmobj->periodic.dof==2)
		{	char pax = fmobj->IsAxisymmetric() ? 'z' : 'y' ;
			sprintf(hline,"%c direction from %g to %g mm",pax,periodic.ymin,periodic.ymax);
			cout << hline << endl;
			if(periodic.fixDv)
			{	sprintf(hline,"   Displacement jump fixed at %g mm for e%c%c = %g",periodic.dv,pax,pax,periodic.dv/(periodic.ymax-periodic.ymin));
				cout << hline << endl;
			}
			if(periodic.fixDvdx)
			{	sprintf(hline,"   Displacement jump slope fixed at %g",periodic.dvdx);
				cout << hline << endl;
			}
		}
		
		// x and y both periodic (not allowed for axisymmetric)
		else
		{	sprintf(hline,"x direction from %g to %g mm",periodic.xmin,periodic.xmax);
			cout << hline << endl;
			if(periodic.fixDu)
			{	sprintf(hline,"   Displacement jump fixed at %g mm for exx = %g",periodic.du,periodic.du/(periodic.xmax-periodic.xmin));
				cout << hline << endl;
			}
			if(periodic.fixDvdx)
			{	sprintf(hline,"   Displacement jump dv fixed at %g",periodic.dvdx);
				cout << hline << endl;
			}
			sprintf(hline,"y direction from %g to %g mm",periodic.ymin,periodic.ymax);
			cout << hline << endl;
			if(periodic.fixDv)
			{	sprintf(hline,"   Displacement jump fixed at %g mm for eyy = %g",periodic.dv,periodic.dv/(periodic.ymax-periodic.ymin));
				cout << hline << endl;
			}
			if(periodic.fixDudy)
			{	sprintf(hline,"   Displacement jump du fixed at %g",periodic.dudy);
				cout << hline << endl;
				
				// if both Dudy and Dvdx then global shear is specified
				if(periodic.fixDvdx)
				{	double b=periodic.dvdx/(periodic.xmax-periodic.xmin);
					double d=periodic.dudy/(periodic.ymax-periodic.ymin);
					sprintf(hline,"    for gxy = %g and rotation = %g",b+d,(d-b)/2);
					cout << hline << endl;
				}
			}
		}
		cout << endl;
	}
}
