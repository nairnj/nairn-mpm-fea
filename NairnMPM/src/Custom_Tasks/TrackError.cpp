/********************************************************************************
    TrackError.cpp
    OSParticulas
    
    Created by Chad Hammerquist on 7/29/2019

	This is a function to track the error on the particle basis over time.
	Currently only the tracking of stress_xx is implemented.
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/TrackError.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Exceptions/CommonException.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"
#include "Read_XML/Expression.hpp"
#include "System/CommonArchiveData.hpp"
#include "System/ArchiveData.hpp"
#include "Global_Quantities/GlobalQuantity.hpp"
#include <fstream>

// number in use
int TrackError::numTETasks = 0;

// Constructors
TrackError::TrackError() : CustomTask()
{
	numExprs = 0;						// number error expression to evaluate
	for(int i=0;i<MAX_COLUMNS;i++)
	{	lp_norm_p[i] = 2;						// default L2 norm
		Cumulative[i] = 0;						// average all time steps (!=0) each step (==0)
		total_error[i] = 0.;
	}
	customArchiveTime = -1.;			// default to use particle time step
	nextCustomArchiveTime = -1.;		// default to start in the beginning
	dummy = new char[1];
	count_n = 0;
	
	// in case more than on tasks
	TrackError::numTETasks++;
	TEnum = TrackError::numTETasks;
}

// Return name of this task
const char *TrackError::TaskName(void) { return "Track Error in Particle Values"; }

// Get input parameters
char *TrackError::InputParam(char *pName,int &input,double &gScaling)
{	
	if (strcmp(pName, "P") == 0)
	{	// order - default is 2 for L2 norm, 0 for RMS
		input = INT_NUM;
		return (char *)&(lp_norm_p[numExprs]);
	}
	else if (strcmp(pName, "cumulative") == 0)
	{	// ==0 to each step, !=0 for cumulative
		input = INT_NUM;
		return (char *)&(Cumulative[numExprs]);
	}
	else if (strcmp(pName, "archiveTime") == 0)
	{	// default to particle archive time, or customize here
		input = DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&customArchiveTime, gScaling, 1.e-3);
	}
	else if (strcmp(pName, "firstArchiveTime") == 0)
	{	// start at beginning or later (only if using custom archive time)
		input = DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&nextCustomArchiveTime, gScaling, 1.e-3);
	}
	else if (strcmp(pName, "function") == 0)
	{	// get an expression
		input = TEXT_PARAMETER;
		return dummy;
	}

	// not found, go to super class
	return CustomTask::InputParam(pName, input, gScaling);
}

// Set up and output
CustomTask *TrackError::Initialize(void)
{
	char fname[300];
	size_t fnsize=300;
	archiver->GetRelativeFilePathNum(fname, fnsize, "%s_Error%d.txt", TEnum);
	cout << "Archive particle error functions to tab-delimited file" << endl;
	cout << "   File: " << fname << endl;

	// time interval
	cout << "   Archive time: ";
	if (customArchiveTime >= 0.)
	{
		cout << customArchiveTime * UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS);
		if (nextCustomArchiveTime < 0.)
		{
			nextCustomArchiveTime = customArchiveTime;
			cout << endl;
		}
		else
		{
			cout << ", starting at " << nextCustomArchiveTime * UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS) << endl;
		}
	}
	else
		cout << "same as particle archives" << endl;


	// print function
	if(numExprs==0)
	{	throw CommonException("No error expressions were provided", "TrackError::Initialize");
	}
	
	// Create and open outfile
	archiver->GetFilePathNum(fname,fnsize,"%s%s_Error%d.txt",TEnum);
	ofstream afile;
	afile.open(fname, ios::out);
	
	cout << "   Error expressions:" << endl;
	for(int i=0;i<numExprs;i++)
	{	if(lp_norm_p[i]<0)
		{	throw CommonException("Invalid P value was used", "TrackError::Initialize");
		}
		
		// expresison number
		cout << "     " << (i+1) << ": ";
		afile << "# Expr_"<< (i+1) << ": ";
		
		// expression style
		switch(lp_norm_p[i])
		{	case 0:
				if(Cumulative[i]==0)
				{	cout << "sqrt[Sum_p ("<<error_expr[i]->GetString()<<")^2/Np]";
					afile << "sqrt[Sum_p ("<<error_expr[i]->GetString()<<")^2/Np]";
				}
				else
				{	cout << "(1/Ns) Sum_s {sqrt[Sum_p ("<<error_expr[i]->GetString()<<")^2/Np]}";
					afile << "(1/Ns) Sum_s {sqrt[Sum_p ("<<error_expr[i]->GetString()<<")^2/Np]}";
				}
				break;
			case 1:
				if(Cumulative[i]==0)
				{	cout << "(1/Np)[Sum_p ("<<error_expr[i]->GetString()<<")]";
					afile << "(1/Np)[Sum_p ("<<error_expr[i]->GetString()<<")]";
				}
				else
				{	cout << "(1/(Ns*Np))[Sum_s Sum_p ("<<error_expr[i]->GetString()<<")]";
					afile << "(1/(Ns*Np))[Sum_s Sum_p ("<<error_expr[i]->GetString()<<")]";
				}
				break;
			default:
				if(Cumulative[i]==0)
				{	cout << "(1/Np)[Sum_p ("<<error_expr[i]->GetString()<<")^"<<lp_norm_p[i]<<"]^(1/"<<lp_norm_p[i]<<")";
					afile << "(1/Np)[Sum_p ("<<error_expr[i]->GetString()<<")]^(1/"<<lp_norm_p[i]<<")";
				}
				else
				{	cout << "(1/(Ns*Np))[Sum_s Sum_p ("<<error_expr[i]->GetString()<<")^"<<lp_norm_p[i]<<"]^(1/"<<lp_norm_p[i]<<")";
					afile << "(1/(Ns*Np))[Sum_s Sum_p ("<<error_expr[i]->GetString()<<")^"<<lp_norm_p[i]<<"]^(1/"<<lp_norm_p[i]<<")";
				}
				break;
		}
		cout << endl;
		afile << endl;
	}
	
	cout << "   where Np = (# particles), Ns = (# steps evaluated)" << endl;
	afile << "# First column = time, Np = (# particles), Ns = (# steps evaluated)" << endl;

	// File Header
	afile << "#setColor";
	for(int i=0;i<numExprs;i++) afile << "\t" << GlobalQuantity::PickColor(i);
	afile << endl;
	afile << "#setName";
	for(int i=0;i<numExprs;i++) afile << "\tExpr_" << (i+1);
	afile << endl;
	afile.close();
	
	return nextTask;
}

// See if need to write to file
// If not writting to file, then only need to do cummulative calculations
CustomTask *TrackError::PrepareForStep(bool &)
{
	if (customArchiveTime >= 0.)
	{
		if (mtime + timestep >= nextCustomArchiveTime)
		{
			doErrorExport = true;
			nextCustomArchiveTime += customArchiveTime;
		}
		else
			doErrorExport = false;
	}
	else
		doErrorExport = archiver->WillArchive(); // archive when particles are archived

	return nextTask;
}

// Evaulate expressions and write to file
CustomTask * TrackError::StepCalculation(void)
{
	// Will it write to a file on this time step
	if(!doErrorExport) return nextTask;
	
	// variables
	unordered_map<string,double> vars;
	double error[10];
	for(int i=0;i<numExprs;i++) error[i] = 0.;

	// loop over non-rigid partiles
	// (I think use of hash map means loop cannot be parallel?)
	int Nfound=0;
	for(int p = 0; p < nmpmsNR; p++)
	{	// the material point
		MPMBase *point = mpm[p];
		if(point->InReservoir()) continue;
		Nfound++;
		
		// current time at end of this time step
		vars["t"] = (mtime+timestep) * UnitsController::Scaling(1.e3);

		// particle positions
		vars["x"] = point->pos.x;
		vars["y"] = point->pos.y;
		vars["z"] = point->pos.z;
		vars["xo"] = point->origpos.x;
		vars["yo"] = point->origpos.y;
		vars["zo"] = point->origpos.z;

		// particle velocity
		vars["vx"] = point->vel.x;
		vars["vy"] = point->vel.y;
		vars["vz"] = point->vel.z;

		// get stress
		double rho0 = theMaterials[point->MatID()]->GetRho(point);
		double rho = rho0 / theMaterials[point->MatID()]->GetCurrentRelativeVolume(point, 0);
		double wt =  rho*UnitsController::Scaling(1.e-6);			// Legacy convert to MPa
		Tensor sp = point->ReadStressTensor();
		vars["sxx"] = wt*sp.xx;
		vars["syy"] = wt*sp.yy;
		vars["szz"] = wt*sp.zz;
		vars["sxy"] = wt*sp.xy;
		vars["sxz"] = wt*sp.xz;
		vars["syz"] = wt*sp.yz;
		
		// Evaluate each error
		for(int i=0;i<numExprs;i++)
		{	double exprVal = error_expr[i]->EvaluateFunction(vars);
			switch(lp_norm_p[i])
			{	case 1:
					error[i] += exprVal;
					break;
				case 0:
					error[i] += exprVal*exprVal;
					break;
				default:
					error[i] += pow(exprVal,(double)lp_norm_p[i]);
					break;
			}
		}
	}

	// convert to final results and/or add to cumulative results
	double Np = (double)Nfound;
	double prevNs = (double)count_n;
	double nextNs = prevNs+1.;
	double rescaleNs = prevNs/nextNs;
	for(int i=0;i<numExprs;i++)
	{	switch(lp_norm_p[i])
		{	case 1:
				// ordinary mean
				if(Cumulative[i]==0)
					total_error[i] = error[i]/Np;
				else
					total_error[i] = rescaleNs*total_error[i] + error[i]/(Np*nextNs);
				break;
			case 0:
				// RMS
				if(Cumulative[i]==0)
					total_error[i] = sqrt(error[i]/Np);
				else
					total_error[i] = rescaleNs*total_error[i] + (sqrt(error[i]/Np))/nextNs;
				break;
			default:
			{	// p norm
				double lpower = 1./(double)lp_norm_p[i];
				if(Cumulative[i]==0)
					total_error[i] = pow(error[i],lpower)/Np;
				else
				{	double lpsum = pow(prevNs*Np*total_error[i],(double)lp_norm_p[i]) + error[i];
					total_error[i] = pow(lpsum,lpower)/(nextNs*Np);
				}
				break;
			}
		}
	}
	// for cumulative terms
	count_n++;
				
	// export to file
	
	// open error file
	char fname[300];
	size_t fnsize=300;
	archiver->GetFilePathNum(fname,fnsize,"%s%s_Error%d.txt",TEnum);
	ofstream afile;
	afile.open(fname, ios::out | ios::app);
	if (!afile.is_open())
		archiver->FileError("Cannot open a Error archive file", fname, "TrackError::StepCalculation");
	
	// Write a row of results
	afile << mtime * UnitsController::Scaling(1.e3);
	for(int i=0;i<numExprs;i++)
		afile << "\t" << total_error[i];
	afile << endl;
	afile.close();

	return nextTask;
}

// set the next function
void TrackError::SetTextParameter(char *fxn, char *ptr)
{
	if(numExprs>=MAX_COLUMNS)
		ThrowSAXException("Too many TrackError functions used");
	
	// Load_Function
	if(ptr == dummy)
	{	if(fxn == NULL)
			ThrowSAXException("TrackError function is missing");
		if(strlen(fxn) == 0)
			ThrowSAXException("TrackError function is empty");

		try
		{	error_expr[numExprs] = Expression::CreateExpression(fxn, "TrackError function is not valid");
			numExprs++;
			
			// if room, make same as this
			if(numExprs<MAX_COLUMNS)
			{	lp_norm_p[numExprs] = lp_norm_p[numExprs-1];
				Cumulative[numExprs] = Cumulative[numExprs-1];
			}
		}
		catch(...)
		{	ThrowSAXException("TrackError function is not valid");
		}
	}
	else
		CustomTask::SetTextParameter(fxn, ptr);
}
