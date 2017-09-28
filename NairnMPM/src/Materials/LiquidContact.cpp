/********************************************************************************
	LiquidContact.cpp
	nairn-mpm-fea
 
	"Slippery" stick where sliding traction is slightly less then
	stick traction depending on viscosity of the liquid. It is intended
	to set shear stress proportation to velocity gradient on the surface.
 
	Created by John Nairn, Feb 9, 2017.
		Copyright (c) 2015 John A. Nairn, All rights reserved.
 ********************************************************************************/

#include "stdafx.h"
#include "Materials/LiquidContact.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "System/UnitsController.hpp"

#pragma mark LiquidContact::Constructors and Destructors

// Constructors
LiquidContact::LiquidContact() {}

LiquidContact::LiquidContact(char *matName) : CoulombFriction(matName)
{
	liquidPhaseID = -1;
	liquidPhase = NULL;
	
	// law is coeff*viscosity*delta V
	// default is 2 mm^-1 or zone of 500 microns
	frictionCoeff = 2.0;
}

#pragma mark LiquidContact::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *LiquidContact::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	// look for solid material
	if(strcmp(xName,"LiquidPhase")==0)
	{	input = PHASE_NUM;
		return (char *)&liquidPhaseID;
	}

	else if(strcmp(xName,"coeff")==0)
	{	input=DOUBLE_NUM;
		return (char *)&frictionCoeff;				// Legacy units 1/mm
	}
	
	// does not all any from MaterialBase
	return (char *)NULL;
}

// Verify input properties do calculations; if problem return string with an error message
// Don't pass on to material base
const char *LiquidContact::VerifyAndLoadProperties(int np)
{
	// check liquid material
	if(liquidPhaseID<1 || liquidPhaseID>nmat)
		return "The liquid phase material was not specied or is invalid";
	liquidPhase = theMaterials[liquidPhaseID-1];
	if(liquidPhase->GetViscosity(-1.)<0.)
		return "The liquid phase material does not have a valid viscosity";
	
	// must be positive
	if(frictionCoeff<0.)
		return "The contact law inverse length must be positive";
	
	// set to be sure parent class friction handled correctly
	frictionStyle = FRICTIONAL;
	
	// no super class tasks allowed
	return NULL;
}

// print contact law details to output window
void LiquidContact::PrintContactLaw(void) const
{
	// stiffness
	cout << "Contact shear stress = k * (viscosity of liquid) * delta(v) at wall" << endl;
	cout << "     k: " << frictionCoeff << " 1/" << UnitsController::Label(CULENGTH_UNITS);
	
	// liquid phase
	cout << "  liquid: ";
	if(liquidPhase!=NULL)
		cout << liquidPhase->name << " (" << liquidPhaseID << ")" << endl;
	else
		cout << "(invalid ID " << liquidPhaseID << ")" << endl;
}

#pragma mark LiquidContact:Step Methods

//#define DEBUG_SECANT

// Return SslideAcDt = Swall * SStickAcDt/(Swall + mred/(viscosity*Ac*deltime))
// Swall is stored in contact variable frictionCoeff (units 1/(length units))
double LiquidContact::GetSslideAcDt(double NAcDt,double SStickAcDt,double mred,
								 double contactArea,bool &inContact,double deltime) const
{
	double x1,y1,x2,y2,xk,yk;
	
	// constants
	double factor = frictionCoeff*contactArea*deltime/mred;			// Sw Ac dt/mired
	double gmaxdot = frictionCoeff*SStickAcDt/mred;					// max final shear rate (for slip)
	
	// ask liquid for solution or brackets
	double viscosity = liquidPhase->BracketContactLawShearRate(gmaxdot,factor,x1,y1,x2,y2);
	
	// secant method if not solved
	if(viscosity<0.)
	{	// numerical solution when liquid has shear rate dependent viscosity
		// false position method - see code example in https://en.wikipedia.org/wiki/False_position_method
		// Using Anderson-Bjork modification (scaling factor below)
		double e = 1.e-5, scale, xkm1=x1;
		int n, m = 20;
		int side = 0;

#ifdef DEBUG_SECANT
		cout << "# Bracket (log(g1(dot),y1,eta1) to (log(g2(dot),y2,eta2) for gmax(dot): (" << log10(gmaxdot*x1) << "," << y1 << "," << liquidPhase->GetViscosity(x1*gmaxdot)
		<< ") to (" << log10(gmaxdot*x2) << "," << y2 << "," << liquidPhase->GetViscosity(x2*gmaxdot)<< ") for " << log10(gmaxdot) << endl;
#endif
		
		// iterate until converged
		for(n=0;n<m;n++)
		{	// get next value
			xk = (x1*y2-x2*y1)/(y2-y1);
			
			// new value
			viscosity = liquidPhase->GetViscosity(xk*gmaxdot);
			yk = xk*(1+factor*viscosity)-1.;
			
#ifdef DEBUG_SECANT
			cout << "#    n=" << n << " (x1,xk,x2)=(" << x1 << "," << xk << "," << x2 << ")"
					<< " (y1,yk,y2)=(" << y1 << "," << yk << "," << y2 << ")" << endl;
#endif
			
			if(yk*y2 > 0.)
			{	// yk and y2 have same sign, copy k to x2 and retain x1
				scale = 1.-yk/y2;
				x2 = xk;
				y2 = yk;
				if(side==-1) y1 = scale > 0. ? y1*scale : 0.5*y1;
				side = -1;
			}
			else if(y1*yk > 0.)
			{	// y1 and yk have same sign, copy k to x1 and retain x2
				scale = 1.- yk/y1;
				x1 = xk;
				y1 = yk;
				if(side==1) y2 = scale > 0. ? y2*scale : 0.5*y2;
				side = 1;
			}
			else
			{	// essentially zero
				break;
			}
			
			if(fabs(xk-xkm1)<e) break;
			xkm1 = xk;
		}
#ifdef DEBUG_SECANT
		cout << "#    done eta = " << viscosity << " at log10(rate) " << log10(xk*gmaxdot) << endl;
#endif
	}
	
	// the final shear rate has the solved viscosity
	return frictionCoeff*SStickAcDt/(frictionCoeff + mred/(viscosity*contactArea*deltime));
}

#pragma mark ContactLaw::Accessors

// return unique, short name for this material
const char *LiquidContact::MaterialType(void) const { return "Liquid/Wall Shear-Rate Dependent Contact"; }

// Give details aobut frictional contact
bool LiquidContact::ContactLawNeedsContactArea(void) const { return true; }
