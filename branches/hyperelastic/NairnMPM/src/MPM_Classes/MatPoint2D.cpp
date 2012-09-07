/********************************************************************************
    MatPoint2D.cpp
    NairnMPM
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "MPM_Classes/MatPoint2D.hpp"
#include "Materials/MaterialBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"

#pragma mark MatPoint2D::Constructors and Destructors

// Constructors
MatPoint2D::MatPoint2D() {}

// Constructors
MatPoint2D::MatPoint2D(int inElemNum,int theMatl,double angin,double thickin) : MPMBase(inElemNum,theMatl,angin)
{	thick=thickin;
}

#pragma mark MatPoint2D:Calculations and Incrementers

// Update Strains for this particle
// Velocities for all fields on present on the nodes
void MatPoint2D::UpdateStrain(double strainTime,int secondPass,int np)
{
	// this particle's material
	MaterialBase *matRef=theMaterials[MatID()];
	
	// exit if rigid
	if(matRef->Rigid()) return;
	
	// make sure have mechanical properties for this material and angle
	matRef->LoadMechanicalProps(this,np);
	
	// get field number
	int matfld=matRef->GetField();
	
	int i,numnds,nds[MaxShapeNds];
    double fn[MaxShapeNds],xDeriv[MaxShapeNds],yDeriv[MaxShapeNds];
	Vector vel;
    double dvxx,dvyy,dvxy,dvyx;
    
	// find shape functions and derviatives
	int iel=ElemID();
	theElements[iel]->GetShapeGradients(&numnds,fn,nds,&ncpos,xDeriv,yDeriv,NULL);
    
    // Find strain rates at particle from current grid velocities
	//   and using the velocity field for that particle and each node and the right material
    dvxx=dvyy=dvxy=dvyx=0.;
    for(i=1;i<=numnds;i++)
	{	vel=nd[nds[i]]->GetVelocity((short)vfld[i],matfld);
        dvxx+=vel.x*xDeriv[i];
        dvyy+=vel.y*yDeriv[i];
        dvxy+=vel.x*yDeriv[i];
        dvyx+=vel.y*xDeriv[i];
    }
	    
    // save velocity gradient (if needed for J integral calculation)
    SetVelocityGradient(dvxx,dvyy,dvxy,dvyx,secondPass);
    
    // convert to strain increments (e.g., now dvxx = dvx/dx * dt = d/dx(du/dt) * dt = d/dt(du/dx) * dt = du/dx)
    dvxx*=strainTime;
    dvyy*=strainTime;
    dvxy*=strainTime;
    dvyx*=strainTime;
    
	// find effective particle transport property - find it from the grid results
	if(transportTasks)
	{	// initialize
		TransportTask *nextTransport=transportTasks;
		while(nextTransport!=NULL)
			nextTransport=nextTransport->ZeroValueExtrap();
			
		// loop over nodes
		for(i=1;i<=numnds;i++)
		{	nextTransport=transportTasks;
			while(nextTransport!=NULL)
				nextTransport=nextTransport->IncrementValueExtrap(nd[nds[i]],fn[i]);
		}
		
		// find change
		nextTransport=transportTasks;
		while(nextTransport!=NULL)
			nextTransport=nextTransport->GetDeltaValue(this);
	}
	
	// If thermal ramp, but no conduction, load temperature change into dTemperature variable
	if(!ConductionTask::active)
	{	ConductionTask::dTemperature=pTemperature-pPreviousTemperature;
		pPreviousTemperature=pTemperature;
	}
	
    // update particle strain and stress using its constituitive law
    matRef->MPMConstLaw(this,dvxx,dvyy,dvxy,dvyx,strainTime,np);
}

// Move position (2D) (in mm)
// external work units g-mm^2/sec^2 (* 10^-9 to get J)
void MatPoint2D::MovePosition(double delTime,Vector *dv)
{	double dx=delTime*dv->x;
	double dy=delTime*dv->y;
	pos.x+=dx;
    pos.y+=dy;
    extWork+=dx*pFext.x+dy*pFext.y;
}

// Move velocity (2D) (in mm/sec) possibly with damping
void MatPoint2D::MoveVelocity(double delTime,double damping,Vector *vstar)
{	vel.x+=delTime*(acc.x-damping*vstar->x);
    vel.y+=delTime*(acc.y-damping*vstar->y);
}

// Scale velocity (2D) (in mm/sec) optionally with rigid materials
void MatPoint2D::SetVelocitySpeed(double speed)
{	double norm=sqrt(vel.x*vel.x+vel.y*vel.y);
	if(DbleEqual(norm,0.)) return;
	vel.x*=speed/norm;
    vel.y*=speed/norm;
}

#pragma mark MatPoint2D::Accessors

// set original position (2D) (in mm)
void MatPoint2D::SetOrigin(Vector *pt)
{	origpos.x=pt->x;
    origpos.y=pt->y;
	origpos.z=0.;
}

// set position (2D) (in mm)
void MatPoint2D::SetPosition(Vector *pt)
{	pos.x=pt->x;
    pos.y=pt->y;
	pos.z=0.;
}

// set position (2D) (in mm/sec)
void MatPoint2D::SetVelocity(Vector *pt)
{	vel.x=pt->x;
    vel.y=pt->y;
	vel.z=0.;
}

// thickness (in mm)
double MatPoint2D::thickness() { return thick; }

// set current volume using current strains - but only 2D strains
void MatPoint2D::SetDilatedVolume(void)
{	double rho=theMaterials[MatID()]->rho*0.001;			// in g/mm^3
	double dilate=(1.+ep.xx)*(1.+ep.yy);
	volume=dilate*mp/rho;									// in mm^3
}

// return internal force as -mp sigma.deriv * 1000. which converts to g mm/sec^2 or micro N
void MatPoint2D::Fint(Vector &fout,double xDeriv,double yDeriv,double zDeriv)
{	fout.x=-mp*(sp.xx*xDeriv+sp.xy*yDeriv)*1000.;
	fout.y=-mp*(sp.xy*xDeriv+sp.yy*yDeriv)*1000.;
	fout.z=0.;
}

// return external force (times a shape function)
void MatPoint2D::Fext(Vector &fout,double fni)
{	fout.x=fni*pFext.x;
	fout.y=fni*pFext.y;
	fout.z=0.;
}

// zero the temperature gradient
void MatPoint2D::AddTemperatureGradient(void)
{	pTemp->DT.x=0.;
    pTemp->DT.y=0.;
	pTemp->DT.z=0.;
}

// add to the temperature gradient
void MatPoint2D::AddTemperatureGradient(Vector *grad)
{	pTemp->DT.x+=grad->x;
    pTemp->DT.y+=grad->y;
}

// return conduction force = - V [D] Grad T . Grad S
double MatPoint2D::FCond(double dshdx,double dshdy,double dshdz)
{
	Tensor *kten=theMaterials[MatID()]->GetkCondTensor();
	return -volume*((kten->xx*pTemp->DT.x + kten->xy*pTemp->DT.y)*dshdx
						+ (kten->xy*pTemp->DT.x + kten->yy*pTemp->DT.y)*dshdy);
}

// zero the concentration gradient
void MatPoint2D::AddConcentrationGradient(void)
{	pDiffusion->Dc.x=0.;
    pDiffusion->Dc.y=0.;
	pDiffusion->Dc.z=0.;
}

// add to the concentration gradient
void MatPoint2D::AddConcentrationGradient(Vector *grad)
{	pDiffusion->Dc.x+=grad->x;
    pDiffusion->Dc.y+=grad->y;
}

// return diffusion force = - V [D] Grad C . Grad S
double MatPoint2D::FDiff(double dshdx,double dshdy,double dshdz)
{
	Tensor *Dten=theMaterials[MatID()]->GetDiffusionTensor();
	return -volume*((Dten->xx*pDiffusion->Dc.x + Dten->xy*pDiffusion->Dc.y)*dshdx
						+ (Dten->xy*pDiffusion->Dc.x + Dten->yy*pDiffusion->Dc.y)*dshdy);
}

// return kinetic energy
double MatPoint2D::KineticEnergy(void)
{	return 0.5*mp*(vel.x*vel.x+vel.y*vel.y);
}
	

