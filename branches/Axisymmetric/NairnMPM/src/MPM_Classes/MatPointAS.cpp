/********************************************************************************
	MatPointAS.hpp
	NairnMPM

	Created by John Nairn on 10/23/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
 
	When MPM is axisymmetric, the x axis becomes radial (or R) direction, the
		y axis because axial direction (or Z direction), and the z direction
		become the tangential direction (or theta direction)
 
	The thickness (in thick) is set to the initial radial position of the
		material point
********************************************************************************/

#include "MPM_Classes/MatPointAS.hpp"
#include "Materials/MaterialBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Exceptions/MPMTermination.hpp"

#pragma mark MatPointAS::Constructors and Destructors

// Constructors
MatPointAS::MatPointAS() {}

// Constructors
MatPointAS::MatPointAS(int inElemNum,int theMatl,double angin,double thickin) : MatPoint2D(inElemNum,theMatl,angin,thickin)
{	
}

#pragma mark MatPointAS:Calculations and Incrementers

// Update Strains for this particle
// Velocities for all fields on present on the nodes
void MatPointAS::UpdateStrain(double strainTime,int secondPass,int np)
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
    double fn[MaxShapeNds],xDeriv[MaxShapeNds],yDeriv[MaxShapeNds],zDeriv[MaxShapeNds];
	Vector vel;
    double dvrr,dvzz,dvrz,dvzr,dvtt;
    
	// find shape functions and derviatives
	int iel=ElemID();
	theElements[iel]->GetShapeGradients(&numnds,fn,nds,&ncpos,xDeriv,yDeriv,zDeriv,this);
    
    // Find strain rates at particle from current grid velocities
	//   and using the velocity field for that particle and each node and the right material
    // In axisymmetric x->r, y->z, and z->hoop
    dvrr=dvzz=dvrz=dvzr=dvtt=0.;
    for(i=1;i<=numnds;i++)
	{	vel=nd[nds[i]]->GetVelocity((short)vfld[i],matfld);
        dvrr += vel.x*xDeriv[i];
        dvzz += vel.y*yDeriv[i];
        dvrz += vel.x*yDeriv[i];
        dvzr += vel.y*xDeriv[i];
        dvtt += vel.x*zDeriv[i];
    }
    
    // save velocity gradient (if needed for J integral calculation)
    SetVelocityGradient(dvrr,dvzz,dvrz,dvzr,secondPass);
    
    // convert to strain increments
    // e.g., now dvrr = dvr/dr * dt = d/dr(du/dt) * dt = d/dt(du/dr) * dt = du/dr)
    dvrr *= strainTime;
    dvzz *= strainTime;
    dvrz *= strainTime;
    dvzr *= strainTime;
    dvtt *= strainTime;
    
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
    matRef->MPMConstLaw(this,dvrr,dvzz,dvrz,dvzr,dvtt,strainTime,np);
}

#pragma mark MatPoint2D::Accessors

// set original position (2D) (in mm) and here set thick to be radial position
void MatPointAS::SetOrigin(Vector *pt)
{	origpos.x=pt->x;
    origpos.y=pt->y;
	origpos.z=0.;
	thick = pt->x;
}

// return internal force as -mp sigma.deriv * 1000. which converts to g mm/sec^2 or micro N
void MatPointAS::Fint(Vector &fout,double xDeriv,double yDeriv,double zDeriv)
{	fout.x=-mp*(sp.xx*xDeriv+sp.xy*yDeriv+sp.zz*zDeriv)*1000.;
	fout.y=-mp*(sp.xy*xDeriv+sp.yy*yDeriv)*1000.;
	fout.z=0.;
}

// To support traction boundary conditions, find the deformed edge, natural coordinates of
// the corners along the edge, elements for those edges, and a normal vector in direction
// of the traction
// return ratio of second nodal weight to first one
double MatPointAS::GetTractionInfo(int face,int dof,int *cElem,Vector *corners,Vector *tscaled,int *numDnds)
{
    *numDnds = 2;
    double faceWt,ratio=1.;
    
    // always UNIFORM_GIMP_AS
	
	// initial vectors only
	double r1x = mpmgrid.gridx*0.25;
	double r2y = mpmgrid.gridy*0.25;
	double rp = pos.x;
        
	Vector c1,c2;
	switch(face)
	{	case 1:
			// lower edge
			c1.x = pos.x-r1x;
			c2.x = pos.x+r1x;
			c1.y = c2.y = pos.y-r2y;
			if(c1.x < 0)
			{	c1.x = 0.;
				r1x = 0.5*c2.x;
				rp = r1x;
			}
			faceWt = r1x*(rp + r1x/3.);				// node 3, node 4 should be minus
			ratio = r1x*(rp - r1x/3.)/faceWt;		// find the ratio
			break;
			
		case 2:
			// right edge
			c1.x = c2.x = pos.x+r1x;
			c1.y = pos.y-r2y;
			c2.y = pos.y+r2y;
			faceWt = r2y*c1.x;
			break;
			
		case 3:
			// top edge
			c1.x = pos.x+r1x;
			c2.x = pos.x-r1x;
			c1.y = c2.y = pos.y+r2y;
			if(c2.x<0.)
			{	c2.x = 0.;
				r1x = 0.5*c1.x;
				rp = r1x;
			}
			faceWt = r1x*(rp - r1x/3.);				// node 1, node 2 should be plus
			ratio = r1x*(rp + r1x/3.)/faceWt;		// find the ratio
			break;
			
		default:
			// left edge
			c1.x = c2.x = pos.x-r1x;
			c1.y = pos.y+r2y;
			c2.y = pos.y-r2y;
			if(c1.x<0.)
			{	c1.x = c2.x= 0.;
				faceWt = 0.;
			}
			else
				faceWt = r2y*c1.x;
			break;
	}
	
	// get elements
	try
	{	cElem[0] = mpmgrid.FindElementFromPoint(&c1)-1;
		theElements[cElem[0]]->GetXiPos(&c1,&corners[0]);
		
		cElem[1] = mpmgrid.FindElementFromPoint(&c2)-1;
		theElements[cElem[1]]->GetXiPos(&c2,&corners[1]);
	}
	catch(...)
	{	throw MPMTermination("A Traction edge node has left the grid.","MatPointAS::GetTractionInfo");
	}
	
    // get traction normal vector by radial integral for first node no the edge
    ZeroVector(tscaled);
	switch(dof)
	{	case 1:
			// normal is x direction
			tscaled->x = faceWt;
			break;
        case 2:
            // normal is y direction
            tscaled->y = faceWt;
		default:
			// normal is z direction (not used here)
            tscaled->z = faceWt;
			break;
	}
	
	return ratio;
}



