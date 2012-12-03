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
#include "Boundary_Conditions/BoundaryCondition.hpp"

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

// get dilated current volume using current deformation gradient
// only used for contact (cracks and multimaterial) and for transport tasks
// when actual is false, get t*Ap, where Ap is deformed particle area
double MatPointAS::GetVolume(bool volumeType)
{	double rho=theMaterials[MatID()]->rho*0.001;			// in g/mm^3
	if(volumeType==DEFORMED_VOLUME)
		return GetRelativeVolume()*mp/rho;					// in mm^3
	
	// get deformed area per unit radial position or per radian
	double pF[3][3];
	GetDeformationGradient(pF);
	
	// note that mp/rho = rho Ap rp0 / rho = Ap rp0
	double areaPerRadian = (pF[0][0]*pF[1][1]-pF[1][0]*pF[0][1])*mp/rho;
	return volumeType==DEFORMED_AREA ? areaPerRadian/origpos.x : areaPerRadian ;
}

// get unscaled volume for use only in contact and imperfect interface calculations
// For axisymmetric, really needs area of the initial particle domain in mm^3
// Calculations will need to multiply by radial position to get local volume
// Here thickness is the original radial position of the particle
double MatPointAS::GetUnscaledVolume(void)
{	double rho=theMaterials[MatID()]->rho*0.001;			// in g/mm^3
	return mp/(rho*thickness());                            // in mm^3 per unit radial position
}


// To support CPDI find nodes in the particle domain, find their elements,
// their natural coordinates, and weighting values for gradient calculations
// Should be done only once per time step and here only for axisynmmetric
// and linear CPDI
void MatPointAS::GetCPDINodesAndWeights(int cpdiType)
{
	// get particle 2D deformation gradient
	double pF[3][3];
	GetDeformationGradient(pF);
	
	// always LINEAR_CPDI_AS
	
	// get polygon vectors - these are from particle to edge
    //      and generalize semi width lp in 1D GIMP
	Vector r1,r2,c;
	r1.x = pF[0][0]*mpmgrid.gridx*0.25;
	r1.y = pF[1][0]*mpmgrid.gridx*0.25;
	r2.x = pF[0][1]*mpmgrid.gridy*0.25;
	r2.y = pF[1][1]*mpmgrid.gridy*0.25;

#ifdef TRUNCATE
	// shrink domain if any have x < 0, but keep particle in the middle
	if(pos.x-fabs(r1.x+r2.x)<0.)
	{	// make pos.x-shrink*fabs(r1.x+r2.x) very small and positive
		double shrink = (pos.x-mpmgrid.gridx*1.e-10)/fabs(r1.x+r2.x);
		r1.x *= shrink;
		r1.y *= shrink;
		r2.x *= shrink;
		r2.y *= shrink;
	}
#endif
    
    // Particle domain area is area of the full parallelogram
    // Assume positive due to orientation of initial vectors, and sign probably does not matter
    double Ap = 4.*(r1.x*r2.y - r1.y*r2.x);
    
	try
	{	// nodes at four courves in ccw direction
		c.x = pos.x-r1.x-r2.x;
		c.y = pos.y-r1.y-r2.y;
		cpdi[0]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[0]->inElem]->GetXiPos(&c,&cpdi[0]->ncpos);
		
		c.x = pos.x+r1.x-r2.x;
		c.y = pos.y+r1.y-r2.y;
		cpdi[1]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[1]->inElem]->GetXiPos(&c,&cpdi[1]->ncpos);
		
		c.x = pos.x+r1.x+r2.x;
		c.y = pos.y+r1.y+r2.y;
		cpdi[2]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[2]->inElem]->GetXiPos(&c,&cpdi[2]->ncpos);
		
		c.x = pos.x-r1.x+r2.x;
		c.y = pos.y-r1.y+r2.y;
		cpdi[3]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[3]->inElem]->GetXiPos(&c,&cpdi[3]->ncpos);
		
		// shape weights
		double rxsum = r1.x+r2.x;
		double rxdiff = r1.x-r2.x;
		cpdi[0]->ws = 0.25 - rxsum/(12.*pos.x);
		cpdi[1]->ws = 0.25 + rxdiff/(12.*pos.x);
		cpdi[2]->ws = 0.25 + rxsum/(12.*pos.x);
		cpdi[3]->ws = 0.25 - rxdiff/(12.*pos.x);
		
		// gradient weighting values
		Ap = 1./Ap;
		double beta = r1.x*r1.y - r2.x*r2.y;
		double rysum = r1.y+r2.y;
		double rydiff = r1.y-r2.y;
		cpdi[0]->wg.x = (rydiff - beta/(3.*pos.x))*Ap;
		cpdi[0]->wg.y = -rxdiff*Ap*4.*cpdi[0]->ws;
		double tipwt =  0.25/pos.x;
		cpdi[0]->wg.z = tipwt;
		
		cpdi[1]->wg.x = (rysum + beta/(3.*pos.x))*Ap;
		cpdi[1]->wg.y = -rxsum*Ap*4.*cpdi[1]->ws;
		cpdi[1]->wg.z = tipwt;
		
		cpdi[2]->wg.x = -(rydiff + beta/(3.*pos.x))*Ap;
		cpdi[2]->wg.y = rxdiff*Ap*4.*cpdi[2]->ws;
		cpdi[2]->wg.z = tipwt;
		
		cpdi[3]->wg.x = -(rysum - beta/(3.*pos.x))*Ap;
		cpdi[3]->wg.y = rxsum*Ap*4.*cpdi[3]->ws;;
		cpdi[3]->wg.z = tipwt;
	}
		
	catch(...)
	{	throw MPMTermination("A CPDI partical domain node has left the grid.","MatPoint2D::GetCPDINodesAndWeights");
	}
    
    // traction BC area saves 1/2 surface area of particle domain on the various edges
    if(faceArea!=NULL)
    {   faceArea->x = sqrt(r1.x*r1.x+r1.y*r1.y)*mpmgrid.GetThickness();			// edges 1 and 3
        faceArea->y = sqrt(r2.x*r2.x+r2.y*r2.y)*mpmgrid.GetThickness();			// edges 2 and 4
    }
}

// To support traction boundary conditions, find the deformed edge, natural coordinates of
// the corners along the edge, elements for those edges, and a normal vector in direction
// of the traction
// return ratio of second nodal weight to first one
double MatPointAS::GetTractionInfo(int face,int dof,int *cElem,Vector *corners,Vector *tscaled,int *numDnds)
{
    *numDnds = 2;
    double faceWt,ratio=1.,r1mag,r2mag;
    double rp=pos.x;
    Vector r1,r2;;
    
    if(ElementBase::useGimp==UNIFORM_GIMP_AS)
    {   r1.x = mpmgrid.gridx*0.25;
        r1.y = 0.;
        r2.x = 0.;
        r2.y = mpmgrid.gridy*0.25;

#ifdef TRUNCATE
        // truncate if extends into r<0
        if(pos.x-r1.x<0.)
        {   r1.x = 0.5*(pos.x+r1.x);
            rp = r1.x;
        }
#endif
        
        // get magnitudes
        r1mag = r1.x;
        r2mag = r2.y;
    }
    
    else
    {   // always LINEAR_CPDI_AS
	
        // get particle 2D deformation gradient
        double pF[3][3];
        GetDeformationGradient(pF);
	
        // get polygon vectors - these are from particle to edge
        //      and generalize semi width lp in 1D GIMP
        r1.x = pF[0][0]*mpmgrid.gridx*0.25;
        r1.y = pF[1][0]*mpmgrid.gridx*0.25;
        r2.x = pF[0][1]*mpmgrid.gridy*0.25;
        r2.y = pF[1][1]*mpmgrid.gridy*0.25;
        
#ifdef TRUNCATE
        // shrink domain if any have x < 0, but keep particle in the middle
        if(pos.x-fabs(r1.x+r2.x)<0.)
        {	// make pos.x-shrink*fabs(r1.x+r2.x) very small and positive
            double shrink = (pos.x-mpmgrid.gridx*1.e-10)/fabs(r1.x+r2.x);
            r1.x *= shrink;
            r1.y *= shrink;
            r2.x *= shrink;
            r2.y *= shrink;
        }
#endif
        
        // get magnitudes
        r1mag = sqrt(r1.x*r1.x + r1.y*r1.y);
        r2mag = sqrt(r2.x*r2.x + r2.y*r2.y);
    }
    
    // Find corners
	Vector c1,c2;
	switch(face)
	{	case 1:
			// lower edge
			c1.x = rp-r1.x-r2.x;         // node 1
			c1.y = pos.y-r1.y-r2.y;
			c2.x = rp+r1.x-r2.x;         // node 2
            c2.y = pos.y+r1.y-r2.y;
			faceWt = r1mag*(rp - r2.x - r1.x/3.);			// node 1, node 2 should be plus
			ratio = r1mag*(rp - r2.x + r1.x/3.)/faceWt;		// find the ratio
			break;
			
		case 2:
			// right edge
			c1.x = rp+r1.x-r2.x;         // node 2
			c1.y = pos.y+r1.y-r2.y;
			c2.x = rp+r1.x+r2.x;         // node 3
			c2.y = pos.y+r1.y+r2.y;
			faceWt = r2mag*(rp + r1.x - r2.x/3.);			// node 2, node 3 should be plus
			ratio = r2mag*(rp + r1.x + r2.x/3.)/faceWt;		// find the ratio
			break;
			
		case 3:
			// top edge
			c1.x = rp+r1.x+r2.x;         // node 3
			c1.y = pos.y+r1.y+r2.y;
			c2.x = rp-r1.x+r2.x;         // node 4
			c2.y = pos.y-r1.y+r2.y;
			faceWt = r1mag*(rp + r2.x + r1.x/3.);			// node 3, node 4 should be minus
			ratio = r1mag*(rp + r2.x - r1.x/3.)/faceWt;		// find the ratio
			break;
			
		default:
			// left edge
			c1.x = rp-r1.x+r2.x;         // node 4
			c1.y = pos.y-r1.y+r2.y;
			c2.x = rp-r1.x-r2.x;         // node 1
			c2.y = pos.y-r1.y-r2.y;
			faceWt = r2mag*(rp - r1.x + r2.x/3.);			// node 4, node 1 should be minus
			ratio = r2mag*(rp - r1.x - r2.x/3.)/faceWt;		// find the ratio
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
	double ex,ey,enorm;
	switch(dof)
	{	case 1:
			// normal is x direction
			tscaled->x = faceWt;
			break;
        case 2:
            // normal is y direction
            tscaled->y = faceWt;
            break;
		case N_DIRECTION:
			// cross product of edge vector with (0,0,1) = (ey, -ex)
			ex = c2.y-c1.y;
			ey = c1.x-c2.x;
		case T1_DIRECTION:
			if(dof==T1_DIRECTION)
			{	ex = c2.x-c1.x;
				ey = c2.y-c1.y;
			}
			// load in direction specified by normalized (ex,ey)
			enorm = sqrt(ex*ex+ey*ey);
			tscaled->x = ex*faceWt/enorm;
			tscaled->y = ey*faceWt/enorm;
			break;
		default:
			// normal is z direction (not used here)
            tscaled->z = faceWt;
			break;
	}
	
	return ratio;
}



