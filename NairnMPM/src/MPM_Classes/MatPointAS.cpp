/********************************************************************************
	MatPointAS.hpp
	nairn-mpm-fea

	Created by John Nairn on 10/23/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
 
	When MPM is axisymmetric, the x axis becomes radial (or R) direction, the
		y axis because axial direction (or Z direction), and the z direction
		become the tangential direction (or theta direction)
 
	The thickness (in thick) is set to the initial radial position of the
		material point
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MatPointAS.hpp"
#include "Materials/MaterialBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Exceptions/CommonException.hpp"

#include "System/MPMPrefix.hpp"

#pragma mark MatPointAS::Constructors and Destructors

// Constructors
MatPointAS::MatPointAS() {}

// Constructors
MatPointAS::MatPointAS(int inElemNum,int theMatl,double angin,double thickin) : MatPoint2D(inElemNum,theMatl,angin,thickin)
{	
}

#pragma mark MatPointAS:Calculations and Incrementers

// Update Strains for this particle
// Velocities for all fields are present on the nodes
// matRef is the material and properties have been loaded, matFld is the material field
// postUpdate true if after particle position updated
void MatPointAS::UpdateStrain(double strainTime,int secondPass,int np,void *props,int matFld,bool postUpdate)
{
#ifdef CONST_ARRAYS
	int ndsArray[MAX_SHAPE_NODES];
	double fn[MAX_SHAPE_NODES],xDeriv[MAX_SHAPE_NODES],yDeriv[MAX_SHAPE_NODES],zDeriv[MAX_SHAPE_NODES];
#else
	int ndsArray[maxShapeNodes];
	double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
#endif
	int i,numnds;
	Vector vel;
    Matrix3 dv;
	Tensor *gStressPtr=NULL;

	// don't need to zero zDeriv because always set in axisymmetric shape functions
    
	// find shape functions and derviatives
	const ElementBase *elemRef = theElements[ElemID()];
	int *nds = ndsArray;
	elemRef->GetShapeGradients(fn,&nds,xDeriv,yDeriv,zDeriv,this);
	numnds = nds[0];
    
    // Find strain rates at particle from current grid velocities
	//   and using the velocity field for that particle and each node and the right material
    // In axisymmetric x->r, y->z, and z->hoop
    for(i=1;i<=numnds;i++)
	{	vel=nd[nds[i]]->GetVelocity((short)vfld[i],matFld);
        dv += Matrix3(vel.x*xDeriv[i],vel.x*yDeriv[i],vel.y*xDeriv[i],vel.y*yDeriv[i],vel.x*zDeriv[i]);
    }
    
    // save velocity gradient (if needed for J integral calculation)
    SetVelocityGradient(dv(0,0),dv(1,0),dv(0,1),dv(1,0),secondPass);
    
	// convert to strain increments
    // e.g., now dvrr = dvr/dr * dt = d/dr(du/dt) * dt = d/dt(du/dr) * dt = du/dr)
    dv.Scale(strainTime);
    
	// pass on to material class to handle
	ResidualStrains res = ScaledResidualStrains(secondPass);
	PerformConstitutiveLaw(dv,strainTime,np,props,&res,gStressPtr);
}

// Pass on to material class
void MatPointAS::PerformConstitutiveLaw(Matrix3 dv,double strainTime,int np,void *props,ResidualStrains *res,Tensor *gStress)
{
    // update particle strain and stress using its constitutive law
	const MaterialBase *matRef = theMaterials[MatID()];
    matRef->MPMConstitutiveLaw(this,dv,strainTime,np,props,res,0,gStress);
}

#pragma mark MatPoint2D::Accessors

// set original position (2D) (in mm) and here set thick to be radial position
void MatPointAS::SetOrigin(Vector *pt)
{	origpos.x=pt->x;
    origpos.y=pt->y;
	origpos.z=0.;
	thick = pt->x;
}

// Find internal force as -mp sigma.deriv in g mm/sec^2 or micro N
// add external force (times a shape function)
// Store in buffer
// (note: stress is specific stress in units N/m^2 mm^3/g which is (g-mm^2/sec^2)/g
void MatPointAS::GetFintPlusFext(Vector *theFrc,double fni,double xDeriv,double yDeriv,double zDeriv)
{
	theFrc->x = -mp*((sp.xx-pressure)*xDeriv+sp.xy*yDeriv+(sp.zz-pressure)*zDeriv) + fni*pFext.x;
	theFrc->y = -mp*(sp.xy*xDeriv+(sp.yy-pressure)*yDeriv) + fni*pFext.y;
	theFrc->z = 0.0;
}

// Get dilated current volume using current deformation gradient
// (only used for crack contact, multimaterial contact, and transport tasks)
// When volumeType is DEFORMED_AREA or DEFORMED_AREA_FOR_GRADIENT, get Area deformed particle area in r-z plane
double MatPointAS::GetVolume(int volumeType)
{	double rho=GetRho();							// in g/mm^3
	if(volumeType==DEFORMED_VOLUME)
		return GetRelativeVolume()*mp/rho;			// in mm^3
	
	// get deformed area in the r-z plane
	// note that mp/rho = rho Ap0 rp0 / rho = Ap0 rp0
	double pF[3][3];
	GetDeformationGradient(pF);
	double areaPerRadian = (pF[0][0]*pF[1][1]-pF[1][0]*pF[0][1])*mp/rho;
	return areaPerRadian/origpos.x;
}

// get unscaled volume for use only in contact and imperfect interface calculations
// For axisymmetric, really needs area of the initial particle domain in mm^3
// Calculations will need to multiply by radial position to get local volume
// Here thickness is the original radial position of the particle
double MatPointAS::GetUnscaledVolume(void)
{	double rho=GetRho();							// in g/mm^3
	return mp/(rho*thickness());                    // in mm^3 per unit radial position
}

// To support CPDI find nodes in the particle domain, find their elements,
// their natural coordinates, and weighting values for gradient calculations
// Should be done only once per time step and here only for axisynmmetric
// and linear CPDI
// throws CommonException() if particle corner has left the grid
void MatPointAS::GetCPDINodesAndWeights(int cpdiType)
{	
	// get polygon vectors - these are from particle to edge
    //      and generalize semi width lp in 1D GIMP
	Vector r1,r2,c;
    GetSemiSideVectors(&r1,&r2,NULL);

	// truncate domain by shrinking if any have x < 0
	if(pos.x-fabs(r1.x)-fabs(r2.x)<=0.)
	{	// make pos.x=shrink*(fabs(r1.x)+fabs(r2.x))
		double shrink = pos.x/(fabs(r1.x)+fabs(r2.x));
		r1.x *= shrink;
		r1.y *= shrink;
		r2.x *= shrink;
		r2.y *= shrink;
	}
	
	// rescale if needed
	ScaleSemiSideVectorsForCPDI(&r1,&r2,NULL);
	
    // Particle domain area is area of the full parallelogram
    // Assume positive due to orientation of initial vectors, and sign probably does not matter
    double Ap = 4.*(r1.x*r2.y - r1.y*r2.x);
    
	try
	{	CPDIDomain **cpdi = GetCPDIInfo();
		
		// nodes at four courves in ccw direction
		c.x = pos.x-r1.x-r2.x;
		c.y = pos.y-r1.y-r2.y;
		cpdi[0]->inElem = mpmgrid.FindElementFromPoint(&c,this)-1;
		theElements[cpdi[0]->inElem]->GetXiPos(&c,&cpdi[0]->ncpos);
		
		c.x = pos.x+r1.x-r2.x;
		c.y = pos.y+r1.y-r2.y;
		cpdi[1]->inElem = mpmgrid.FindElementFromPoint(&c,this)-1;
		theElements[cpdi[1]->inElem]->GetXiPos(&c,&cpdi[1]->ncpos);
		
		c.x = pos.x+r1.x+r2.x;
		c.y = pos.y+r1.y+r2.y;
		cpdi[2]->inElem = mpmgrid.FindElementFromPoint(&c,this)-1;
		theElements[cpdi[2]->inElem]->GetXiPos(&c,&cpdi[2]->ncpos);
		
		c.x = pos.x-r1.x+r2.x;
		c.y = pos.y-r1.y+r2.y;
		cpdi[3]->inElem = mpmgrid.FindElementFromPoint(&c,this)-1;
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
		
    catch(CommonException& err)
    {   char msg[200];
		size_t msgSize=200;
		snprintf(msg,msgSize,"A CPDI particle domain node has left the grid: %s",err.Message());
        throw CommonException(msg,"MatPointAS::GetCPDINodesAndWeights");
    }
    
    // traction BC area saves 1/2 surface area of particle domain on the various edges
    if(faceArea!=NULL)
    {   faceArea->x = sqrt(r1.x*r1.x+r1.y*r1.y)*mpmgrid.GetThickness();			// edges 1 and 3
        faceArea->y = sqrt(r2.x*r2.x+r2.y*r2.y)*mpmgrid.GetThickness();			// edges 2 and 4
    }
}
