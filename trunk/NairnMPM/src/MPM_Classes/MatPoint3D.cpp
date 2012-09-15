/********************************************************************************
    MatPoint3D.cpp
    NairnMPM
    
    Created by John Nairn on 7/21/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "MPM_Classes/MatPoint3D.hpp"
#include "Materials/MaterialBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"

#pragma mark MatPoint3D::Constructors and Destructors

// Constructors
MatPoint3D::MatPoint3D()
{
}

// Constructors
MatPoint3D::MatPoint3D(int inElemNum,int theMatl,double angin) : MPMBase(inElemNum,theMatl,angin)
{
}

#pragma mark MatPoint2D:Calculations and Incrementers

// Update Strains for this particle
// Velocities for all fields on present on the nodes
void MatPoint3D::UpdateStrain(double strainTime,int secondPass,int np)
{
	// material for this particle
	MaterialBase *matRef=theMaterials[MatID()];
	
	// exit if rigid
	if(matRef->Rigid()) return;
	
	// make sure mechanical properties for this material and angle
	matRef->LoadMechanicalProps(this,np);
	
	// material field
	int matfld=matRef->GetField();
	
	int i,numnds,nds[MaxShapeNds];
    double fn[MaxShapeNds],xDeriv[MaxShapeNds],yDeriv[MaxShapeNds],zDeriv[MaxShapeNds];
	Vector vel;
    double dvxx,dvyy,dvzz,dvxy,dvyx,dvxz,dvzx,dvyz,dvzy;
    
	// find shape functions and derviatives
	int iel=ElemID();
	theElements[iel]->GetShapeGradients(&numnds,fn,nds,&ncpos,xDeriv,yDeriv,zDeriv,this);
    
    // Find strain rates at particle from current grid velocities
	//   and using the velocity field for that particle with each node
    dvxx=dvyy=dvxy=dvyx=dvzz=dvxz=dvzx=dvyz=dvzy=0.;
    for(i=1;i<=numnds;i++)
	{	vel=nd[nds[i]]->GetVelocity((short)vfld[i],matfld);
        dvxx+=vel.x*xDeriv[i];
        dvyy+=vel.y*yDeriv[i];
        dvzz+=vel.z*zDeriv[i];
        dvxy+=vel.x*yDeriv[i];
        dvyx+=vel.y*xDeriv[i];
        dvxz+=vel.x*zDeriv[i];
        dvzx+=vel.z*xDeriv[i];
        dvyz+=vel.y*zDeriv[i];
        dvzy+=vel.z*yDeriv[i];
    }
	    
    // convert to strain increments
    dvxx*=strainTime;
    dvyy*=strainTime;
	dvzz*=strainTime;
    dvxy*=strainTime;
    dvyx*=strainTime;
    dvxz*=strainTime;
    dvzx*=strainTime;
    dvyz*=strainTime;
    dvzy*=strainTime;
    
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
    matRef->MPMConstLaw(this,dvxx,dvyy,dvzz,dvxy,dvyx,dvxz,dvzx,dvyz,dvzy,strainTime,np);
}

// Move position (3D) (in mm)
// external work units g-mm^2/sec^2 (* 10^-9 to get J)
void MatPoint3D::MovePosition(double delTime,Vector *dv)
{	double dx=delTime*dv->x;
	double dy=delTime*dv->y;
	double dz=delTime*dv->z;
	pos.x+=dx;
    pos.y+=dy;
    pos.z+=dz;
    extWork+=dx*pFext.x+dy*pFext.y+dz*pFext.z;
}

// Move velocity (3D) (in mm/sec) possible with damping
void MatPoint3D::MoveVelocity(double delTime,double damping,Vector *vstar)
{	vel.x+=delTime*(acc.x-damping*vstar->x);
    vel.y+=delTime*(acc.y-damping*vstar->y);
    vel.z+=delTime*(acc.z-damping*vstar->z);
}

// Scale velocity (3D) (in mm/sec) optionally with rigid materials
void MatPoint3D::SetVelocitySpeed(double speed)
{	double norm=sqrt(vel.x*vel.x+vel.y*vel.y+vel.z*vel.z);
	if(DbleEqual(norm,0.)) return;
	vel.x*=speed/norm;
    vel.y*=speed/norm;
    vel.z*=speed/norm;
}

#pragma mark MatPoint3D::Accessors

// set original position (3D) (in mm)
void MatPoint3D::SetOrigin(Vector *pt) { origpos=*pt; }

// set position (3D) (in mm)
void MatPoint3D::SetPosition(Vector *pt) { pos=*pt; }

// set position (3D) (in mm/sec)
void MatPoint3D::SetVelocity(Vector *v) { vel=*v; }

// no thickness
double MatPoint3D::thickness() { return -1.; }

// set current volume using current strains
void MatPoint3D::SetDilatedVolume(void)
{	double rho=theMaterials[MatID()]->rho*0.001;			// in g/mm^3
	double dilate=(1.+ep.xx)*(1.+ep.yy)*(1+ep.zz);
	volume=dilate*mp/rho;									// in mm^3
}

// calculate internal force as -mp sigma.deriv * 1000.
void MatPoint3D::Fint(Vector &fout,double xDeriv,double yDeriv,double zDeriv)
{	fout.x=-mp*(sp.xx*xDeriv+sp.xy*yDeriv+sp.xz*zDeriv)*1000.;
	fout.y=-mp*(sp.xy*xDeriv+sp.yy*yDeriv+sp.yz*zDeriv)*1000.;
	fout.z=-mp*(sp.xz*xDeriv+sp.yz*yDeriv+sp.zz*zDeriv)*1000.;
}

// external force (times a shape function)
void MatPoint3D::Fext(Vector &fout,double fni)
{	fout.x=fni*pFext.x;
	fout.y=fni*pFext.y;
	fout.z=fni*pFext.z;
}

// zero the temperature gradient
void MatPoint3D::AddTemperatureGradient(void)
{	pTemp->DT.x=0.;
    pTemp->DT.y=0.;
    pTemp->DT.z=0.;
}

// add to the concentration gradient
void MatPoint3D::AddTemperatureGradient(Vector *grad)
{	pTemp->DT.x+=grad->x;
    pTemp->DT.y+=grad->y;
    pTemp->DT.z+=grad->z;
}

// return conduction force = - V [D] Grad T . Grad S
double MatPoint3D::FCond(double dshdx,double dshdy,double dshdz)
{
	Tensor *kten=theMaterials[MatID()]->GetkCondTensor();
	return -volume*((kten->xx*pTemp->DT.x + kten->xy*pTemp->DT.y + kten->xz*pTemp->DT.z)*dshdx
						+ (kten->xy*pTemp->DT.x + kten->yy*pTemp->DT.y + kten->yz*pTemp->DT.z)*dshdy
						+ (kten->xz*pTemp->DT.x + kten->yz*pTemp->DT.y + kten->zz*pTemp->DT.z)*dshdz);
}

// zero the concentration gradient
void MatPoint3D::AddConcentrationGradient(void)
{	pDiffusion->Dc.x=0.;
    pDiffusion->Dc.y=0.;
    pDiffusion->Dc.z=0.;
}

// add to the concentration gradient
void MatPoint3D::AddConcentrationGradient(Vector *grad)
{	pDiffusion->Dc.x+=grad->x;
    pDiffusion->Dc.y+=grad->y;
    pDiffusion->Dc.z+=grad->z;
}

// return diffusion force = - V [D] Grad C . Grad S
double MatPoint3D::FDiff(double dshdx,double dshdy,double dshdz)
{
	Tensor *Dten=theMaterials[MatID()]->GetDiffusionTensor();
	return -volume*((Dten->xx*pDiffusion->Dc.x + Dten->xy*pDiffusion->Dc.y + Dten->xz*pDiffusion->Dc.z)*dshdx
						+ (Dten->xy*pDiffusion->Dc.x + Dten->yy*pDiffusion->Dc.y + Dten->yz*pDiffusion->Dc.z)*dshdy
						+ (Dten->xz*pDiffusion->Dc.x + Dten->yz*pDiffusion->Dc.y + Dten->zz*pDiffusion->Dc.z)*dshdz);
}

// return kinetic energy
double MatPoint3D::KineticEnergy(void)
{	return 0.5*mp*(vel.x*vel.x+vel.y*vel.y+vel.z*vel.z);
}

// get deformation gradient, which is stored in strain and rotation tensors
void MatPoint3D::GetDeformationGradient(double F[][3])
{
	// current deformation gradient in 2D
	F[0][0] = 1. + ep.xx + eplast.xx;
	F[1][1] = 1. + ep.yy + eplast.yy;
	F[2][2] = 1. + ep.zz + eplast.zz;
    double exy = ep.xy + eplast.xy;
	F[0][1] = 0.5*(exy - wrot.xy);
	F[1][0] = 0.5*(exy + wrot.xy);
    double exz = ep.xz + eplast.xz;
	F[0][2] = 0.5*(exz - wrot.xz);
	F[2][0] = 0.5*(exz + wrot.xz);
    double eyz = ep.yz + eplast.yz;
	F[1][2] = 0.5*(eyz - wrot.yz);
	F[2][1] = 0.5*(eyz + wrot.yz);
}

// get relative volume from det J for large deformation material laws
double MatPoint3D::GetRelativeVolume(void)
{   double pF[3][3];
    GetDeformationGradient(pF);
    return pF[0][0]*(pF[1][1]*pF[2][2]-pF[2][1]*pF[1][2])
                - pF[0][1]*(pF[1][0]*pF[2][2]-pF[2][0]*pF[1][2])
                + pF[0][2]*(pF[1][0]*pF[2][1]-pF[2][0]*pF[1][1]);
}

// to support CPDI return nodes for corners (or for 9 nodes) and weights for shape functions
//	and shape function gradients
void MatPoint3D::GetCPDINodesAndWeights(int cpdiType)
{
	throw "CPDI nodes not yet available for 3D";
}

// To support traction boundary conditions, find the deformed edge, natural coordinates of
// the corners around the face, elements for those faces, and a normal vector in direction
// of the traction. Input vectors need to be length 4
void MatPoint3D::GetTractionInfo(int face,int dof,int *cElem,Vector *corners,Vector *tnorm,int *numDnds)
{
	throw "Traction boundary conditions not yet available for 3D";
}

