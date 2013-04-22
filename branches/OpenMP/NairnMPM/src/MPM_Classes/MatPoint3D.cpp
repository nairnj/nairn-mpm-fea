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
#include "NairnMPM_Class/MeshInfo.hpp"

// NEWINCLUDE
#include "Exceptions/CommonException.hpp"

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
// Velocities for all fields are present on the nodes
// matRef is the material and properties have been loaded, matFld is the material field
void MatPoint3D::UpdateStrain(double strainTime,int secondPass,int np,void *props,int matFld)
{
	int i,numnds,nds[maxShapeNodes];
    double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
	Vector vel;
    Matrix3 dv;
    
	// find shape functions and derviatives
	const ElementBase *elemRef = theElements[ElemID()];
	elemRef->GetShapeGradients(&numnds,fn,nds,&ncpos,xDeriv,yDeriv,zDeriv,this);
    
    // Find strain rates at particle from current grid velocities
	//   and using the velocity field for that particle with each node
    for(i=1;i<=numnds;i++)
	{	vel=nd[nds[i]]->GetVelocity((short)vfld[i],matFld);
        dv += Matrix3(vel.x,vel.y,vel.z,xDeriv[i],yDeriv[i],zDeriv[i]);
    }
	    
    // convert to strain increments
    dv.Scale(strainTime);
    
	// find effective particle transport properties from grid results
	ResidualStrains res;
	res.dT = 0;
	res.dC = 0.;
	if(!ConductionTask::active)
	{	res.dT = pTemperature-pPreviousTemperature;
		pPreviousTemperature = pTemperature;
	}
	else
	{	for(i=1;i<=numnds;i++)
		res.dT += conduction->IncrementValueExtrap(nd[nds[i]],fn[i]);
		res.dT = conduction->GetDeltaValue(this,res.dT);
	}
	if(DiffusionTask::active)
	{	for(i=1;i<=numnds;i++)
		res.dC += diffusion->IncrementValueExtrap(nd[nds[i]],fn[i]);
		res.dC = diffusion->GetDeltaValue(this,res.dC);
	}

    // update particle strain and stress using its constituitive law
    const MaterialBase *matRef = theMaterials[MatID()];
    matRef->MPMConstitutiveLaw(this,dv,strainTime,np,props,&res);
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

// calculate internal force as -mp sigma.deriv * 1000.
// add external force (times a shape function)
// store in buffer
// (note: stress is specific stress in units N/m^2 cm^3/g, Multiply by 1000 to make it mm/sec^2)
void MatPoint3D::GetFintPlusFext(Vector *theFrc,double fni,double xDeriv,double yDeriv,double zDeriv)
{	
	double mpug = mp*1000.;
	theFrc->x = -mpug*((sp.xx-pressure)*xDeriv+sp.xy*yDeriv+sp.xz*zDeriv) + fni*pFext.x;
	theFrc->y = -mpug*(sp.xy*xDeriv+(sp.yy-pressure)*yDeriv+sp.yz*zDeriv) + fni*pFext.y;
	theFrc->z = -mpug*(sp.xz*xDeriv+sp.yz*yDeriv+(sp.zz-pressure)*zDeriv) + fni*pFext.z;
}

// add to the concentration gradient (non-rigid particles only)
void MatPoint3D::AddTemperatureGradient(Vector *grad)
{	pTemp->DT.x+=grad->x;
    pTemp->DT.y+=grad->y;
    pTemp->DT.z+=grad->z;
}

// return conduction force = - V [D] Grad T . Grad S
double MatPoint3D::FCond(double dshdx,double dshdy,double dshdz,TransportProperties *t)
{
	Tensor *kten = &(t->kCondTensor);
	return -mp*GetRelativeVolume()*((kten->xx*pTemp->DT.x + kten->xy*pTemp->DT.y + kten->xz*pTemp->DT.z)*dshdx
						+ (kten->xy*pTemp->DT.x + kten->yy*pTemp->DT.y + kten->yz*pTemp->DT.z)*dshdy
						+ (kten->xz*pTemp->DT.x + kten->yz*pTemp->DT.y + kten->zz*pTemp->DT.z)*dshdz);
}

// add to the concentration gradient
void MatPoint3D::AddConcentrationGradient(Vector *grad)
{	pDiffusion->Dc.x+=grad->x;
    pDiffusion->Dc.y+=grad->y;
    pDiffusion->Dc.z+=grad->z;
}

// return diffusion force = - V [D] Grad C . Grad S
double MatPoint3D::FDiff(double dshdx,double dshdy,double dshdz,TransportProperties *t)
{
	Tensor *Dten = &(t->diffusionTensor);
	return -GetVolume(DEFORMED_VOLUME)*((Dten->xx*pDiffusion->Dc.x + Dten->xy*pDiffusion->Dc.y + Dten->xz*pDiffusion->Dc.z)*dshdx
						+ (Dten->xy*pDiffusion->Dc.x + Dten->yy*pDiffusion->Dc.y + Dten->yz*pDiffusion->Dc.z)*dshdy
						+ (Dten->xz*pDiffusion->Dc.x + Dten->yz*pDiffusion->Dc.y + Dten->zz*pDiffusion->Dc.z)*dshdz);
}

// return kinetic energy (g mm^2/sec^2) = nanoJ
double MatPoint3D::KineticEnergy(void)
{	return 0.5*mp*(vel.x*vel.x+vel.y*vel.y+vel.z*vel.z);
}

// get deformation gradient, which is stored in strain and rotation tensors
Matrix3 MatPoint3D::GetDeformationGradientMatrix(void)
{	double F[3][3];
	GetDeformationGradient(F);
	Matrix3 Fm(F[0][0],F[0][1],F[0][2],F[1][0],F[1][1],F[1][2],F[2][0],F[2][1],F[2][2]);
	return Fm;
}

// get the symmetric elastic Left-Cauchy tensor in a Matrix3
Matrix3 MatPoint3D::GetElasticLeftCauchyMatrix(void)
{   return Matrix3(eplast.xx,eplast.xy,eplast.xz,eplast.xy,eplast.yy,eplast.yz,
                    eplast.xz,eplast.yz,eplast.zz);
}

// get deformation gradient, which is stored in strain and rotation tensors
void MatPoint3D::GetDeformationGradient(double F[][3])
{
	// current deformation gradient in 3D
    if(theMaterials[MatID()]->PartitionsElasticAndPlasticStrain())
    {   F[0][0] = 1. + ep.xx + eplast.xx;
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
    else
    {   F[0][0] = 1. + ep.xx;
        F[1][1] = 1. + ep.yy;
        F[2][2] = 1. + ep.zz;
        F[0][1] = 0.5*(ep.xy - wrot.xy);
        F[1][0] = 0.5*(ep.xy + wrot.xy);
        F[0][2] = 0.5*(ep.xz - wrot.xz);
        F[2][0] = 0.5*(ep.xz + wrot.xz);
        F[1][2] = 0.5*(ep.yz - wrot.yz);
        F[2][1] = 0.5*(ep.yz + wrot.yz);
    }
}

// get relative volume from det J for large deformation material laws
double MatPoint3D::GetRelativeVolume(void)
{   double pF[3][3];
    GetDeformationGradient(pF);
    return pF[0][0]*(pF[1][1]*pF[2][2]-pF[2][1]*pF[1][2])
                - pF[0][1]*(pF[1][0]*pF[2][2]-pF[2][0]*pF[1][2])
                + pF[0][2]*(pF[1][0]*pF[2][1]-pF[2][0]*pF[1][1]);
}

// get dilated current volume using current deformation gradient
// only used for contact (cracks and multimaterial) and for transport tasks
double MatPoint3D::GetVolume(bool volumeType)
{	double rho=theMaterials[MatID()]->rho*0.001;				// in g/mm^3
	return GetRelativeVolume()*mp/rho;							// in mm^3
}

// to support CPDI return nodes for corners (or for 9 nodes) and weights for shape functions
//	and shape function gradients
// throws CommonException() if particle corner has left the grid
void MatPoint3D::GetCPDINodesAndWeights(int cpdiType)
{
	// get particle 2D deformation gradient
	double pF[3][3];
	GetDeformationGradient(pF);
	
	// get polygon vectors - these are from particle to edge
    //      and generalize semi width lp in 1D GIMP
	Vector r1,r2,r3,c;
	r1.x = pF[0][0]*mpmgrid.partx;
	r1.y = pF[1][0]*mpmgrid.partx;
	r1.z = pF[2][0]*mpmgrid.partx;
	r2.x = pF[0][1]*mpmgrid.party;
	r2.y = pF[1][1]*mpmgrid.party;
	r2.z = pF[2][1]*mpmgrid.party;
	r3.x = pF[0][2]*mpmgrid.partz;
	r3.y = pF[1][2]*mpmgrid.partz;
	r3.z = pF[2][2]*mpmgrid.partz;
	
    // Particle domain volume is 8 * volume of the parallelepiped defined by r1, r2, and r3
	// V = 8 * (r1 . (r2 X r3))
    // Assume positive due to orientation of initial vectors, and sign probably does not matter
    double Vp = 8.*(r1.x * (r2.y*r3.z - r2.z*r3.y)
                    + r1.y * (r2.z*r3.x - r2.x*r3.z)
                     + r1.z * (r2.x*r3.y - r2.y*r3.x) );
	
	// always LINEAR_CPDI
    
	try
	{	// nodes at 8 node
		c.x = pos.x-r1.x-r2.x-r3.x;
		c.y = pos.y-r1.y-r2.y-r3.y;
		c.z = pos.z-r1.z-r2.z-r3.z;
		cpdi[0]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[0]->inElem]->GetXiPos(&c,&cpdi[0]->ncpos);
		
		c.x = pos.x+r1.x-r2.x-r3.x;
		c.y = pos.y+r1.y-r2.y-r3.y;
		c.z = pos.z+r1.z-r2.z-r3.z;
		cpdi[1]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[1]->inElem]->GetXiPos(&c,&cpdi[1]->ncpos);
		
		c.x = pos.x+r1.x+r2.x-r3.x;
		c.y = pos.y+r1.y+r2.y-r3.y;
		c.z = pos.z+r1.z+r2.z-r3.z;
		cpdi[2]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[2]->inElem]->GetXiPos(&c,&cpdi[2]->ncpos);
		
		c.x = pos.x-r1.x+r2.x-r3.x;
		c.y = pos.y-r1.y+r2.y-r3.y;
		c.z = pos.z-r1.z+r2.z-r3.z;
		cpdi[3]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[3]->inElem]->GetXiPos(&c,&cpdi[3]->ncpos);
		
		c.x = pos.x-r1.x-r2.x+r3.x;
		c.y = pos.y-r1.y-r2.y+r3.y;
		c.z = pos.z-r1.z-r2.z+r3.z;
		cpdi[4]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[4]->inElem]->GetXiPos(&c,&cpdi[4]->ncpos);
		
		c.x = pos.x+r1.x-r2.x+r3.x;
		c.y = pos.y+r1.y-r2.y+r3.y;
		c.z = pos.z+r1.z-r2.z+r3.z;
		cpdi[5]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[5]->inElem]->GetXiPos(&c,&cpdi[5]->ncpos);
		
		c.x = pos.x+r1.x+r2.x+r3.x;
		c.y = pos.y+r1.y+r2.y+r3.y;
		c.z = pos.z+r1.z+r2.z+r3.z;
		cpdi[6]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[6]->inElem]->GetXiPos(&c,&cpdi[6]->ncpos);
		
		c.x = pos.x-r1.x+r2.x+r3.x;
		c.y = pos.y-r1.y+r2.y+r3.y;
		c.z = pos.z-r1.z+r2.z+r3.z;
		cpdi[7]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
		theElements[cpdi[7]->inElem]->GetXiPos(&c,&cpdi[7]->ncpos);

		// gradient weighting values
		Vp = 1./Vp;
		cpdi[0]->wg.x = (r1.z*r2.y - r1.y*r2.z - r1.z*r3.y + r2.z*r3.y + r1.y*r3.z - r2.y*r3.z)*Vp;
		cpdi[0]->wg.y = (-(r1.z*r2.x) + r1.x*r2.z + r1.z*r3.x - r2.z*r3.x - r1.x*r3.z + r2.x*r3.z)*Vp;
		cpdi[0]->wg.z = (r1.y*r2.x - r1.x*r2.y - r1.y*r3.x + r2.y*r3.x + r1.x*r3.y - r2.x*r3.y)*Vp;
		cpdi[1]->wg.x = (r1.z*r2.y - r1.y*r2.z - r1.z*r3.y - r2.z*r3.y + r1.y*r3.z + r2.y*r3.z)*Vp;
		cpdi[1]->wg.y = (-(r1.z*r2.x) + r1.x*r2.z + r1.z*r3.x + r2.z*r3.x - r1.x*r3.z - r2.x*r3.z)*Vp;
		cpdi[1]->wg.z = (r1.y*r2.x - r1.x*r2.y - r1.y*r3.x - r2.y*r3.x + r1.x*r3.y + r2.x*r3.y)*Vp;
		cpdi[2]->wg.x = (r1.z*r2.y - r1.y*r2.z + r1.z*r3.y - r2.z*r3.y - r1.y*r3.z + r2.y*r3.z)*Vp;
		cpdi[2]->wg.y = (-(r1.z*r2.x) + r1.x*r2.z - r1.z*r3.x + r2.z*r3.x + r1.x*r3.z - r2.x*r3.z)*Vp;
		cpdi[2]->wg.z = (r1.y*r2.x - r1.x*r2.y + r1.y*r3.x - r2.y*r3.x - r1.x*r3.y + r2.x*r3.y)*Vp;
		cpdi[3]->wg.x = (r1.z*r2.y - r1.y*r2.z + r1.z*r3.y + r2.z*r3.y - r1.y*r3.z - r2.y*r3.z)*Vp;
		cpdi[3]->wg.y = (-(r1.z*r2.x) + r1.x*r2.z - r1.z*r3.x - r2.z*r3.x + r1.x*r3.z + r2.x*r3.z)*Vp;
		cpdi[3]->wg.z = (r1.y*r2.x - r1.x*r2.y + r1.y*r3.x + r2.y*r3.x - r1.x*r3.y - r2.x*r3.y)*Vp;
		cpdi[4]->wg.x = (-(r1.z*r2.y) + r1.y*r2.z - r1.z*r3.y + r2.z*r3.y + r1.y*r3.z - r2.y*r3.z)*Vp;
		cpdi[4]->wg.y = (r1.z*r2.x - r1.x*r2.z + r1.z*r3.x - r2.z*r3.x - r1.x*r3.z + r2.x*r3.z)*Vp;
		cpdi[4]->wg.z = (-(r1.y*r2.x) + r1.x*r2.y - r1.y*r3.x + r2.y*r3.x + r1.x*r3.y - r2.x*r3.y)*Vp;
		cpdi[5]->wg.x = (-(r1.z*r2.y) + r1.y*r2.z - r1.z*r3.y - r2.z*r3.y + r1.y*r3.z + r2.y*r3.z)*Vp;
		cpdi[5]->wg.y = (r1.z*r2.x - r1.x*r2.z + r1.z*r3.x + r2.z*r3.x - r1.x*r3.z - r2.x*r3.z)*Vp;
		cpdi[5]->wg.z = (-(r1.y*r2.x) + r1.x*r2.y - r1.y*r3.x - r2.y*r3.x + r1.x*r3.y + r2.x*r3.y)*Vp;
		cpdi[6]->wg.x = (-(r1.z*r2.y) + r1.y*r2.z + r1.z*r3.y - r2.z*r3.y - r1.y*r3.z + r2.y*r3.z)*Vp;
		cpdi[6]->wg.y = (r1.z*r2.x - r1.x*r2.z - r1.z*r3.x + r2.z*r3.x + r1.x*r3.z - r2.x*r3.z)*Vp;
		cpdi[6]->wg.z = (-(r1.y*r2.x) + r1.x*r2.y + r1.y*r3.x - r2.y*r3.x - r1.x*r3.y + r2.x*r3.y)*Vp;
		cpdi[7]->wg.x = (-(r1.z*r2.y) + r1.y*r2.z + r1.z*r3.y + r2.z*r3.y - r1.y*r3.z - r2.y*r3.z)*Vp;
		cpdi[7]->wg.y = (r1.z*r2.x - r1.x*r2.z - r1.z*r3.x - r2.z*r3.x + r1.x*r3.z + r2.x*r3.z)*Vp;
		cpdi[7]->wg.z = (-(r1.y*r2.x) + r1.x*r2.y + r1.y*r3.x + r2.y*r3.x - r1.x*r3.y - r2.x*r3.y)*Vp;
	}
	catch(...)
	{	throw CommonException("A CPDI particle domain node has left the grid.","MatPoint3D::GetCPDINodesAndWeights");
	}
    
    // traction BC area saves 1/4 the total face area
    if(faceArea!=NULL)
	{	// edges 1 and 3 = |r1 X r3|
		c.x = r1.y*r3.z - r1.z*r3.y;
		c.y = r1.z*r3.x - r1.x*r3.z;
		c.z = r1.x*r3.y - r1.y*r3.x;
		faceArea->x = sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
		// edges 2 and 4 = |r2 X r3|
		c.x = r2.y*r3.z - r2.z*r3.y;
		c.y = r2.z*r3.x - r2.x*r3.z;
		c.z = r2.x*r3.y - r2.y*r3.x;
		faceArea->y = sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
		// top and bottom = |r1 X r2|
		c.x = r1.y*r2.z - r1.z*r2.y;
		c.y = r1.z*r2.x - r1.x*r2.z;
		c.z = r1.x*r2.y - r1.y*r2.x;
		faceArea->z = sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
    }
}

// To support traction boundary conditions, find the deformed edge, natural coordinates of
// the corners around the face, elements for those faces, and a normal vector in direction
// of the traction. Input vectors need to be length 4
// throws CommonException() in traction edge as left the grid
double MatPoint3D::GetTractionInfo(int face,int dof,int *cElem,Vector *corners,Vector *tscaled,int *numDnds)
{
    *numDnds = 4;
    double faceWt;
	
	// always UNIFORM_GIMP or LINEAR_CPDI
    
    // which GIMP method (cannot be used in POINT_GIMP)
    if(ElementBase::useGimp==UNIFORM_GIMP)
    {   // initial vectors only
        double r1x = mpmgrid.partx;
        double r2y = mpmgrid.party;
		double r3z = mpmgrid.partz;
        
        Vector c1,c2,c3,c4;
        switch(face)
        {	case 1:
                // lower face n = (0,-1,0)
                c1.x = c3.x = pos.x-r1x;
                c2.x = c4.x = pos.x+r1x;
                c1.y = c2.y = c3.y = c4.y = pos.y-r2y;
				c1.z = c2.z = pos.z-r3z;
				c3.z = c4.z = pos.z+r3z;
                faceWt = r1x*r3z;
                break;
                
            case 2:
                // right face n = (1,0,0)
                c1.x = c3.x = c2.x = c4.x = pos.x+r1x;
                c1.y = c3.y = pos.y-r2y;
                c2.y = c4.y = pos.y+r2y;
				c1.z = c2.z = pos.z-r3z;
				c3.z = c4.z = pos.z+r3z;
                faceWt = r2y*r3z;
                break;
                
            case 3:
                // top face n = (0,1,0)
                c1.x = c3.x = pos.x+r1x;
                c2.x = c4.x = pos.x-r1x;
                c1.y = c2.y = c3.y = c4.y = pos.y+r2y;
				c1.z = c2.z = pos.z-r3z;
				c3.z = c4.z = pos.z+r3z;
                faceWt = r1x*r3z;
                break;
                
            case 4:
                // left face n = (-1,0,0)
                c1.x = c3.x = c2.x = c4.x = pos.x-r1x;
                c1.y = c3.y = pos.y+r2y;
                c2.y = c4.y = pos.y-r2y;
				c1.z = c2.z = pos.z-r3z;
				c3.z = c4.z = pos.z+r3z;
                faceWt = r2y*r3z;
                break;
			
			case 5:
                // bottom face n = (0,0,-1)
                c1.x = c2.x = pos.x-r1x;
                c3.x = c4.x = pos.x+r1x;
                c1.y = c3.y = pos.y-r2y;
                c2.y = c4.y = pos.y+r2y;
				c1.z = c2.z = c3.z = c4.z = pos.z-r3z;
                faceWt = r2y*r1x;
                break;
			
			default:
                // top face n = (0,0,1)
                c1.x = c2.x = pos.x-r1x;
                c3.x = c4.x = pos.x+r1x;
                c1.y = c3.y = pos.y-r2y;
                c2.y = c4.y = pos.y+r2y;
				c1.z = c2.z = c3.z = c4.z = pos.z+r3z;
                faceWt = r2y*r1x;
                break;
       }
        
        // get elements
        try
        {	cElem[0] = mpmgrid.FindElementFromPoint(&c1)-1;
            theElements[cElem[0]]->GetXiPos(&c1,&corners[0]);
            
            cElem[1] = mpmgrid.FindElementFromPoint(&c2)-1;
            theElements[cElem[1]]->GetXiPos(&c2,&corners[1]);
			
            cElem[2] = mpmgrid.FindElementFromPoint(&c3)-1;
            theElements[cElem[2]]->GetXiPos(&c3,&corners[2]);
			
            cElem[3] = mpmgrid.FindElementFromPoint(&c4)-1;
            theElements[cElem[3]]->GetXiPos(&c4,&corners[3]);
        }
        catch(...)
        {	throw CommonException("A Traction edge node has left the grid.","MatPoint2D::GetTractionInfo");
        }
    }
    else
    {   // get deformed corners, but get element and natural coordinates
        //  from CPDI info because corners have moved by here for any
        //  simulations that update strains between initial extrapolation
        //  and the grid forces calculation
        int d1,d2,d3,d4;
        switch(face)
        {	case 1:
                // lower face n = (0,-1,0)
                d1=0;
                d2=1;
				d3=4;
				d4=5;
                break;
                
            case 2:
                // right face n = (1,0,0)
                d1=1;
                d2=2;
				d3=5;
				d4=6;
                break;
                
            case 3:
                // top face n = (0,1,0)
                d1=2;
                d2=3;
				d3=6;
				d4=7;
                break;
                
            case 4:
                // left face n = (-1,0,0)
                d1=3;
                d2=0;
				d3=7;
				d4=4;
                break;
			
			case 5:
				// bottom face n = (0,0,-1)
                d1=0;
                d2=1;
				d3=2;
				d4=3;
                break;
			
			default:
				// top face n = (0,0,1)
                d1=4;
                d2=5;
				d3=6;
				d4=7;
                break;
        }
        
        // copy for initial state at start of time step
        cElem[0] = cpdi[d1]->inElem;
		CopyVector(&corners[0], &cpdi[d1]->ncpos);
        cElem[1] = cpdi[d2]->inElem;
		CopyVector(&corners[1], &cpdi[d2]->ncpos);
        cElem[2] = cpdi[d3]->inElem;
		CopyVector(&corners[2], &cpdi[d3]->ncpos);
        cElem[3] = cpdi[d4]->inElem;
		CopyVector(&corners[3], &cpdi[d4]->ncpos);
        
        // get weighting factor as 1/4 of face area
        // 1/4 is the get average of the four nodes
        if(face==1 || face==3)
            faceWt = faceArea->x;
        else if(face==2 || face==4)
            faceWt = faceArea->y;
		else
			faceWt = faceArea->z;
        
    }
	
    // get traction normal vector
    ZeroVector(tscaled);
	switch(dof)
	{	case 1:
			// normal is x direction
			tscaled->x = faceWt;
			break;
        case 2:
            // normal is y direction
            tscaled->y = faceWt;
            break;
		default:
			// normal is z direction (not used here)
            tscaled->z = faceWt;
			break;
	}
	
	return 1.;
}

