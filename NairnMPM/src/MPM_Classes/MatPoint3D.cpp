/********************************************************************************
    MatPoint3D.cpp
    nairn-mpm-fea
    
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
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"

// locations of the 8 corners relative to the center
// when adding semi-side vectors
static double r1s[8]={-1.,1.,1.,-1.,-1.,1.,1.,-1.};
static double r2s[8]={-1.,-1.,1.,1.,-1.,-1.,1.,1.};
static double r3s[8]={-1.,-1.,-1.,-1.,1.,1.,1.,1.};

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
	
	// don't need to zero zDeriv because always set in 3D elements
    
	// find shape functions and derviatives
	const ElementBase *elemRef = theElements[ElemID()];
	elemRef->GetShapeGradients(&numnds,fn,nds,xDeriv,yDeriv,zDeriv,this);
    
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

	// pass on to material class to handle
	PerformConstitutiveLaw(dv,strainTime,np,props,&res);
}

// Pass on to material class
void MatPoint3D::PerformConstitutiveLaw(Matrix3 dv,double strainTime,int np,void *props,ResidualStrains *res)
{
    // update particle strain and stress using its constitutive law
	const MaterialBase *matRef = theMaterials[MatID()];
    matRef->MPMConstitutiveLaw(this,dv,strainTime,np,props,res,0);
}

// Move position (3D) (in mm) possibly with particle damping and accWt = dt/2
// vstar is velocity extrapolated from grid to particle and end of time step with grid damping
//        vstar = v(g->p)(n+1) - alpha_g(t)*v(g->p)(n)*dt
// acc is acceleration with grid damping
//        acc = a(g->p)(n) - alpha_g(t)*v(g->p)(n)
// Update with addition of particle damping is
//        dx = (vstar - alpha_g(t)*vp(n)*dt)*dt - (1/2)(acc-alpha_g(t)*vp(n))*dt^2
//           = (vstar - (dt/2)*(acc+alpha_g(t)*vp(n)))*dt
// Must be called BEFORE velocity update, because is needs vp(n) at start of timestep
void MatPoint3D::MovePosition(double delTime,Vector *vstar,double accWt,double particleAlpha)
{	
	double dx = delTime*(vstar->x - accWt*(acc.x + particleAlpha*vel.x));
	double dy = delTime*(vstar->y - accWt*(acc.y + particleAlpha*vel.y));
	double dz = delTime*(vstar->z - accWt*(acc.z + particleAlpha*vel.z));
	pos.x += dx;
    pos.y += dy;
    pos.z += dz;
}

// Move velocity (3D) (in mm/sec) possibly with particle damping
// acc is acceleration with grid damping
//        acc = a(g->p)(n) - alpha_g(t)*v(g->p)(n)
// Update with addition of particle damping is
//        v = vp + (acc-alpha_g(t)*vp(n))*dt
//          = vp*(1-alpha_g(t)*dt) + acc*dt
void MatPoint3D::MoveVelocity(double delTime,double particleAlpha)
{
	vel.x = vel.x*(1. - particleAlpha*delTime) + delTime*acc.x;
    vel.y = vel.y*(1. - particleAlpha*delTime) + delTime*acc.y;
    vel.z = vel.z*(1. - particleAlpha*delTime) + delTime*acc.z;
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

// calculate internal force as -mp sigma.deriv
// add external force (times a shape function)
// store in buffer
// (note: stress is specific stress in units N/m^2 mm^3/g which is (g-mm^2/sec^2)/g
void MatPoint3D::GetFintPlusFext(Vector *theFrc,double fni,double xDeriv,double yDeriv,double zDeriv)
{	
	theFrc->x = -mp*((sp.xx-pressure)*xDeriv+sp.xy*yDeriv+sp.xz*zDeriv) + fni*pFext.x;
	theFrc->y = -mp*(sp.xy*xDeriv+(sp.yy-pressure)*yDeriv+sp.yz*zDeriv) + fni*pFext.y;
	theFrc->z = -mp*(sp.xz*xDeriv+sp.yz*yDeriv+(sp.zz-pressure)*zDeriv) + fni*pFext.z;
}

// add to the concentration gradient (non-rigid particles only)
void MatPoint3D::AddTemperatureGradient(int offset,Vector *grad)
{	pTemp[offset]+=grad->x;
    pTemp[offset+1]+=grad->y;
    pTemp[offset+2]+=grad->z;
}

// return conduction force = - mp (Vp/V0) [k/rho0] Grad T . Grad S (units nJ/sec)
// and k/rho0 is stored in k in units (nJ mm^2/(sec-K-g))
//  (non-rigid particles only)
double MatPoint3D::FCond(int offset,double dshdx,double dshdy,double dshdz,TransportProperties *t)
{
	Tensor *kten = &(t->kCondTensor);
	return -mp*GetRelativeVolume()*((kten->xx*pTemp[offset] + kten->xy*pTemp[offset+1] + kten->xz*pTemp[offset+2])*dshdx
						+ (kten->xy*pTemp[offset] + kten->yy*pTemp[offset+1] + kten->yz*pTemp[offset+2])*dshdy
						+ (kten->xz*pTemp[offset] + kten->yz*pTemp[offset+1] + kten->zz*pTemp[offset+2])*dshdz);
}

// add to the concentration gradient
void MatPoint3D::AddConcentrationGradient(Vector *grad)
{	pDiffusion[gGRADx]+=grad->x;
    pDiffusion[gGRADy]+=grad->y;
    pDiffusion[gGRADz]+=grad->z;
}

// return diffusion force = - V [D] Grad C . Grad S
double MatPoint3D::FDiff(double dshdx,double dshdy,double dshdz,TransportProperties *t)
{
	Tensor *Dten = &(t->diffusionTensor);
	return -GetVolume(DEFORMED_VOLUME)*((Dten->xx*pDiffusion[gGRADx] + Dten->xy*pDiffusion[gGRADy] + Dten->xz*pDiffusion[gGRADz])*dshdx
						+ (Dten->xy*pDiffusion[gGRADx] + Dten->yy*pDiffusion[gGRADy] + Dten->yz*pDiffusion[gGRADz])*dshdy
						+ (Dten->xz*pDiffusion[gGRADx] + Dten->yz*pDiffusion[gGRADy] + Dten->zz*pDiffusion[gGRADz])*dshdz);
}

// return kinetic energy (g mm^2/sec^2) = nanoJ
double MatPoint3D::KineticEnergy(void)
{	return 0.5*mp*(vel.x*vel.x+vel.y*vel.y+vel.z*vel.z);
}

// get deformation gradient, which is stored in strain and rotation tensors
Matrix3 MatPoint3D::GetDeformationGradientMatrix(void) const
{	double F[3][3];
	GetDeformationGradient(F);
	Matrix3 Fm(F[0][0],F[0][1],F[0][2],F[1][0],F[1][1],F[1][2],F[2][0],F[2][1],F[2][2]);
	return Fm;
}

//sset deformation gradient, which is stored in strain and rotation tensors
void MatPoint3D::SetDeformationGradientMatrix(Matrix3 F)
{	
	// Normal strains
	ep.xx = F(0,0) - 1.;
	ep.yy = F(1,1) - 1.;
	ep.zz = F(2,2) - 1.;
	
	// Shear Strains
	ep.xy = F(1,0) + F(0,1);
	ep.xz = F(2,0) + F(0,2);
	ep.yz = F(2,1) + F(1,2);
	
	// Rotation Strains
	wrot.xy = F(1,0) - F(0,1);
	wrot.xz = F(2,0) - F(0,2);
	wrot.yz = F(2,1) - F(1,2);
}

// get displacement gradient from grad u = F - I
Matrix3 MatPoint3D::GetDisplacementGradientMatrix(void) const
{	double F[3][3];
	GetDeformationGradient(F);
	Matrix3 Fm(F[0][0],F[0][1],F[0][2],F[1][0],F[1][1],F[1][2],F[2][0],F[2][1],F[2][2]);
	
	Fm(0,0) -= 1.;
	Fm(1,1) -= 1.;
	Fm(2,2) -= 1.;
	
	return Fm;
}

// get the symmetric elastic Left-Cauchy tensor in a Matrix3
Matrix3 MatPoint3D::GetElasticLeftCauchyMatrix(void)
{   return Matrix3(eplast.xx,eplast.xy,eplast.xz,eplast.xy,eplast.yy,eplast.yz,
                    eplast.xz,eplast.yz,eplast.zz);
}

// get deformation gradient, which is stored in strain and rotation tensors
void MatPoint3D::GetDeformationGradient(double F[][3]) const
{
	// current deformation gradient in 3D
	F[0][0] = 1. + ep.xx;
	F[1][1] = 1. + ep.yy;
	F[2][2] = 1. + ep.zz;
	F[0][1] = 0.5*(ep.xy - wrot.xy);
	F[1][0] = 0.5*(ep.xy + wrot.xy);
	F[0][2] = 0.5*(ep.xz - wrot.xz);
	F[2][0] = 0.5*(ep.xz + wrot.xz);
	F[1][2] = 0.5*(ep.yz - wrot.yz);
	F[2][1] = 0.5*(ep.yz + wrot.yz);
}

// Get R.sqrt(B).RT-I = V-I for elastic biot strain rotated into current particle orientation
// Assume that elastic B matrix is in the alt strain tensor
Matrix3 MatPoint3D::GetElasticBiotStrain(void)
{	// Get Sqrt(B)
	Tensor *B = GetAltStrainTensor();
	Matrix3 Be = Matrix3(B->xx,B->xy,B->xz,B->xy,B->yy,B->yz,B->xy,B->xz,B->zz);
	Vector lam;
	Matrix3 Q = Be.Eigenvectors(lam);
	Matrix3 LamI = Matrix3(sqrt(lam.x),0.,0.,sqrt(lam.y),sqrt(lam.z));
	Matrix3 U = LamI.RMRT(Q);
	
	// Get rotation matrix
	Matrix3 F = GetDeformationGradientMatrix();
	Matrix3 R;
	Matrix3 Utot = F.RightDecompose(&R,NULL);
	
	// Rotate U to V, subtract I and return the result
	Matrix3 V = U.RMRT(R);
	V(0,0) -= 1.;
	V(1,1) -= 1.;
	V(2,2) -= 1.;
	return V;
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
// only used for crack contact, multimaterial contact, and transport tasks
double MatPoint3D::GetVolume(int volumeType)
{	double rho=GetRho();						// in g/mm^3
	return GetRelativeVolume()*mp/rho;			// in mm^3
}

// Get vectors from particle to edge
// r1 = F.(psz.x,0,0), r2 = F.(0,psz.y,0), and r3 = F.(0,0,psz.z)
// where psz is particle size
void MatPoint3D::GetSemiSideVectors(Vector *r1,Vector *r2,Vector *r3) const
{
	// get particle 2D deformation gradient
	double pF[3][3];
	GetDeformationGradient(pF);
	
	// get polygon vectors - these are from particle to edge
    //      and generalize semi width lp in 1D GIMP
    Vector psz = mpmgrid.GetParticleSize();
	r1->x = pF[0][0]*psz.x;
	r1->y = pF[1][0]*psz.x;
	r1->z = pF[2][0]*psz.x;
    
	r2->x = pF[0][1]*psz.y;
	r2->y = pF[1][1]*psz.y;
	r2->z = pF[2][1]*psz.y;
    
	r3->x = pF[0][2]*psz.z;
	r3->y = pF[1][2]*psz.z;
	r3->z = pF[2][2]*psz.z;
}

// Get undeformed size (only called by GIMP traction)
// This uses mpmgrid so membrane particles must override
void MatPoint3D::GetUndeformedSemiSides(double *r1x,double *r2y,double *r3z) const
{   Vector psz = mpmgrid.GetParticleSize();
    *r1x = psz.x;
    *r2y = psz.y;
    *r3z = psz.z;
}

// to support CPDI return nodes for corners (or for 9 nodes) and weights for shape functions
//	and shape function gradients
// throws CommonException() if particle corner has left the grid
void MatPoint3D::GetCPDINodesAndWeights(int cpdiType)
{   
	// get polygon vectors - these are from particle to edge
    //      and generalize semi width lp in 1D GIMP
	Vector c,r1,r2,r3;
    GetSemiSideVectors(&r1,&r2,&r3);
	
    // Particle domain volume is 8 * volume of the parallelepiped defined by r1, r2, and r3
	// V = 8 * (r1 . (r2 X r3))
    // Assume positive due to orientation of initial vectors, and sign probably does not matter
    double Vp = 8.*(r1.x * (r2.y*r3.z - r2.z*r3.y)
                    + r1.y * (r2.z*r3.x - r2.x*r3.z)
                     + r1.z * (r2.x*r3.y - r2.y*r3.x) );
    
	// always LINEAR_CPDI
    
	try
	{	// find all 8 corner nodes
        for(int i=0;i<8;i++)
        {   c.x = pos.x + r1s[i]*r1.x + r2s[i]*r2.x + r3s[i]*r3.x;
            c.y = pos.y + r1s[i]*r1.y + r2s[i]*r2.y + r3s[i]*r3.y;
            c.z = pos.z + r1s[i]*r1.z + r2s[i]*r2.z + r3s[i]*r3.z;
            cpdi[i]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
            theElements[cpdi[i]->inElem]->GetXiPos(&c,&cpdi[i]->ncpos);
        }

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
    catch(CommonException err)
    {   char msg[200];
        sprintf(msg,"A CPDI particle domain node has left the grid: %s",err.Message());
        throw CommonException(msg,"MatPoint3D::GetCPDINodesAndWeights");
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
	Vector e1,e2;
	
	// always UNIFORM_GIMP or LINEAR_CPDI
    
    // which GIMP method (cannot be used in POINT_GIMP)
    if(ElementBase::useGimp==UNIFORM_GIMP)
    {   // initial vectors only
        double r1x,r2y,r3z;
        GetUndeformedSemiSides(&r1x,&r2y,&r3z);
        
		// edges are c1 to c2 to c4 to c3
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
        catch(CommonException err)
        {   char msg[200];
            sprintf(msg,"A Traction edge node has left the grid: %s",err.Message());
            throw CommonException(msg,"MatPoint3D::GetTractionInfo");
        }
		
		if(dof==N_DIRECTION)
		{	e1 = c2;
			SubVector(&e1,&c1);
			e2 = c3;
			SubVector(&e2,&c1);
		}
    }
    else
    {   // get deformed corners, but get element and natural coordinates
        //  from CPDI info because corners have moved by here for any
        //  simulations that update strains between initial extrapolation
        //  and the grid forces calculation
		
		// edges are d1 to d2 to d4 to d3
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
        
		// get edge vector
		if(dof==N_DIRECTION)
		{	theElements[cElem[1]]->GetPosition(&corners[1],&e1);
			theElements[cElem[2]]->GetPosition(&corners[2],&e2);
			Vector e0;
			theElements[cElem[0]]->GetPosition(&corners[0],&e0);
			SubVector(&e1,&e0);
			SubVector(&e2,&e0);
		}
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
		case N_DIRECTION:
		{	Vector cp;
			CrossProduct(&cp,&e1,&e2);
			double enorm = sqrt(DotVectors(&cp,&cp));
			tscaled->x = cp.x*faceWt/enorm;
			tscaled->y = cp.y*faceWt/enorm;
			tscaled->z = cp.z*faceWt/enorm;
			break;
		}
		default:
			// normal is z direction (not used here)
            tscaled->z = faceWt;
			break;
	}
	
	// always 1 in 3D (used in AS)
	return 1.;
}

// Get Rotation matrix for initial material orientation (anistropic only)
Matrix3 MatPoint3D::GetInitialRotation(void)
{	double z = GetAnglez0InRadians();
	double y = GetAngley0InRadians();
	double x = GetAnglex0InRadians();
	double cx = cos(x);
	double sx = sin(x);
	double cy = cos(y);
	double sy = sin(y);
	double cz = cos(z);
	double sz = sin(z);
	return Matrix3( cy*cz,  cz*sx*sy+cx*sz,  -cx*cz*sy+sx*sz,
				   -cy*sz,  cx*cz-sx*sy*sz,   cz*sx+cx*sy*sz,
				   sy,     -cy*sx,            cx*cy );
	//return Matrix3( cy*cz,         -cy*sz,           sy,
	//			    cz*sx*sy+cx*sz, cx*cz-sx*sy*sz, -cy*sx,
	//			   -cx*cz*sy+sx*sz, cz*sx+cx*sy*sz,  cx*cy );
}

