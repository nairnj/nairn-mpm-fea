/********************************************************************************
    MatPoint3D.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 7/21/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/NairnMPM.hpp"
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
// postUpdate true means after particle position updated
void MatPoint3D::UpdateStrain(double strainTime,int secondPass,int np,void *props,int matFld,bool postUpdate)
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

	// don't need to zero zDeriv because always set in 3D elements
	
	// find shape functions and derviatives
	const ElementBase *elemRef = theElements[ElemID()];
	int *nds = ndsArray;
	elemRef->GetShapeGradients(fn,&nds,xDeriv,yDeriv,zDeriv,this);
	numnds = nds[0];
	
    // Find strain rates at particle from current grid velocities
	//   and using the velocity field for that particle with each node
    for(i=1;i<=numnds;i++)
	{	vel=nd[nds[i]]->GetVelocity((short)vfld[i],matFld);
		Vector grad = MakeVector(xDeriv[i],yDeriv[i],zDeriv[i]);
		dv += Matrix3(&vel,&grad);
    }
	
	// convert to strain increments
    dv.Scale(strainTime);
    
	// pass on to material class to handle
	ResidualStrains res = ScaledResidualStrains(secondPass);
	PerformConstitutiveLaw(dv,strainTime,np,props,&res,gStressPtr);
}

// Pass on to material class
void MatPoint3D::PerformConstitutiveLaw(Matrix3 dv,double strainTime,int np,void *props,ResidualStrains *res,Tensor *gStress)
{
    // update particle strain and stress using its constitutive law
	const MaterialBase *matRef = theMaterials[MatID()];
    matRef->MPMConstitutiveLaw(this,dv,strainTime,np,props,res,0,gStress);
}

// Move position and velocity (2D)
void MatPoint3D::MoveParticle(GridToParticleExtrap *gp)
{
	// get vm = S(v-a dt) for FLIP and Sv^+ FMPM (and PIC)
	Vector vm;
	if(gp->m>0)
	{	// FMPM damps with Sk^+(k)
		vm = gp->Svtilde;
	}
	else if(gp->m>-2)
	{	// FLIP and XPIC(1) damps with Sv (which is initial lumped velocity)
		vm = gp->Svtilde;
		AddScaledVector(&vm,&gp->Sacc,-timestep);
	}
	else
	{	// XPIC(k>1) damps by initial lumped velocity, and Svtilde holds Sv(k)
		vm = gp->Svlumped;
		AddScaledVector(&vm,&gp->Sacc,-timestep);
	}
	
	// find Adamp0
	Vector Adamp0;
	Adamp0 = SetScaledVector(&vel,gp->particleAlpha);
	AddScaledVector(&Adamp0,&vm,gp->gridAlpha);
	
	Vector delV;
	if(gp->m>0)
	{	// save initial particle velocity
		Vector delXRate = vel;;
		
		// FMPM(m) update
		vel = vm;
		AddScaledVector(&vel,&Adamp0,-timestep);
		
		// find change in velocity
		delV = vel;
		SubVector(&delV,&delXRate);
		
		// finish midpoint position update
		// 2*del X = (V(n+1)+V(n))
        AddVector(&delXRate,&vel);
		AddScaledVector(&pos,&delXRate,0.5*timestep);
	}
	else if(gp->m==0)
	{	// FLIP update velocity change
		delV.x = (gp->Sacc.x - Adamp0.x)*timestep;
		delV.y = (gp->Sacc.y - Adamp0.y)*timestep;
		delV.z = (gp->Sacc.z - Adamp0.z)*timestep;
		
		// velocity
		AddVector(&vel,&delV);
		
		// position update
		Vector delXRate;
		delXRate.x = vm.x + 0.5*delV.x;
		delXRate.y = vm.y + 0.5*delV.y;
		delXRate.z = vm.z + 0.5*delV.z;
		AddScaledVector(&pos,&delXRate,timestep);
	}
	else
	{	// XPIC update
		
		// save initial particle velocity
		Vector delXRate = vel;
		
		// XPIC(k) velocity update
		// For XPIC(1) Svtilde has Sv+, for XPIC(k>1) Svtilde has S(v(k)+a*dt)
		vel.x = gp->Svtilde.x - Adamp0.x*timestep;
		vel.y = gp->Svtilde.y - Adamp0.y*timestep;
		vel.z = gp->Svtilde.z - Adamp0.z*timestep;

		// find change in velocity
		delV.x = vel.x - delXRate.x;
		delV.y = vel.y - delXRate.y;
		delV.z = vel.z - delXRate.z;

		// position update
		delXRate.x = vm.x + 0.5*delV.x;
		delXRate.y = vm.y + 0.5*delV.y;
		delXRate.z = vm.z + 0.5*delV.z;
		AddScaledVector(&pos,&delXRate,timestep);
	}

	// J Integral needs effective particle acceleration
	acc = MakeVector(delV.x/timestep,delV.y/timestep,delV.z/timestep);
		
#ifdef CHECK_NAN
	if(pos.x!=pos.x || pos.y!=pos.y || pos.z!=pos.z || vel.x!=vel.x || vel.y!=vel.y || vel.z!=vel.z)
	{
#pragma omp critical (output)
		{	cout << "\n# MatPoint2D::MoveParticle: bad pos or vel" << endl;
			Describe();
		}
	}
#endif
}

// Move rigid particle by new current velocity
void MatPoint3D::MovePosition(double delTime)
{	pos.x += delTime*vel.x;
    pos.y += delTime*vel.y;
	pos.z += delTime*vel.z;
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

// particle semi size in actual units
Vector MatPoint3D::GetParticleSize(void) const
{	Vector part = theElements[ElemID()]->GetDeltaBox();
	part.x *= 0.5*mpm_lp.x;
	part.y *= 0.5*mpm_lp.y;
	part.z *= 0.5*mpm_lp.z;
	return part;
}

// get minimum particle side
double MatPoint3D::GetMinParticleLength(void) const
{	Vector part = GetParticleSize();
	double minPart = part.x;
	if (part.y < minPart) minPart = part.y;
	if (part.z < minPart) minPart = part.z;
	return 2.*minPart;
}

// calculate internal force as -mp sigma.deriv
// add external force (times a shape function)
// store in buffer
// (note: stress is specific stress in units N/m^2 mm^3/g which is (g-mm^2/sec^2)/g
void MatPoint3D::GetFintPlusFext(Vector *theFrc,double fni,double xDeriv,double yDeriv,double zDeriv)
{	
	theFrc->x = -mp*((sp.xx-pressure)*xDeriv+sp.xy*yDeriv+sp.xz*zDeriv) + fni*pFext.x;
	theFrc->y = -mp*(sp.xy*xDeriv+(sp.yy-pressure)*yDeriv+sp.yz*zDeriv) + fni*pFext.y;
	theFrc->z = -mp*(sp.xz*xDeriv+sp.yz*yDeriv+(sp.zz-pressure)*zDeriv) + fni*pFext.z;
	
#ifdef CHECK_NAN
	if(theFrc->x!=theFrc->x || theFrc->y!=theFrc->y || theFrc->z!=theFrc->z)
	{
#pragma omp critical (output)
		{	cout << "\n# MatPoint3D::GetFintPlusFext: bad nodal fint+fext";
			PrintVector(" = ",theFrc);
			cout << endl;
			cout << "# sp = (" << sp.xx << "," << sp.yy << "," << sp.zz << "," << sp.yz
				<< "," << sp.xz << "," << sp.xy << ") P = " << pressure << endl;
			PrintVector("# pFext = ",&pFext);
			cout << endl;
			cout << "# Sip = " << fni << " Gip = (" << xDeriv << "," << yDeriv << "," << zDeriv << ")" << endl;
		}
	}
#endif
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
void MatPoint3D::AddConcentrationGradient(int dnum,Vector *grad)
{	pDiff[dnum]->grad.x += grad->x;
	pDiff[dnum]->grad.y += grad->y;
	pDiff[dnum]->grad.z += grad->z;
}

// return diffusion force = - V [D] Grad C . Grad S
double MatPoint3D::FDiff(double dshdx,double dshdy,double dshdz,TransportProperties *t,int number)
{
	Tensor *Dten = &(t->diffusionTensor);
    Vector pgrad = pDiff[number]->grad;
	return -GetVolume(DEFORMED_VOLUME)*((Dten->xx*pgrad.x + Dten->xy*pgrad.y + Dten->xz*pgrad.z)*dshdx
                                      + (Dten->xy*pgrad.x + Dten->yy*pgrad.y + Dten->yz*pgrad.z)*dshdy
                                      + (Dten->xz*pgrad.x + Dten->yz*pgrad.y + Dten->zz*pgrad.z)*dshdz);
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

// get displacement gradient in current configuration
// Uses grad u = RFR^T - I = RV - I
Matrix3 MatPoint3D::GetDisplacementGradientMatrix(void) const
{	Matrix3 F = GetDeformationGradientMatrix();
	Matrix3 Rn;
	Matrix3 V = F.LeftDecompose(&Rn,NULL);
	V = Rn*V;
	V(0,0) -= 1.;
	V(1,1) -= 1.;
	V(2,2) -= 1.;
	return V;
}

// get displacement gradient from grad u = F - I
Matrix3 MatPoint3D::GetDisplacementGradientForJ(const MaterialBase *matref)
{
	throw "Displacement gradient for 3D J calculations not available";
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

// Get sqrt(B)-I = V-I for elastic biot strain rotated into current particle orientation
// Assume that elastic B matrix is in the alt strain tensor
Matrix3 MatPoint3D::GetElasticBiotStrain(void)
{	// Get Sqrt(B)
	Tensor *B = GetAltStrainTensor();
	Matrix3 Be = Matrix3(B->xx,B->xy,B->xz,B->xy,B->yy,B->yz,B->xy,B->xz,B->zz);
	Vector lam;
	Matrix3 Q = Be.Eigenvectors(lam);
	Matrix3 LamI = Matrix3(sqrt(lam.x),0.,0.,sqrt(lam.y),sqrt(lam.z));
	Matrix3 V = LamI.RMRT(Q);
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
    Vector psz = GetParticleSize();
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

// rescale vectors if needed in CPDI calcualtions
void MatPoint3D::ScaleSemiSideVectorsForCPDI(Vector *r1,Vector *r2,Vector *r3) const
{	// skip if not activated
	if(ElementBase::rcrit<0.) return;
	
	// convert rcrit to dimensions based smallest dimension this element
	Vector cellSize = theElements[ElemID()]->GetDeltaBox();
	double rcrit = ElementBase::rcrit*fmin(cellSize.x,fmin(cellSize.y,cellSize.z));
	
	// check r1+r2+r3
	bool rescale = false;
	Vector la = MakeVector(r1->x+r2->x+r3->x, r1->y+r2->y+r3->y, r1->z+r2->z+r3->z);
	double lamag = sqrt(la.x*la.x+la.y*la.y+la.z*la.z);
	if(lamag>rcrit)
	{	ScaleVector(&la,rcrit/lamag);
		rescale = true;
	}
	
	// check r1-r2+r3
	Vector lb = MakeVector(r1->x-r2->x+r3->x, r1->y-r2->y+r3->y, r1->z-r2->z+r3->z);
	double lbmag = sqrt(lb.x*lb.x+lb.y*lb.y+lb.z*lb.z);
	if(lbmag>rcrit)
	{	ScaleVector(&lb,rcrit/lbmag);
		rescale = true;
	}
	
	// check -r1+r2+r3
	Vector lc = MakeVector(-r1->x+r2->x+r3->x, -r1->y+r2->y+r3->y, -r1->z+r2->z+r3->z);
	double lcmag = sqrt(lc.x*lc.x+lc.y*lc.y+lc.z*lc.z);
	if(lcmag>rcrit)
	{	ScaleVector(&lc,rcrit/lcmag);
		rescale = true;
	}
	
	// check -r1-r2+r3
	Vector ld = MakeVector(-r1->x-r2->x+r3->x, -r1->y-r2->y+r3->y, -r1->z-r2->z+r3->z);
	double ldmag = sqrt(ld.x*ld.x+ld.y*ld.y+ld.z*ld.z);
	if(ldmag>rcrit)
	{	ScaleVector(&ld,rcrit/ldmag);
		rescale = true;
	}
	
	// exit if no scale
	if(!rescale) return;
	
	// redefine r1 and r2 and r3
	r1->x = 0.25*(la.x+lb.x-lc.x-ld.x);
	r1->y = 0.25*(la.y+lb.y-lc.y-ld.y);
	r1->z = 0.25*(la.z+lb.z-lc.z-ld.z);
	
	r2->x = 0.25*(la.x-lb.x+lc.x-ld.x);
	r2->y = 0.25*(la.y-lb.y+lc.y-ld.y);
	r2->z = 0.25*(la.z-lb.z+lc.z-ld.z);
	
	r3->x = 0.25*(la.x+lb.x+lc.x+ld.x);
	r3->y = 0.25*(la.y+lb.y+lc.y+ld.y);
	r3->z = 0.25*(la.z+lb.z+lc.z+ld.z);
}

#define ELLIPSE_METHOD_3D
// Get distance from particle center to surface of deformed particle
// along the provided unit vector
double MatPoint3D::GetDeformedRadius(Vector *norm) const
{
#ifdef ELLIPSE_METHOD_3D
	// unnormalized normal in undeformed config (n0x,n0y) = F^{-1}n
	Matrix3 F = GetDeformationGradientMatrix();
	Matrix3 Finv = F.Inverse();
	Vector n0 = Finv.Times(norm);
	
	// return radius to deformed inscribed ellipsoid
	// see JANOSU-014-82
	Vector psz = GetParticleSize();
	double rd = sqrt(n0.x*n0.x/(psz.x*psz.x) + n0.y*n0.y/(psz.y*psz.y) + n0.z*n0.z/(psz.z*psz.z));
	return 1./rd;
#else
	// See JANOSU-13-74 to 77
	Vector a,b,c;
	GetSemiSideVectors(&a,&b,&c);
	Vector aXb,aXc,bXc;
	CrossProduct(&aXb,&a,&b);
	CrossProduct(&aXc,&a,&c);
	CrossProduct(&bXc,&b,&c);
	double nab = fabs(DotVectors(norm,&aXb));
	double nac = fabs(DotVectors(norm,&aXc));
	double nbc = fabs(DotVectors(norm,&bXc));
	return fabs(DotVectors(&a,&bXc))/fmax(nab,fmax(nac,nbc));
#endif
}

// Get undeformed size (only called by GIMP traction and Custom thermal ramp)
// This uses mpmgrid so membrane particles must override
// Assumes undeformation particle aligned with x-y-z axes
void MatPoint3D::GetUndeformedSemiSides(double *r1x,double *r2y,double *r3z) const
{   Vector psz = GetParticleSize();
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
	
	// rescale if needed
	ScaleSemiSideVectorsForCPDI(&r1,&r2,&r3);
	
    // Particle domain volume is 8 * volume of the parallelepiped defined by r1, r2, and r3
	// V = 8 * (r1 . (r2 X r3))
    // Assume positive due to orientation of initial vectors, and sign probably does not matter
    double Vp = 8.*(r1.x * (r2.y*r3.z - r2.z*r3.y)
                    + r1.y * (r2.z*r3.x - r2.x*r3.z)
                     + r1.z * (r2.x*r3.y - r2.y*r3.x) );
    
	// always LINEAR_CPDI or BSPLINE_CPDI
    
	try
	{	CPDIDomain **cpdi = GetCPDIInfo();
		
		// find all 8 corner nodes
        for(int i=0;i<8;i++)
        {   c.x = pos.x + r1s[i]*r1.x + r2s[i]*r2.x + r3s[i]*r3.x;
            c.y = pos.y + r1s[i]*r1.y + r2s[i]*r2.y + r3s[i]*r3.y;
            c.z = pos.z + r1s[i]*r1.z + r2s[i]*r2.z + r3s[i]*r3.z;
            cpdi[i]->inElem = mpmgrid.FindElementFromPoint(&c,this)-1;
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
    catch(CommonException& err)
    {   char msg[200];
		size_t msgSize=200;
		snprintf(msg,msgSize,"A CPDI particle domain node has left the grid: %s",err.Message());
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

// Given face, find dimensionless coordinates of corners (4 in 3D), elements for those corners (4 in 3D),
//     and deformed particle radii (3 in 3D)
// Return number of corners in numDnds
// Function result is face weight (which is Area/4 in 3D)
// On input cElem, corners, and radii are arrays of length equal to number of corners (4 in 3D)
Vector MatPoint3D::GetSurfaceInfo(int face,int dof,int *cElem,Vector *corners,Vector *radii,int *numDnds,double *redge)
{
	*numDnds = 4;
	Vector c1, c2, c3, c4, c5;		// four corners of face, other direction node 1
	Vector r1,r2,r3;
	
	// get polygon vectors - these are from particle to edge
	GetSemiSideVectors(&r1,&r2,&r3);
	double uGIMPsize = -1.;
#ifndef TRACTION_ALWAYS_DEFORMED
    // this option uses undeformed edge for all GIMP except finite GIMP
    if(!ElementBase::UsingCPDIMethod() && ElementBase::useGimp!=FINITE_GIMP)
    {	double rx,ry,rz;
        
        // get deformed area
        Vector Ap;
        switch(face)
        {	case 1:
            case 3:
                CrossProduct(&Ap,&r1,&r3);
                break;
            case 4:
            case 2:
                CrossProduct(&Ap,&r2,&r3);
                break;
            case 5:
            default:
                CrossProduct(&Ap,&r1,&r2);
                break;

        }
        uGIMPsize = sqrt(DotVectors(&Ap,&Ap));
        
        // for extrapolations, however, switch to uGIMP nodes
        GetUndeformedSemiSides(&rx,&ry,&rz);
        ZeroVector(&r1);
        ZeroVector(&r2);
        ZeroVector(&r3);
        r1.x = rx;
        r2.y = ry;
        r3.z = rz;
	}
#endif
	
	switch(face)
	{	case 1:
			// lower face n = (0,-1,0) nodes 1,2,6,5
			c1.x = pos.x-r1.x-r2.x-r3.x;	// 1
			c1.y = pos.y-r1.y-r2.y-r3.y;
			c1.z = pos.z-r1.z-r2.z-r3.z;
			c2.x = pos.x+r1.x-r2.x-r3.x;	// 2
			c2.y = pos.y+r1.y-r2.y-r3.y;
			c2.z = pos.z+r1.z-r2.z-r3.z;
			c3.x = pos.x+r1.x-r2.x+r3.x;	// 6
			c3.y = pos.y+r1.y-r2.y+r3.y;
			c3.z = pos.z+r1.z-r2.z+r3.z;
			c4.x = pos.x-r1.x-r2.x+r3.x;	// 5
			c4.y = pos.y-r1.y-r2.y+r3.y;
			c4.z = pos.z-r1.z-r2.z+r3.z;
			// fifth node 4
			c5.x = pos.x-r1.x+r2.x-r3.x;	// 4
			c5.y = pos.y-r1.y+r2.y-r3.y;
			c5.y = pos.z-r1.z+r2.z-r3.z;
			break;
			
		case 2:
			// right face n = (1,0,0) nodes 2,3,7,6
			c1.x = pos.x+r1.x-r2.x-r3.x;	// 2
			c1.y = pos.y+r1.y-r2.y-r3.y;
			c1.z = pos.z+r1.z-r2.z-r3.z;
			c2.x = pos.x+r1.x+r2.x-r3.x;	// 3
			c2.y = pos.y+r1.y+r2.y-r3.y;
			c2.z = pos.z+r1.z+r2.z-r3.z;
			c3.x = pos.x+r1.x+r2.x+r3.x;	// 7
			c3.y = pos.y+r1.y+r2.y+r3.y;
			c3.z = pos.z+r1.z+r2.z+r3.z;
			c4.x = pos.x+r1.x-r2.x+r3.x;	// 6
			c4.y = pos.y+r1.y-r2.y+r3.y;
			c4.z = pos.z+r1.z-r2.z+r3.z;
			// fifth node 1
			c5.x = pos.x-r1.x-r2.x-r3.x;	// 1
			c5.y = pos.y-r1.y-r2.y-r3.y;
			c5.z = pos.z-r1.z-r2.z-r3.z;
			break;
			
		case 3:
			// top face n = (0,1,0) nodes 3,4,8,7
			c1.x = pos.x+r1.x+r2.x-r3.x;	// 3
			c1.y = pos.y+r1.y+r2.y-r3.y;
			c1.z = pos.z+r1.z+r2.z-r3.z;
			c2.x = pos.x-r1.x+r2.x-r3.x;	// 4
			c2.y = pos.y-r1.y+r2.y-r3.y;
			c2.z = pos.z-r1.z+r2.z-r3.z;
			c3.x = pos.x-r1.x+r2.x+r3.x;	// 8
			c3.y = pos.y-r1.y+r2.y+r3.y;
			c3.z = pos.z-r1.z+r2.z+r3.z;
			c4.x = pos.x+r1.x+r2.x+r3.x;	// 7
			c4.y = pos.y+r1.y+r2.y+r3.y;
			c4.z = pos.z+r1.z+r2.z+r3.z;
			// fifth node 2
			c5.x = pos.x+r1.x-r2.x-r3.x;	// 2
			c5.y = pos.y+r1.y-r2.y-r3.y;
			c5.z = pos.z+r1.z-r2.z-r3.z;
			break;
			
		case 4:
			// left face n = (-1,0,0) nodes 4,1,5,8
			c1.x = pos.x-r1.x+r2.x-r3.x;	// 4
			c1.y = pos.y-r1.y+r2.y-r3.y;
			c1.z = pos.z-r1.z+r2.z-r3.z;
			c2.x = pos.x-r1.x-r2.x-r3.x;	// 1
			c2.y = pos.y-r1.y-r2.y-r3.y;
			c2.z = pos.z-r1.z-r2.z-r3.z;
			c3.x = pos.x-r1.x-r2.x+r3.x;	// 5
			c3.y = pos.y-r1.y-r2.y+r3.y;
			c3.z = pos.z-r1.z-r2.z+r3.z;
			c4.x = pos.x-r1.x+r2.x+r3.x;	// 8
			c4.y = pos.y-r1.y+r2.y+r3.y;
			c4.z = pos.z-r1.z+r2.z+r3.z;
			// fifth node 3
			c5.x = pos.x+r1.x+r2.x-r3.x;	// 3
			c5.y = pos.y+r1.y+r2.y-r3.y;
			c5.z = pos.z+r1.z+r2.z-r3.z;
			break;
			
		case 5:
			// bottom face n = (0,0,-1) nodes 2,1,4,3
			c1.x = pos.x+r1.x-r2.x-r3.x;	// 2
			c1.y = pos.y+r1.y-r2.y-r3.y;
			c1.z = pos.z+r1.z-r2.z-r3.z;
			c2.x = pos.x-r1.x-r2.x-r3.x;	// 1
			c2.y = pos.y-r1.y-r2.y-r3.y;
			c2.z = pos.z-r1.z-r2.z-r3.z;
			c3.x = pos.x-r1.x+r2.x-r3.x;	// 4
			c3.y = pos.y-r1.y+r2.y-r3.y;
			c3.z = pos.z-r1.z+r2.z-r3.z;
			c4.x = pos.x+r1.x+r2.x-r3.x;	// 3
			c4.y = pos.y+r1.y+r2.y-r3.y;
			c4.z = pos.z+r1.z+r2.z-r3.z;
			// fifth node 6
			c5.x = pos.x+r1.x-r2.x+r3.x;	// 6
			c5.y = pos.y+r1.y-r2.y+r3.y;
			c5.z = pos.z+r1.z-r2.z+r3.z;
			break;
			
		default:
			// top face n = (0,0,1) nodes 5,6,7,8
			c1.x = pos.x-r1.x-r2.x+r3.x;	// 5
			c1.y = pos.y-r1.y-r2.y+r3.y;
			c1.z = pos.z-r1.z-r2.z+r3.z;
			c2.x = pos.x+r1.x-r2.x+r3.x;	// 6
			c2.y = pos.y+r1.y-r2.y+r3.y;
			c2.z = pos.z+r1.z-r2.z+r3.z;
			c3.x = pos.x+r1.x+r2.x+r3.x;	// 7
			c3.y = pos.y+r1.y+r2.y+r3.y;
			c3.z = pos.z+r1.z+r2.z+r3.z;
			c4.x = pos.x-r1.x+r2.x+r3.x;	// 8
			c4.y = pos.y-r1.y+r2.y+r3.y;
			c4.z = pos.z-r1.z+r2.z+r3.z;
			// fifth node 1
			c5.x = pos.x-r1.x-r2.x-r3.x;	// 1
			c5.y = pos.y-r1.y-r2.y-r3.y;
			c5.z = pos.z-r1.z-r2.z-r3.z;
			break;
	}
	
	// get elements and xipos for the two corners
	try
	{	cElem[0] = mpmgrid.FindElementFromPoint(&c1,this)-1;
		theElements[cElem[0]]->GetXiPos(&c1,&corners[0]);
		
		cElem[1] = mpmgrid.FindElementFromPoint(&c2,this)-1;
		theElements[cElem[1]]->GetXiPos(&c2,&corners[1]);
		
		cElem[2] = mpmgrid.FindElementFromPoint(&c3,this)-1;
		theElements[cElem[2]]->GetXiPos(&c3,&corners[2]);
		
		cElem[3] = mpmgrid.FindElementFromPoint(&c4,this)-1;
		theElements[cElem[3]]->GetXiPos(&c4,&corners[3]);
	}
	catch(CommonException& err)
	{   char msg[200];
		size_t msgSize=200;
		snprintf(msg,msgSize,"A Traction edge node has left the grid: %s",err.Message());
		throw CommonException(msg,"MatPoint3D::GetSurfaceInfo");
	}
	
	// get relevant edge radii for this face
	radii[0] = MakeVector(0.5*(c2.x-c1.x),0.5*(c2.y-c1.y),0.5*(c2.z-c1.z));
	radii[1] = MakeVector(0.5*(c4.x-c1.x),0.5*(c4.y-c1.y),0.5*(c4.z-c1.z));
	radii[2] = MakeVector(0.5*(c5.x-c1.x),0.5*(c5.y-c1.y),0.5*(c5.z-c1.z));
	
	// Face Area/4 = mag of cross product vector
	Vector Ap;
	CrossProduct(&Ap,&radii[0],&radii[1]);
	double faceWt = uGIMPsize<0. ? sqrt(DotVectors(&Ap,&Ap)) : uGIMPsize ;

	// get traction normal vector times 1/4 area
	Vector wtNorm;
	ZeroVector(&wtNorm);
	switch(dof)
	{	case 1:
			// normal is x direction
			wtNorm.x = faceWt;
			break;
		case 2:
			// normal is y direction
			wtNorm.y = faceWt;
			break;
		case N_DIRECTION:
			// vector Ap is normal to the face with magnitude area/4
			if(uGIMPsize>0.)
				ScaleVector(&Ap,uGIMPsize/sqrt(DotVectors(&Ap,&Ap)));
			wtNorm = Ap;
			break;
		default:
			// normal is z direction
			wtNorm.z = faceWt;
			break;
	}

	// Return hat n * Area/4 as vector
	return wtNorm;
}

// Get Rotation matrix for initial material orientation (anistropic only)
// Matrix in Rz.Ry.Rx with
//    Rz=((cz,sz,0),(-sz,cz,0),(0,0,1))
//    Ry=((cy,0,-sy),(0,1,0),(sy,0,cy))
//    Rx=((1,0,0),(0,cx,sx),(0,-sx,cx))
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
	
	// Return Rz.Ry.Rx
	return Matrix3( cy*cz,  cz*sx*sy+cx*sz,  -cx*cz*sy+sx*sz,
				   -cy*sz,  cx*cz-sx*sy*sz,   cz*sx+cx*sy*sz,
				   sy,     -cy*sx,            cx*cy );
	
	// This is (Rz.Ry.Rx)^T = RzT.RyT.RxT
	//return Matrix3( cy*cz,         -cy*sz,           sy,
	//			    cz*sx*sy+cx*sz, cx*cz-sx*sy*sz, -cy*sx,
	//			   -cx*cz*sy+sx*sz, cz*sx+cx*sy*sz,  cx*cy );
}

