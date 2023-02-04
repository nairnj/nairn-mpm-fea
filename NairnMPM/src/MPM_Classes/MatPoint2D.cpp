/********************************************************************************
    MatPoint2D.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MatPoint2D.hpp"
#include "Materials/MaterialBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark MatPoint2D::Constructors and Destructors

// Constructors
MatPoint2D::MatPoint2D() {}

// Constructors
MatPoint2D::MatPoint2D(int inElemNum,int theMatl,double angin,double thickin) : MPMBase(inElemNum,theMatl,angin)
{	thick=thickin;
}

#pragma mark MatPoint2D:Calculations and Incrementers

// Update Strains for this particle
// Velocities for all fields are present on the nodes
// matRef is the material and properties have been loaded, matFld is the material field
// postUpdate true means called after updating particle position
void MatPoint2D::UpdateStrain(double strainTime,int secondPass,int np,void *props,int matFld,bool postUpdate)
{
#ifdef CONST_ARRAYS
	int ndsArray[MAX_SHAPE_NODES];
	double fn[MAX_SHAPE_NODES],xDeriv[MAX_SHAPE_NODES],yDeriv[MAX_SHAPE_NODES],zDeriv[MAX_SHAPE_NODES];
#else
	int ndsArray[maxShapeNodes];
	double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
#endif
	Vector vel;
    Matrix3 dv;
	Tensor *gStressPtr=NULL;

	// don't need to zero zDeriv for 2D planar because never used in this function
    // Note that both plane strain and plain stress will get dv(2,2)=0. In plane stress,
    //    constitutive law should calculate ezz strain.
    
	// find shape functions and derviatives
	const ElementBase *elemRef = theElements[ElemID()];
	int *nds = ndsArray;
	elemRef->GetShapeGradients(fn,&nds,xDeriv,yDeriv,zDeriv,this);
	int numnds = nds[0];

    // Find strain rates at particle from current grid velocities
	//   and using the velocity field for that particle and each node and the right material
    for(int i=1;i<=numnds;i++)
	{	vel = nd[nds[i]]->GetVelocity((short)vfld[i],matFld);
        dv += Matrix3(vel.x*xDeriv[i],vel.x*yDeriv[i],vel.y*xDeriv[i],vel.y*yDeriv[i],0.);
    }

    // save velocity gradient (if needed for J integral calculation)
    SetVelocityGradient(dv(0,0),dv(1,1),dv(0,1),dv(1,0),secondPass);

    // convert to strain increments (e.g., now dvxx = dvx/dx * dt = d/dx(du/dt) * dt = d/dt(du/dx) * dt = du/dx)
    dv.Scale(strainTime);
	
	// pass on to material class to handle
	ResidualStrains res = ScaledResidualStrains(secondPass);
	PerformConstitutiveLaw(dv,strainTime,np,props,&res,gStressPtr);
}

// Pass on to material class
void MatPoint2D::PerformConstitutiveLaw(Matrix3 dv,double strainTime,int np,void *props,ResidualStrains *res,Tensor *gStress)
{
    // update particle strain and stress using its constitutive law
	const MaterialBase *matRef = theMaterials[MatID()];
    matRef->MPMConstitutiveLaw(this,dv,strainTime,np,props,res,0,gStress);
}

// Move position and velocity (2D)
void MatPoint2D::MoveParticle(GridToParticleExtrap *gp)
{
	// get vm = S(v-a dt) for FLIP and Sv^+ FMPM (and PIC)
	Vector vm;
	if(gp->m>0)
	{	// FMPM damps with Sk^+(k)
		vm.x = gp->Svtilde.x;
		vm.y = gp->Svtilde.y;
	}
	else if(gp->m>-2)
	{	// FLIP and XPIC(1) damps with Sv (which is initial lumped velocity)
		vm.x = gp->Svtilde.x - gp->Sacc.x*timestep;
		vm.y = gp->Svtilde.y - gp->Sacc.y*timestep;
	}
	else
	{	// XPIC(k>1) damps by initial lumped velocity, and Svtilde holds Sv(k)
		vm.x = gp->Svlumped.x - gp->Sacc.x*timestep;
		vm.y = gp->Svlumped.y - gp->Sacc.y*timestep;
	}

	// find Adamp0
	Vector Adamp0;
	Adamp0.x = gp->gridAlpha*vm.x + gp->particleAlpha*vel.x;
	Adamp0.y = gp->gridAlpha*vm.y + gp->particleAlpha*vel.y;
	
	Vector delV;
	if(gp->m>0)
	{	// FMPM update
		
		// save initial particle velocity
		Vector delXRate;
		delXRate.x = vel.x;
		delXRate.y = vel.y;
		
		// FMPM(m) update
		vel.x = vm.x - Adamp0.x*timestep;
		vel.y = vm.y - Adamp0.y*timestep;
		
		// find change in velocity
		delV.x = vel.x - delXRate.x;
		delV.y = vel.y - delXRate.y;
		
		// position update
		// del X = 0.5*(V(n+1)+V(n))
        delXRate.x += vel.x;
        delXRate.y += vel.y;
		pos.x += 0.5*delXRate.x*timestep;
		pos.y += 0.5*delXRate.y*timestep;
	}
	else if(gp->m==0)
	{	// FLIP update velcity change
		delV.x = (gp->Sacc.x - Adamp0.x)*timestep;
		delV.y = (gp->Sacc.y - Adamp0.y)*timestep;
		
		// velocity
		vel.x += delV.x;
		vel.y += delV.y;
		
		// position update
		Vector delXRate;
		delXRate.x = vm.x + 0.5*delV.x;
		delXRate.y = vm.y + 0.5*delV.y;
		pos.x += delXRate.x*timestep;
		pos.y += delXRate.y*timestep;
	}
	else
	{	// XPIC update
		
		// save initial particle velocity
		Vector delXRate;
		delXRate.x = vel.x;
		delXRate.y = vel.y;
		
		// XPIC(k) velocity update
		// For XPIC(1) Svtilde has Sv+, for XPIC(k>1) Svtilde has S(v(k)+a*dt)
		vel.x = gp->Svtilde.x - Adamp0.x*timestep;
		vel.y = gp->Svtilde.y - Adamp0.y*timestep;
		
		// find change in velocity
		delV.x = vel.x - delXRate.x;
		delV.y = vel.y - delXRate.y;
		
		// position update
		delXRate.x = vm.x + 0.5*delV.x;
		delXRate.y = vm.y + 0.5*delV.y;
		pos.x += delXRate.x*timestep;
		pos.y += delXRate.y*timestep;
	}
	
	// J Integral needs effective particle acceleration
	acc = MakeVector(delV.x/timestep,delV.y/timestep,0.);
		
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
void MatPoint2D::MovePosition(double delTime)
{	pos.x += delTime*vel.x;
    pos.y += delTime*vel.y;
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

// set velocity (2D) (in mm/sec)
void MatPoint2D::SetVelocity(Vector *pt)
{	vel.x=pt->x;
    vel.y=pt->y;
	vel.z=0.;
}

// thickness (in mm)
double MatPoint2D::thickness() { return thick; }

// Find internal force as -mp sigma.deriv in g mm/sec^2 or micro N
// add external force (times a shape function)
// Store in buffer
// (note: stress is specific stress in units N/m^2 mm^3/g which is (g-mm^2/sec^2)/g
void MatPoint2D::GetFintPlusFext(Vector *theFrc,double fni,double xDeriv,double yDeriv,double zDeriv)
{	
	theFrc->x = -mp*((sp.xx-pressure)*xDeriv+sp.xy*yDeriv) + fni*pFext.x;
	theFrc->y = -mp*(sp.xy*xDeriv+(sp.yy-pressure)*yDeriv) + fni*pFext.y;
	theFrc->z = 0.0;
	
#ifdef CHECK_NAN
	if(theFrc->x!=theFrc->x || theFrc->y!=theFrc->y)
	{
#pragma omp critical (output)
		{	cout << "\n# MatPoint2D::GetFintPlusFext: bad nodal fint+fext";
			PrintVector(" = ",theFrc);
			cout << endl;
			cout << "# sp = (" << sp.xx << "," << sp.yy << "," << sp.xy << ") P = " << pressure << endl;
			PrintVector("# pFext = ",&pFext);
			cout << endl;
			cout << "# Sip = " << fni << " Gip = (" << xDeriv << "," << yDeriv << "," << zDeriv << ")" << endl;
		}
	}
#endif
}

// add to the temperature gradient (non-rigid particles only)
void MatPoint2D::AddTemperatureGradient(int offset,Vector *grad)
{	pTemp[offset]+=grad->x;
    pTemp[offset+1]+=grad->y;
}

// return conduction force = -Vp [k/rho0] Grad T . Grad S = - mp (Vp/V0) [k/rho0] Grad T . Grad S (units nJ/sec)
// and k/rho0 is stored in k in units (nJ mm^2/(sec-K-g))
//  (non-rigid particles only)
double MatPoint2D::FCond(int offset,double dshdx,double dshdy,double dshdz,TransportProperties *t)
{
	Tensor *kten = &(t->kCondTensor);
	return -mp*GetRelativeVolume()*((kten->xx*pTemp[offset] + kten->xy*pTemp[offset+1])*dshdx
						+ (kten->xy*pTemp[offset] + kten->yy*pTemp[offset+1])*dshdy);
}

// add to the concentration gradient (1/mm) (non-rigid particles only)
void MatPoint2D::AddConcentrationGradient(int dnum,Vector *grad)
{	pDiff[dnum]->grad.x += grad->x;
	pDiff[dnum]->grad.y += grad->y;
}

// return diffusion force = - V [D] Grad C . Grad S in (mm^3) (mm^2/sec) (1/mm) (1/mm) = mm^3/sec
// (non-rigid particles only)
double MatPoint2D::FDiff(double dshdx,double dshdy,double dshdz,TransportProperties *t,int number)
{
	Tensor *Dten = &(t->diffusionTensor);
    Vector pgrad = pDiff[number]->grad;
	return -GetVolume(DEFORMED_VOLUME)*((Dten->xx*pgrad.x + Dten->xy*pgrad.y)*dshdx
                                      + (Dten->xy*pgrad.x + Dten->yy*pgrad.y)*dshdy);
}

// return kinetic energy (g mm^2/sec^2) = nanoJ
double MatPoint2D::KineticEnergy(void)
{	return 0.5*mp*(vel.x*vel.x+vel.y*vel.y);
}

// get deformation gradient, which is stored in strain and rotation tensors
Matrix3 MatPoint2D::GetDeformationGradientMatrix(void) const
{	double F[3][3];
	GetDeformationGradient(F);
	Matrix3 Fm(F[0][0],F[0][1],F[1][0],F[1][1],F[2][2]);
	return Fm;
}

// get deformation gradient, which is stored in strain and rotation tensors
void MatPoint2D::SetDeformationGradientMatrix(Matrix3 F)
{	
	// Normal strains
	ep.xx = F(0,0) - 1.;
	ep.yy = F(1,1) - 1.;
	ep.zz = F(2,2) - 1.;
	
	// Shear Strain
	ep.xy = F(1,0) + F(0,1);
	
	// Rotation Strain
	wrot.xy = F(1,0) - F(0,1);
}

// get displacement gradient in current configuration
// Uses grad u = RFR^T - I = RV - I
Matrix3 MatPoint2D::GetDisplacementGradientMatrix(void) const
{	Matrix3 F = GetDeformationGradientMatrix();
	Matrix3 Rn;
	Matrix3 V = F.LeftDecompose(&Rn,NULL);
	V = Rn*V;
	V(0,0) -= 1.;
	V(1,1) -= 1.;
	V(2,2) -= 1.;
	return V;
}

// get displacement gradient in current conficguration using
//     grad u = R(V-ec)-I
// where ec is cracking strain in softening materials
// this method us current for 2D only
Matrix3 MatPoint2D::GetDisplacementGradientForJ(const MaterialBase *matref)
{
	Matrix3 F = GetDeformationGradientMatrix();
	Matrix3 Rn;
	Matrix3 V = F.LeftDecompose(&Rn,NULL);
	
	// if needed, subtract cracking strain
	Tensor ecrack;
	if(matref->GetCrackingStrain(this,&ecrack,true,&Rn))
	{	V(0,0) -= ecrack.xx;
		V(0,1) -= 0.5*ecrack.xy;
		V(1,0) -= 0.5*ecrack.xy;
		V(1,1) -= ecrack.yy;
	}
	
	// subtract identity
	V = Rn*V;
	V(0,0) -= 1.;
	V(1,1) -= 1.;
	return V;
}

// get the symmetric elastic Left-Cauchy tensor in a Matrix3
Matrix3 MatPoint2D::GetElasticLeftCauchyMatrix(void)
{   return Matrix3(eplast.xx,eplast.xy,eplast.xy,eplast.yy,eplast.zz);
}

// get deformation gradient, which is stored in strain and rotation tensors
void MatPoint2D::GetDeformationGradient(double F[][3]) const
{
	// current deformation gradient in 2D
	F[0][0] = 1. + ep.xx;
	F[0][1] = 0.5*(ep.xy - wrot.xy);
	F[1][0] = 0.5*(ep.xy + wrot.xy);
	F[1][1] = 1. + ep.yy;
	F[2][2] = 1. + ep.zz;
}

// Get sqrt(B)-I = V-I for elastic biot strain rotated into current particle orientation
// Assume that elastic B matrix is in the alt strain tensor
Matrix3 MatPoint2D::GetElasticBiotStrain(void)
{	// Get Sqrt(B)
	Tensor *B = GetAltStrainTensor();
	Matrix3 Be = Matrix3(B->xx,B->xy,B->xy,B->yy,B->zz);
	Vector lam = Be.Eigenvalues();
	Matrix3 Q = Be.Eigenvectors(lam);
	Matrix3 LamI = Matrix3(sqrt(lam.x),0.,0.,sqrt(lam.y),sqrt(lam.z));
	Matrix3 V = LamI.RMRT(Q);
	V(0,0) -= 1.;
	V(1,1) -= 1.;
	V(2,2) -= 1.;
	return V;
}

// get relative volume from det J for large deformation material laws
double MatPoint2D::GetRelativeVolume(void)
{
	double pF[3][3];
    GetDeformationGradient(pF);
    return pF[2][2]*(pF[0][0]*pF[1][1]-pF[1][0]*pF[0][1]);
}

// Get dilated current volume using current deformation gradient
// (only used for crack contact, multimaterial contact, and transport tasks)
// when volumeType is DEFORMED_AREA or DEFORMED_AREA_FOR_GRADIENT, get t0*Ap, where
//    Ap is deformed particle area and t0 in initial thickness
double MatPoint2D::GetVolume(int volumeType)
{	double rho=GetRho();								// in g/mm^3
	if(volumeType==DEFORMED_VOLUME)
		return GetRelativeVolume()*mp/rho;				// in mm^3
	
	// get thickness times initial area (for contact and for gradient)
	// Note that mp/rho = rho Ap0 t0/rho = Ap0 * t0
	double pF[3][3];
	GetDeformationGradient(pF);
	return (pF[0][0]*pF[1][1]-pF[1][0]*pF[0][1])*mp/rho;
}

// Get vectors from particle to edge
// r1 = F.(psz.x,0) and r2 = F.(0,psz.y)
// where psz is semiparticle size
// This uses mpmgrid so membrane particles must override
void MatPoint2D::GetSemiSideVectors(Vector *r1,Vector *r2,Vector *r3) const
{
	// get particle 2D deformation gradient
	double pF[3][3];
	GetDeformationGradient(pF);
	
	// get polygon vectors - these are from particle to edge
    //      and generalize semi width lp in 1D GIMP
    Vector psz = GetParticleSize();
	r1->x = pF[0][0]*psz.x;
	r1->y = pF[1][0]*psz.x;
	r1->z = 0.;
    
	r2->x = pF[0][1]*psz.y;
	r2->y = pF[1][1]*psz.y;
	r2->z = 0.;
}

// rescale vectors if needed in CPDI calculations
void MatPoint2D::ScaleSemiSideVectorsForCPDI(Vector *r1,Vector *r2,Vector *r3) const
{	// skip if not activated
	if(ElementBase::rcrit<0.) return;
	
	// convert rcrit to dimensions based smallest dimension this element
	Vector cellSize = theElements[ElemID()]->GetDeltaBox();
	double rcrit = ElementBase::rcrit*fmin(cellSize.x,cellSize.y);
	
	bool rescale = false;
	
	// check r1+r2
	Vector la = MakeVector(r1->x+r2->x, r1->y+r2->y, 0.);
	double lamag = sqrt(la.x*la.x+la.y*la.y);
	if(lamag>rcrit)
	{	ScaleVector(&la,rcrit/lamag);
		rescale = true;
	}
	
	// check r1-r2
	Vector lb = MakeVector(r1->x-r2->x, r1->y-r2->y, 0.);
	double lbmag = sqrt(lb.x*lb.x+lb.y*lb.y);
	if(lbmag>rcrit)
	{	ScaleVector(&lb,rcrit/lbmag);
		rescale = true;
	}
	
	// exit if no scale
	if(!rescale) return;
	
	// redefine r1 and r2
	r1->x = 0.5*(la.x+lb.x);
	r1->y = 0.5*(la.y+lb.y);
	
	r2->x = 0.5*(la.x-lb.x);
	r2->y = 0.5*(la.y-lb.y);
}

#define ELLIPSE_METHOD_2D
// Get distance from particle center to surface of deformed particle
// along the provided unit vector
double MatPoint2D::GetDeformedRadius(Vector *norm) const
{
#ifdef ELLIPSE_METHOD_2D
	// unnormalized normal in undeformed config (n0x,n0y) = F^{-1}n
	double F[3][3];
	GetDeformationGradient(F);
	double rdet2D = 1./(F[0][0]*F[1][1]-F[0][1]*F[1][0]);
	double n0x =  (F[1][1]*norm->x - F[1][0]*norm->y)*rdet2D;
	double n0y = (-F[0][1]*norm->x + F[0][0]*norm->y)*rdet2D;
	
	// return radius to deformed inscribed ellipse
	// see JANOSU-014-82
	Vector psz = GetParticleSize();
	double rd = sqrt(n0x*n0x/(psz.x*psz.x) + n0y*n0y/(psz.y*psz.y));
	return 1./rd;
#else
	Vector r1,r2,r3;
	GetSemiSideVectors(&r1,&r2,&r3);
	double amag = fabs(norm->x*r1.y-norm->y*r1.x);
	double bmag = fabs(norm->x*r2.y-norm->y*r2.x);
	return fabs(r1.x*r2.y-r1.y*r2.x)/fmax(amag,bmag);
#endif
}

// Get undeformed size (only called by GIMP traction and Custom thermal ramp)
// This uses mpmgrid so membrane particles must override
// Assumes undeformation particle aligned with x-y-z axes
void MatPoint2D::GetUndeformedSemiSides(double *r1x,double *r2y,double *r3z) const
{   Vector psz = GetParticleSize();
    *r1x = psz.x;
    *r2y = psz.y;
}

// To support CPDI find nodes in the particle domain, find their element,
// their natural coordinates, and weighting values for gradient calculations
// Should be done only once per time step
// throws CommonException() if particle corner has left the grid
void MatPoint2D::GetCPDINodesAndWeights(int cpdiType)
{
	// get polygon vectors - these are from particle to edge
    //      and generalize semi width lp in 1D GIMP
	Vector r1,r2,c;
    GetSemiSideVectors(&r1,&r2,NULL);
	
	// rescale if needed
	ScaleSemiSideVectorsForCPDI(&r1,&r2,NULL);
	
    // Particle domain area is area of the full parallelogram
    // Assume positive due to orientation of initial vectors, and sign probably does not matter
    double Ap = 4.*(r1.x*r2.y - r1.y*r2.x);
    
	try
	{	CPDIDomain **cpdi = GetCPDIInfo();
		
		if(cpdiType != QUADRATIC_CPDI)
		{	// nodes at four corners in ccw direction
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
			
			// gradient weighting values
			// wg[i] = (1/Ap)*(xi[i](r2.y,-r2.x) + eta[i](-r1.y,r1.x))
			// where xi = {-1,1,1,-1} and eta = {-1,-1,1,1}
			Ap = 1./Ap;
			cpdi[0]->wg.x = (r1.y-r2.y)*Ap;
			cpdi[1]->wg.x = (r1.y+r2.y)*Ap;
			cpdi[2]->wg.x = (-r1.y+r2.y)*Ap;
			cpdi[3]->wg.x = (-r1.y-r2.y)*Ap;
			
			cpdi[0]->wg.y = (-r1.x+r2.x)*Ap;
			cpdi[1]->wg.y = (-r1.x-r2.x)*Ap;
			cpdi[2]->wg.y = (r1.x-r2.x)*Ap;
			cpdi[3]->wg.y = (r1.x+r2.x)*Ap;
		}
		
		else
		{	// nodes at four corners in ccw direction
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
			
			// nodes at four edges in ccw direction
			c.x = pos.x-r2.x;
			c.y = pos.y-r2.y;
			cpdi[4]->inElem = mpmgrid.FindElementFromPoint(&c,this)-1;
			theElements[cpdi[4]->inElem]->GetXiPos(&c,&cpdi[4]->ncpos);
			
			c.x = pos.x+r1.x;
			c.y = pos.y+r1.y;
			cpdi[5]->inElem = mpmgrid.FindElementFromPoint(&c,this)-1;
			theElements[cpdi[5]->inElem]->GetXiPos(&c,&cpdi[5]->ncpos);
			
			c.x = pos.x+r2.x;
			c.y = pos.y+r2.y;
			cpdi[6]->inElem = mpmgrid.FindElementFromPoint(&c,this)-1;
			theElements[cpdi[6]->inElem]->GetXiPos(&c,&cpdi[6]->ncpos);
			
			c.x = pos.x-r1.x;
			c.y = pos.y-r1.y;
			cpdi[7]->inElem = mpmgrid.FindElementFromPoint(&c,this)-1;
			theElements[cpdi[7]->inElem]->GetXiPos(&c,&cpdi[7]->ncpos);
			
			// node on material point
			cpdi[8]->inElem = ElemID();
			theElements[cpdi[8]->inElem]->GetXiPos(&pos,&cpdi[8]->ncpos);
			
			// gradient weighting values - use linear weights
			Ap = 1./(3.*Ap);
			cpdi[0]->wg.x = (r1.y-r2.y)*Ap;
			cpdi[0]->wg.y = (-r1.x+r2.x)*Ap;
			cpdi[1]->wg.x = (r1.y+r2.y)*Ap;
			cpdi[1]->wg.y = (-r1.x-r2.x)*Ap;
			cpdi[2]->wg.x = (-r1.y+r2.y)*Ap;
			cpdi[2]->wg.y = (r1.x-r2.x)*Ap;
			cpdi[3]->wg.x = (-r1.y-r2.y)*Ap;
			cpdi[3]->wg.y = (r1.x+r2.x)*Ap;
			cpdi[4]->wg.x = 4.*r1.y*Ap;
			cpdi[4]->wg.y = -4.*r1.x*Ap;
			cpdi[5]->wg.x = 4.*r2.y*Ap;
			cpdi[5]->wg.y = -4.*r2.x*Ap;
			cpdi[6]->wg.x = -4.*r1.y*Ap;
			cpdi[6]->wg.y = 4.*r1.x*Ap;
			cpdi[7]->wg.x = -4.*r2.y*Ap;
			cpdi[7]->wg.y = 4.*r2.x*Ap;
			cpdi[8]->wg.x = 0.;
			cpdi[8]->wg.y = 0.;
			
			/*
			int i;
			for(i=4;i<9;i++)
			{	cpdi[i]->wg.x = 0.;
				cpdi[i]->wg.y = 0.;
			}
			*/
		}
	}
    catch(CommonException& err)
    {   char msg[200];
		size_t msgSize=200;
        snprintf(msg,msgSize,"A CPDI particle domain node has left the grid: %s",err.Message());
        throw CommonException(msg,"MatPoint2D::GetCPDINodesAndWeights");
    }
    
    // traction BC area saves 1/2 surface area of particle domain on the various edges
    if(faceArea!=NULL)
    {   faceArea->x = sqrt(r1.x*r1.x+r1.y*r1.y)*mpmgrid.GetThickness();			// edges 1 and 3
        faceArea->y = sqrt(r2.x*r2.x+r2.y*r2.y)*mpmgrid.GetThickness();			// edges 2 and 4
    }
}

// Called once per time step in the initialization task by GetShapeFunctionData()
// throws CommonException() if particle corner has left the grid
void MatPoint2D::GetFiniteGIMP_Integrals(void)
{
	// get polygon vectors - these are from particle to edge
	//      and generalize semi width lp in 1D GIMP
	Vector r1,r2,minc,maxc;
	GetSemiSideVectors(&r1,&r2,NULL);
	Vector particle[4];  // particle corners
	
	try
	{	FiniteGIMPInfo *fgimp = GetFiniteGIMPInfo();
		
		// nodes at four corners in ccw direction
		double rx = fabs(r1.x)+fabs(r2.x);
		double ry = fabs(r1.y)+fabs(r2.y);
		minc.x = pos.x-rx;
		maxc.x = pos.x+rx;
		minc.y = pos.y-ry;
		maxc.y = pos.y+ry;		
		
		// zero-based coordinates to extreme elements
		int lx,ly,lz,ux,uy,uz;
		mpmgrid.FindElementCoordinatesFromPoint(&minc,lx,ly,lz);
		mpmgrid.FindElementCoordinatesFromPoint(&maxc,ux,uy,uz);
				
		/*if(check to see if we have enough nodes){
			char msg[200];
            size_t msgSize=200;
			snprintf(msg,msgSize,"Found more nodes than maximum allowed. Might need to change maximum: %s",err.Message());
			throw CommonException(msg,"MatPoint2D::GetFiniteGIMP_Integrals");
		}*/ 
		
		// Find the corners of the particle 
		particle[0].x = pos.x-r1.x-r2.x;  
		particle[0].y = pos.y-r1.y-r2.y;  // corner1
		particle[1].x = pos.x+r1.x-r2.x;
		particle[1].y = pos.y+r1.y-r2.y;  //corner 2
		particle[2].x = pos.x+r1.x+r2.x;
		particle[2].y = pos.y+r1.y+r2.y;  //corner 3 
		particle[3].x = pos.x-r1.x+r2.x; 
		particle[3].y = pos.y-r1.y+r2.y;  //corner 4
		
		// r1.x*r2.y - r1.y*r2.x = Ap/4 (Ap is parallelgram area). Get 1./(4*Ap)
		double DivideByVolume = 1./(16.*fabs(r1.x*r2.y - r1.y*r2.x));

		// lower right at (row,col)=(lx,ly) and upper right at (row,col)=(ux,uy)
		// Search all elements in the window for intersection the paterial point
		int foundElement = 0;
		int i_1;
        double interp;
		for(int col=lx;col<=ux;col++)
		{	
			for(int row=ly;row<=uy;row++)
			{	
				// zero based element number
				int iel = row*(mpmgrid.yplane-1) + col;
				Vector Polygon[12];    // polygon for intersection
				Vector Polygon2[12];    // polygon for intersection
				ElementBase *elem = theElements[iel];	
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// first clip on right side
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
				bool previous_vertex_inside = (particle[3].x <= elem->xmax);
				bool current_vertex_inside=false;
				int nSides = 0; // number of corners in Polygon
				for(int i=0;i<4;i++){
					current_vertex_inside = (particle[i].x <= elem->xmax);
					if(current_vertex_inside && previous_vertex_inside){
						Polygon[nSides].x = particle[i].x;
						Polygon[nSides].y = particle[i].y;
						nSides++;  // add this point to the list
					}else if(!current_vertex_inside && previous_vertex_inside){
						// do intersection add this point to list
						i_1 = (i==0)?3:i-1;  // because c++ modulus doesn't work on negative numbers
						interp = ((elem->xmax-particle[i].x)/(particle[i_1].x-particle[i].x));
						
						Polygon[nSides].y = interp*particle[i_1].y+(1.0-interp)*particle[i].y;
						Polygon[nSides].x = elem->xmax;
						nSides++;  // add this point to the list
					}else if(current_vertex_inside && !previous_vertex_inside){
						// do intersection and add this point to list
						i_1 = (i==0)?3:i-1; // because c++ modulus doesn't work on negative numbers
						interp = ((elem->xmax-particle[i].x)/(particle[i_1].x-particle[i].x));
						
						Polygon[nSides].y = interp*particle[i_1].y+(1.0-interp)*particle[i].y;
						Polygon[nSides].x = elem->xmax;
						nSides++;  // add this point to the list
						
						Polygon[nSides].x = particle[i].x;
						Polygon[nSides].y = particle[i].y;
						nSides++;  // add this point to the list
					}
					previous_vertex_inside = current_vertex_inside;
				}
				if(nSides<3) continue; // no overlap, skip this  element
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// now clip on left side
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
				previous_vertex_inside = (Polygon[nSides-1].x >= elem->xmin);
				current_vertex_inside=false;
				int nSides2 =0; // number of corners in Polygon2
				for(int i=0;i<nSides;i++){
					current_vertex_inside = (Polygon[i].x >= elem->xmin);
					if(current_vertex_inside && previous_vertex_inside){
						Polygon2[nSides2].x = Polygon[i].x;
						Polygon2[nSides2].y = Polygon[i].y;
						nSides2++;  // add this point to the list
					}else if(!current_vertex_inside && previous_vertex_inside){
						// do intersection and add this point to list
						i_1 = (i==0)?(nSides-1):i-1;  // because c++ modulus doesn't work on negative numbers
						interp = ((elem->xmin-Polygon[i].x)/(Polygon[i_1].x-Polygon[i].x));
						
						Polygon2[nSides2].y = interp*Polygon[i_1].y+(1.0-interp)*Polygon[i].y;
						Polygon2[nSides2].x = elem->xmin;
						nSides2++;  // add this point to the list
					}else if(current_vertex_inside && !previous_vertex_inside){
						// do intersection and add this point to list
						i_1 = (i==0)?(nSides-1):i-1;  // because c++ modulus doesn't work on negative numbers
						interp = ((elem->xmin-Polygon[i].x)/(Polygon[i_1].x-Polygon[i].x));
						
						Polygon2[nSides2].y = interp*Polygon[i_1].y+(1.0-interp)*Polygon[i].y;
						Polygon2[nSides2].x = elem->xmin;
						nSides2++;  // add this point to the list
						
						Polygon2[nSides2].x = Polygon[i].x;
						Polygon2[nSides2].y = Polygon[i].y;
						nSides2++;  // add this point to the list
					}
					previous_vertex_inside = current_vertex_inside;
				}
				if(nSides2<3) continue; // no overlap, skip this  element
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// now clip on bottom
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
				previous_vertex_inside = (Polygon2[nSides2-1].y >= elem->ymin);
				current_vertex_inside=false;
				nSides =0;
				for(int i=0;i<nSides2;i++){
					current_vertex_inside = (Polygon2[i].y >= elem->ymin);
					if(current_vertex_inside && previous_vertex_inside){
						Polygon[nSides].x = Polygon2[i].x;
						Polygon[nSides].y = Polygon2[i].y;
						nSides++;  // add this point to the list
					}else if(!current_vertex_inside && previous_vertex_inside){
						// do intersection and add this point to list
						i_1 = (i==0)?(nSides2-1):i-1;  // because c++ modulus doesn't work on negative numbers
						interp = ((elem->ymin-Polygon2[i].y)/(Polygon2[i_1].y-Polygon2[i].y));
						
						Polygon[nSides].x = interp*Polygon2[i_1].x+(1.0-interp)*Polygon2[i].x;
						Polygon[nSides].y = elem->ymin;
						nSides++;  // add this point to the list
					}else if(current_vertex_inside && !previous_vertex_inside){
						// do intersection and add this point to list
						i_1 = (i==0)?(nSides2-1):i-1;  // because c++ modulus doesn't work on negative numbers
						interp = ((elem->ymin-Polygon2[i].y)/(Polygon2[i_1].y-Polygon2[i].y));
						
						Polygon[nSides].x = interp*Polygon2[i_1].x+(1.0-interp)*Polygon2[i].x;
						Polygon[nSides].y = elem->ymin;
						nSides++;  // add this point to the list
						
						// also current point to list
						Polygon[nSides].x = Polygon2[i].x;
						Polygon[nSides].y = Polygon2[i].y;
						nSides++;  // add this point to the list
					}
					previous_vertex_inside = current_vertex_inside;
				}
				if(nSides<3) continue; // no overlap, skip this element
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// finally clip on the top
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
				previous_vertex_inside = (Polygon[nSides-1].y <= elem->ymax);
				current_vertex_inside=false;
				nSides2 = 0;
				for(int i=0;i<nSides;i++){
					current_vertex_inside = (Polygon[i].y <= elem->ymax);
					if(current_vertex_inside && previous_vertex_inside){
						Polygon2[nSides2].x = Polygon[i].x;
						Polygon2[nSides2].y = Polygon[i].y;
						nSides2++;  // add this point to the list
					}else if(!current_vertex_inside && previous_vertex_inside){
						// do intersection and add this point to list
						i_1 = (i==0)?(nSides-1):i-1;  // because c++ modulus doesn't work on negative numbers
						interp = ((elem->ymax-Polygon[i].y)/(Polygon[i_1].y-Polygon[i].y));
						
						Polygon2[nSides2].x = interp*Polygon[i_1].x+(1.0-interp)*Polygon[i].x;
						Polygon2[nSides2].y = elem->ymax;
						nSides2++;  // add this point to the list
					}else if(current_vertex_inside && !previous_vertex_inside){
						// do intersection and add this point to list
						i_1 = (i==0)?(nSides-1):i-1;  // because c++ modulus doesn't work on negative numbers
						interp = ((elem->ymax-Polygon[i].y)/(Polygon[i_1].y-Polygon[i].y));
						
						Polygon2[nSides2].x = interp*Polygon[i_1].x+(1.0-interp)*Polygon[i].x;
						Polygon2[nSides2].y = elem->ymax;
						nSides2++;  // add this point to the list
						
						// also add current point to list
						Polygon2[nSides2].x = Polygon[i].x;
						Polygon2[nSides2].y = Polygon[i].y;
						nSides2++;  // add this point to the list
					}
					previous_vertex_inside = current_vertex_inside;
				}
				
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// If we found enough corners, calculate moments
				// Here Polygon2 of nSides2 is intersection of element with the material point
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				nSides=nSides2; // for my sanity
				if(nSides>2){
					
					for(int i=0;i<nSides;i++){
						elem->GetXiPos(&Polygon2[i],&Polygon[i]);// transform to element coordinate system
					}

					// Find Moments or integral over area of (1, x, y, xy)
					double moment_0 = 0.;
					double moment_x = 0.;
					double moment_y = 0.;
					double moment_xy = 0.;
					for(int i=0;i<nSides;i++){
						int k = (i+1)%nSides;
						double cross_product = Polygon[i].x*Polygon[k].y-Polygon[k].x*Polygon[i].y;
						moment_0 += cross_product;
						moment_x += cross_product*(Polygon[i].x+Polygon[k].x);
						moment_y += cross_product*(Polygon[i].y+Polygon[k].y);
						moment_xy += cross_product*(2.0*(Polygon[i].x*Polygon[i].y+Polygon[k].x*Polygon[k].y)+ Polygon[i].x*Polygon[k].y+Polygon[k].x*Polygon[i].y);
					}
					
					//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					// Output stuff
					//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					if(moment_0>1.0e-16){	//	if the area is zero relative to 1, then didn't find an element	
						if(foundElement >= maxElementIntersections){
							char msg[200];
							size_t msgSize=200;
							snprintf(msg,msgSize,"Found more elements than allocated. Might need to change maximum number of intersected elements");
							throw CommonException(msg,"MatPoint2D::GetFiniteGIMP_Integrals");
						}
					
						// scale by (det J)/Ap = (dx*dy/4)/Ap = (dx*dy) * (1/(4*Ap)) (latter from above)
						//double Scale4Grid = elem->GetDeltaX()*elem->GetDeltaY()*DivideByVolume;
						double Scale4Grid = elem->GetDeltaX()*elem->GetDeltaY()*DivideByVolume;
						fgimp->InTheseElements[foundElement+1] = iel;
						fgimp->moment_0[foundElement] = Scale4Grid*moment_0/8.0;  // divide by 2 for moment, by 4 for shape function
						fgimp->moment_x[foundElement] = Scale4Grid*moment_x/24.0;  // divide by 6 for moment, by 4 for shape function
						fgimp->moment_y[foundElement] = Scale4Grid*moment_y/24.0;  // divide by 6 for moment, by 4 for shape function
						fgimp->moment_xy[foundElement] = Scale4Grid*moment_xy/96.0;   // divide by 24 for moment, by 2 for shape function
						foundElement++; // found another element
					}
				}
			}
		}
		// we found this many elements
		fgimp->InTheseElements[0] = foundElement;
	}
	
	catch(CommonException& err)
	{   char msg[200];
		size_t msgSize=200;
		snprintf(msg,msgSize,"A particle corner node has left the grid: %s",err.Message());
		throw CommonException(msg,"MatPoint2D::GetFiniteGIMP_Integrals");
	}
	
}

// To support exact tractions get a list of elements containining the edge, the coordinates
// of the endpoints of intersection of edge with that element, and a scaling factor normal
// to traction and proportional to length of intersected segment.
// On input numDnds is maximum number of elements that can be found
// throws CommonException() if edge has left the grid
void MatPoint2D::GetExactTractionInfo(int face,int dof,int *cElem,Vector *corners,Vector *tscaled,int *numDnds) const
{
	int maxFoundElements = *numDnds;
	
	// get polygon vectors - these are from particle to edge
	//      and generalize semi width lp in 1D GIMP
	Vector r1,r2,c1,c2;
	GetSemiSideVectors(&r1,&r2,NULL);
	
	// find end points of edge 1, 2, 3, or 4 and c1 and c2
	switch(face)
	{	case 1:
			// lower edge 1 to 2
			c1.x = pos.x-r1.x-r2.x;
			c1.y = pos.y-r1.y-r2.y;
			c2.x = pos.x+r1.x-r2.x;
			c2.y = pos.y+r1.y-r2.y;
			break;
			
		case 2:
			// right edge 2 to 3
			c1.x = pos.x+r1.x-r2.x;
			c1.y = pos.y+r1.y-r2.y;
			c2.x = pos.x+r1.x+r2.x;
			c2.y = pos.y+r1.y+r2.y;
			break;
			
		case 3:
			// top edge 3 to 4
			c1.x = pos.x+r1.x+r2.x;
			c1.y = pos.y+r1.y+r2.y;
			c2.x = pos.x-r1.x+r2.x;
			c2.y = pos.y-r1.y+r2.y;
			break;
			
		default:
			// left edge 4 to 1
			c1.x = pos.x-r1.x+r2.x;
			c1.y = pos.y-r1.y+r2.y;
			c2.x = pos.x-r1.x-r2.x;
			c2.y = pos.y-r1.y-r2.y;
			break;
	}
	c1.z = 0.;
	c2.z = 0.;
	
	// Find range of elements with this edge
	int lx,ly,lz,ux,uy,uz;
	try
	{	mpmgrid.FindElementCoordinatesFromPoint(&c1,lx,ly,lz);
		mpmgrid.FindElementCoordinatesFromPoint(&c2,ux,uy,uz);
	}
	catch(CommonException& err)
	{   char msg[200];
		size_t msgSize=200;
		snprintf(msg,msgSize,"A Traction edge node has left the grid: %s",err.Message());
		throw CommonException(msg,"MatPoint2D::GetExactTractionInfo");
	}
	
	// swap if needed
	if(uy<ly)
	{	int temp = uy;
		uy = ly;
		ly = temp;
	}
	if(ux<lx)
	{	int temp = ux;
		ux = lx;
		lx = temp;
	}
	
	// get edge vector
	double ex = c2.x-c1.x;
	double ey = c2.y-c1.y;
	double p[4],q[4];
	p[0] = -ex;
	p[1] = ex;
	p[2] = -ey;
	p[3] = ey;
	
	// get normal vector in direction of the traction divided by 4 (n/4)
	Vector tnorm;
	double enorm;
	switch(dof)
	{	case 1:
			// normal is x direction
			tnorm = MakeVector(0.25,0.,0.);
			break;
		case 2:
			// normal is y direction
			tnorm = MakeVector(0.,0.25,0.);
			break;
		case N_DIRECTION:
			// cross product of edge vector with (0,0,1) = (ey, -ex)
			enorm = ey;
			ey = -ex;
			ex = enorm;
		case T1_DIRECTION:
			// load in direction specified by normalized (ex,ey)
			enorm = 4.*sqrt(ex*ex+ey*ey);
			tnorm = MakeVector(ex/enorm,ey/enorm,0.);
			break;
		default:
			// normal is z direction (not used in 2D)
			break;
	}
	
	// loop over possible elements
	int iel,iel0,foundElement=0,foundCorner=0;
	for(int row=ly;row<=uy;row++)
	{	// element at start of row
		iel0 = row*(mpmgrid.yplane-1);
		
		for(int col=lx;col<=ux;col++)
		{	// zero based element number
			iel = iel0 + col;
			ElementBase *elem = theElements[iel];
			
			// Liang-Barsky algorithm terms
			q[0] = c1.x - elem->xmin;
			q[1] = elem->xmax - c1.x;
			q[2] = c1.y - elem->ymin;
			q[3] = elem->ymax - c1.y;
			
			// check for vertical line in element, but not on the right edge
			double t1 = 0., t2 = 1.;
			if(p[0]==0.)
			{	// keep only if intersects this element
				if(c1.x<elem->xmin || c1.x>elem->xmax) continue;
				if(p[2]<0.)
				{	t1 = fmax(t1,q[2]/p[2]);
					t2 = fmin(t2,q[3]/p[3]);
				}
				else
				{	t1 = fmax(t1,q[3]/p[3]);
					t2 = fmin(t2,q[2]/p[2]);
				}
			}
			
			// check for horizontal line in element, but not on the top edge
			else if(p[2]==0.)
			{	// keep only if intersects this element
				if(c1.y<elem->ymin || c1.y>elem->ymax) continue;
				if(p[0]<0.)
				{	t1 = fmax(t1,q[0]/p[0]);
					t2 = fmin(t2,q[1]/p[1]);
				}
				else
				{	t1 = fmax(t1,q[1]/p[1]);
					t2 = fmin(t2,q[0]/p[0]);
				}
			}
			
			// non-vertical and non-horizontal lines
			else
			{	for(int k=0;k<4;k++)
				{	if(p[k]<0.)
						t1 = fmax(t1,q[k]/p[k]);
					else
						t2 = fmin(t2,q[k]/p[k]);
				}
			}
			
			// line is from t1 to t2
			if(t1>=t2) continue;				// no intersection with this element
			
			// Is there room?
			if(foundElement>=maxFoundElements)
			{   char msg[200];
				size_t msgSize=200;
				snprintf(msg,msgSize,"A Traction edge intersects more than %d elements",maxFoundElements);
				throw CommonException(msg,"MatPoint2D::GetExactTractionInfo");
			}
			
			// store element and dimensioned end points
			cElem[foundElement] = iel;
			corners[foundCorner] = MakeVector(c1.x + t1*p[1],c1.y + t1*p[3],0.);
			Vector x2 = MakeVector(c1.x + t2*p[1],c1.y + t2*p[3],0.);
			corners[foundCorner+1] = x2;
			SubVector(&x2,&corners[foundCorner]);
			
			// Store weighting factor
			// B vec n |x2-x1| / 4 for planar and vec n |x2-x1| / 4 for axisymmetric
			// Calling code will need to add the B factor if not axisymmetric
			CopyScaleVector(&tscaled[foundElement],&tnorm,sqrt(DotVectors2D(&x2,&x2)));
			
			// increment number of found elements
			foundElement++;
			foundCorner+=2;
		}
	}
	
	// store number of elements found
	*numDnds = foundElement;
}

// Given face, find dimensionless coordinates of corners (2 in 2D), elements for those corners (2 in 2D),
//     and deformed particle radii (2 in 2D)
// Return number of corners in numDnds
// Function result is vector in BC direction * face weight (which is B*|r1| = Area/2 in 2D)
// On input cElem, corners, and radii are areas of length equal to number of corners (2 in 2D)
Vector MatPoint2D::GetSurfaceInfo(int face,int dof,int *cElem,Vector *corners,Vector *radii,int *numDnds,double *redge)
{
	*numDnds = 2;
	Vector c1, c2, c3;		// edge from 1 to 2, other direction from 1 to 3
	Vector r1,r2;
	
	// get polygon vectors - these are from particle to edge
	GetSemiSideVectors(&r1,&r2,NULL);
	double uGIMPsize = -1.;
#ifndef TRACTION_ALWAYS_DEFORMED
    // this option uses undeformed edge for all GIMP except finite GIMP
    if(!ElementBase::UsingCPDIMethod() && ElementBase::useGimp!=FINITE_GIMP)
    {	double rx,ry;
        switch(face)
        {	case 1:
            case 3:
                uGIMPsize = sqrt(DotVectors(&r1,&r1));
                break;
                
            case 2:
            default:
                uGIMPsize = sqrt(DotVectors(&r2,&r2));
                break;
        }
        GetUndeformedSemiSides(&rx,&ry,NULL);
        ZeroVector(&r1);
        ZeroVector(&r2);
        r1.x = rx;
        r2.y = ry;
    }
#endif
	
	// handle axisymmetric
	if(fmobj->IsAxisymmetric())
	{	// truncate by shrinking domain if any have x < 0, but keep particle in the middle
		if(pos.x-fabs(r1.x+r2.x)<0.)
		{	// make pos.x-shrink*fabs(r1.x+r2.x) = eps (small and positive)
			double shrink = (pos.x-theElements[ElemID()]->GetDeltaX()*1.e-10)/fabs(r1.x+r2.x);
			r1.x *= shrink;
			r1.y *= shrink;
			r2.x *= shrink;
			r2.y *= shrink;
		}
	}
	
	// Find positions along edge from c1 o c2
	// The third point is node before c1
	switch(face)
	{	case 1:
			// lower edge nodes 1 to 2
			c1.x = pos.x-r1.x-r2.x;
			c1.y = pos.y-r1.y-r2.y;
			c2.x = pos.x+r1.x-r2.x;
			c2.y = pos.y+r1.y-r2.y;
			// third is node 4
			c3.x = pos.x-r1.x+r2.x;
			c3.y = pos.y-r1.y+r2.y;
			break;
			
		case 2:
			// right edge nodes 2 to 3
			c1.x = pos.x+r1.x-r2.x;
			c1.y = pos.y+r1.y-r2.y;
			c2.x = pos.x+r1.x+r2.x;
			c2.y = pos.y+r1.y+r2.y;
			// third is node 1
			c3.x = pos.x-r1.x-r2.x;
			c3.y = pos.y-r1.y-r2.y;
			break;
			
		case 3:
			// top edge nodes 3 to 4
			c1.x = pos.x+r1.x+r2.x;
			c1.y = pos.y+r1.y+r2.y;
			c2.x = pos.x-r1.x+r2.x;
			c2.y = pos.y-r1.y+r2.y;
			// third is node 2
			c3.x = pos.x+r1.x-r2.x;
			c3.y = pos.y+r1.y-r2.y;
			break;
			
		default:
			// left edge nodes 4 to 1
			c1.x = pos.x-r1.x+r2.x;
			c1.y = pos.y-r1.y+r2.y;
			c2.x = pos.x-r1.x-r2.x;
			c2.y = pos.y-r1.y-r2.y;
			// third is node 3
			c3.x = pos.x+r1.x+r2.x;
			c3.y = pos.y+r1.y+r2.y;
			break;
	}
	c1.z=c2.z=c3.z=0.;
	
	// get elements and xipos for the two corners
	try
	{	cElem[0] = mpmgrid.FindElementFromPoint(&c1,this)-1;
		theElements[cElem[0]]->GetXiPos(&c1,&corners[0]);
		
		cElem[1] = mpmgrid.FindElementFromPoint(&c2,this)-1;
		theElements[cElem[1]]->GetXiPos(&c2,&corners[1]);
	}
	catch(CommonException& err)
	{   char msg[200];
		size_t msgSize=200;
		snprintf(msg,msgSize,"A Traction edge node has left the grid: %s",err.Message());
		throw CommonException(msg,"MatPoint2D::GetSurfaceInfo");
	}
		
	// get relevant edge radii for this edge
	radii[0].x = 0.5*(c2.x-c1.x);
	radii[0].y = 0.5*(c2.y-c1.y);
	radii[1].x = 0.5*(c3.x-c1.x);
	radii[1].y = 0.5*(c3.y-c1.y);
	
	// get face weighting
 	// if planar then 1/2 the face area = |r|B
	// if axisymmetric then |r| and load redge
	double faceWt = uGIMPsize<0. ? sqrt(DotVectors(&radii[0],&radii[0])) : uGIMPsize ;
	if(fmobj->IsAxisymmetric())
		*redge = pos.x - radii[1].x;
	else
		faceWt *= thickness();
	
	// get traction normal vector multiplied by 1/2 the face area = |r|B
	double enorm,ex,ey;
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
			// cross product of edge vector with (0,0,1) = (r1y, -r1x)
			ex = radii[0].y;
			ey = -radii[0].x;
			enorm = sqrt(ex*ex+ey*ey);
			wtNorm.x = ex*faceWt/enorm;
			wtNorm.y = ey*faceWt/enorm;
			break;
		case T1_DIRECTION:
			// load in direction specified by normalized (ex,ey)
			ex = radii[0].x;
			ey = radii[0].y;
			enorm = sqrt(ex*ex+ey*ey);
			wtNorm.x = ex*faceWt/enorm;
			wtNorm.y = ey*faceWt/enorm;
			break;
		default:
			// normal is z direction (not used in 2D or axisymmetric)
			wtNorm.z = faceWt;
			break;
	}
	
	// Return hat n * Area/2 as vector
	return wtNorm;
}

// Get Rotation matrix for initial material orientation (anistropic only)
Matrix3 MatPoint2D::GetInitialRotation(void)
{	double theta = GetAnglez0InRadians();
	double cs = cos(theta);
	double sn = sin(theta);
	return Matrix3(cs,sn,-sn,cs,1.);
}

