/********************************************************************************
    MPMBase.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "MPM_Classes/MPMBase.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Boundary_Conditions/MatPtHeatFluxBC.hpp"
#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "Materials/MaterialBase.hpp"
#include "Materials/RigidMaterial.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/UnitsController.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Global_Quantities/BodyForce.hpp"

// globals
MPMBase **mpm;		// list of material points
int nmpmsNR=0;		// number for last non-rigid material point
// MODEL_RIGID_BODIES: when not defined, nmpmsRB will always equal to nmpmsNR
int nmpmsRB=0;		// number for last rigid block material
int nmpmsRC=0;		// number for last rigid contact material
int nmpms=0;		// number of material points

#pragma mark MPMBase::Constructors and Destructors

// Constructors
MPMBase::MPMBase()
{
}

// throws std::bad_alloc
MPMBase::MPMBase(int elem,int theMatl,double angin)
{
    inElem=elem;
	mp=-1.;						// calculated in PreliminaryParticleCalcs, unless set in input file
    matnum=theMatl;
	SetAnglez0InDegrees(angin);
	SetAngley0InDegrees(0.0);
	SetAnglex0InDegrees(0.0);
    
	// space to hold velocity fields
	// these only change if there are cracks
	vfld = new char[maxShapeNodes];
	
	// for J integral (on non-rigid only)
	// allocated by CalcJKTask if activated
	velGrad=NULL;
	
	// material history data (as needed by material class)
	matData=NULL;
	
	// custom taks history
	customMatData=NULL;
	
	// diffusion gradients (c units/mm) (if needed and zeroed each time step)
	// solvent or poroelasticity task always uses pDiff[0], extra tasks will vary
	pDiff=NULL;
	
	// concentration (c units)
	SetConcentration(0.,false);
    
	// temperature (degrees)
	SetTemperature(0.,0.);
	
	// temperature gradient (degrees/mm) (if needed and zeroed each time step)
	pTemp=NULL;

    // CPDI and GIMP data (recalculated each time step)
    cpdi_or_gimp = NULL;
    
    // This is set when using flux BCs in GetCPDINodesAndWeights()
    // But it not appear to be used by any code?
    faceArea = NULL;

	// rotation matrix (when tracked)
	Rtot = NULL;
	
	// set all to zero
	ResetMaterialPoint();
}

// Reset material point when deleted
// or initialize when created
void MPMBase::ResetMaterialPoint(void)
{
	for(int i=1;i<maxShapeNodes;i++)
		vfld[i]=NO_CRACK;
	
	// zero stresses and strains (including deformation gradient)
	ZeroTensor(&sp);
	pressure = 0.;
	ZeroTensor(&ep);
	ZeroTensorAntisym(&wrot);
	ZeroTensor(&eplast);
	ZeroVector(&acc);
	
	// zero increment initial residual strains
	dTrans.dT = 0.;
	dTrans.dC = 0.;
	dTrans.doopse = 0.;		// for generalized plane stress or strain
	oopIncrement = 0.;		// out-of-plane increment (stress or strain)
	
	// zero energies and buffers
	plastEnergy=0.;
	prev_dTad=0.;
	buffer_dTad=0.;
	workEnergy=0.;
	resEnergy=0.;
	heatEnergy=0.;
	entropy=0;

	// PS - when point created, velocity and position and ext force should be set too
	ZeroVector(&vel);
	
	// counts crossing and sign is whether or not left the grid
	elementCrossings=0;
}

// allocation velGrad tensor data if need in this calculations (non rigid only)
// throws std::bad_alloc
void MPMBase::AllocateJStructures(void)
{	velGrad=new Tensor;
}

// allocate structures when needed for particle CPDI or GIMP domains
// throws std::bad_alloc
bool MPMBase::AllocateCPDIorGIMPStructures(int gimpType,bool isThreeD)
{
    int cpdiSize=0;
    
	// If CPDI find size of CPDI datastructuress, otherwise pass on to GIMP structur allocation
    if(gimpType==LINEAR_CPDI || gimpType==LINEAR_CPDI_AS || gimpType==BSPLINE_CPDI)
        cpdiSize = isThreeD ? 8 : 4 ;
    else if(gimpType==QUADRATIC_CPDI)
        cpdiSize = 9;
	else if(gimpType==FINITE_GIMP)
		return AllocateFiniteGIMPStructures(isThreeD);
	else
		return true;
	
    // create memory for cpdiSize pointers
    CPDIDomain **cpdi = new (nothrow) CPDIDomain *[cpdiSize];
    if(cpdi == NULL) return false;
	
    // create each one
    int i;
    for(i=0;i<cpdiSize;i++)
    {   cpdi[i] = new CPDIDomain;
        cpdi[i]->wg.z = 0.;             // set zero once for 2D calculations
        cpdi[i]->ncpos.z = 0.;          // set zero once for 2D calculations
		
		// weights constant except for axisymmetric CPDI
		if(gimpType==LINEAR_CPDI || gimpType==BSPLINE_CPDI)
			cpdi[i]->ws = isThreeD ? 0.125 : 0.25 ;
		else if(gimpType==QUADRATIC_CPDI)
		{	if(i<4)
				cpdi[i]->ws = 1./36.;
			else if(i<8)
				cpdi[i]->ws = 1./9.;
			else
				cpdi[i]->ws = 4./9.;
		}
    }
    
    // save face areas (or lengths in 2D) if traction, heat flux,
    // or any type of diffusion flux BC is being used
    if(firstTractionPt!=NULL || firstFluxPt!=NULL || firstHeatFluxPt!=NULL)
		faceArea = new Vector;
	
	// load to generic variable
	cpdi_or_gimp = (char *)cpdi;
    
    return true;
}

// Destructor (and it is virtual)
MPMBase::~MPMBase() { }

#pragma mark MPMBase::Methods

// hold or reverse the direction (should only be done for rigid material particles)
// holdFirst == true, store velocity in acc and zero the velocity
// holdFirst == false, if holding, reverse using stored velocity otherwise just reverse
void MPMBase::ReverseParticle(bool holdFirst,bool holding)
{
	if(holdFirst)
    {   acc = vel;
        ZeroVector(&vel);
    }
    else if(holding)
    {   CopyScaleVector(&vel,&acc,-1.);
        ZeroVector(&acc);
    }
    else
    {	vel.x=-vel.x;
        vel.y=-vel.y;
        vel.z=-vel.z;
    }
}

// stop the particle (should only be done for rigid material particles)
void MPMBase::StopParticle(void)
{	vel.x=0.;
    vel.y=0.;
	vel.z=0.;
}

// Subclass must override to support exact tractions
void MPMBase::GetExactTractionInfo(int face,int dof,int *cElem,Vector *corners,Vector *tscaled,int *numDnds) const {}

#pragma mark FINITE GIMP SUPPORT

// allocate structures when needed for finite GIMP calculations
// throws std::bad_alloc
bool MPMBase::AllocateFiniteGIMPStructures(bool isThreeD)
{
    // create memory for cpdiSize pointers
    FiniteGIMPInfo *finfo = new FiniteGIMPInfo;
    if(finfo == NULL) return false;
    
    // save for intersected element
    //        InTheseElements: [0] for total number, [1] to [maxElementIntersections] for elements
    //        Rest: [0] to [maxElementIntersections-1] for element data
    // this is not final number of nodes with non-zero shape function, but number of
    // elements that might intersect a particle
    finfo->InTheseElements =  new int [maxElementIntersections+1];
    finfo->moment_0 =  new double [maxElementIntersections];
    finfo->moment_x =  new double [maxElementIntersections];
    finfo->moment_y =  new double [maxElementIntersections];
    finfo->moment_xy =  new double [maxElementIntersections];

    // load to generic variable
    cpdi_or_gimp = (char *)finfo;
    
    return true;
}

// Subclass must override to support finite GIMP
void MPMBase::GetFiniteGIMP_Integrals(void) {}

// Get finite gimp info accessor
FiniteGIMPInfo *MPMBase::GetFiniteGIMPInfo(void) { return (FiniteGIMPInfo *)cpdi_or_gimp; }

#pragma mark SPIN MOMENTUM FEATURES

// Get spatial velocity gradient on the particles
// For future version that might track velocity gradients on particles
Matrix3 MPMBase::GetParticleGradVp(bool spatial)
{
	// this creates zerod matrix in contructor
	Matrix3 spatialGradVp;
	
	// return the matrix
	return spatialGradVp;
}

// Get angular momentum on the particles
#ifndef OSPARTICULAS
// For future version that might track velocity gradients on particles
#endif
Vector MPMBase::GetParticleAngMomentum(void)
{
	Vector partLp = MakeVector(0.,0.,0.);
	return partLp;
}

#pragma mark MPMBase::Accessors

// scale residual strains for current update method
// secondPass implies in USAVG_METHOD mode, therefore simply set to dTrans unless using USAVG+/-
ResidualStrains MPMBase::ScaledResidualStrains(int secondPass)
{
	ResidualStrains res = dTrans;
	if(secondPass)
	{	res.dT *= (1.0-fractionUSF);
		res.dC *= (1.0-fractionUSF);
		res.doopse *= (1.0-fractionUSF);
	}
	else if(fmobj->mpmApproach==USAVG_METHOD)
	{	res.dT *= fractionUSF;
		res.dC *= fractionUSF;
		res.doopse *= fractionUSF;
	}
	return res;
}

// set mass in PreliminaryParticleCalcs, but only if input file did not set it first
// When tracking spin, make sure particle momentum is initialized
// throws std::bad_alloc
void MPMBase::InitializeMass(double rho,double volPerParticle,bool trackSpin)
{	if(mp<0.) mp = rho*volPerParticle;
}

// set or average velocity gradient (2D) (only needed for second J Term)
// for USF and SZS, sets gradient on its only pass
// for USAVG, stores in on first pass and then averages on second (J calc follows the second pass)
void MPMBase::SetVelocityGradient(double dvxx,double dvyy,double dvxy,double dvyx,int secondPass)
{
    if(JTerms!=2) return;
    if(!secondPass)
    {	velGrad->xx = dvxx;
        velGrad->yy = dvyy;
        velGrad->xy = dvxy;
        velGrad->zz = dvyx;						// yx stored in zz
    }
    else
    {	velGrad->xx = (velGrad->xx+dvxx)*0.5;
        velGrad->yy = (velGrad->yy+dvyy)*0.5;
        velGrad->xy = (velGrad->xy+dvxy)*0.5;
        velGrad->zz = (velGrad->zz+dvyx)*0.5;		// yx stored in zz
    }
}

// Set concentration on the particle. This is called
//  1. When material point is created (with pConc=0). This call will count diffusion
//      tasks (once) and create the pDiff list for paticle data.
//  2. When material point added to the MpsController. Now has initial concentration
//      which might differ from zero. Sets pDiff[0] only
//  3. When particle deleted, reservoir will reset to provided concentration and
//      set all other diffusion to zero (will need to change this if ever
//      have other diffusion tasks that needs non-zero reset)
void MPMBase::SetConcentration(double pConc,bool fromReservoir)
{	// Make sure diffusion counted, then exit if none
	DiffusionTask::CountDiffusionTasks();
	if(numDiffusion<=0) return;
	
	// create all pDiff if needed on the material point
	if(pDiff==NULL)
	{	pDiff = new DiffusionInfo *[numDiffusion];
		for(int i=0;i<numDiffusion;i++)
		{	pDiff[i] = new DiffusionInfo;
			pDiff[i]->grad = MakeVector(0.,0.,0.);
			pDiff[i]->bufferSource = 0.;
		}
	}
	
	// we only set moisture values now (always in [0] if there)
	// or it might
    int diffNum = 0;
	if(diffusion!=NULL)
	{	pDiff[0]->conc = pConc;
		pDiff[0]->prevConc = fmax(diffusion->reference,0.);
        diffNum = 1;
		
	}
    
    // particle deletion should zero the rest too (only when called by reservoir)
    if(numDiffusion>diffNum && fromReservoir)
    {   for(int i=diffNum;i<numDiffusion;i++)
        {   pDiff[i]->conc = 0.;
            pDiff[i]->prevConc = 0.;
        }
    }
}

// set temperature (only and for all particles, during initialization)
void MPMBase::SetTemperature(double pTempSet,double pRefTemp)
{	pTemperature=pTempSet;
	pPreviousTemperature=pRefTemp;
}

// zero the temperature gradient (non-rigid particles only)
void MPMBase::AddTemperatureGradient(int offset)
{   pTemp[offset]=0.;
    pTemp[offset+1]=0.;
	pTemp[offset+2]=0.;
}

// zero the concentration gradient
void MPMBase::AddConcentrationGradient(int dnum)
{	pDiff[dnum]->grad.x=0.;
    pDiff[dnum]->grad.y=0.;
    pDiff[dnum]->grad.z=0.;
}

// material ID (convert to zero based)
int MPMBase::MatID(void) const { return matnum-1; }			// zero based material array in data storage
int MPMBase::ArchiveMatID(void) const { return matnum; }		// one based for archiving

// element ID (convert to zero based)
bool MPMBase::InReservoir(void) const { return inElem==1; }
int MPMBase::ElemID(void) const { return inElem-1; }				 // zero based element array in data storage
void MPMBase::ChangeElemID(int newElem,bool adjust)
{	// adjust dimensionless size if needed
	if(adjust)
	{	Vector prev = theElements[ElemID()]->GetDeltaBox();
		Vector next = theElements[newElem]->GetDeltaBox();
		mpm_lp.x *= prev.x/next.x;
		mpm_lp.y *= prev.y/next.y;
		mpm_lp.z *= prev.z/next.z;			// always 1 in 2D
	}
	
	// reset element
	inElem=newElem+1;                   // set using zero basis
	IncrementElementCrossings();		// count crossing
}
int MPMBase::ArchiveElemID(void) { return inElem; }			// one based for archiving

// return current element crossings for archiving and reset to zero
int MPMBase::GetElementCrossings(void) { return elementCrossings>=0 ? elementCrossings : -elementCrossings; }
void MPMBase::SetElementCrossings(int ec) { elementCrossings = ec; }
void MPMBase::IncrementElementCrossings(void)
{	elementCrossings = elementCrossings>=0 ? elementCrossings+1 : elementCrossings-1;
}
// negative element crossings means the particle previously left the grid
bool MPMBase::HasLeftTheGridBefore(void) { return elementCrossings<0; }
void MPMBase::SetHasLeftTheGridBefore(bool setting) { elementCrossings = setting ? -GetElementCrossings() : GetElementCrossings(); }

// get unscaled volume for use only in contact and imperfect interface calculations
// return result in mm^3
double MPMBase::GetUnscaledVolume(void)
{	return mp/GetRho();
}

// Get membrane particle size as fraction of cell size in each direction
// Called to archive membrane geometry
void MPMBase::GetInitialSize(Vector &lp) const { lp = mpm_lp; }

// Set particle size as fraction of cell size in each direction in current element
// It's radius is -1 to 1 natural coordinates
void MPMBase::SetDimensionlessSize(Vector *lp)
{	mpm_lp = *lp;
	mpmgrid.TrackMinParticleSize(GetParticleSize());
}
void MPMBase::SetDimensionlessByPts(int pointsPerSide)
{	// linear number of points
	double lp = 1./(double)pointsPerSide;
	if(fmobj->IsThreeD())
		mpm_lp = MakeVector(lp,lp,lp);
	else
		mpm_lp = MakeVector(lp,lp,1.);
	mpmgrid.TrackMinParticleSize(GetParticleSize());
}

// Get particle size as fraction of cell size in each direction
// Called by GIMP shape function code
void MPMBase::GetDimensionlessSize(Vector &lp) const { lp = mpm_lp; }

// Particle semi-size in actual units (3D overrides to add z component)
// note that in 2D, part.z will be 1
Vector MPMBase::GetParticleSize(void) const
{	Vector part = theElements[inElem-1]->GetDeltaBox();
	part.x *= 0.5*mpm_lp.x;
	part.y *= 0.5*mpm_lp.y;
	return part;
}
double MPMBase::GetParticleXSize(void) const { return 0.5*mpm_lp.x*theElements[inElem-1]->GetDeltaX(); }
double MPMBase::GetParticleYSize(void) const { return 0.5*mpm_lp.y*theElements[inElem-1]->GetDeltaY(); }
double MPMBase::GetParticleZSize(void) const { return 0.5*mpm_lp.z*theElements[inElem-1]->GetDeltaZ(); }

// get shortest side (3D overrides to check z component) 
double MPMBase::GetMinParticleLength(void) const
{	Vector part = GetParticleSize();
	double minPart = part.x;
	if (part.y < minPart) minPart = part.y;
	return 2.*minPart;
}

// get rotation angle for 2D calculations
// Use polar decomposition to get sin(theta) and cos(theta) for ccw rotation
//		from initial material position to current position
void MPMBase::Get2DSinCos(double *sintheta,double *costheta)
{	// particle rotation
	double F11plusF22 = 2. + ep.xx + ep.yy;
	double denom = 1./sqrt(F11plusF22*F11plusF22+wrot.xy*wrot.xy);
	double sinthetap = wrot.xy*denom;
	double costhetap = F11plusF22*denom;
	
	// total rotation
	// sin(anglez0-thetap) = cos(anglez0)*sin(thetap) - sin(anglez0)*cos(thetap)
	// cos(anglez0-thetap) = cos(anglez0)*cos(thetap) + sin(snglez0)*sin(thetap)
	double c0 = cos(anglez0);
	double s0 = sin(anglez0);
	*sintheta = c0*sinthetap - s0*costhetap;
	*costheta = c0*costhetap + s0*sinthetap;
}

// Rotation matrix (for materials that track it)
// Currently only tracked in in Transisotropic (and subclasses) when in 3D and then
//    it is rotation all the way to material axes
// Tracking activated by call InitRtot in SetInitialParticleState()
// Rtot is set when
// 1. Small Rotation - when call FillElasticProperties3D() just before constitutive law
//		it will have Rnm1 or rotation before F undated
// 2. Large Rotation - when call constitutive law - it stores Rn
// It is called for use
// 1. Anisotropic plasticity
// 2. GetTransportProps() in Tranisotropic (and subclasses)
Matrix3 *MPMBase::GetRtotPtr(void) { return Rtot; }
void MPMBase::SetRtot(Matrix3 newR) { Rtot->set(newR); }

// throws std::bad_alloc
void MPMBase::InitRtot(Matrix3 newR)
{	if(Rtot==NULL) Rtot = new Matrix3();
	Rtot->set(newR);
}

// decompose current F to get rotation
// apply initial rotation to account for aniostropy
Matrix3 MPMBase::GetRtot(void)
{	
	// get deformation gradient
	Matrix3 pF = GetDeformationGradientMatrix();
	
	// Decompose for R (U is not needed)
    Matrix3 R;
	//Matrix3 U = pF.RightDecompose(&R,NULL);
    pF.RightDecompose(&R,NULL);
	
	// apply initial rotation total rotation from initial material axes
	Matrix3 R0 = GetInitialRotation();
	Matrix3 Rtot = R*R0;
	
	// return it
	return Rtot;
}

// anglez0 is initial z cw orientation angle from global to material (2D and 3D z,y,x scheme)
// In 2D and small rotation, 0.5*wrot.xy is ccw from initial to current axes, thus
//		anglez0-0.5*wrot.xy is cw from current to material or ccw from material to current
double MPMBase::GetRotationZ(void) { return anglez0-0.5*wrot.xy; }
double MPMBase::GetParticleRotationZ(void) { return -0.5*wrot.xy; }		// rotation in simulation, ignoring initial angle
double MPMBase::GetRotationZInDegrees(void) { return 180.*(anglez0-0.5*wrot.xy)/PI_CONSTANT; }
void MPMBase::SetAnglez0(double angle) { anglez0=angle; }
void MPMBase::SetAnglez0InDegrees(double angle) { anglez0=PI_CONSTANT*angle/180.; }
double MPMBase::GetAnglez0InDegrees(void) { return 180.*anglez0/PI_CONSTANT; }
double MPMBase::GetAnglez0InRadians(void) { return anglez0; }

// angley0 is initial y cw orientation angle (3D z,y,x scheme)
double MPMBase::GetRotationY(void) { return angley0+0.5*wrot.xz; }
double MPMBase::GetParticleRotationY(void) { return +0.5*wrot.xz; }		// rotation in simulation, ignoring initial angle
double MPMBase::GetRotationYInDegrees(void) { return 180.*(angley0+0.5*wrot.xz)/PI_CONSTANT; }
void MPMBase::SetAngley0(double angle) { angley0=angle; }
void MPMBase::SetAngley0InDegrees(double angle) { angley0=PI_CONSTANT*angle/180.; }
double MPMBase::GetAngley0InDegrees(void) { return 180.*angley0/PI_CONSTANT; }
double MPMBase::GetAngley0InRadians(void) { return angley0; }

// anglex0 is initial x cw orientation angle (3D z,y,x scheme)
double MPMBase::GetRotationX(void) { return anglex0-0.5*wrot.yz; }
double MPMBase::GetParticleRotationX(void) { return -0.5*wrot.yz; }		// rotation in simulation, ignoring initial angle
double MPMBase::GetRotationXInDegrees(void) { return 180.*(anglex0-0.5*wrot.yz)/PI_CONSTANT; }
void MPMBase::SetAnglex0(double angle) { anglex0=angle; }
void MPMBase::SetAnglex0InDegrees(double angle) { anglex0=PI_CONSTANT*angle/180.; }
double MPMBase::GetAnglex0InDegrees(void) { return 180.*anglex0/PI_CONSTANT; }
double MPMBase::GetAnglex0InRadians(void) { return anglex0; }

// 2D and 3D rotational strain tensor
void MPMBase::IncrementRotationStrain(double rotXY) { wrot.xy+=rotXY; }
void MPMBase::IncrementRotationStrain(double rotXY,double rotXZ,double rotYZ)
{	wrot.xy+=rotXY;
	wrot.xz+=rotXZ;
	wrot.yz+=rotYZ;
}
TensorAntisym *MPMBase::GetRotationStrainTensor(void) { return &wrot; }

// For 2D or axisymmetric, increment Fzz to Fzz(n) = (1+dezz)*Fzz(n-1)
// In the strain tensor, we store Fzz(n)-1 and Fzz(n-1) = 1 + ep.zz
//		Fzz(n)-1 = (1+dezz)*(1 + ep.zz)-1 = ep.zz + dezz*(1+ep.zz)
void MPMBase::IncrementDeformationGradientZZ(double dezz)
{	ep.zz += dezz*(1.+ep.zz);
}

// Get V-I which is Biot strain rotated into current particle orientation
// It is based on total strain
Matrix3 MPMBase::GetBiotStrain(void) const
{
	Matrix3 F = GetDeformationGradientMatrix();
	Matrix3 V = F.LeftDecompose(NULL,NULL);
	V(0,0) -= 1.;
	V(1,1) -= 1.;
	V(2,2) -= 1.;
	return V;
}

// 2D and 3D material points override to get initial rotation matrix
Matrix3 MPMBase::GetInitialRotation(void)
{	return Matrix3(1.,0.,0.,1.,1.);
}

// energies per unit mass
double MPMBase::GetPlastEnergy(void) { return plastEnergy; }
void MPMBase::AddPlastEnergy(double energyInc) { plastEnergy+=energyInc; }

// return buffer_dTad and clear it too
double MPMBase::GetClear_dTad(void)
{	double dTad = buffer_dTad;
    buffer_dTad = 0.;
    return dTad;
}

// Only called in IncrementHeatEnergy()
// a material should never call this directly
void MPMBase::Add_dTad(double dTInc) { buffer_dTad+=dTInc; }

// buffer diffusion source term (meaning depends on the task)
void MPMBase::AddParticleDiffusionSource(int diffNumber,double dpud)
{	pDiff[diffNumber]->bufferSource+=dpud;
}

// return diffusion buffer source and clear it too
// Nver called in transport only mode
bool MPMBase::GetClearParticleDiffusionSource(int diffNumber,double &dpud)
{	dpud = pDiff[diffNumber]->bufferSource;
	pDiff[diffNumber]->bufferSource = 0.;
	return dpud!=0.;
}

double MPMBase::GetWorkEnergy(void) { return workEnergy; }
void MPMBase::SetWorkEnergy(double energyTot) { workEnergy=energyTot; }
void MPMBase::AddWorkEnergy(double energyInc) { workEnergy+=energyInc; }

void MPMBase::AddWorkEnergyAndResidualEnergy(double energyInc,double resInc)
{	workEnergy+=energyInc;
	resEnergy+=resInc;
}

double MPMBase::GetResidualEnergy(void) { return resEnergy; }
void MPMBase::SetResidualEnergy(double resInc) { resEnergy=resInc; }
void MPMBase::AddResidualEnergy(double resInc) { resEnergy+=resInc; }

double MPMBase::GetStrainEnergy(void) { return workEnergy - resEnergy; }

double MPMBase::GetHeatEnergy(void) { return heatEnergy; }
void MPMBase::SetHeatEnergy(double energyTot) { heatEnergy = energyTot; }
void MPMBase::AddHeatEnergy(double energyInc) { heatEnergy += energyInc; }

void MPMBase::SetEntropy(double entropyTot) { entropy = entropyTot; }
double MPMBase::GetEntropy(void) { return entropy; }
// isothermal entropy change due to heat dq and constant temperature
void MPMBase::AddEntropy(double dq,double Tgp)
{	entropy += dq/Tgp;
}
// adiabatic entropy change from T1 to T2 at constant heat capacity Cv
void MPMBase::AddEntropy(double Cv,double T1,double T2)
{	entropy += Cv*log(T2/T1);
}


double MPMBase::GetInternalEnergy(void) { return workEnergy + heatEnergy; }

// pointers to variables
Vector *MPMBase::GetPFext(void) { return &pFext; }
Vector *MPMBase::GetNcpos(void) { return &ncpos; }
CPDIDomain **MPMBase::GetCPDIInfo(void) { return (CPDIDomain **)cpdi_or_gimp; }
Vector *MPMBase::GetAcc(void) { return &acc; }
Tensor *MPMBase::GetVelGrad(void) { return velGrad; }
Tensor *MPMBase::GetStrainTensor(void) { return &ep; }

// material classes only should call GetStressTensor and can change the stress
// all others should use ReadStressTensor() because that gives materials a chance to
//  customize the way stresses are stored, such as to separate pressure and deviatoric stress
Tensor *MPMBase::GetStressTensor(void) { return &sp; }
Tensor MPMBase::ReadStressTensor(void) { return theMaterials[MatID()]->GetStress(&sp,pressure,this); }
void MPMBase::StoreStressTensor(Tensor *newsp) { theMaterials[MatID()]->SetStress(newsp,this); }
void MPMBase::StoreThicknessStressIncrement(double dszz)
{	theMaterials[MatID()]->IncrementThicknessStress(dszz,this);
	oopIncrement = dszz;
}
void MPMBase::StoreThicknessStrainIncrement(double dezz)
{	IncrementDeformationGradientZZ(dezz);
	oopIncrement = dezz;
}
void MPMBase::IncrementPressure(double dP) { pressure += dP; }
void MPMBase::SetPressure(double P) { pressure = P; }
double MPMBase::GetPressure(void) { return pressure; }

// These methods return pointer to the same symmetric tensor. A material class must
// decide how it will use the eplast tensor variable and use only that one purpose.
// The only reason for two names is to make the code more readable
Tensor *MPMBase::GetAltStrainTensor(void) { return &eplast; }

// history data (offset is in bytes and not in number of doubles)
char *MPMBase::GetHistoryPtr(int offset) { return matData+offset; }
void MPMBase::SetHistoryPtr(char *createMatData)
{	if(matData!=NULL) delete [] matData;
	matData=createMatData;
}
double MPMBase::GetHistoryDble(int index,int offset)
{	// assumes matData is array of doubles and gets one at index (0 based)
	// offset allows relocatable history and only used by materials with child materials
	double *h=(double *)(matData+offset);
	return h[index];
}
void MPMBase::SetHistoryDble(int index,double history,int offset)
{	// assumes matData is array of doubles and sets one at index (0 based)
	// offset allows relocatable history and only used by materials with child materials
	double *h=(double *)(matData+offset);
	h[index]=history;
}

// custom history data
char *MPMBase::GetCustomHistoryPtr(void) { return customMatData; }
void MPMBase::SetCustomHistoryPtr(char *createMatData)
{	if(customMatData!=NULL) delete [] customMatData;
	customMatData=createMatData;
}
double MPMBase::GetCustomHistoryDble(int index)
{	// assumes customMatData is array of doubles and gets one at index (0 based)
	// call responsible to make sure valid index
	double *h=(double *)(customMatData+index);
	return *h;
}
void MPMBase::SetCustomHistoryDble(int index,double history)
{	// assumes customMatData is array of doubles and sets one at index (0 based)
	// caller responsible to make sure valid index
	double *h=(double *)(customMatData+index);
	*h=history;
}

// get density from the material
double MPMBase::GetRho(void)
{	return theMaterials[matnum-1]->GetRho(this);
}

// get saturation concentration from the material
double MPMBase::GetConcSaturation()
{	return theMaterials[matnum-1]->GetMaterialConcentrationSaturation(this);
}

// set relative toughness
void MPMBase::SetRelativeStrength(double rel)
{   theMaterials[matnum-1]->SetRelativeStrength(this,rel);
}

// set relativee toughness
void MPMBase::SetRelativeToughness(double rel)
{   theMaterials[matnum-1]->SetRelativeToughness(this,rel);
}

// Decribe for debugging use or output on some errors
void MPMBase::Describe(void)
{   cout << "# Step=" << fmobj->mstep;
    cout << " pt: #=" << num << " pos=(" << pos.x << "," << pos.y << "," << pos.z << ") mass=" << mp <<
                " matl=" << matnum << " elem=" << inElem << endl;
    cout << "#     vel=(" << vel.x << "," << vel.y << "," << vel.z << ") " << UnitsController::Label(CUVELOCITY_UNITS);
	cout << " origpos=(" << origpos.x << "," << origpos.y << "," << origpos.z << ")" << endl;
    Matrix3 pF = GetDeformationGradientMatrix();
    cout << "#       F=" << pF << ", |F|=" << pF.determinant() << endl;
    double rho0=GetRho();
    double rho = rho0*UnitsController::Scaling(1.e-6)/theMaterials[MatID()]->GetCurrentRelativeVolume(this,0);
    cout << "#       P(tracked)=" << pressure*rho << " " << UnitsController::Label(PRESSURE_UNITS);
    cout << " P(calc)=" << (-(sp.xx+sp.yy+sp.zz)*rho/3.) << " " << UnitsController::Label(PRESSURE_UNITS) << endl;
    cout << "# sigmaii=(" << sp.xx*rho << "," << sp.yy*rho << "," << sp.zz*rho << ") " << UnitsController::Label(PRESSURE_UNITS) << endl;
    cout << "#   tauij=(" << sp.xy*rho << "," << sp.xz*rho << "," << sp.yz*rho << ") " << UnitsController::Label(PRESSURE_UNITS) << endl;
	cout << "#       T= " << pTemperature << " prev T=" << pPreviousTemperature << endl;
    for(int i=0;i<numDiffusion;i++)
		cout << "  D[" << i << "] c= " << pDiff[i]->conc << " prev c=" << pDiff[i]->prevConc << endl;
    
    // nonzero history data
    char *hist = GetHistoryPtr(0);
    if(hist!=NULL)
    {   cout << "#    nonzero history: ";
        int numRow=0;
        for(int i=1;i<=20;i++)
        {   double histVal = theMaterials[MatID()]->GetHistory(i,hist);
            if(histVal!=0.)
            {   if(numRow==4)
                {   cout << "\n#                     ";
                    numRow = 0;
                }
                cout << i << "=" << histVal << " ";
                numRow++;
            }
        }
        cout << endl;
    }
}

// track index into mpm[] array. Pass zero-based address, but stored as 1-based address
void MPMBase::SetNum(int zeroBasedNumber) { num = zeroBasedNumber+1; }
int MPMBase::GetNum(void) const { return num; }

