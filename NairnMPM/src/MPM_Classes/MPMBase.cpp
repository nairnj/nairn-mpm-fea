/********************************************************************************
    MPMBase.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "MPM_Classes/MPMBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Boundary_Conditions/MatPtHeatFluxBC.hpp"
#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/UnitsController.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"

// globals
MPMBase **mpm;		// list of material points
int nmpms=0;		// number of material points
int nmpmsNR=0;		// number for last non-rigid material point
int nmpmsRC=0;		// number for last rigid contact material

#pragma mark MPMBase::Constructors and Destructors

// Constructors
MPMBase::MPMBase()
{
}

// throws std::bad_alloc
MPMBase::MPMBase(int elem,int theMatl,double angin)
{
    int i;
    
    inElem=elem;
	mp=-1.;						// calculated in PreliminaryCalcs, unless set in input file
    matnum=theMatl;
	SetAnglez0InDegrees(angin);
	SetAngley0InDegrees(0.0);
	SetAnglex0InDegrees(0.0);
    
	// space to hold velocity fields
	// these only change if there are cracks
	vfld = new char[maxShapeNodes];
    for(i=1;i<maxShapeNodes;i++)
        vfld[i]=NO_CRACK;
	
    // zero stresses and strains
	ZeroTensor(&sp);
    pressure = 0.;
	ZeroTensor(&ep);
	ZeroTensor(&eplast);
	ZeroTensorAntisym(&wrot);
	ZeroVector(&acc);
    
    // zero energies
    plastEnergy=0.;
#ifndef NEW_HEAT_METHOD
    prev_dTad=0.;
#endif
    buffer_dTad=0.;
    workEnergy=0.;
    resEnergy=0.;
    heatEnergy=0.;
    entropy=0;
	
	// for J integral if needed (on non-rigid only)
	velGrad=NULL;
	
	// material data is needed
	matData=NULL;
	
	// concentration (c units) and gradient (c units/mm)
	pConcentration=0.;
	pDiffusion=NULL;
    
	// temperature (degrees) and gradient (degrees/mm)
	SetTemperature(0.,0.);
	pTemp=NULL;

    // CPDI and GIMP data
    cpdi_or_gimp = NULL;
    faceArea = NULL;

    // PS - when point created, velocity and position and ext force should be set too
	ZeroVector(&vel);
	
	// counts crossing and sign is whether or not left the grid
	elementCrossings=0;
	
	// rotation matrix (when tracked)
	Rtot = NULL;
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
    
	// If CPDU find size of CPDI datastructuress, other pass on to GIMP structur allocation
    if(gimpType==LINEAR_CPDI || gimpType==LINEAR_CPDI_AS)
        cpdiSize = isThreeD ? 8 : 4 ;
    else if(gimpType==QUADRATIC_CPDI)
        cpdiSize = 9;
    else
	{
#ifdef LOAD_GIMP_INFO
        return AllocateGIMPStructures(gimpType,isThreeD);
#else
		return true;
#endif
	}
    
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
		if(gimpType==LINEAR_CPDI)
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
    
    // save face areas (or lengths in 2D)
    if(firstTractionPt!=NULL || firstFluxPt!=NULL || firstHeatFluxPt!=NULL)
		faceArea = new Vector;
	
	// load to generic variable
	cpdi_or_gimp = (char *)cpdi;
    
    return true;
}

#ifdef LOAD_GIMP_INFO
// allocate structures when needed for particle domains (CPDI or GIMP)
// throws std::bad_alloc
bool MPMBase::AllocateGIMPStructures(int gimpType,bool isThreeD)
{
	// not needed for point GIMP
	if(gimpType==POINT_GIMP) return true;
	
	// maximum nodes
	int maxnds = isThreeD ? 27 : 9 ;
	
	// gimp info
	GIMPNodes *gimp = new GIMPNodes;
    if(gimp == NULL) return false;
	
	// arrays
	gimp->nds = new (nothrow) int[maxnds+1];
	if(gimp->nds == NULL) return false;
	
	gimp->ndIDs = new (nothrow) unsigned char[maxnds];
	if(gimp->ndIDs == NULL) return false;

	// load to generic variable
	cpdi_or_gimp = (char *)gimp;
    
	return true;
}
#endif

// Destructor (and it is virtual)
MPMBase::~MPMBase() { }

#pragma mark MPMBase::Methods

// hold or reverse the direction (should only be done for rigid material particles)
// holdFirst == true, store velocity in acc and zero the velocity
// holdFirst == false, if holding, reverse using stored velocity otherwise just reverse
void MPMBase::ReverseParticle(bool holdFirst,bool holding)
{   if(holdFirst)
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

#pragma mark MPMBase::Accessors

// set mass in PreliminaryCalcs, but only if input file did not set it first
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

// set concentration (only, and for all particles, during initialization)
void MPMBase::SetConcentration(double pConc,double pRefConc)
{	pConcentration=pConc;
	pPreviousConcentration=pRefConc;
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
void MPMBase::AddConcentrationGradient(void)
{	pDiffusion[gGRADx]=0.;
    pDiffusion[gGRADy]=0.;
    pDiffusion[gGRADz]=0.;
}

// material ID (convert to zero based)
int MPMBase::MatID(void) const { return matnum-1; }			// zero based material array in data storage
int MPMBase::ArchiveMatID(void) const { return matnum; }		// one based for archiving

// element ID (convert to zero based)
int MPMBase::ElemID(void) const { return inElem-1; }					// zero based element array in data storage
void MPMBase::ChangeElemID(int newElem,bool adjust)
{	inElem=newElem+1;                   // set using zero basis
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
// It is radius in -1 to 1 natural coordinates
void MPMBase::SetDimensionlessSize(Vector *lp) { mpm_lp = *lp; }
void MPMBase::SetDimensionlessByPts(int pointsPerCell)
{	double lp;
    
	// surface length per particle on edge
	switch(pointsPerCell)
	{	case 1:
			lp = 1.0;
			break;
        case 9:
        case 27:
            // 9 is always 2D and 27 is always 3D
            lp = 1./3.;
            break;
        case 16:
            lp = 0.25;
            break;
        case 25:
            lp = 0.20;
            break;
		case 4:
		case 8:
		default:
			// 4 is always 2D and 8 is always 3D
			lp = 0.5;
			break;
	}
	mpm_lp = MakeVector(lp,lp,lp);
}

// Get particle size as fraction of cell size in each direction
// Called by GIMP shape function code
void MPMBase::GetDimensionlessSize(Vector &lp) const { lp = mpm_lp; }

// Particle semi-size in actual units (3D overrides to add z component)
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
// All hypoelastic materials when in large rotation mode and tracking is activated
// 3D, small rotation, and anisotropic tracks rotation always. It is calculated when fill elastic
//    properties and then stored for later use by anisotropic plastic or transport properties
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

// For 2D or axisymmetric, iuncrement Fzz to Fzz(n) = (1+dezz)*Fzz(n-1)
// In the strain tensor, we store Fzz(n)-1 and Fzz(n-1) = 1 + ep.zz
//		Fzz(n)-1 = (1+dezz)*(1 + ep.zz) = ep.zz + dezz*(1+ep.zz)
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

// energies
double MPMBase::GetPlastEnergy(void) { return plastEnergy; }
void MPMBase::AddPlastEnergy(double energyInc) { plastEnergy+=energyInc; }

#ifdef NEW_HEAT_METHOD

// return buffer_dTad and clear it too
double MPMBase::GetClear_dTad(void)
{	double dTad = buffer_dTad;
    buffer_dTad = 0.;
    return dTad;
}

#else

// Get only in UpdateParticles task to add adiabiatic temperature rise to the particles
// Buffer amount for next heat increment and clear it too
double MPMBase::GetBufferClear_dTad(void)
{	prev_dTad = buffer_dTad;
    buffer_dTad = 0.;
    return prev_dTad;
}

// Only called in IncrementHeatEnergy()
double MPMBase::GetClearPrevious_dTad(void)
{	double dTad = prev_dTad;
    prev_dTad = 0.;
    return dTad;
}

#endif

// Only called in IncrementHeatEnergy()
// a material should never call this directly
void MPMBase::Add_dTad(double dTInc) { buffer_dTad+=dTInc; }

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
void MPMBase::AddEntropy(double entropyInc) { entropy += entropyInc; }

double MPMBase::GetInternalEnergy(void) { return workEnergy + heatEnergy; }

// pointers to variables
Vector *MPMBase::GetPFext(void) { return &pFext; }
Vector *MPMBase::GetNcpos(void) { return &ncpos; }
CPDIDomain **MPMBase::GetCPDIInfo(void) { return (CPDIDomain **)cpdi_or_gimp; }
#ifdef LOAD_GIMP_INFO
GIMPNodes *MPMBase::GetGIMPInfo(void) { return (GIMPNodes *)cpdi_or_gimp; }
#endif
Vector *MPMBase::GetAcc(void) { return &acc; }
Tensor *MPMBase::GetVelGrad(void) { return velGrad; }
Tensor *MPMBase::GetStrainTensor(void) { return &ep; }

// material classes only should call GetStressTensor and can change the stress
// all others should use ReadStressTensor() because that gives materials a chance to
//  customize the way stresses are stored, such as to separate pressure and deviatoric stress
Tensor *MPMBase::GetStressTensor(void) { return &sp; }
Tensor MPMBase::ReadStressTensor(void) { return theMaterials[MatID()]->GetStress(&sp,pressure,this); }
void MPMBase::IncrementPressure(double dP) { pressure += dP; }
void MPMBase::SetPressure(double P) { pressure = P; }
double MPMBase::GetPressure(void) { return pressure; }

// These methods return pointer to the same symmetric tensor. A material class must
// decide how it will use the eplast tensor variable and use only that one purpose.
// The only reason for two names is to make the code more readable
Tensor *MPMBase::GetAltStrainTensor(void) { return &eplast; }

// history data
char *MPMBase::GetHistoryPtr(int offset) { return matData+offset; }
void MPMBase::SetHistoryPtr(char *createMatData)
{	if(matData!=NULL) delete [] matData;
	matData=createMatData;
}
double MPMBase::GetHistoryDble(int index,int offset)
{	// assumes matData is array of doubles and gets one at index (0 based)
	// offset allows relocatable history and only used my materials with child materials
	double *h=(double *)(matData+offset);
	return h[index];
}
void MPMBase::SetHistoryDble(int index,double history,int offset)
{	// assumes matData is array of doubles and sets one at index (0 based)
	// offset allows relocatable history and only used my materials with child materials
	double *h=(double *)(matData+offset);
	h[index]=history;
}

// get density from the material
double MPMBase::GetRho(void)
{	return theMaterials[matnum-1]->GetRho(this);
}

// get density from the material
double MPMBase::GetConcSaturation(void)
{	return theMaterials[matnum-1]->GetConcSaturation(this);
}

// Decribe for debugging use or output on some errors
void MPMBase::Describe(void)
{	cout << "# pt: pos=(" << pos.x << "," << pos.y << "," << pos.z << ") mass=" << mp << 
                " matl=" << matnum << " elem=" << inElem << endl;
    cout << "#     vel=(" << vel.x << "," << vel.y << "," << vel.z << ") " << UnitsController::Label(CUVELOCITY_UNITS) << endl;
    Matrix3 pF = GetDeformationGradientMatrix();
    cout << "#       F=" << pF << ", |F|=" << pF.determinant() << endl;
    double rho0=GetRho();
    double rho = rho0*UnitsController::Scaling(1.e-6)/theMaterials[MatID()]->GetCurrentRelativeVolume(this,0);
    cout << "#       P= " << pressure*rho << " " << UnitsController::Label(PRESSURE_UNITS) << endl;
    cout << "# sigmaii=(" << sp.xx*rho << "," << sp.yy*rho << "," << sp.zz*rho << ") " << UnitsController::Label(PRESSURE_UNITS) << endl;
    cout << "#   tauij=(" << sp.xy*rho << "," << sp.xz*rho << "," << sp.yz*rho << ") " << UnitsController::Label(PRESSURE_UNITS) << endl;
	cout << "#       T= " << pTemperature << " prev T=" << pPreviousTemperature << endl;
}



