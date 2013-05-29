/********************************************************************************
    MatPoint2D.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

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
void MatPoint2D::UpdateStrain(double strainTime,int secondPass,int np,void *props,int matFld)
{
	int i,numnds,nds[maxShapeNodes];
    double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
	Vector vel;
    Matrix3 dv;
	
	// don't need to zero zDeriv for 2D planar because never used in this function
    
	// find shape functions and derviatives
	const ElementBase *elemRef = theElements[ElemID()];
	elemRef->GetShapeGradients(&numnds,fn,nds,xDeriv,yDeriv,zDeriv,this);
    
    // Find strain rates at particle from current grid velocities
	//   and using the velocity field for that particle and each node and the right material
    for(i=1;i<=numnds;i++)
	{	vel = nd[nds[i]]->GetVelocity((short)vfld[i],matFld);
        dv += Matrix3(vel.x*xDeriv[i],vel.x*yDeriv[i],vel.y*xDeriv[i],vel.y*yDeriv[i],0.);
    }
	    
    // save velocity gradient (if needed for J integral calculation)
    SetVelocityGradient(dv(0,0),dv(1,1),dv(0,1),dv(1,0),secondPass);
    
    // convert to strain increments (e.g., now dvxx = dvx/dx * dt = d/dx(du/dt) * dt = d/dt(du/dx) * dt = du/dx)
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

// Move position (2D) (in mm)
// external work units g-mm^2/sec^2 (* 10^-9 to get J)
void MatPoint2D::MovePosition(double delTime,Vector *dv)
{	double dx=delTime*dv->x;
	double dy=delTime*dv->y;
	pos.x+=dx;
    pos.y+=dy;
    extWork+=dx*pFext.x+dy*pFext.y;
}

// Move velocity (2D) (in mm/sec) possibly with damping
void MatPoint2D::MoveVelocity(double delTime,double damping,Vector *vstar)
{	vel.x+=delTime*(acc.x-damping*vstar->x);
    vel.y+=delTime*(acc.y-damping*vstar->y);
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

// set position (2D) (in mm/sec)
void MatPoint2D::SetVelocity(Vector *pt)
{	vel.x=pt->x;
    vel.y=pt->y;
	vel.z=0.;
}

// thickness (in mm)
double MatPoint2D::thickness() { return thick; }

// Find internal force as -mp sigma.deriv * 1000. which converts to g mm/sec^2 or micro N
// add external force (times a shape function)
// Store in buffer
// (note: stress is specific stress in units N/m^2 cm^3/g, Multiply by 1000 to make it mm/sec^2)
void MatPoint2D::GetFintPlusFext(Vector *theFrc,double fni,double xDeriv,double yDeriv,double zDeriv)
{	
	double mpug = mp*1000.;
	theFrc->x = -mpug*((sp.xx-pressure)*xDeriv+sp.xy*yDeriv) + fni*pFext.x;
	theFrc->y = -mpug*(sp.xy*xDeriv+(sp.yy-pressure)*yDeriv) + fni*pFext.y;
	theFrc->z = 0.0;
}

// add to the temperature gradient (non-rigid particles only)
void MatPoint2D::AddTemperatureGradient(Vector *grad)
{	pTemp->DT.x+=grad->x;
    pTemp->DT.y+=grad->y;
}

// return conduction force = - mp (Vp/V0) [k/rho0] Grad T . Grad S (units N-mm/sec)
// and k/rho0 is stored in k in units (N mm^3/(sec-K-g))
//  (non-rigid particles only)
double MatPoint2D::FCond(double dshdx,double dshdy,double dshdz,TransportProperties *t)
{
	Tensor *kten = &(t->kCondTensor);
	return -mp*GetRelativeVolume()*((kten->xx*pTemp->DT.x + kten->xy*pTemp->DT.y)*dshdx
						+ (kten->xy*pTemp->DT.x + kten->yy*pTemp->DT.y)*dshdy);
}

// add to the concentration gradient (1/mm) (non-rigid particles only)
void MatPoint2D::AddConcentrationGradient(Vector *grad)
{	pDiffusion->Dc.x+=grad->x;
    pDiffusion->Dc.y+=grad->y;
}

// return diffusion force = - V [D] Grad C . Grad S in (mm^3) (mm^2/sec) (1/mm) (1/mm) = mm^3/sec
// (non-rigid particles only)
double MatPoint2D::FDiff(double dshdx,double dshdy,double dshdz,TransportProperties *t)
{
	Tensor *Dten = &(t->diffusionTensor);
	return -GetVolume(DEFORMED_VOLUME)*((Dten->xx*pDiffusion->Dc.x + Dten->xy*pDiffusion->Dc.y)*dshdx
						+ (Dten->xy*pDiffusion->Dc.x + Dten->yy*pDiffusion->Dc.y)*dshdy);
}

// return kinetic energy (g mm^2/sec^2) = nanoJ
double MatPoint2D::KineticEnergy(void)
{	return 0.5*mp*(vel.x*vel.x+vel.y*vel.y);
}

// get deformation gradient, which is stored in strain and rotation tensors
Matrix3 MatPoint2D::GetDeformationGradientMatrix(void)
{	double F[3][3];
	GetDeformationGradient(F);
	Matrix3 Fm(F[0][0],F[0][1],F[1][0],F[1][1],F[2][2]);
	return Fm;
}

// get the symmetric elastic Left-Cauchy tensor in a Matrix3
Matrix3 MatPoint2D::GetElasticLeftCauchyMatrix(void)
{   return Matrix3(eplast.xx,eplast.xy,eplast.xy,eplast.yy,eplast.zz);
}

// get deformation gradient, which is stored in strain and rotation tensors
void MatPoint2D::GetDeformationGradient(double F[][3])
{
	// current deformation gradient in 2D
    if(theMaterials[MatID()]->PartitionsElasticAndPlasticStrain())
    {   F[0][0] = 1. + ep.xx + eplast.xx;
        double exy = ep.xy + eplast.xy;
        F[0][1] = 0.5*(exy - wrot.xy);
        F[1][0] = 0.5*(exy + wrot.xy);
        F[1][1] = 1. + ep.yy + eplast.yy;
        F[2][2] = 1. + ep.zz + eplast.zz;
    }
    else
    {   F[0][0] = 1. + ep.xx;
        F[0][1] = 0.5*(ep.xy - wrot.xy);
        F[1][0] = 0.5*(ep.xy + wrot.xy);
        F[1][1] = 1. + ep.yy;
        F[2][2] = 1. + ep.zz;
    }
}

// get relative volume from det J for large deformation material laws
double MatPoint2D::GetRelativeVolume(void)
{   double pF[3][3];
    GetDeformationGradient(pF);
    return pF[2][2]*(pF[0][0]*pF[1][1]-pF[1][0]*pF[0][1]);
}

// get dilated current volume using current deformation gradient
// only used for contact (cracks and multimaterial) and for transport tasks
// when actual is false, get t*Ap, where Ap is deformed particle area
double MatPoint2D::GetVolume(bool volumeType)
{	double rho=theMaterials[MatID()]->rho*0.001;		// in g/mm^3
	if(volumeType==DEFORMED_VOLUME)
		return GetRelativeVolume()*mp/rho;						// in mm^3
	
	// get thickness times area (for contact and for gradient)
	double pF[3][3];
	GetDeformationGradient(pF);
	return (pF[0][0]*pF[1][1]-pF[1][0]*pF[0][1])*mp/rho;
}

// To support CPDI find nodes in the particle domain, find their element,
// their natural coordinates, and weighting values for gradient calculations
// Should be done only once per time step
// throws CommonException() if particle corner has left the grid
void MatPoint2D::GetCPDINodesAndWeights(int cpdiType)
{
	// get particle 2D deformation gradient
	double pF[3][3];
	GetDeformationGradient(pF);
	
	// get polygon vectors - these are from particle to edge
    //      and generalize semi width lp in 1D GIMP
	Vector r1,r2,c;
	r1.x = pF[0][0]*mpmgrid.partx;
	r1.y = pF[1][0]*mpmgrid.partx;
	r2.x = pF[0][1]*mpmgrid.party;
	r2.y = pF[1][1]*mpmgrid.party;
	
    // Particle domain area is area of the full parallelogram
    // Assume positive due to orientation of initial vectors, and sign probably does not matter
    double Ap = 4.*(r1.x*r2.y - r1.y*r2.x);
    
	try
	{	if(cpdiType == LINEAR_CPDI)
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
			
			// gradient weighting values
			Ap = 1./Ap;
			cpdi[0]->wg.x = (r1.y-r2.y)*Ap;
			cpdi[0]->wg.y = (-r1.x+r2.x)*Ap;
			cpdi[1]->wg.x = (r1.y+r2.y)*Ap;
			cpdi[1]->wg.y = (-r1.x-r2.x)*Ap;
			cpdi[2]->wg.x = (-r1.y+r2.y)*Ap;
			cpdi[2]->wg.y = (r1.x-r2.x)*Ap;
			cpdi[3]->wg.x = (-r1.y-r2.y)*Ap;
			cpdi[3]->wg.y = (r1.x+r2.x)*Ap;
		}
		
		else
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
			
			// nodes at four edges in ccw direction
			c.x = pos.x-r2.x;
			c.y = pos.y-r2.y;
			cpdi[4]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
			theElements[cpdi[4]->inElem]->GetXiPos(&c,&cpdi[4]->ncpos);
			
			c.x = pos.x+r1.x;
			c.y = pos.y+r1.y;
			cpdi[5]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
			theElements[cpdi[5]->inElem]->GetXiPos(&c,&cpdi[5]->ncpos);
			
			c.x = pos.x+r2.x;
			c.y = pos.y+r2.y;
			cpdi[6]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
			theElements[cpdi[6]->inElem]->GetXiPos(&c,&cpdi[6]->ncpos);
			
			c.x = pos.x-r1.x;
			c.y = pos.y-r1.y;
			cpdi[7]->inElem = mpmgrid.FindElementFromPoint(&c)-1;
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
    catch(CommonException err)
    {   char msg[200];
        sprintf(msg,"A CPDI particle domain node has left the grid: %s",err.Message());
        throw CommonException(msg,"MatPoint2D::GetCPDINodesAndWeights");
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
// throws CommonException() in edge has left the grid
double MatPoint2D::GetTractionInfo(int face,int dof,int *cElem,Vector *corners,Vector *tscaled,int *numDnds)
{
    *numDnds = 2;
    double faceWt,ex,ey,enorm;
	Vector c1,c2;
	
	// always UNIFORM_GIMP or LINEAR_CPDI
    
    // which GIMP method (cannot be used in POINT_GIMP)
    if(ElementBase::useGimp==UNIFORM_GIMP)
    {   // initial vectors only
        double r1x = mpmgrid.partx;
        double r2y = mpmgrid.party;

        switch(face)
        {	case 1:
                // lower edge
                c1.x = pos.x-r1x;
                c2.x = pos.x+r1x;
                c1.y = c2.y = pos.y-r2y;
                faceWt = r1x*mpmgrid.GetThickness();
                break;
                
            case 2:
                // right edgt
                c1.x = c2.x = pos.x+r1x;
                c1.y = pos.y-r2y;
                c2.y = pos.y+r2y;
                faceWt = r2y*mpmgrid.GetThickness();
                break;
                
            case 3:
                // top edge
                c1.x = pos.x+r1x;
                c2.x = pos.x-r1x;
                c1.y = c2.y = pos.y+r2y;
                faceWt = r1x*mpmgrid.GetThickness();
                break;
                
            default:
                // left edge
                c1.x = c2.x = pos.x-r1x;
                c1.y = pos.y+r2y;
                c2.y = pos.y-r2y;
                faceWt = r2y*mpmgrid.GetThickness();
                break;
        }
        
        // get elements
        try
        {	cElem[0] = mpmgrid.FindElementFromPoint(&c1)-1;
            theElements[cElem[0]]->GetXiPos(&c1,&corners[0]);
            
            cElem[1] = mpmgrid.FindElementFromPoint(&c2)-1;
            theElements[cElem[1]]->GetXiPos(&c2,&corners[1]);
        }
        catch(CommonException err)
        {   char msg[200];
            sprintf(msg,"A Traction edge node has left the grid: %s",err.Message());
            throw CommonException(msg,"MatPoint2D::GetTractionInfo");
        }
		
		// get edge vector
		ex = c2.x-c1.x;
		ey = c2.y-c1.y;
    }
    else
    {   // get deformed corners, but get from element and natural coordinates
        //  from CPDI info because corners may have moved by here for any
        //  simulations that update strains between initial extrapolation
        //  and the grid forces calculation
        int d1,d2;
        switch(face)
        {	case 1:
                // lower edge
                d1=0;
                d2=1;
                break;
                
            case 2:
                // right edge
                d1=1;
                d2=2;
                break;
                
            case 3:
                // top edge
                d1=2;
                d2=3;
                break;
                
            default:
                // left edge
                d1=3;
                d2=0;
                break;
        }
        
        // copy for initial state at start of time step
        cElem[0] = cpdi[d1]->inElem;
        corners[0].x = cpdi[d1]->ncpos.x;
        corners[0].y = cpdi[d1]->ncpos.y;
        corners[0].z = 0.;
        cElem[1] = cpdi[d2]->inElem;
        corners[1].x = cpdi[d2]->ncpos.x;
        corners[1].y = cpdi[d2]->ncpos.y;
        corners[1].z = 0.;
        
        // get weighting factor as 1/2 of face area
        // the 1/2 is weighting factor to average the two nodes
        if(face==1 || face==3)
            faceWt = faceArea->x;
        else
            faceWt = faceArea->y;
		
		// get edge vector
		if(dof==N_DIRECTION || dof==T1_DIRECTION)
		{	theElements[cElem[0]]->GetPosition(&corners[0],&c1);
			theElements[cElem[1]]->GetPosition(&corners[1],&c2);
			ex = c2.x-c1.x;
			ey = c2.y-c1.y;
		}
        
    }
	
    // get traction normal vector by 1/2 the face area
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
			// cross product of edge vector with (0,0,1) = (ey, -ex)
			enorm = ey;
			ey = -ex;
			ex = enorm;
		case T1_DIRECTION:
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
	
	// always 1 in 2D planar (used in AS)
	return 1.;
}
	
