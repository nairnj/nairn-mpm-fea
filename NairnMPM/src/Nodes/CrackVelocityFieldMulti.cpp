/********************************************************************************
	CrackVelocityFieldMulti.cpp
	NairnMPM

	Created by John Nairn on 21 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "Exceptions/CommonException.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Materials/MaterialBase.hpp"
//#include "MPM_Classes/MPMBase.hpp"

#pragma mark INITIALIZATION

// constructor
CrackVelocityFieldMulti::CrackVelocityFieldMulti(short theLoc,int cnum) : CrackVelocityField(theLoc,cnum)
{	numberMaterials=0;
}

// zero all active material velocity fiels
void CrackVelocityFieldMulti::ZeroMatFields(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->Zero();
	}
	numberMaterials=0;
}

#pragma mark TASK 1 METHODS

// Add to mass gradient
void CrackVelocityFieldMulti::AddMassGradient(int matfld,double mp,double dNdx,double dNdy,double dNdz,MPMBase *mptr)
{	//Tensor *ep=mptr->GetStrainTensor();
	//TensorAntisym *w=mptr->GetRotationStrainTensor();
	//double dvdx=(ep->xy+w->xy)/2.;
	//double dudy=(ep->xy-w->xy)/2.;
	mvf[matfld]->massGrad->x+=mp*dNdx;
	mvf[matfld]->massGrad->y+=mp*dNdy;
	mvf[matfld]->massGrad->z+=mp*dNdz;
}

// Delete empty velocity fields, count number of materials, and return total mass
double CrackVelocityFieldMulti::GetTotalMassTask1(void)
{	int i;
	double mass=0.;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	numberMaterials++;
			mass+=mvf[i]->mass;
		}
	}
	return mass;
}

#pragma mark TASK 3 METHODS

// Add to fint spread out over the materials to each has same extra accerations = f/M
void CrackVelocityFieldMulti::AddFintSpreadTask3(Vector *f)
{	int i;
	
	// special case for only one material
	if(numberMaterials==1)
	{	for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveField(mvf[i]))
			{	AddVector(&mvf[i]->fint,f);
				break;
			}
		}
	}
	
	// more than one material
	else
	{	double totMass=GetTotalMass();
		for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveField(mvf[i]))
				AddScaledVector(&mvf[i]->fint,f,mvf[i]->mass/totMass);
		}
	}
}

// Add to fext spread out over the materials to each has same extra accerations = f/M
void CrackVelocityFieldMulti::AddFextSpreadTask3(Vector *f)
{	int i;
	
	// special case for only one material
	if(numberMaterials==1)
	{	for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveField(mvf[i]))
			{	AddVector(&mvf[i]->fext,f);
				break;
			}
		}
	}
	
	// more than one material
	else
	{	double totMass=GetTotalMass();
		for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveField(mvf[i]))
				AddScaledVector(&mvf[i]->fext,f,mvf[i]->mass/totMass);
		}
	}
}

// Calculate total force at a node from current values
// 		now m*a in g mm/sec^2
void CrackVelocityFieldMulti::CalcFtotTask3(double extDamping)
{	int i;
    for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->CalcFtotTask3(extDamping);
	}
}

#pragma mark TASK 4 METHODS

// update momenta for this MPM step
//  pk(i+1) = pk(i) + ftot * dt
void CrackVelocityFieldMulti::UpdateMomentaTask4(double timestep)
{	// update momenta
	int i;
    for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			AddScaledVector(&mvf[i]->pk,&mvf[i]->ftot,timestep);
    }
}

#pragma mark TASK 6 METHODS

// zero momentum and displacement at a node for new calculations
void CrackVelocityFieldMulti::RezeroNodeTask6(void)
{	int i;
    for(i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveField(mvf[i]))
		{	ZeroVector(&mvf[i]->pk);
			ZeroVector(&mvf[i]->disp);
			ZeroVector(mvf[i]->massGrad);
		}
    }
}

#pragma mark MATERIAL CONTACT

// Called in multimaterial mode to check contact at nodes with multiple materials
void CrackVelocityFieldMulti::MaterialContact(int nodenum,int vfld,bool postUpdate,double deltime)
{
	// exit if no contact
	if(numberMaterials<=1) return;
	
	// get center of mass results
	int i,j;
	Vector Pc,dispc;
	ZeroVector(&Pc);
	double Mc=0.;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	AddVector(&Pc,&mvf[i]->pk);
			Mc+=mvf[i]->mass;
		}
	}
	if(contact.displacementCheck)
	{	dispc=GetCMDisplacement();
		ScaleVector(&dispc,1./Mc);
	}
	
	// loop over each material
	for(i=0;i<maxMaterialFields;i++)
    {	if(!MatVelocityField::ActiveField(mvf[i])) continue;
		
		// some variables
		Vector norm,delPi;
		double rho=theMaterials[MaterialBase::fieldMatIDs[i]]->rho*0.001;	// in g/mm^3
		double dotn,massi=mvf[i]->mass,massRatio=massi/Mc;
		
		// First determine contact law from other material with most volume
		double maxOtherMaterialVolume=0.,rhoj,rhopaired=rho;
		int ipaired=0;
		for(j=0;j<maxMaterialFields;j++)
		{	if(j==i || !MatVelocityField::ActiveField(mvf[j])) continue;
			rhoj=theMaterials[MaterialBase::fieldMatIDs[j]]->rho*0.001;	// in g/mm^3
			double matUnscaledVolume=mvf[j]->mass/rhoj;
			if(matUnscaledVolume>maxOtherMaterialVolume)
			{	maxOtherMaterialVolume=matUnscaledVolume;
				ipaired=j;
				rhopaired=rhoj;
			}
		}
		// problem if ipaired not found, but it will be found
		int maxContactLaw=contact.GetMaterialContactLaw(i,ipaired);
		double maxFriction=contact.GetMaterialFriction(i,ipaired);
		
		if(maxContactLaw!=NOCONTACT)
		{	// check nodal volume
			if(unscaledVolume/mpmgrid.GetCellVolume()<contact.materialContactVmin) continue;
		
			// ignore very small mass nodes
			if(massRatio<1.e-6 || massRatio>0.999999) continue;
			
			// find -mi(vi-vc) = (ma/mc)pc-pi or momentum change to match ctr of mass momentum
			CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
			AddScaledVector(&delPi,&Pc,massRatio);

			// Get normal vector by various options
			switch(contact.materialNormalMethod)
			{	case MAXIMUM_VOLUME_GRADIENT:
				{	// Use mat with largest magnitude volume gradient
					Vector normi,normj;
					nd[nodenum]->GetMassGradient(vfld,i,&normi,1.);
					nd[nodenum]->GetMassGradient(vfld,ipaired,&normj,-1.);
					double magi=DotVectors(&normi,&normi);
					double magj=DotVectors(&normj,&normj);
					if(magi/rho >= magj/rhopaired)
						CopyScaleVector(&norm,&normi,1./sqrt(magi));
					else
						CopyScaleVector(&norm,&normj,1./sqrt(magj));
					break;
				}
				case MAXIMUM_VOLUME:
					// Use mat with most volume
					if(massi/rho >= maxOtherMaterialVolume)
						nd[nodenum]->GetMassGradient(vfld,i,&norm,1.);
					else
						nd[nodenum]->GetMassGradient(vfld,ipaired,&norm,-1.);
					CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
					break;
				/*
				case EACH_MATERIALS_MASS_GRADIENT:
					// Use each mat as is
					nd[nodenum]->GetMassGradient(vfld,i,&norm,1.);
					CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
					break;
				
				case AVERAGE_MAT_VOLUME_GRADIENTS:
					// Take an average of the two volume gradients
					break;
				*/
				default:
					break;
			}			
			
			// get approach direction momentum form delPi.n (actual (vc-vi).n = delPi/mi)
			dotn=DotVectors(&delPi,&norm);
			
			// With this check, any movement apart will be taken as noncontact
			// Also, frictional contact assumes dotn<0
			if(dotn>=0.) continue;
			
			// displacement check
			if(contact.displacementCheck)
			{	// get other mass and ignore if very small mass in other materials
				double scaleDisp=(Mc-massi)/Mc;
				if(scaleDisp<1.e-6) continue;
				
				// scale displacements
				Vector dispcScaled,dispi;
				CopyScaleVector(&dispcScaled,&dispc,1./scaleDisp);
				scaleDisp*=massi;
				CopyScaleVector(&dispi,&mvf[i]->disp,1./scaleDisp);
				double dvel = dotn/scaleDisp;
				if(contact.MaterialContact(&dispi,&dispcScaled,&norm,dvel,postUpdate,deltime)==SEPARATED) continue;
			}
		}
		else
		{	// for no contact rule only get single velocity field conditions
			CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
			AddScaledVector(&delPi,&Pc,massRatio);
		}
		
		// the material is in contact
		Vector tang;
		double dott;
		
		switch(maxContactLaw)
		{	case STICK:
			case NOCONTACT:
				break;
				
			case FRICTIONLESS:
				CopyScaleVector(&delPi,&norm,dotn);
				break;
				
			case FRICTIONAL:
				CopyVector(&tang,&delPi);
				AddScaledVector(&tang,&norm,-dotn);
				dott=sqrt(DotVectors(&tang,&tang));
				if(!DbleEqual(dott,0.))
				{	ScaleVector(&tang,1./dott);
					dott=DotVectors(&delPi,&tang);
					if(dott<0.)
					{	ScaleVector(&tang,-1.);
						dott=-dott;
					}
					if(dott>-maxFriction*dotn)
					{	AddScaledVector(&norm,&tang,-maxFriction);
						CopyScaleVector(&delPi,&norm,dotn);
					}
				}
				break;
				
			default:
				break;
		}
		
		// on post update contact, do not change nodes with boundary conditions
		unsigned char fixedDirection=nd[nodenum]->fixedDirection;
		if(postUpdate && fixedDirection)
		{	if(fixedDirection&X_DIRECTION) delPi.x=0.;
			if(fixedDirection&Y_DIRECTION) delPi.y=0.;
			if(fixedDirection&Z_DIRECTION) delPi.z=0.;
		}
		
		// change momenta
		mvf[i]->ChangeMomentum(&delPi,postUpdate,deltime);
		
		// special case two materials for efficiency (or if both will find normal the same way)
		if(numberMaterials==2)
		{	mvf[ipaired]->ChangeMomentum(ScaleVector(&delPi,-1.),postUpdate,deltime);
			break;
		}
	}
	
	/* This code looks at each pair of materials for contact. This might be good, but it is hard
	 to know. For instance if the same material changes its momentum twice, the second change
	 might cause new contact at the first one.
	 
	 An alternative used by Unitah is to change all vs center of mass, but this has issues too and
	 also would make it hard to do displacement-based contact.
	 
	 Both methods are the same if only two material interact at the point which will be the
	 most common contact situation
	 */
	
	/*
	// Loop over each pair
	int i,j;
	for(i=0;i<maxMaterialFields-1;i++)
    {	if(!MatVelocityField::ActiveField(mvf[i])) continue;

		for(j=i+1;j<maxMaterialFields;j++)
		{	if(!MatVelocityField::ActiveField(mvf[j])) continue;
			
			// handle contact between material i and j
			Vector norm,norma,normb,delPa,*pka,*pkb;
			double dotn=0.,mnode,massa,massb;
			
			// First determine if there is contact
			if(contact.GetMaterialContactLaw(i,j)!=NOCONTACT)
			{	// ignore edge nodes
				if(unscaledVolume/mpmgrid.GetCellVolume()<contact.materialContactVmin) continue;
				
				// ignore imbalanced nodes
				massa=mvf[i]->mass;
				massb=mvf[j]->mass;
				mnode=1./(massa+massb);
				double mfraction=massa*mnode;
				if(mfraction<1.e-6 || 1.-mfraction<1.e-6) continue;
				
				// Find -ma(va-vc) which is parallel to (vb-va)
				pka=&mvf[i]->pk;
				pkb=&mvf[j]->pk;
				CopyScaleVector(&delPa,pkb,massa*mnode);
				AddScaledVector(&delPa,pka,-massb*mnode);
				
				// average the two normals from the mass gradient and normalize
				nd[nodenum]->GetMassGradient(vfld,i,&norma,1.);
				nd[nodenum]->GetMassGradient(vfld,j,&normb,1.);
				CopyVector(&norm,&norma);
				SubVector(&norm,&normb);
				ScaleVector(&norm,1./sqrt(DotVectors(&norm,&norm)));
				
				// get approach direction momentum (actual (vb-va).n = dotn*(ma+mb)/(ma*mb)
				dotn=DotVectors(&delPa,&norm);
				
				// With this check, any movement apart will be taken as noncontact
				// Also, frictional contact assumes dotn<0
				if(dotn>=0.) continue;
				
				// displacement check
				if(contact.displacementCheck)
				{	Vector dispa,dispb;
					CopyScaleVector(&dispa,&mvf[i]->disp,1./massa);
					CopyScaleVector(&dispb,&mvf[j]->disp,1./massb);
					double dvel = DbleEqual(massa,0.) || DbleEqual(massb,0.) ? 0 : (massa+massb)*dotn/(massa*massb);
					if(contact.MaterialContact(&dispa,&dispb,&norm,dvel,postUpdate,deltime)==SEPARATED) continue;
				}
			}
			else
			{	// for no contact rule only get stick conditions
				massa=mvf[i]->mass;
				massb=mvf[j]->mass;
				pka=&mvf[i]->pk;
				pkb=&mvf[j]->pk;
				mnode=1./(massa+massb);
				mnode=1./(massa+massb);
				CopyScaleVector(&delPa,pkb,massa*mnode);
				AddScaledVector(&delPa,pka,-massb*mnode);
			}
			
			// the two materials are in contact
			Vector delPb,tang;
			double dott,mu;
			
			switch(contact.GetMaterialContactLaw(i,j))
			{	case STICK:
				case NOCONTACT:
					break;
					
				case FRICTIONLESS:
					CopyScaleVector(&delPa,&norm,dotn);
					break;
					
				case FRICTIONAL:
					CopyVector(&tang,&delPa);
					AddScaledVector(&tang,&norm,-dotn);
					dott=sqrt(DotVectors(&tang,&tang));
					if(!DbleEqual(dott,0.))
					{	ScaleVector(&tang,1./dott);
						dott=DotVectors(&delPa,&tang);
						if(dott<0.)
						{	ScaleVector(&tang,-1.);
							dott=-dott;
						}
						mu=-contact.GetMaterialFriction(i,j);
						if(dott>mu*dotn)
						{	AddScaledVector(&norm,&tang,mu);
							CopyScaleVector(&delPa,&norm,dotn);
						}
					}
					break;
					
				default:
					break;
			}
			CopyScaleVector(&delPb,&delPa,-1.);
				
			// on post update contact, do not change nodes with boundary conditions
			unsigned char fixedDirection=nd[nodenum]->fixedDirection;
			if(postUpdate && fixedDirection)
			{	if(fixedDirection&X_DIRECTION) delPa.x=delPb.x=0.;
				if(fixedDirection&Y_DIRECTION) delPa.y=delPb.y=0.;
				if(fixedDirection&Z_DIRECTION) delPa.z=delPb.z=0.;
			}
			
			// change momenta
			mvf[i]->ChangeMomentum(&delPa,postUpdate,deltime);
			mvf[j]->ChangeMomentum(&delPb,postUpdate,deltime);
		}
	}
	 */
}

// retrieve mass gradient
void CrackVelocityFieldMulti::GetMassGradient(int matfld,Vector *grad,double scale) { CopyScaleVector(grad,mvf[matfld]->massGrad,scale); }

#pragma mark VELOCITY METHODS

// Calculate velocity at a node from current momentum and mass matrix in all velocity fields
void CrackVelocityFieldMulti::CalcVelocityForStrainUpdate(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->CalcVelocityForStrainUpdate();
	}
}

#pragma mark BOUNDARY CONDITIONS

// zero x moment and velocity
void CrackVelocityFieldMulti::SetXMomVel(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	mvf[i]->pk.x=0.;
			mvf[i]->vk.x=0.;
		}
	}
}

// zero y moment and velocity
void CrackVelocityFieldMulti::SetYMomVel(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	mvf[i]->pk.y=0.;
			mvf[i]->vk.y=0.;
		}
	}
}

// zero z moment and velocity
void CrackVelocityFieldMulti::SetZMomVel(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	mvf[i]->pk.z=0.;
			mvf[i]->vk.z=0.;
		}
	}
}

// zero momentum in direction (cos(angle), -sin(angle)) or vector rotated from postive x axis
// by clockwise angle. The desired vector is (p.t)t where t is unit vector normal to
// skew direction and here t = (sin(angle), cos(angle))
void CrackVelocityFieldMulti::SetSkewMomVel(double angle)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	double c=cos(angle),s=sin(angle);
			Vector *npk=&mvf[i]->pk;
			double momx=npk->x*s*s + npk->y*c*s;
			double momy=npk->x*c*s + npk->y*c*c;
			npk->x=momx;
			npk->y=momy;
			mvf[i]->vk.x=momx/mvf[i]->mass;
			mvf[i]->vk.y=momy/mvf[i]->mass;
		}
	}
}

// add to x momentum and velocity from BCs
void CrackVelocityFieldMulti::AddXMomVel(double vx)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	mvf[i]->pk.x+=mvf[i]->mass*vx;
			mvf[i]->vk.x+=vx;
		}
	}
}

// add to y momentum and velocity from BCs
void CrackVelocityFieldMulti::AddYMomVel(double vy)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	mvf[i]->pk.y+=mvf[i]->mass*vy;
			mvf[i]->vk.y+=vy;
		}
	}
}

// add to z momentum and velocity from BCs
void CrackVelocityFieldMulti::AddZMomVel(double vz)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	mvf[i]->pk.z+=mvf[i]->mass*vz;
			mvf[i]->vk.z+=vz;
		}
	}
}

// Add velocity in direction (cos(angle), -sin(angle)) or vector rotated from postive x axis
// by clockwise angle.
void CrackVelocityFieldMulti::AddSkewMomVel(double vel,double angle)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	double velx=cos(angle)*vel;
			double vely=-sin(angle)*vel;
			mvf[i]->pk.x+=mvf[i]->mass*velx;
			mvf[i]->vk.x+=velx;
			mvf[i]->pk.y+=mvf[i]->mass*vely;
			mvf[i]->vk.y+=vely;
		}
	}
}

// set x force to -p(interpolated)/time such that updated momentum
//    of pk.x + deltime*ftot.x will be zero
void CrackVelocityFieldMulti::SetXFtot(double deltime)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->ftot.x=-mvf[i]->pk.x/deltime;
	}
}

// set y force to -p(interpolated)/time such that updated momentum
//    of pk.y + deltime*ftot.y will be zero
void CrackVelocityFieldMulti::SetYFtot(double deltime)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->ftot.y=-mvf[i]->pk.y/deltime;
	}
}

// set z force to -p(interpolated)/time such that updated momentum
//    of pk.z + deltime*ftot.z will be zero
void CrackVelocityFieldMulti::SetZFtot(double deltime)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->ftot.z=-mvf[i]->pk.z/deltime;
	}
}

// Change current ftot such that updated momentum of (pk.x + deltime*ftot.x, pk.y + deltime*ftot.y) will have zero
// momentum in the (cos(theta), -sin(angle)) direction (or direction clockwise from positive x axis by angle).
// Superpose force to induce zero momentum in the skew direction:
//			f dt = (p.t)t - p(interpolated)
// where t = (sin(angle),cos(angle)) is tangential vector, with the existing component of total force
// in the t direction or (f.t)t
void CrackVelocityFieldMulti::SetSkewFtot(double deltime,double angle)
{	
	 int i;
	 for(i=0;i<maxCrackFields;i++)
	 {   if(MatVelocityField::ActiveField(mvf[i]))
		{	double c=cos(angle),s=sin(angle);
			double dfxdt=-mvf[i]->pk.x*c*c + mvf[i]->pk.y*c*s;	// to get zero skew momentum
			double dfydt=mvf[i]->pk.x*c*s - mvf[i]->pk.y*s*s;
			double fx=mvf[i]->ftot.x*s*s + mvf[i]->ftot.y*c*s;	// f normal to skew direction
			double fy=mvf[i]->ftot.x*c*s + mvf[i]->ftot.y*c*c;
			mvf[i]->ftot.x=fx + dfxdt/deltime;
			mvf[i]->ftot.y=fy + dfydt/deltime;
		}
	 }
}

// add to x force such that updated momentum will be mass*velocity
void CrackVelocityFieldMulti::AddXFtot(double deltime,double velx)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->ftot.x+=mvf[i]->mass*velx/deltime;
	}
}

// add to y force such that updated momentum will be mass*velocity
void CrackVelocityFieldMulti::AddYFtot(double deltime,double vely)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->ftot.y+=mvf[i]->mass*vely/deltime;
	}
}

// add to z force such that updated momentum will be mass*velocity
void CrackVelocityFieldMulti::AddZFtot(double deltime,double velz)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->ftot.z+=mvf[i]->mass*velz/deltime;
	}
}

// set force in the skew direction (cos(angle),-sin(angle)), or direction rotated from postive
// x axis by clockwise angle.
void CrackVelocityFieldMulti::AddSkewFtot(double deltime,double vel,double angle)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	double velx=cos(angle)*vel;
			double vely=-sin(angle)*vel;
			mvf[i]->ftot.x+=mvf[i]->mass*velx/deltime;
			mvf[i]->ftot.y+=mvf[i]->mass*vely/deltime;
		}
	}
}

#pragma mark ACCESSORS

// location for crack in this field
// total mass all velocity fields
double CrackVelocityFieldMulti::GetTotalMass(void)
{	int i;
	double mass=0;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mass+=mvf[i]->mass;
	}
	return mass;
}

// get one mass
double CrackVelocityFieldMulti::GetMass(int matfld)
{	if(MatVelocityField::ActiveField(mvf[matfld]))
		return mvf[matfld]->mass;
	else
		return 0.;
}

// get center of mass momentum for all material fields in this crack velocity field
Vector CrackVelocityFieldMulti::GetCMatMomentum(void)
{	Vector pk;
	ZeroVector(&pk);
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			AddVector(&pk,&mvf[i]->pk);
	}
	return pk;
}

// get center of mass (displacement*mass so displacement is vector/total mass)
Vector CrackVelocityFieldMulti::GetCMDisplacement(void)
{	Vector dk;
	ZeroVector(&dk);
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			AddVector(&dk,&mvf[i]->disp);
	}
	return dk;
}

// get center of mass momentum for all material fields in this crack velocity field
Vector CrackVelocityFieldMulti::GetCMatFtot(void)
{	Vector fk;
	ZeroVector(&fk);
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			AddVector(&fk,&mvf[i]->ftot);
	}
	return fk;
}

/* in response to crack contact, change moment by changing velocity of all materials
	the same amount
 
   Change velocity by dP/M, where M is total mass
   Material i velocity becomes vi = pi/mi + dP/M
   Material i momentum change is mi vi = pi + mi dP/M
*/
void CrackVelocityFieldMulti::ChangeMomentum(Vector *delP,bool postUpdate,double deltime)
{
	int i;
	
	// special case for only one material
	if(numberMaterials==1)
	{	for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveField(mvf[i]))
			{	mvf[i]->ChangeMomentum(delP,postUpdate,deltime);
				break;
			}
		}
	}
	
	// more than one material
	else
	{	Vector partialDelP;
		double totMass=GetTotalMass();
		for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveField(mvf[i]))
				mvf[i]->ChangeMomentum(CopyScaleVector(&partialDelP, delP, mvf[i]->mass/totMass),postUpdate,deltime);
		}
	}
}

// copy all material velocity fields for boundary conditions methods, returning new offset into the save array
int CrackVelocityFieldMulti::CopyFieldMomenta(Vector *holdPk,int offset)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	holdPk[offset]=mvf[i]->pk;
			offset++;
		}
	}
	return offset;
}

// paste all material velocity fields back for boundary conditions methods, returning new offset into the saved array
int CrackVelocityFieldMulti::PasteFieldMomenta(Vector *holdPk,int offset)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	mvf[i]->pk=holdPk[offset];
			offset++;
		}
	}
	return offset;
}
	