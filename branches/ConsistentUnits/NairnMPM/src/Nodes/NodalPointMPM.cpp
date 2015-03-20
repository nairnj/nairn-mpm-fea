/********************************************************************************
    More NodalPoint.cpp for MPM only
    nairn-mpm-fea
    
    Created by John Nairn on Mar 17 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	NodalPoint: Methods for MPM only
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Nodes/NodalPoint.hpp"
#include "System/ArchiveData.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Exceptions/CommonException.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackNode.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "Nodes/MatVelocityField.hpp"
#include "Nodes/MaterialInterfaceNode.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"

// class statics
double NodalPoint::interfaceEnergy=0.;

#pragma mark INITIALIZATION

// MPM Destructor
NodalPoint::~NodalPoint()
{	// Delete crack velocity fields
	if(cvf!=NULL)
	{	for(int i=0;i<maxCrackFields;i++)
        {   if(cvf[i]!=NULL)
				delete cvf[i];
		}
		free(cvf);
	}
}

// Called by PreliminaryCalcs() just before the MPM analysis starts
// Can allocate things that were not known while reading the XML file
void NodalPoint::PrepareForFields()
{
	// need maxCrackFields Crack Velocity fields
    //    (1 if no cracks or MAX_FIELDS_FOR_CRACKS/MAX_FIELDS_FOR_ONE_CRACK if any cracks)
	cvf=(CrackVelocityField **)malloc(sizeof(CrackVelocityField *)*maxCrackFields);
	if(cvf==NULL) throw CommonException("Memory error allocating crack velocity field pointers.",
										"NodalPoint::PrepareForFields");

	// cvf[0] is always created, others created as needed and left in place once needed
	cvf[0]=CrackVelocityField::CreateCrackVelocityField(0,0,0);
	if(cvf[0]==NULL) throw CommonException("Memory error allocating crack velocity field 0.",
										   "NodalPoint::PrepareForFields");
	int i;
	for(i=1;i<maxCrackFields;i++) cvf[i]=NULL;
}


// zero data and reduce to one field at start of a step
void NodalPoint::InitializeForTimeStep(void)
{	// Zero the crack velocity fields
	for(int i=0;i<maxCrackFields;i++)
	{	if(cvf[i]!=NULL)
			cvf[i]->Zero(0,0,TRUE);
	}
	
	// for conduction and diffusion
	gVolume=0.;
	gConcentration=0.;
	gTemperature=0.;
	gMpCp=0.;
	fdiff=0.;
	fcond=0.;
}

#pragma mark TASK 0 METHODS

// When there are cracks, call this method to allocate crack and material velocity fields
// that are needed on this time step. Called once in initialization task
// throws CommonException() on meory error
// SCWarning - may look at [2] and [3]
short NodalPoint::AddCrackVelocityField(int matfld,CrackField *cfld)
{
	// find and return velocity field, allocate memory if needed
	short vfld=0;
	
	// If there are two cracks, try to use field [3], but if problems, set cfld[1].loc to NO_CRACK
	// and try to use field [0], [1], or [2] below
	if(cfld[1].loc != NO_CRACK)
	{	if(!CrackVelocityField::ActiveCrackField(cvf[3]))
		{	// a new field - put [0] into first crack and [1] into second crack (number and orientation)
			if(cvf[3]==NULL)
			{	cvf[3]=CrackVelocityField::CreateCrackVelocityField(3,cfld[0].loc,cfld[0].crackNum);
				if(cvf[3]==NULL)
                    throw CommonException("Memory error allocating crack velocity field #3","NodalPoint::AddCrackVelocityField");
			}
			else
				cvf[3]->SetLocationAndCrack(cfld[0].loc,cfld[0].crackNum,FIRST_CRACK);
			cvf[3]->SetLocationAndCrack(cfld[1].loc,cfld[1].crackNum,SECOND_CRACK);
			
			// add the crack normals to selected field (and done)
			vfld = 3;
			cvf[3]->AddNormals(&cfld[0].norm,FIRST_CRACK);
			cvf[3]->AddNormals(&cfld[1].norm,SECOND_CRACK);
		}
		else
		{	// make sure same two cracks
			int c1=-1,c2=-1;
			if(cfld[0].crackNum==cvf[3]->crackNumber(FIRST_CRACK))
			{	// crack in [0] matches first crack in [3]
				if(cfld[1].crackNum==cvf[3]->crackNumber(SECOND_CRACK))
				{	// crack in [1] matches second crack in [3] (handle in (c1>=0) block below)
					c1=0;
					c2=1;
				}
				else
				{	// crack [1] is a third crack at this node.Try again with single crack that matches [0]
					cfld[1].crackNum = NO_CRACK;
					
					// warn
					warnings.Issue(CrackHeader::warnThreeCracks,21);
				}
			}
			else if(cfld[0].crackNum==cvf[3]->crackNumber(SECOND_CRACK))
			{	// crack in [0] matches second crack in [3]
				if(cfld[1].crackNum==cvf[3]->crackNumber(FIRST_CRACK))
				{	// crack in [1] matches first crack in [3] (handle in (c1>=0) block below)
					c1=1;
					c2=0;
				}
				else
				{	// crack [1] is a third crack at this node. Try again with single crack that matches [0]
					cfld[1].loc = NO_CRACK;
					
					// warn
					warnings.Issue(CrackHeader::warnThreeCracks,22);
				}
			}
			else if(cfld[1].crackNum==cvf[3]->crackNumber(FIRST_CRACK) || cfld[1].crackNum==cvf[3]->crackNumber(SECOND_CRACK))
			{	// crack in [1] matches first or second crack in [3], but crack in [0] is unknown
				// move [1] to [0] and try agin with single crack that matches
				cfld[0].crackNum = cfld[1].crackNum;
				cfld[0].loc = cfld[1].loc;
				cfld[0].norm = cfld[1].norm;
				cfld[1].loc = NO_CRACK;
				
				// warn
				warnings.Issue(CrackHeader::warnThreeCracks,23);
			}
			else
			{	// neither crack matches a crack in [3], which is pretty bad, just switch to field zero (and done)
				// add both normals to same crack (but may not be needed in field [0])
				vfld = 0;
				cvf[0]->AddNormals(&cfld[0].norm,FIRST_CRACK);
				cvf[0]->AddNormals(&cfld[1].norm,FIRST_CRACK);
				
				// warn
				warnings.Issue(CrackHeader::warnThreeCracks,24);
			}
			
			// continue to this block only if the two crack match the two cracks in [3]
			if(c1>=0)
			{	// Now that the two cracks match the two cracks in [3], check their sides
				if(cfld[c1].loc==cvf[3]->location(FIRST_CRACK))
				{	// first crack is correct orientation
					if(cfld[c2].loc=cvf[3]->location(SECOND_CRACK))
					{	// they both match so add the normals (and done)
						vfld = 3;
						cvf[3]->AddNormals(&cfld[c1].norm,FIRST_CRACK);
						cvf[3]->AddNormals(&cfld[c2].norm,SECOND_CRACK);
					}
					else
					{	// c2 crack is wrong or node may be on c2 crack, might be better to use [1] or [2]
						// (i.e., to change side of c2 to same side and try again with a single cracks
						// move c1 to [0] crack info if needed
						if(c1!=0)
						{	cfld[0].crackNum = cfld[1].crackNum;
							cfld[0].loc = cfld[1].loc;
							cfld[0].norm = cfld[1].norm;
						}
						cfld[1].loc = NO_CRACK;
						
						// warn
						if(warnings.Issue(CrackHeader::warnNodeOnCrack,21)==GAVE_WARNING) Describe();
					}
				}
				else if(cfld[c2].loc==cvf[3]->location(SECOND_CRACK))
				{	// c2 crack correct, but c1 is wrong or node may be on c1 crack, might be better to use [1] pr [2]
					// (i.e., change side of c1 to same side and try again with a single crack)
					// switch c2 to [0] crack info
					if(c2!=0)
					{	cfld[0].crackNum = cfld[1].crackNum;
						cfld[0].loc = cfld[1].loc;
						cfld[0].norm = cfld[1].norm;
					}
					cfld[1].loc = NO_CRACK;
					
					// warn
					if(warnings.Issue(CrackHeader::warnNodeOnCrack,22)==GAVE_WARNING) Describe();
				}
				else
				{	// both cracks wrong, not sure what to do so just revert to [0] (and done)
					// add both normals to only crack in [0] (but may not be needed in field [0])
					vfld=0;
					cvf[0]->AddNormals(&cfld[0].norm,FIRST_CRACK);
					cvf[0]->AddNormals(&cfld[1].norm,FIRST_CRACK);
					
					// warn
					if(warnings.Issue(CrackHeader::warnNodeOnCrack,23)==GAVE_WARNING) Describe();
				}
			}
		}
	}
	
	// If only one crack (or two cracks are trying again with a single crack), continue the calculations
	if(cfld[1].loc==NO_CRACK)
	{	switch(cfld[0].loc)
		{	case NO_CRACK:
				vfld = 0;
				break;
				
			case ABOVE_CRACK:
			case BELOW_CRACK:
				// if cvf[1] is empty, then use it now
				if(!CrackVelocityField::ActiveCrackField(cvf[1]))
                {   if(cvf[1]==NULL)
                    {	cvf[1]=CrackVelocityField::CreateCrackVelocityField(1,cfld[0].loc,cfld[0].crackNum);
                        if(cvf[1]==NULL) 
                            throw CommonException("Memory error allocating crack velocity field #1","NodalPoint::AddCrackVelocityField"); 
                    }
                    else
                        cvf[1]->SetLocationAndCrack(cfld[0].loc,cfld[0].crackNum,FIRST_CRACK);
					vfld=1;
				}
				
                // if cvf[1] is not empty, but refers to same crack, the two options are:
                //  a. same location so use it OR b. different location which means node on crack warning and switch to [0]
				else if(cvf[1]->crackNumber(FIRST_CRACK)==cfld[0].crackNum)
				{	if(cvf[1]->location(FIRST_CRACK)==cfld[0].loc)
					{	// found another point for [1]
						vfld=1;
					}
					else
					{	// it can only be field [0]
						vfld=0;

						// Here means both above and below same crack for field [1], which can only happen if a
                        // node is on a crack
                        if(warnings.Issue(CrackHeader::warnNodeOnCrack,11)==GAVE_WARNING) Describe();
                    
                        // tell field [0] it has crack from field 1 (but info currently not used)
                        if(cvf[0]->location(FIRST_CRACK)==NO_CRACK)
                            cvf[0]->SetLocationAndCrack(cfld[0].loc,cfld[0].crackNum,FIRST_CRACK);
					}
				}
				
				// Here means found a new crack, which hopefully will be appropriate for [2]
				// Here means cvf[1]->crackNumber(FIRST_CRACK)!=cfld[0].crackNum
                // Note that this code should never be reached if only have one crack
				else
				{	// if [2] is empty, then use it now
					if(!CrackVelocityField::ActiveCrackField(cvf[2]))
					{   if(cvf[2]==NULL)
                        {	cvf[2]=CrackVelocityField::CreateCrackVelocityField(2,cfld[0].loc,cfld[0].crackNum);
                            if(cvf[2]==NULL)
                                throw CommonException("Memory error allocating crack velocity field #2","NodalPoint::AddCrackVelocityField"); 
                        }
                        else
                            cvf[2]->SetLocationAndCrack(cfld[0].loc,cfld[0].crackNum,FIRST_CRACK);
						vfld=2;
					}
                    
                    // if [2] is not empty, but refers to same crack, the two options are:
                    //  a. same location so use it OR b. different location which means node on crack warning
					else if(cvf[2]->crackNumber(FIRST_CRACK)==cfld[0].crackNum)
					{	if(cvf[2]->location(FIRST_CRACK)==cfld[0].loc)
						{	// found another point for [2]
							vfld=2;
						}
						else
						{	// it can only be field 0
							vfld=0;

							// Here means both above and below crack for field [2], which can only happen if a
                            // node is on a crack
                            if(warnings.Issue(CrackHeader::warnNodeOnCrack,12)==GAVE_WARNING) Describe();
                                
                            // tell field [0] it has crack from field 2 (but info currently not used)
                            if(cvf[0]->location(SECOND_CRACK)==NO_CRACK)
                                cvf[0]->SetLocationAndCrack(cfld[0].loc,cfld[0].crackNum,SECOND_CRACK);
						}
					}
                    
                    // here means crack differs from both [1] and [2] so it is a third crack on this node
					else
					{	// found a third crack at this node, try to use [0]
						vfld=0;
						
						// warn
						warnings.Issue(CrackHeader::warnThreeCracks,11);
 					}
				}
				
				// add the crack normals to selected field
                cvf[vfld]->AddNormals(&cfld[0].norm,FIRST_CRACK);
				
				break;
			default:
				break;
		}
	}
	
	
	// return crack velocity field that was used
	return vfld;
}

// When there are cracks or multimedia, call this method to allocate material velocity fields
// that are needed on this time step. Called onlyb in initialization task
void NodalPoint::AddMatVelocityField(short vfld,int matfld)
{	cvf[vfld]->AddMatVelocityField(matfld);
}

// precheck to avoid too much critical code in parallel block
bool NodalPoint::NeedsMatVelocityField(short vfld,int matfld) const
{   return cvf[vfld]->NeedsMatVelocityField(matfld);
}

// When has crack and multimaterial velocity fields, make sure ghost node has copy of needed fields
void NodalPoint::CopyFieldInitialization(NodalPoint *ghost)
{	ghost->UseTheseFields(cvf);
}

// Create fields on this this node that match the supplied fields
void NodalPoint::UseTheseFields(CrackVelocityField **rcvf)
{	
	for(int i=0;i<maxCrackFields;i++)
	{	if(rcvf[i]==NULL) continue;
		
		if(cvf[i]==NULL)
		{	// create one in ghost node if not there already
			cvf[i] = CrackVelocityField::CreateCrackVelocityField(i,0,0);
			if(cvf[i]==NULL) throw CommonException("Memory error allocating crack velocity field 3.",
												   "NodalPoint::UseTheseFields");
		}
		
		// make these match
		cvf[i]->MatchRealFields(rcvf[i]);
	}
}

#pragma mark TASK 1 AND 6 METHODS

// In mass and momentum task
// 1. Add momentum
// 2. If cracks or multimaterials, add displacements and volume
// 3. If multimaterials, add volume gradient
// 4. Add mass
// 5. Transport tasks (for non-rigid only)
void NodalPoint::AddMassMomentum(MPMBase *mptr,short vfld,int matfld,double shape,double dNdx,double dNdy,double dNdz,
								 int numPts,bool nonRigid)
{
	// add momentum
	double mp = mptr->mp;
	double fnmp = shape*mp;
	Vector wtvel;
    
    // crash with KERN_INVALID_ADDRESS at 0x0000000000000000 might mean NULL field encountered
    /*
    if(cvf[vfld]==0)
    {   cout << "\n NULL field for vfld=" << vfld << endl;
        Describe();
        mptr->Describe();
        throw CommonException("NULL crack velocity field","NodalPoint::AddMassMomentum()");
    }
     */
    
	cvf[vfld]->AddMomentumTask1(matfld,CopyScaleVector(&wtvel,&mptr->vel,fnmp),&mptr->vel,numPts);
	
	// crack contact calculations
	// (only if cracks or multimaterial mode, i.e., contact is being done)
	if(firstCrack!=NULL || fmobj->multiMaterialMode)
	{	// displacement or position for contact calculations
		if(mpmgrid.GetContactByDisplacements())
		{	// extrapolate displacements
			Vector pdisp=mptr->pos;
			cvf[vfld]->AddDisplacement(matfld,fnmp,SubVector(&pdisp,&mptr->origpos));
		}
		else
		{	// extraplate positions
			cvf[vfld]->AddDisplacement(matfld,fnmp,&mptr->pos);
		}
		
		// add dilated volume, only used by transport tasks, contact, and imperfect interfaces
		if(nonRigid)
			cvf[vfld]->AddVolume(matfld,shape*mptr->GetVolume(DEFORMED_AREA));
		else
			cvf[vfld]->AddVolume(matfld,shape*mptr->GetUnscaledVolume());

		// material contact calculations (only if multimaterial mode)
		if(fmobj->multiMaterialMode)
			cvf[vfld]->AddVolumeGradient(matfld,mptr,dNdx,dNdy,dNdz);
	}
	
	// more for non-rigid contact materials
	if(nonRigid)
	{	// add to lumped mass matrix
		cvf[vfld]->AddMass(matfld,fnmp);
        
        // transport calculations
        TransportTask *nextTransport=transportTasks;
        while(nextTransport!=NULL)
            nextTransport=nextTransport->Task1Extrapolation(this,mptr,shape);
	}
	else
	{	// for rigid particles, let the crack velocity field know
		cvf[vfld]->AddMassTask1(matfld,fnmp,1);
	}
}

// copy ghost node mass and momentum to real node
void NodalPoint::CopyMassAndMomentum(NodalPoint *real)
{	for(int vfld=0;vfld<maxCrackFields;vfld++)
    {	if(CrackVelocityField::ActiveField(cvf[vfld]))
            cvf[vfld]->CopyMassAndMomentum(real,vfld);
    }
}

// In 2nd mass and momentum extrapolation when update strains last (and only for non-rigid particles)
// 1. Add momentum
// 2. If cracks or multimaterials, add displacements and volume
// 3. If multimaterials, add volume gradient
void NodalPoint::AddMassMomentumLast(MPMBase *mptr,short vfld,int matfld,double shape,double dNdx,double dNdy,double dNdz)
{
	// add momentum
	double mp = mptr->mp;
	double fnmp = shape*mp;
	cvf[vfld]->AddMomentumTask6(matfld,fnmp,&mptr->vel);
	
	// crack contact calculations
	// (only if cracks or multimaterial mode, i.e., contact is being done)
	if(firstCrack!=NULL || fmobj->multiMaterialMode)
	{	// displacement or position for contact calculations
		if(mpmgrid.GetContactByDisplacements())
		{	// extrapolate displacements
			Vector pdisp=mptr->pos;
			cvf[vfld]->AddDisplacement(matfld,fnmp,SubVector(&pdisp,&mptr->origpos));
		}
		else
		{	// extraplate positions
			cvf[vfld]->AddDisplacement(matfld,fnmp,&mptr->pos);
		}
		
		// add dilated volume, only used by transport tasks, contact, and imperfect interfaces
        cvf[vfld]->AddVolume(matfld,shape*mptr->GetVolume(DEFORMED_AREA));
        
		// material contact calculations (only if multimaterial mode)
		if(fmobj->multiMaterialMode)
			cvf[vfld]->AddVolumeGradient(matfld,mptr,dNdx,dNdy,dNdz);
	}
}

// copy ghost node mass an momentum to real node
void NodalPoint::CopyMassAndMomentumLast(NodalPoint *real)
{	for(int vfld=0;vfld<maxCrackFields;vfld++)
    {	if(CrackVelocityField::ActiveField(cvf[vfld]))
            cvf[vfld]->CopyMassAndMomentumLast(real,vfld);
    }
}

// Add to momentum vector to selected field (in g mm/sec)(both 2D and 3D)
void NodalPoint::AddMomentumTask1(short vfld,int matfld,double wt,Vector *vel,int numPts)
{	Vector wtvel;
	cvf[vfld]->AddMomentumTask1(matfld,CopyScaleVector(&wtvel,vel,wt),vel,numPts);
}

// Add mass for selected field
void NodalPoint::AddMass(short vfld,int matfld,double mnode) { cvf[vfld]->AddMass(matfld,mnode); }

// for rigid particles, adding mass is counting number of rigid particles
void NodalPoint::AddMassTask1(short vfld,int matfld,double mnode,int numPts) { cvf[vfld]->AddMassTask1(matfld,mnode,numPts); }

// Copy volume gradient when copying from ghost to real node
void NodalPoint::CopyVolumeGradient(short vfld,int matfld,Vector *grad)
{	cvf[vfld]->CopyVolumeGradient(matfld,grad);
}

// Calculate total mass and count number of materials on this node
// Changes only this->numberMaterials and then only if in multimaterial mode
// Calculations might need to exclude nodes whose mass is too small in crack calculations
void NodalPoint::CalcTotalMassAndCount(void)
{	int i;
	nodalMass = 0.;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			nodalMass += cvf[i]->GetTotalMassAndCount();
	}
}

#ifdef COMBINE_RIGID_MATERIALS

// When has rigid particles, multimaterial mode, and cracks, copy all rigid particles in
//   field [0] to other fields
// This is only called with option to copy rigid materials to all crack fields
void NodalPoint::CopyRigidParticleField(void)
{
    // if field [0] has no rigid particles, then nothing to copy
    int rigidFieldNum;
    MatVelocityField *theRigidField = ((CrackVelocityFieldMulti *)cvf[0])->GetRigidMaterialField(&rigidFieldNum);
    if(theRigidField==NULL) return;
    
	// transfer the rigid field to all other crack fields
	for(int i=1;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(cvf[i]))
		{	// copy rigid material from field 0 to field i
			((CrackVelocityFieldMulti *)cvf[i])->CopyRigidFrom(theRigidField,rigidFieldNum);
		}
	}
}

#endif

#pragma mark TASK 3 METHODS

// Add to internal force
void NodalPoint::AddFtotTask3(short vfld,int matfld,Vector *f)
{	if(cvf[vfld]==NULL) throw CommonException("NULL crack velocity field in grid forces test","NodalPoint::AddFtotTask3");
	cvf[vfld]->AddFtotTask3(matfld,f);
}

// copy ghost node forces to real node (nonrigid only)
void NodalPoint::CopyGridForces(NodalPoint *real)
{	for(int vfld=0;vfld<maxCrackFields;vfld++)
	{	if(CrackVelocityField::ActiveNonrigidField(cvf[vfld]))
			cvf[vfld]->CopyGridForces(real,vfld);
	}
}

// Add to internal force spread out over materials for same acceleration on each
// Only called by AddTractionForce() and CrackInterfaceForce()
void NodalPoint::AddFtotSpreadTask3(short vfld,Vector f) { cvf[vfld]->AddFtotSpreadTask3(&f); }

// Add to traction force (g-mm/sec^2)
void NodalPoint::AddTractionTask3(MPMBase *mpmptr,short vfld,int matfld,Vector *f)
{	if(CrackVelocityField::ActiveField(cvf[vfld]))
	{	cvf[vfld]->AddFtotTask3(matfld,f);
	}
}

#pragma mark TASK 4 METHODS

// update momenta for this MPM step
void NodalPoint::UpdateMomentaOnNode(double timestep)
{	// update momenta in all crack velocity fields
    for(int i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(cvf[i]))
			cvf[i]->UpdateMomentaOnField(timestep);
    }
}

#pragma mark TASK 5 METHODS

// Increment velocity and acceleration for this material point using one velocity field
void NodalPoint::IncrementDelvaTask5(short vfld,int matfld,double fi,Vector *delv,Vector *dela) const
{   cvf[vfld]->IncrementDelvaTask5(matfld,fi,delv,dela);
}

#pragma mark TASK 6 METHODS

// zero momentum at a node for new calculations
void NodalPoint::RezeroNodeTask6(double deltaTime)
{	int i;
    for(i=0;i<maxCrackFields;i++)
    {	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->RezeroNodeTask6(deltaTime);
    }
}

// Add to momentum on second pass for selected field that must be there
void NodalPoint::AddMomentumTask6(short vfld,int matfld,double wt,Vector *vel)
{	cvf[vfld]->AddMomentumTask6(matfld,wt,vel);
}

#pragma mark TASK 7 CRACK CALCULATIONS FOR J AND K

// Initialize fields for grid extrapolations for strains, etc.
void NodalPoint::ZeroDisp(void)
{	
	int i;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(cvf[i]))
			cvf[i]->CreateStrainField();
	}
}

// Find field [0] from velocity field on same side of cracks as this node when
//	  do J integral calculations
// But when a contour cross a crack, phantom nodes are inserted on the crack plane
//    and those nodes use field[0] for contour integration leading up to that node but
//    field [1] for the segment after that node
// If crackNum is not 0, then phantom is crossing a different crack crackNum. A phantom
//    node is placed on that crack and use field [0] for integration up that node. For
//    integration stating on that node, use field [1] or [2], whichever is for that crack.
int NodalPoint::GetFieldForCrack(bool phantomNode,bool firstNode,DispField **dfld,int crackNum)
{
	// If phantom node and it is at the start of a segment, get field [1]
	// or get [1] or [2] depending on crackNum (if it is not zero)
	if(phantomNode && firstNode)
	{	bool active1 = CrackVelocityField::ActiveNonrigidField(cvf[1]);
		
		if(crackNum==0)
		{	// always get [1] here when do J contour
			if(!active1)
			{	*dfld = NULL;
				return 0;
			}
			*dfld = cvf[1]->df;
			return cvf[1]->GetNumberPoints();
		}
		else
		{	// Find [1] or [2] on opposite side of tcrack crackNum when doing interpolation
			
			// First check [1] if active. If it is correct crack, then use it
			if(active1)
			{	if(cvf[1]->crackNumber(FIRST_CRACK)==crackNum)
                {	*dfld = cvf[1]->df;
                    return cvf[1]->GetNumberPoints();
                }
			}
			
			// If [1] fails use [2] instead
            // SCWarning - if [1] is empty it looks at [2], but should return NULL if only one crack
			// This does not check crackNum, but unlikley to have crackNum match neither one
			bool active2 = CrackVelocityField::ActiveNonrigidField(cvf[2]);
			if(!active2)
			{	*dfld = NULL;
				return 0;
			}
			*dfld = cvf[2]->df;
			return cvf[2]->GetNumberPoints();
		}
	}
	
	// For real nodes or phantom nodes at the end of a segment, get field [0]
	bool active0 = CrackVelocityField::ActiveNonrigidField(cvf[0]);
	if(!active0)
	{	*dfld = NULL;
		return 0;
	}
	*dfld = cvf[0]->df;
	return cvf[0]->GetNumberPoints();
}

// Initialize fields on a ghost node for grid extrapolations for strains, etc.
void NodalPoint::ZeroDisp(NodalPoint *real)
{
	int i;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(real->cvf[i]))
            cvf[i]->CreateStrainField();
	}
}


// delete any strain fields that were created
void NodalPoint::DeleteDisp(void)
{
	int i;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(cvf[i]))
			cvf[i]->DeleteStrainField();
	}
}

// delete any strain fields that were created in a ghost node
void NodalPoint::DeleteDisp(NodalPoint *real)
{
	int i;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(real->cvf[i]))
            cvf[i]->DeleteStrainField();
	}
}

// Add to displacement gradient and track material type
void NodalPoint::AddUGradient(short vfld,double wt,double dudx,double dudy,double dvdx,double dvdy,int activeMatField,double mp)
{	DispField *df = cvf[vfld]->df;
	df->du.x += wt*dudx;
	df->du.y += wt*dudy;
	df->dv.x += wt*dvdx;
	df->dv.y += wt*dvdy;
    df->mass += wt;
	
	// if more than one material get shape function extrapolation to each node
	if(numActiveMaterials>1)
		df->matWeight[activeMatField] += wt/mp;
}

// add to nodal stress (in-plane only)
// wt includes rho thus stress is N/m^2
// Must scale by 1e-6 to get N/mm^2
void NodalPoint::AddStress(short vfld,double wt,Tensor *stress)
{	DispField *df = cvf[vfld]->df;
	df->stress.xx += wt*stress->xx;
	df->stress.yy += wt*stress->yy;
	df->stress.xy += wt*stress->xy;
}

// GRID_JTERMS
// Add to velocity to get kinetic energy on the grid
// wt includes sqrt(density (g/mm^3)), if axisymmtric include r for mm/radian
// kinetic energy is twice actual (to save divide by two) and units are
//		final vx^2 = g/mm^3 mm^2/sec^2 = N/m^2 = J/m^3
//		to get N/mm^2 and account for 1/2, multiply by 0.5*1e-6
void NodalPoint::AddGridVelocity(short vfld,double wt,double vx,double vy)
{	DispField *df = cvf[vfld]->df;
	df->vx += wt*vx;
	df->vy += wt*vy;
}

// Add to kinetic energy and strain energy
// wt includes density (g/mm^3), if axisymmtric include r for mm/radian
// kinetic energy is twice actual (to save divide by two) and units are
//		g/mm^3 mm^2/sec^2 = N/m^2 = J/m^3
//		to get N/mm^2 and account for 1/2, multiply by 0.5*1e-6
// work has units N/m^2 = J/m^3
void NodalPoint::AddEnergy(short vfld,double wt,double vx,double vy,double work)
{	DispField *df = cvf[vfld]->df;
	df->kinetic += wt*(vx*vx+vy*vy);
	df->work += wt*work;
}

// Copy J integral terms from ghost node to real node
void NodalPoint::CopyUGradientStressEnergy(NodalPoint *real)
{	for(int vfld=0;vfld<maxCrackFields;vfld++)
    {   CrackVelocityField *rcvf = real->cvf[vfld];
		if(CrackVelocityField::ActiveNonrigidField(rcvf))
        {   DispField *gdf = cvf[vfld]->df;
			DispField *rdf = rcvf->df;
			rdf->du.x += gdf->du.x;
			rdf->du.y += gdf->du.y;
			rdf->dv.x += gdf->dv.x;
			rdf->dv.y += gdf->dv.y;
			rdf->stress.xx += gdf->stress.xx;
			rdf->stress.yy += gdf->stress.yy;
			rdf->stress.xy += gdf->stress.xy;
			rdf->kinetic += gdf->kinetic;
			rdf->work += gdf->work;
            rdf->mass += gdf->mass;
            rdf->vx += gdf->vx;
            rdf->vy += gdf->vy;
            
            // if more than one material get shape function extrapolation to each node
            if(numActiveMaterials>1)
            {   for(int i=0;i<numActiveMaterials;i++)
                {   rdf->matWeight[i] += gdf->matWeight[i];
                }
            }
		}
	}
}

// add material weights to an array
// called by crack segment when finding crack tip materials
// Future might want to keep above and below separate to find interface crack
// Maybe could do with old GetFieldForCrack()
void NodalPoint::AddMatWeights(double wt,double *matWeight)
{	int i,j;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(cvf[i]))
		{	DispField *df = cvf[i]->df;
			if(df!=NULL)
			{	for(j=0;j<numActiveMaterials;j++)
					matWeight[j] += wt*df->matWeight[j];
			}
		}
	}
}

// Finish extrapolated strain field terms
//	Divide by mass to get final results
//	Scale to compatible units
void NodalPoint::CalcStrainField(void)
{
	int j;
	double mnode;
	
	// do all strain fields
	for(j=0;j<maxCrackFields;j++)
	{	if(!CrackVelocityField::ActiveNonrigidField(cvf[j])) continue;
		DispField *df=cvf[j]->df;
		if(df==NULL) continue;
		
		mnode=1./df->mass;
		df->du.x *= mnode;				// no units
		df->du.y *= mnode;
		df->dv.x *= mnode;
		df->dv.y *= mnode;
		
		// GRID_JTERMS
		if(JGridEnergy)
        {	df->vx *= mnode;            // sqrt(nJ/mm^3) = sqrt(uN/mm^2)
			df->vy *= mnode;
			
			df->stress.xx *= mnode;				// uN/mm^2
			df->stress.yy *= mnode;
			df->stress.xy *= mnode;
			
			// find grid energy from nodal extrapolations
			df->kinetic = 0.5*(df->vx*df->vx + df->vy*df->vy);					// uN/mm^2
			df->work = 0.5*(df->stress.xx*df->du.x + df->stress.yy*df->dv.y
							+ df->stress.xy*(df->du.y+df->dv.x));				// uN/mm^2
		}
		
		else
		{	df->stress.xx *= mnode;						// uN/mm^2
			df->stress.yy *= mnode;
			df->stress.xy *= mnode;
			
			// find energy by extrapolating particle energies
			df->kinetic *= mnode*.5;					// uN/mm^2
			df->work *= mnode;							// uN/mm^2
		}
    }
}

// Interpolate two nodes (near crack plane). This method is only called in J calculation for the
//		phantom nodes placed on the crack plane of cracks that cross the countour
// Interpolate [0] from node1 and [1] or [2] that crosses CrackNum from node2 to phantom [0]
// Interpolate [0] from node2 and [1] or [2] that crosses CrackNum from nod1 to phantom [1]
// Symbolically gets [0] = (1-fract)*n1[0] + fract*n2[i]
//                   [1] = (1-fract)*n1[i] + fract*n2[0]
// where [i] is [1] or [2] for field on opposite side of crack crackNum
void NodalPoint::Interpolate(NodalPoint *n1,NodalPoint *n2,double fract,int crackNum)
{
	// need strain field in first crack velocity field and entire second
	// crack velocity field for this phantom node (it may be zero)
	cvf[0]->CreateStrainField();
	cvf[1]=CrackVelocityField::CreateCrackVelocityField(1,BELOW_CRACK,crackNum);
	if(cvf[1]==NULL) throw CommonException("Memory error allocating crack velocity field 1.",
										   "NodalPoint::Interpolate");
	cvf[1]->CreateStrainField();
	
	// fetch [0] from node 1 and [1] or [2] from node 2
	DispField *a1fld,*a2fld;
	int a1 = n1->GetFieldForCrack(false,false,&a1fld,0);				// gets [0]
	int a2 = n2->GetFieldForCrack(true,true,&a2fld,crackNum);		    // gets [1] or [2]
	AverageStrain(cvf[0]->df,a1fld,a2fld,fract);
	cvf[0]->SetNumberPoints(a1+a2);

	// fetch [1] or [2] from node 1 and [0] from node 2
	a1 = n1->GetFieldForCrack(true,true,&a1fld,crackNum);			// gets [1] or [2]
	a2 = n2->GetFieldForCrack(false,false,&a2fld,0);		    	// gets [0]
	AverageStrain(cvf[1]->df,a1fld,a2fld,fract);
	cvf[1]->SetNumberPoints(a1+a2);
}

// interpolate between two fields and store in destination field
// fract is fracture of distance from first field to the point
//     thus want (1-fract)*src1 + fract*src2
void NodalPoint::AverageStrain(DispField *dest,DispField *src1,DispField *src2,double fract)
{	if(src1!=NULL && src2!=NULL)
    {   // mass and distance weighted average
        double f1 = (1.-fract);
        double f2 = fract;
        dest->du.x =      f1*src1->du.x +      f2*src2->du.x;
        dest->du.y =      f1*src1->du.y +      f2*src2->du.y;
        dest->dv.x =      f1*src1->dv.x +      f2*src2->dv.x;
        dest->dv.y =      f1*src1->dv.y +      f2*src2->dv.y;
        dest->kinetic =   f1*src1->kinetic +   f2*src2->kinetic;
        dest->work =      f1*src1->work +      f2*src2->work;
        dest->stress.xx = f1*src1->stress.xx + f2*src2->stress.xx;
        dest->stress.yy = f1*src1->stress.yy + f2*src2->stress.yy;
        dest->stress.zz = f1*src1->stress.zz + f2*src2->stress.zz;
        dest->stress.xy = f1*src1->stress.xy + f2*src2->stress.xy;
        dest->mass =      f1*src1->mass +      f2*src2->mass;
    }
    else if(src1!=NULL)
    {   dest->du.x =      (1.-fract)*src1->du.x;
        dest->du.y =      (1.-fract)*src1->du.y;
        dest->dv.x =      (1.-fract)*src1->dv.x;
        dest->dv.y =      (1.-fract)*src1->dv.y;
        dest->kinetic =   (1.-fract)*src1->kinetic;
        dest->work =      (1.-fract)*src1->work;
        dest->stress.xx = (1.-fract)*src1->stress.xx;
        dest->stress.yy = (1.-fract)*src1->stress.yy;
        dest->stress.zz = (1.-fract)*src1->stress.zz;
        dest->stress.xy = (1.-fract)*src1->stress.xy;
        dest->mass =      (1.-fract)*src1->mass;
    }
    else if(src2!=NULL)
    {   dest->du.x =      fract*src2->du.x;
        dest->du.y =      fract*src2->du.y;
        dest->dv.x =      fract*src2->dv.x;
        dest->dv.y =      fract*src2->dv.y;
        dest->kinetic =   fract*src2->kinetic;
        dest->work =      fract*src2->work;
        dest->stress.xx = fract*src2->stress.xx;
        dest->stress.yy = fract*src2->stress.yy;
        dest->stress.zz = fract*src2->stress.zz;
        dest->stress.xy = fract*src2->stress.xy;
        dest->mass =      fract*src2->mass;
    }
    else
    {	// outside the grid
        dest->du.x =      0.;
        dest->du.y =      0.;
        dest->dv.x =      0.;
        dest->dv.y =      0.;
        dest->kinetic =   0.;
        dest->work =      0.;
        dest->stress.xx = 0.;
        dest->stress.yy = 0.;
        dest->stress.zz = 0.;
        dest->stress.xy = 0.;
        dest->mass =      0.;
    }
}


#pragma mark TASK 8 METHODS

// Increment velocity for crack surface
short NodalPoint::IncrementDelvSideTask8(short side,int crackNumber,double fi,Vector *delv,double *surfaceMass,CrackSegment *seg) const
{
	// get velocity field to use for surface particle to node pair
	short vfld = GetFieldForSurfaceParticle(side,crackNumber,seg);
	
	// Exit if no field
	if(vfld<0) return false;
	if(!CrackVelocityField::ActiveNonrigidField(cvf[vfld])) return false;

	// increment the velocity if enough mass
	double fieldMass = GetNodalMass(true);
	if(cvf[vfld]->IncrementDelvTask8(fi,delv,&fieldMass))
	{	*surfaceMass+=fi;
		return true;
	}
	return false;
}

// Find velocity field to use when for crack surface particle at (x1,y1) to this node
// The side is a side of crack crackNumber. Look for that crack inthe fields first, otherwise
//		do a crack crossing calculation
// SCWarning - may look at [2] and [3] - perhaps two sections - one for one crack and one for more
short NodalPoint::GetFieldForSurfaceParticle(short side,int crackNumber,CrackSegment *seg) const
{
	int otherSide;
	short vfld=-1;
	double x1=seg->surfx[side-1];
	double y1=seg->surfy[side-1];
	
	// 1. Has field [1] with crackNumber
	//      Node with: [1], [1]&[2], [0]&[1], [1]&[3], [0]&[1]&[2], [0]&[1]&[3], [1]&[2]&[3], or [0]&[1]&[2]&[3]
    //      This block always finds vfld 0, 1, or 3
	if(CrackVelocityField::ActiveNonrigidField(cvf[1],crackNumber))
	{	if(side==cvf[1]->location(FIRST_CRACK))
		{	vfld=1;
			
			// maybe switch [1] to [3]
			if(CrackVelocityField::ActiveNonrigidField(cvf[3]))
			{	// if line crosses second crack in [3], switch to [3]
				int otherCrack=cvf[3]->OppositeCrackTo(crackNumber,side,&otherSide);
				if(otherCrack>0)
				{	if(SurfaceCrossesOneCrack(x1,y1,x,y,otherCrack)!=NO_CRACK)
						vfld=3;
				}
			}
			else if(SurfaceCrossesOtherCrack(x1,y1,x,y,crackNumber)!=NO_CRACK)
			{	// it crosses any other crack, need [3], but not there
				vfld=-1;
			}
		}
		else
		{	if(CrackVelocityField::ActiveNonrigidField(cvf[1]))
			   vfld=0;
			
			// maybe switch [0] to [2]
			if(CrackVelocityField::ActiveNonrigidField(cvf[2]))
			{	// if line crosses the crack found in [2], then switch to [2]
				if(SurfaceCrossesOneCrack(x1,y1,x,y,cvf[2]->crackNumber(FIRST_CRACK))!=NO_CRACK)
				   vfld=2;
			}
			else if(SurfaceCrossesOtherCrack(x1,y1,x,y,crackNumber)!=NO_CRACK)
			{	// it crosses any other, needs [2], but not there
				vfld=-1;
			}
		}
	}
	
	// 2. Has field [2] with crackNumber (may have [0], [1], and [3] as well)
    //      This block always finds vfld (even if not there)
	else if(CrackVelocityField::ActiveNonrigidField(cvf[2],crackNumber))
	{	if(side==cvf[2]->location(FIRST_CRACK))
		{	vfld=2;
	
			// maybe switch [2] to [3]
			if(CrackVelocityField::ActiveNonrigidField(cvf[3]))
			{	// if line crosses second crack in [3], switch to [3]
				int otherCrack=cvf[3]->OppositeCrackTo(crackNumber,side,&otherSide);
				if(otherCrack>0)
				{	if(SurfaceCrossesOneCrack(x1,y1,x,y,otherCrack)!=NO_CRACK)
						vfld=3;
				}
			}
			else if(SurfaceCrossesOtherCrack(x1,y1,x,y,crackNumber)!=NO_CRACK)
			{	// it crosses any other crack, need [3], but not there
				vfld=-1;
			}
		}
		else
		{	if(CrackVelocityField::ActiveNonrigidField(cvf[1]))
				vfld=0;
			
			// maybe switch [0] to [1]
			if(CrackVelocityField::ActiveNonrigidField(cvf[1]))
			{	// if line crosses the crack in [1], then switch to [1]
				if(SurfaceCrossesOneCrack(x1,y1,x,y,cvf[1]->crackNumber(FIRST_CRACK))!=NO_CRACK)
					vfld=1;
			}
			else if(SurfaceCrossesOtherCrack(x1,y1,x,y,crackNumber)!=NO_CRACK)
			{	// it crosses another crack, needs [1], but not there
				vfld=-1;
			}
		}
	}
	
	// 3. If not found in [1] or [2]
	else
	{	// Look in [3], but only use if it has this crack and a line cross the other crack too
		//		Most likely node with [3], [0]&[3], [1]&[3], [0]&[1]&[3] with crackNumber corresponding to missing [1] or [2]
		if(CrackVelocityField::ActiveNonrigidField(cvf[3]))
		{	// verify has correct field and retreive other crack
			int otherCrack=cvf[3]->OppositeCrackTo(crackNumber,side,&otherSide);
			if(otherCrack>0)
			{	// if crack found, then can use [3] if crosses the other crack in [3] too
				if(SurfaceCrossesOneCrack(x1,y1,x,y,otherCrack)!=NO_CRACK)
					vfld=3;
			}
		}
	
		// if still not found, slightly move (may not be needed) and then look for cracks other
		// then the one associated with this surface
		if(vfld<0)
		{	Vector moved=seg->SlightlyMovedIfNotMovedYet(side);
			CrackField cfld[2];
			SurfaceCrossesCracks(moved.x,moved.y,x,y,cfld);
			if(cfld[0].loc==NO_CRACK)
			{   // no crack found so use field [0]
				if(CrackVelocityField::ActiveNonrigidField(cvf[0]))
					vfld=0;
			}
			else if(cfld[1].loc==NO_CRACK)
			{	// found one crack, but only
				if(CrackVelocityField::ActiveNonrigidField(cvf[1]))
				{	if(cfld[0].crackNum==cvf[1]->crackNumber(FIRST_CRACK) && cfld[0].loc==cvf[1]->location(FIRST_CRACK))
						vfld=1;
				}
				else if(vfld<0 && CrackVelocityField::ActiveNonrigidField(cvf[2]))
				{	if(cfld[0].crackNum==cvf[2]->crackNumber(FIRST_CRACK) && cfld[0].loc==cvf[2]->location(FIRST_CRACK))
						vfld=2;
				}
			}
			else if(CrackVelocityField::ActiveNonrigidField(cvf[3]))
			{	// found two cracks, but only use if same two cracks and sides that are in [3]
				if(cfld[0].crackNum==cvf[3]->crackNumber(FIRST_CRACK) && cfld[1].crackNum==cvf[3]->crackNumber(SECOND_CRACK))
				{	if(cfld[0].loc==cvf[3]->location(FIRST_CRACK) && cfld[1].loc==cvf[3]->location(SECOND_CRACK))
						vfld=3;
				}
				else if(cfld[1].crackNum==cvf[3]->crackNumber(FIRST_CRACK) && cfld[0].crackNum==cvf[3]->crackNumber(SECOND_CRACK))
				{	if(cfld[1].loc==cvf[3]->location(FIRST_CRACK) && cfld[0].loc==cvf[3]->location(SECOND_CRACK))
						vfld=3;
				}
			}
		}
	}
	
	// return result. Positive result guarantees field that sees cracks, otherwise no field
	return vfld;
}
	
// Determine if line from crack surface particle (1) to node (2) crosses 
// one or more cracks
void NodalPoint::SurfaceCrossesCracks(double x1,double y1,double x2,double y2,CrackField *cfld) const
{
	CrackHeader *nextCrack=firstCrack;
	int cfound=0;
    short vfld;
	Vector norm;
	
	cfld[0].loc=NO_CRACK;			// NO_CRACK, ABOVE_CRACK, or BELOW_CRACK
	cfld[1].loc=NO_CRACK;
	
	while(nextCrack!=NULL)
	{	vfld=nextCrack->CrackCross(x1,y1,x2,y2,&norm);
		if(vfld!=NO_CRACK)
		{	cfld[cfound].loc=vfld;
			cfld[cfound].crackNum=nextCrack->GetNumber();
			cfound++;
			if(cfound>1) break;			// stop if found two, if there are more then two, physics may be off
		}
		nextCrack=(CrackHeader *)nextCrack->GetNextObject();
	}
}

// Determine if line from crack surface particle (1) to node (2) crosses specific crack
// one or more cracks
int NodalPoint::SurfaceCrossesOneCrack(double x1,double y1,double x2,double y2,int cnum) const
{
	CrackHeader *nextCrack=firstCrack;
	Vector norm;
	
	while(nextCrack!=NULL)
	{	if(nextCrack->GetNumber()==cnum)
			return nextCrack->CrackCross(x1,y1,x2,y2,&norm);
		nextCrack=(CrackHeader *)nextCrack->GetNextObject();
	}
	return NO_CRACK;
}

// Determine if line from crack surface particle (x1,y1) to node (x2,y22) crosses any crack other than the
// one provided
int NodalPoint::SurfaceCrossesOtherCrack(double x1,double y1,double x2,double y2,int cnum) const
{
	CrackHeader *nextCrack=firstCrack;
	Vector norm;
	
	while(nextCrack!=NULL)
	{	if(nextCrack->GetNumber()!=cnum)
		{	int cross = nextCrack->CrackCross(x1,y1,x2,y2,&norm);
			if(cross!=NO_CRACK) return cross;
		}
		nextCrack=(CrackHeader *)nextCrack->GetNextObject();
	}
	return NO_CRACK;
}

// Calculate CM velocity at a node and store in nv[0] velocity
//   (only used when contact.GetMoveOnlySurfaces() is FALSE, i.e., when crack particles
//		move in CM velocity field,  and then only when crack plane particles are about to move)
void NodalPoint::CalcCMVelocityTask8(void)
{
	Vector nodePk;
	ZeroVector(&nodePk);
	int hasParticles=0;
	double foundMass = 0.;
	for(int i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(cvf[i]))
		{	if(cvf[i]->CollectMomentaTask8(&nodePk,&foundMass))
				hasParticles = 1;
		}
	}
	
	// store in field zero
	if(hasParticles) ScaleVector(&nodePk,1./foundMass);
	cvf[0]->SetCMVelocityTask8(&nodePk,hasParticles);
}

// Get velocity for center of mass
//   (assumes recently called CalcCMVelocityTask8() for this node)
//   (only used when contact.GetMoveOnlySurfaces() is FALSE and crack plane
//       particles are moving)
bool NodalPoint::GetCMVelocityTask8(Vector *vk) const
{	return cvf[0]->GetCMVelocityTask8(vk);
}


#pragma mark INCREMENTERS

// Add displacements to selected field
void NodalPoint::AddDisplacement(short vfld,int matfld,double wt,Vector *pdisp)
{	cvf[vfld]->AddDisplacement(matfld,wt,pdisp);
}

// Add volume to selected field
void NodalPoint::AddVolume(short vfld,int matfld,double wtVol)
{	cvf[vfld]->AddVolume(matfld,wtVol);
}

#pragma mark VELOCITY FIELDS

// Calculate velocity at a node from current momentum and mass matrix
void NodalPoint::CalcVelocityForStrainUpdate(void)
{	// get velocity for all crack and material velocity fields
	int i;
    for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->CalcVelocityForStrainUpdate();
    }
}

// Get velocity for selected field
Vector NodalPoint::GetVelocity(short vfld,int matfld)
{   return cvf[vfld]->GetVelocity(matfld);
}

// Add total kinetic energy and mass on the grid for all velocity fields to supplied variables
// in nanoJ (mass in g, vel in mm/sec)
void NodalPoint::AddKineticEnergyAndMass(double &kineticEnergy,double &totalMass)
{	int i;
    for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->AddKineticEnergyAndMass(kineticEnergy,totalMass);
	}
}


// Get contact force for selected field
Vector NodalPoint::GetContactForce(short vfld,int matfld)
{	return cvf[vfld]->GetContactForce(matfld);
}

// Some all forces from rigid material velocity fields
Vector NodalPoint::GetTotalContactForce(bool clearForces)
{	Vector fcontact;
	ZeroVector(&fcontact);
    for(int i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->SumAndClearRigidContactForces(&fcontact,clearForces);
	}
	double scale = -UnitsController::Scaling(1.e-6)/timestep;
	fcontact.x*=scale;
	fcontact.y*=scale;
	fcontact.z*=scale;
	return fcontact;
}

#pragma mark MATERIAL CONTACT

// Called in multimaterial mode to check contact at nodes with multiple materials
// On first call in time step, first and last on pointers to MaterialInterfaceNode * because those
//		objects are created for later interface calculations
// postUpdate is TRUE when called between momentum update and particle update and otherwise is FALSE
// throws CommonException() on memory error or rigid material problem
void NodalPoint::MaterialContactOnNode(double deltime,int callType,MaterialInterfaceNode **first,MaterialInterfaceNode **last)
{
	// check each crack velocity field on this node
	for(int i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->MaterialContactOnCVF(this,i,deltime,callType,first,last);
	}
}

// retrieve volume gradient for matnum (1 based) in crack field only (or zero if
// not there or not tracked)
void NodalPoint::GetMatVolumeGradient(int matnum,Vector *grad) const
{
	// Default to zero and exit if not in multimaterial mode
	ZeroVector(grad);
	if(!fmobj->multiMaterialMode) return;
	
	// Only there if field zero is active
    if(CrackVelocityField::ActiveField(cvf[0]))
    {   // convert to field number
		int matfld = theMaterials[matnum-1]->GetField();
        if(cvf[0]->HasVolumeGradient(matfld))
            cvf[0]->GetVolumeGradient(matfld,this,grad,1.);
    }
}

// This node is known to have imperfect interface with forces in cvf[vfld] from material mati
// to matipaired (second material with the max volume)
void NodalPoint::MaterialInterfaceForce(MaterialInterfaceNode *mmnode)
{	
    // recall interface response from material interface node
    Vector fImpInt;
    double energy = mmnode->GetInterfaceTraction(&fImpInt);
    
    // add total force (in g mm/sec^2) to material field
    int vfld,mati,matipaired;
    mmnode->GetFieldInfo(&vfld, &mati, &matipaired);
    cvf[vfld]->AddFtotTask3(mati,&fImpInt);
    
    // add negative force to paired material (in a pair)
    if(matipaired>=0)
    {   ScaleVector(&fImpInt, -1.);
        cvf[vfld]->AddFtotTask3(matipaired,&fImpInt);
    }
    
    // add interface energy in units g-mm^2/sec^2 (multiply by 1e-9 to get J - kg-m^2/sec^2)
    // but only once for pair of interface nodes
    interfaceEnergy+=energy;
}

#pragma mark CRACK SURFACE CONTACT

// Look for crack contact and adjust accordingly.
// first and last are only non-NULL in mass and momentum task and this method will create
//      a crack node for any nodes that need to do crack contact
// poastUpdate is only TRUE when called in the momentum update (and the expectation is that
//      the forces at that time are the contact forces and used for friction)
void NodalPoint::CrackContact(bool postUpdate,double deltime,CrackNode **first,CrackNode **last)
{	// Nothing to do if not near a crack contact surface: Possible fields are
	//  1. Those with no contacts: [0], [1], [3], [0]&[3], [1]&[2]
	//  2. Those with contacts: [0]&[1], [1]&[3], [0]&[1]&[2], [0]&[1]&[3], [1]&[2]&[3], and [0]&[1]&[2]&[3]
	//  3. Never occurs [2], [0]&[2], [2]&[3], [0]&[2]&[3]
	
	// exit if no contact (note thgat mirrored materials only in nonrigid fields that have other nonrigid material too,
	//    thus sufficient to look for any non-rigid particles for active crack fields)
	// SCWarning - has2 and has3 should be set to false in single crack mode
	bool has1=CrackVelocityField::ActiveNonrigidField(cvf[1]);
	bool has2=CrackVelocityField::ActiveNonrigidField(cvf[2]);
	if(!has1 && !has2) return;	// True for [0], [3], and [0]&[3]
	bool has0=CrackVelocityField::ActiveNonrigidField(cvf[0]);
	bool has3=CrackVelocityField::ActiveNonrigidField(cvf[3]);
	if(!has0 && !has3) return;	// True for [1] and [1]&[2]
	
	// store references to this node for future use
	if(first!=NULL)
	{	*last = new CrackNode(this,*last);
		if(*last==NULL) throw CommonException("Memory error allocating storage for a crack node.",
															   "NodalPoint::CrackContact");
		if(*first==NULL) *first = *last;
	}
	
	// between [0] and [1]
	int cnum,cabove;
	if(has0 && has1)
	{	cnum=cvf[1]->crackNumber(FIRST_CRACK);
		if(contact.HasContact(cnum))
		{	cabove=(cvf[1]->location(FIRST_CRACK)==ABOVE_CRACK) ? 1 : 0;
			AdjustContact(cabove,1-cabove,&(cvf[1]->norm[FIRST_CRACK]),cnum,postUpdate,deltime);
		}
	}
	
	// between [0] & [2]
	if(has2 && has0)
	{	cnum=cvf[2]->crackNumber(FIRST_CRACK);
		if(contact.HasContact(cnum))
		{	cabove=(cvf[2]->location(FIRST_CRACK)==ABOVE_CRACK) ? 2 : 0;
			AdjustContact(cabove,2-cabove,&(cvf[2]->norm[FIRST_CRACK]),cnum,postUpdate,deltime);
		}
	}
	
	// with [3]
	if(has3)
	{	// between [1] & [3]
		if(has1) CrackContactThree(1,postUpdate,deltime);

		// between [2] & [3]
		if(has2) CrackContactThree(2,postUpdate,deltime);
	}
	
}

// Contact between field [single] and field [3] when both fields are present
// We expect [3] to be on same side as the crack for [single], thus do contact
//    between [single] and [3] for the crack in [3] that does not match the
//    crack in [single]. In other words, [single] acts like field [0] relative
//    to field [3] crossing that other crack
// Caller must verify that [3] is not empty
void NodalPoint::CrackContactThree(int single,bool postUpdate,double deltime)
{
	// get common crack
	int cnum=cvf[single]->crackNumber(FIRST_CRACK);
	
	// We expect [3] to be on opposite side of crack that DOES NOT match [single]
	// find the other crack
	int otherCrack;
	if(cnum == cvf[3]->crackNumber(FIRST_CRACK))
		otherCrack = SECOND_CRACK;
	else if(cnum == cvf[3]->crackNumber(SECOND_CRACK))
		otherCrack = FIRST_CRACK;
	else
	{	// Field [3] for two other cracks, so skip this calculation
		return;
	}
	
	// skip if otherCrack not doing contact
	if(!contact.HasContact(otherCrack)) return;

	// get above field
	int cabove,cbelow;
	if(cvf[3]->location(otherCrack)==ABOVE_CRACK)
	{	cabove = 3;
		cbelow = single;
	}
	else
	{	cabove=single;
		cbelow=3;
	}
	
	// adjust contact
	AdjustContact(cabove,cbelow,&(cvf[3]->norm[otherCrack]),otherCrack,postUpdate,deltime);
}

// Look for crack contact and adjust accordingly - a for field above and b for field below and both
// fields must be verified as present (1 or more points)
void NodalPoint::AdjustContact(short a,short b,Vector *norm,int crackNumber,bool postUpdate,double deltime)
{
	// see if in contact and get change in momentum
    Vector delP;
	int inContact;
	bool changeMomentum = contact.GetDeltaMomentum(this,&delP,cvf[a],cvf[b],norm,crackNumber,postUpdate,deltime,&inContact);
	
	// exit if no momentum change
	if(!changeMomentum) return;
	
	// on post update contact, do not change nodes with boundary conditions
	if(postUpdate && (fixedDirection&XYZ_SKEWED_DIRECTION))
	{	if(fixedDirection&X_DIRECTION) delP.x=0.;
		if(fixedDirection&Y_DIRECTION) delP.y=0.;
		if(fixedDirection&Z_DIRECTION) delP.z=0.;
	}
	
    // change momenta
	cvf[a]->ChangeCrackMomentum(&delP,postUpdate,deltime);
	Vector delPb;
    cvf[b]->ChangeCrackMomentum(CopyScaleVector(&delPb,&delP,-1.),postUpdate,deltime);
}

// Look for crack contact and adjust accordingly
void NodalPoint::CrackInterfaceForce(void)
{	// Nothing to do if not near a crack contact surface: Possible fields are
	//  1. Those with no contacts: [0], [1], [3], [0]&[3], [1]&[2]
	//  2. Those with contacts: [0]&[1], [1]&[3], [0]&[1]&[2], [0]&[1]&[3], [1]&[2]&[3], and [0]&[1]&[2]&[3]
	//  3. Never occurs [2], [0]&[2], [2]&[3], [0]&[2]&[3]
	
	// skip those with no contact
	// SCWarning - has2 and has3 should be set to false in single crack mode
	bool has1=CrackVelocityField::ActiveNonrigidField(cvf[1]);
	bool has2=CrackVelocityField::ActiveNonrigidField(cvf[2]);
	if(!has1 && !has2) return;		// True for [0], [3], and [0]&[3]
	bool has0=CrackVelocityField::ActiveNonrigidField(cvf[0]);
	bool has3=CrackVelocityField::ActiveNonrigidField(cvf[3]);
	if(!has0 && !has3) return;	// True for [1] and [1]&[2]
    
	
	// between [0] and [1] across first crack
	int cnum,cabove;
	if(has0 && has1)
	{	cnum=cvf[1]->crackNumber(FIRST_CRACK);
		if(contact.IsImperfect(cnum))
		{	cabove=(cvf[1]->location(FIRST_CRACK)==ABOVE_CRACK) ? 1 : 0;
			AddInterfaceForceOnCrack(cabove,1-cabove,&(cvf[1]->norm[FIRST_CRACK]),cnum);
		}
	}
			
	// between [0] and [2] across second crack
	if(has2 && has0)
	{	cnum=cvf[2]->crackNumber(FIRST_CRACK);
		if(contact.IsImperfect(cnum))
		{	cabove=(cvf[2]->location(FIRST_CRACK)==ABOVE_CRACK) ? 2 : 0;
			AddInterfaceForceOnCrack(cabove,2-cabove,&(cvf[2]->norm[FIRST_CRACK]),cnum);
		}
	}
	
	// with [3]
	if(has3)
	{	// between [1] & [3]
		if(has1) InterfaceForceInteractingCracks(1);
	
		// between [2] & [3]
		if(has2) InterfaceForceInteractingCracks(2);
	}
}

// Interface force between field [single] and field [3] when both fields are present
// We expect [3] to be on same side as the crack for [single], thus do contact
//    between [single] and [3] for the crack in [3] that does not match the
//    crack in [single]. In other words, [single] acts like field [0] relative
//    to field [3] crossing that other crack
// Caller must verify that [3] is not emptuy
void NodalPoint::InterfaceForceInteractingCracks(int single)
{
	// get common crack
	int cnum=cvf[single]->crackNumber(FIRST_CRACK);
	
	// We expect [3] to be on opposite side of crack that DOES NOT match [single]
	// find the other crack
	int otherCrack;
	if(cnum == cvf[3]->crackNumber(FIRST_CRACK))
		otherCrack = SECOND_CRACK;
	else if(cnum == cvf[3]->crackNumber(SECOND_CRACK))
		otherCrack = FIRST_CRACK;
	else
	{	// Field [3] for two other cracks, so skip this calculation
		return;
	}
	
	// skip if otherCrack not an interface
	if(!contact.IsImperfect(otherCrack)) return;
	
	// get above field
	int cabove,cbelow;
	if(cvf[3]->location(otherCrack)==ABOVE_CRACK)
	{	cabove = 3;
		cbelow = single;
	}
	else
	{	cabove=single;
		cbelow=3;
	}
	
	// do interface force
	AddInterfaceForceOnCrack(cabove,cbelow,&(cvf[3]->norm[otherCrack]),otherCrack);
}

// Look for cracks as imperfect interfaces and adjust accordingly - a for field above and b for field below
// fields must be verified as present (1 or more points)
void NodalPoint::AddInterfaceForceOnCrack(short a,short b,Vector *norm,int crackNumber)
{	Vector fImpInt;
	
    // Use contact laws to change momenta - returns TRUE or FALSE if adjustment was made
	double rawEnergy;
	if(!contact.GetInterfaceForceOnCrack(this,&fImpInt,cvf[a],cvf[b],norm,crackNumber,&rawEnergy,x))
		return;
	
	// add total force (in g mm/sec^2)
	AddFtotSpreadTask3(a,MakeVector(fImpInt.x,fImpInt.y,0.));
	AddFtotSpreadTask3(b,MakeVector(-fImpInt.x,-fImpInt.y,0.));
	
	// add interface energy in units g-mm^2/sec^2 (multiply by 1e-9 to get J - kg-m^2/sec^2)
	interfaceEnergy += rawEnergy;
}

#pragma mark ACCESSORS

// number of particles for this node
int NodalPoint::NumberParticles(void)
{	int totalParticles=0;
	for(int i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			totalParticles+=cvf[i]->GetNumberPoints();
	}
	return totalParticles;
}

// number of crack fields with non rigid particles
int NodalPoint::NumberNonrigidCracks(void)
{	int totalCracks=0;
	for(int i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(cvf[i]))
		totalCracks++;
	}
	return totalCracks;
}

// Look for presence of non rigid point on this node
bool NodalPoint::NodeHasNonrigidParticles(void) const
{	
	for(int i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
		{	if(cvf[i]->HasPointsNonrigid())
				return true;
		}
	}
	return false;
}

// describe velocity field
void NodalPoint::Describe(void) const
{	cout << "#  node=" << num << " pt=(" << x << "," << y << "," << z << ") mass=" << nodalMass << endl;
	
    bool active = false;
	for(int i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
        {   if(!active)
            {   cout << "#  active crack velocity fields:" << endl;
                active=true;
            }
			cout << "#  " << i << ". ";
			cvf[i]->Describe();
		}
	}
    if(!active) cout << "#  no active crack velocity fields (might be in initialization)" << endl;
}
	
// total nodal mass
// argument not used here, but is in OSParticulas
double NodalPoint::GetNodalMass(bool requireCracks) const { return nodalMass; }

#pragma mark BOUNDARY CONDITION METHODS

// set one component of velocity and momentum to zero
void NodalPoint::SetMomVel(Vector *norm)
{
	for(int i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
            cvf[i]->SetMomVel(norm);
	}
}

// Add one component of velocity and momentum at a node (assumes mass already set)
void NodalPoint::AddMomVel(Vector *norm,double vel)
{	
	for(int i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
            cvf[i]->AddMomVel(norm,vel);
	}
}

// Reflect one component of velocity and momentum from a node
void NodalPoint::ReflectMomVel(Vector *norm,NodalPoint *ndptr)
{
	// only field zero, which assumes no cracks near the symmetry plane
	if(CrackVelocityField::ActiveField(cvf[0]))
	{	if(CrackVelocityField::ActiveField(ndptr->cvf[0]))
			cvf[0]->ReflectMomVel(norm,ndptr->cvf[0]);
	}
}


// set force in direction norm to -p(interpolated)/time such that updated momentum
//    of pk.i + deltime*ftot.i will be zero in that direction
void NodalPoint::SetFtotDirection(Vector *norm,double deltime,Vector *freaction)
{	
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
            cvf[i]->SetFtotDirection(norm,deltime,freaction);
	}
}

// set one component of force such that updated momentum will be mass*velocity
void NodalPoint::AddFtotDirection(Vector *norm,double deltime,double vel,Vector *freaction)
{	
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
            cvf[i]->AddFtotDirection(norm,deltime,vel,freaction);
	}
}

// set one component of force such that updated momentum will match reflected node
void NodalPoint::ReflectFtotDirection(Vector *norm,double deltime,NodalPoint *ndptr,Vector *freaction)
{
	// only field zero, which assumes no cracks near the symmetry plane
	if(CrackVelocityField::ActiveField(cvf[0]))
	{	if(CrackVelocityField::ActiveField(ndptr->cvf[0]))
			cvf[0]->ReflectFtotDirection(norm,deltime,ndptr->cvf[0],freaction);
	}
}

// Mark a direction as fixed by velocity BC
// Assume 1 means x, 2 means y, 4 means z
void NodalPoint::SetFixedDirection(int dir)
{	fixedDirection|=dir;
}

// Unmark a direction as fixed by velocity BC
// Assume 1 means x, 2 means y, 4 means z
void NodalPoint::UnsetFixedDirection(int dir)
{	if(fixedDirection&dir)
		fixedDirection^=dir;
}

#pragma mark CLASS METHODS

// zero all velocity fields at start of time step
void NodalPoint::PreliminaryCalcs(void)
{	int i;
    for(i=1;i<=nnodes;i++)
        nd[i]->PrepareForFields();
}

// Find Grid point velocities
void NodalPoint::GetGridVelocitiesForStrainUpdate(void)
{	int i;
    for(i=1;i<=nnodes;i++)
        nd[i]->CalcVelocityForStrainUpdate();
}

// Find Grid CM velocities (only for cracks when when contact.GetMoveOnlySurfaces() is FALSE)
void NodalPoint::GetGridCMVelocitiesTask8(void)
{	
#pragma omp parallel for
    for(int i=1;i<=nnodes;i++)
        nd[i]->CalcCMVelocityTask8();
}


