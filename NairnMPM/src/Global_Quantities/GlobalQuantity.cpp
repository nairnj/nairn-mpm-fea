/********************************************************************************
    GlobalQuantity.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Mon Jan 12 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved. 
	
	Adding Global Quantity
		1. Add tag in GlobalQuantity.hpp
		2. Respond to string in GlobalQuantity(char *,int) initializer
		3. Add calculations in AppendQuantity()
********************************************************************************/

#include "stdafx.h"
#include "Global_Quantities/GlobalQuantity.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/Reservoir.hpp"
#include "Materials/MaterialBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "System/ArchiveData.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "System/UnitsController.hpp"
#include "Elements/ElementBase.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"

// Single global contact law object
GlobalQuantity *firstGlobal=NULL;
static GlobalQuantity *lastGlobal=NULL;
static int numGlobal=0;
static char quote='"';

/*******************************************************************
	GlobalQuantity: Constructors and Destructors
*******************************************************************/

// Constructors
GlobalQuantity::GlobalQuantity()
{
}

// throws std::bad_alloc
GlobalQuantity::GlobalQuantity(char *quant,int whichOne)
{
	char nameStr[200];
	size_t nameSize=200;
	whichMat=whichOne;
	
	quantity = DecodeGlobalQuantity(quant,&subcode,&whichMat);
	
	// set name
	if(whichMat!=0)
		snprintf(nameStr,nameSize,"%s mat %d",quant,whichMat);
	else
		strcpy(nameStr,quant);
 	FinishNewQuantity(nameStr);
}

// throws std::bad_alloc
GlobalQuantity::GlobalQuantity(char *quant,Vector *ptLoc)
{
	char nameStr[200];
	size_t nameSize=200;
	
	quantity = DecodeGlobalQuantity(quant,&subcode,&whichMat);
	
	// save room for particle number
    snprintf(nameStr,nameSize,"%s pt ",quant);
	FinishNewQuantity(nameStr);
	
	// store point and any material allowed
	ptPos = ptLoc;
	whichMat = 0;
}

// Find closet point to tracer particle coordinates
// Not that ptNum is zero based number
GlobalQuantity *GlobalQuantity::FindTracerParticle(int p,Vector *ploc)
{
	// skip if not a point
	if(ptPos==NULL) return nextGlobal;
	
	// find distance
	double dx = ploc->x-ptPos->x;
	double dy = ploc->y-ptPos->y;
	double dz = fmobj->IsThreeD() ? ploc->z-ptPos->z : 0.;
	double dist = sqrt(dx*dx+dy*dy+dz*dz);
	if(dist<minDist)
	{	ptNum = p;
		minDist = dist;
	}
	
	return nextGlobal;
}

// finish global quantity settings
void GlobalQuantity::FinishNewQuantity(char *nameStr)
{
	// store quantity name
	name=new char[strlen(nameStr)+1];
	strcpy(name,nameStr);

	// set color ID
	colorID=numGlobal % 10;

	// this object is currently the last one
	SetNextGlobal(NULL);

	// adjust previous global quantity or set firstGlobal if this is the first one
	if(lastGlobal!=NULL)
		lastGlobal->SetNextGlobal(this);
	else
		firstGlobal=this;
	lastGlobal=this;

	// count the number of objects
	numGlobal++;
	
	// assume not a point
	ptPos = NULL;
	ptNum = -1;
	minDist = 1.e99;
}

// decode quant it to quantity ID and subcode (used for history variables)
// mcode can change whichMat if needed
int GlobalQuantity::DecodeGlobalQuantity(const char *quant,int *hcode,int *mcode)
{
	int theQuant;
	
	// set quantity and subcode
	*hcode = 0;
	if(strcmp(quant,"sxx")==0 || strcmp(quant,"sRR")==0)
		theQuant=AVG_SXX;
	else if(strcmp(quant,"syy")==0 || strcmp(quant,"sZZ")==0)
		theQuant=AVG_SYY;
	else if(strcmp(quant,"sxy")==0 || strcmp(quant,"sRZ")==0)
		theQuant=AVG_SXY;
	else if(strcmp(quant,"szz")==0 || strcmp(quant,"sTT")==0)
		theQuant=AVG_SZZ;
	else if(strcmp(quant,"sxz")==0)
		theQuant=AVG_SXZ;
	else if(strcmp(quant,"syz")==0)
		theQuant=AVG_SYZ;
	
	else if(strcmp(quant,"exx")==0 || strcmp(quant,"eRR")==0)
		theQuant=AVG_EXX;
	else if(strcmp(quant,"eyy")==0 || strcmp(quant,"eZZ")==0)
		theQuant=AVG_EYY;
	else if(strcmp(quant,"exy")==0 || strcmp(quant,"eRZ")==0)
		theQuant=AVG_EXY;
	else if(strcmp(quant,"ezz")==0 || strcmp(quant,"eTT")==0)
		theQuant=AVG_EZZ;
	else if(strcmp(quant,"exz")==0)
		theQuant=AVG_EXZ;
	else if(strcmp(quant,"eyz")==0)
		theQuant=AVG_EYZ;
	
	else if(strcmp(quant,"Fxx")==0 || strcmp(quant,"FRR")==0)
		theQuant=AVG_FXX;
	else if(strcmp(quant,"Fxy")==0 || strcmp(quant,"FRZ")==0)
		theQuant=AVG_FXY;
	else if(strcmp(quant,"Fxz")==0)
		theQuant=AVG_FXZ;
	else if(strcmp(quant,"Fyx")==0 || strcmp(quant,"FZR")==0)
		theQuant=AVG_FYX;
	else if(strcmp(quant,"Fyy")==0 || strcmp(quant,"FZZ")==0)
		theQuant=AVG_FYY;
	else if(strcmp(quant,"Fyz")==0)
		theQuant=AVG_FYZ;
	else if(strcmp(quant,"Fzx")==0)
		theQuant=AVG_FZX;
	else if(strcmp(quant,"Fzy")==0)
		theQuant=AVG_FZY;
	else if(strcmp(quant,"Fzz")==0 || strcmp(quant,"FTT")==0)
		theQuant=AVG_FZZ;
	
	else if(strcmp(quant,"exxe")==0 || strcmp(quant,"eRRe")==0)
		theQuant=AVG_EXXE;
	else if(strcmp(quant,"eyye")==0 || strcmp(quant,"eZZe")==0)
		theQuant=AVG_EYYE;
	else if(strcmp(quant,"exye")==0 || strcmp(quant,"eRZe")==0)
		theQuant=AVG_EXYE;
	else if(strcmp(quant,"ezze")==0 || strcmp(quant,"eTTe")==0)
		theQuant=AVG_EZZE;
	else if(strcmp(quant,"exze")==0)
		theQuant=AVG_EXZE;
	else if(strcmp(quant,"eyze")==0)
		theQuant=AVG_EYZE;
	
	else if(strcmp(quant,"exxp")==0 || strcmp(quant,"eRRp")==0)
		theQuant=AVG_EXXP;
	else if(strcmp(quant,"eyyp")==0 || strcmp(quant,"eZZp")==0)
		theQuant=AVG_EYYP;
	else if(strcmp(quant,"exyp")==0 || strcmp(quant,"eRZp")==0)
		theQuant=AVG_EXYP;
	else if(strcmp(quant,"ezzp")==0 || strcmp(quant,"eTTp")==0)
		theQuant=AVG_EZZP;
	else if(strcmp(quant,"exzp")==0)
		theQuant=AVG_EXZP;
	else if(strcmp(quant,"eyzp")==0)
		theQuant=AVG_EYZP;
	
	else if(strcmp(quant,"Kinetic Energy")==0)
		theQuant=KINE_ENERGY;
	else if(strcmp(quant,"Grid Kinetic Energy")==0)
		theQuant=GRID_KINE_ENERGY;
	else if(strcmp(quant,"Strain Energy")==0)
		theQuant=STRAIN_ENERGY;
	else if(strcmp(quant,"Heat Energy")==0 || strcmp(quant,"Thermal Energy")==0)
		theQuant=HEAT_ENERGY;
	else if(strcmp(quant,"Entropy")==0)
		theQuant=ENTROPY_ENERGY;
	else if(strcmp(quant,"Internal Energy")==0)
		theQuant=INTERNAL_ENERGY;
	else if(strcmp(quant,"Helmholz Energy")==0)
		theQuant=HELMHOLZ_ENERGY;
	else if(strcmp(quant,"Interface Energy")==0)
		theQuant=INTERFACE_ENERGY;
	else if(strcmp(quant,"Friction Work")==0)
		theQuant=FRICTION_WORK;
	else if(strcmp(quant,"Work Energy")==0)
		theQuant=WORK_ENERGY;
	else if(strcmp(quant,"Plastic Energy")==0)
		theQuant=PLAS_ENERGY;
	
	else if(strcmp(quant,"velx")==0 || strcmp(quant,"velR")==0)
		theQuant=AVG_VELX;
	else if(strcmp(quant,"vely")==0 || strcmp(quant,"velZ")==0)
		theQuant=AVG_VELY;
	else if(strcmp(quant,"velz")==0)
		theQuant=AVG_VELZ;
	
	else if(strcmp(quant,"dispx")==0 || strcmp(quant,"dispR")==0)
		theQuant=AVG_DISPX;
	else if(strcmp(quant,"dispy")==0 || strcmp(quant,"dispZ")==0)
		theQuant=AVG_DISPY;
	else if(strcmp(quant,"dispz")==0)
		theQuant=AVG_DISPZ;
	
	else if(strcmp(quant,"temp")==0 || strcmp(quant,"Temp")==0 || strcmp(quant,"Temperature")==0 || strcmp(quant,"temperature")==0)
		theQuant=AVG_TEMP;
	else if(strcmp(quant,"concentration")==0 || strcmp(quant,"porepressure")==0)
		theQuant=WTFRACT_CONC;
    else if(strcmp(quant,"concFlux")==0)
        theQuant=TOT_FLUX;
	else if(strcmp(quant,"Step number")==0)
		theQuant=STEP_NUMBER;
	else if(strcmp(quant,"CPU time")==0)
		theQuant=CPU_TIME;
	else if(strcmp(quant,"Elapsed time")==0)
		theQuant=ELAPSED_TIME;
	else if(strcmp(quant,"alpha")==0)
		theQuant=GRID_ALPHA;
	else if(strcmp(quant,"palpha")==0)
		theQuant=PARTICLE_ALPHA;
	
	else if(strcmp(quant,"contactx")==0 || strcmp(quant,"contactR")==0)
		theQuant=TOT_FCONX;
	else if(strcmp(quant,"contacty")==0 || strcmp(quant,"contactZ")==0)
		theQuant=TOT_FCONY;
	else if(strcmp(quant,"contactz")==0)
		theQuant=TOT_FCONZ;
	
	else if(strcmp(quant,"reactionx")==0 || strcmp(quant,"reactionR")==0)
		theQuant=TOT_REACTX;
	else if(strcmp(quant,"reactiony")==0 || strcmp(quant,"reactionZ")==0)
		theQuant=TOT_REACTY;
	else if(strcmp(quant,"reactionz")==0)
		theQuant=TOT_REACTZ;
	
	else if(strcmp(quant,"px")==0 || strcmp(quant,"pR")==0)
		theQuant=LINMOMX;
	else if(strcmp(quant,"py")==0 || strcmp(quant,"pZ")==0)
		theQuant=LINMOMY;
	else if(strcmp(quant,"pZ")==0)
		theQuant=LINMOMZ;
	
	else if(strcmp(quant,"Lx")==0)
		theQuant=ANGMOMX;
	else if(strcmp(quant,"Ly")==0)
		theQuant=ANGMOMY;
	else if(strcmp(quant,"Lz")==0)
		theQuant=ANGMOMZ;

	// only nonzero if add particle spin
	else if(strcmp(quant,"Lpx")==0)
		theQuant=LPMOMX;
	else if(strcmp(quant,"Lpy")==0)
		theQuant=LPMOMY;
	else if(strcmp(quant,"Lpz")==0)
		theQuant=LPMOMZ;
	else if(strcmp(quant,"wpx")==0)
		theQuant=ANGVELX;
	else if(strcmp(quant,"wpy")==0)
		theQuant=ANGVELY;
	else if(strcmp(quant,"wpz")==0)
		theQuant=ANGVELZ;
	
	else if(strcmp(quant,"heatWatts")==0)
		theQuant=TOT_REACTQ;
	
	else if(strcmp(quant,"Decohesion")==0)
		theQuant=DECOHESION;
	
	else if(strcmp(quant,"ReservoirSize")==0)
		theQuant=RESERVOIR_SIZE;
	
	else
	{	theQuant=UNKNOWN_QUANTITY;
		
		// possibly a history variable which must be "history n"
		if(strlen(quant)>7)
		{	char nameStr[200];
			strcpy(nameStr,quant);
			nameStr[7]=0;
			if(strcmp(nameStr,"history")==0)
			{	sscanf(quant,"%*s %d",hcode);
				theQuant=HISTORY_VARIABLE;
			}
		}
        
        // possibly a concentration flux which must be "concFlux n"
        if(theQuant==UNKNOWN_QUANTITY && strlen(quant)>8)
        {   char nameStr[200];
            strcpy(nameStr,quant);
            nameStr[8] = 0;
            if(strcmp(nameStr,"concFlux")==0)
            {   sscanf(quant,"%*s %d",hcode);
                theQuant = TOT_FLUX;
             }
        }
        
        // possibly a crack length which must be "crack length n"
        if(theQuant==UNKNOWN_QUANTITY && strlen(quant)>12)
        {   char nameStr[200];
            strcpy(nameStr,quant);
            nameStr[12] = 0;
            if(strcmp(nameStr,"crack length")==0)
            {   sscanf(quant,"%*s %*s %d",hcode);
                theQuant = CRACK_LENGTH;
                *mcode = 0;         // no material allowed
            }
        }
        
        // possibly a concentration which must be "concentration n"
        if(theQuant==UNKNOWN_QUANTITY && strlen(quant)>13)
        {   char nameStr[200];
            strcpy(nameStr,quant);
            nameStr[13] = 0;
            if(strcmp(nameStr,"concentration")==0)
            {   sscanf(quant,"%*s %d",hcode);
                theQuant = WTFRACT_CONC;
            }
        }
        
       // possibly a debonded crack length which must be "debonded crack length n"
        if(theQuant==UNKNOWN_QUANTITY && strlen(quant)>21)
        {   char nameStr[200];
            strcpy(nameStr,quant);
            nameStr[21] = 0;
            if(strcmp(nameStr,"debonded crack length")==0)
            {   sscanf(quant,"%*s %*s %*s %d",hcode);
                theQuant = DEBONDED_CRACK_LENGTH;
                *mcode = 0;         // no material allowed
            }
        }
	}
	
	return theQuant;
}

#pragma mark GlobalQuantity: Methods

// append quantity
GlobalQuantity *GlobalQuantity::AppendName(char *fline)
{
	char nameStr[100];
	size_t nameSize=100;
	
	if(quantity==UNKNOWN_QUANTITY || quantity==DECOHESION)
		return nextGlobal;
	
	snprintf(nameStr,nameSize,"\t%c%s%c",quote,name,quote);
	strcat(fline,nameStr);
	return nextGlobal;
}

// append quantity
GlobalQuantity *GlobalQuantity::AppendQuantity(vector<double> &toArchive)
{
	int p;
	double value=0.,rho0,Vtot=0.,mp,Jp,Vp;
	int matid,qid=0,p0=0,pend=nmpms;
	bool singleParticle = false;
	if(ptNum>=0)
	{	p0 = ptNum;
		pend = ptNum+1;
		singleParticle = true;
	}
	
	switch(quantity)
	{   // stresses (MPa in Legacy)
		case AVG_SZZ:
			qid=ZZ;
		case AVG_SXZ:
			if(quantity==AVG_SXZ) qid=XZ;
		case AVG_SYZ:
			if(quantity==AVG_SYZ) qid=YZ;
	    case AVG_SXX:
			if(quantity==AVG_SXX) qid=XX;
		case AVG_SYY:
			if(quantity==AVG_SYY) qid=YY;
		case AVG_SXY:
			if(quantity==AVG_SXY) qid=XY;
			
			// Volume weighted average is Sum (Vp rhop stressp) / Sum Vp = Sum (mp stressp) / Sum Vp
			// where Vp = J mp/rho0
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid=mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					mp = mpm[p]->mp;
                    Tensor sp = mpm[p]->ReadStressTensor();
 					value += mp*Tensor_i(&sp,qid);
					Vtot += Jp*mp/rho0;
				}
			}
			if(Vtot>0.) value /= Vtot;
			value *= UnitsController::Scaling(1.e-6);
			break;
		
		// Elastic strain (% in Legacy Units)
		// New method small strain = Biot strain - archived plastic strain
		// New method hyperelastic = Biot strain from elastic B in plastic strain
		// Membranes = 0
		case AVG_EZZE:
			qid=ZZ;
	    case AVG_EXZE:
			if(quantity==AVG_EXZE) qid=XZ;
	    case AVG_EYZE:
			if(quantity==AVG_EYZE) qid=YZ;
	    case AVG_EXXE:
			if(quantity==AVG_EXXE) qid=XX;
		case AVG_EYYE:
			if(quantity==AVG_EYYE) qid=YY;
		case AVG_EXYE:
			if(quantity==AVG_EXYE) qid=XY;
			
			// Volume weighted average is Sum (Vp strainp) / Sum Vp 
			// where Vp = J mp/rho0
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid=mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					if(theMaterials[matid]->AltStrainContains()==ENG_BIOT_PLASTIC_STRAIN)
					{	// elastic strain = total strain minus plastic strain
						Matrix3 biot = mpm[p]->GetBiotStrain();
						Tensor *eplast=mpm[p]->GetAltStrainTensor();
						value += Vp*(biot.get(qid,2.) - Tensor_i(eplast,qid));
					}
					else if(theMaterials[matid]->AltStrainContains()==LEFT_CAUCHY_ELASTIC_B_STRAIN)
					{	// get elastic strain from elastic B in alt strain
						Matrix3 biot = mpm[p]->GetElasticBiotStrain();
						value += Vp*biot.get(qid,2.);
					}
					else
					{	// elastic strain = total strain for these materials
						Matrix3 biot = mpm[p]->GetBiotStrain();
						value += Vp*biot.get(qid,2.);
					}
					Vtot += Vp;
				}
			}
			if(Vtot>0.) value /= Vtot;
			value *= UnitsController::Scaling(100.);
			break;

		// plastic strain (% in Legacy)
		// New method small strain = archived plastic strain
		// New method hyperelastic = Biot strain from F - Biot strain from elastic B in plastic strain
		// Membrane = 0
		case AVG_EZZP:
			qid=ZZ;
		case AVG_EXZP:
			if(quantity==AVG_EXZP) qid=XZ;
		case AVG_EYZP:
			if(quantity==AVG_EYZP) qid=YZ;
		case AVG_EXXP:
			if(quantity==AVG_EXXP) qid=XX;
		case AVG_EYYP:
			if(quantity==AVG_EYYP) qid=YY;
		case AVG_EXYP:
			if(quantity==AVG_EXYP) qid=XY;
			
			// Volume weighted average is Sum (Vp strainp) / Sum Vp
			// where Vp = J mp/rho0
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid=mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					if(theMaterials[matid]->AltStrainContains()==ENG_BIOT_PLASTIC_STRAIN)
					{	// plastic strain all ready
						Tensor *eplast=mpm[p]->GetAltStrainTensor();
						value += Vp*Tensor_i(eplast,qid);
					}
					else if(theMaterials[matid]->AltStrainContains()==LEFT_CAUCHY_ELASTIC_B_STRAIN)
					{	// plastic strain = total strain minus elastic strain (from B)
						Matrix3 biotTot = mpm[p]->GetBiotStrain();
						Matrix3 biotElastic = mpm[p]->GetElasticBiotStrain();
						value += Vp*(biotTot.get(qid,2.) - biotElastic.get(qid,2.)) ;
					}
					// others have zero plastic strain
					Vtot += Vp;;
				}
			}
			if(Vtot>0.) value /= Vtot;
			value *= UnitsController::Scaling(100.);
			break;

		// total strain (% in Legacy)
		// New method small strain = Biot strain from F
		// New method hyperelastic = Biot strain from F
		case AVG_EZZ:
			qid=ZZ;
		case AVG_EXZ:
			if(quantity==AVG_EXZ) qid=XZ;
		case AVG_EYZ:
			if(quantity==AVG_EYZ) qid=YZ;
		case AVG_EXX:
			if(quantity==AVG_EXX) qid=XX;
		case AVG_EYY:
			if(quantity==AVG_EYY) qid=YY;
		case AVG_EXY:
			if(quantity==AVG_EXY) qid=XY;
			
			// Volume weighted average is Sum (Vp strainp) / Sum Vp
			// where Vp = J mp/rho0
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid=mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					Matrix3 biot = mpm[p]->GetBiotStrain();
					value += Vp*biot.get(qid,2.);
					Vtot += Vp;
				}
			}
			if(Vtot>0.) value /= Vtot;
			value *= UnitsController::Scaling(100.);
			break;
		
		// Deformation gradient
		case AVG_FZZ:
			qid=ZZ;
		case AVG_FXZ:
			if(quantity==AVG_FXZ) qid=XZ;
		case AVG_FYZ:
			if(quantity==AVG_FYZ) qid=YZ;
		case AVG_FXX:
			if(quantity==AVG_FXX) qid=XX;
		case AVG_FYY:
			if(quantity==AVG_FYY) qid=YY;
		case AVG_FXY:
			if(quantity==AVG_FXY) qid=XY;
		case AVG_FYX:
			if(quantity==AVG_FYX) qid=YX;
		case AVG_FZY:
			if(quantity==AVG_FZY) qid=ZY;
		case AVG_FZX:
			if(quantity==AVG_FZX) qid=ZX;
			
			// Volume weighted average is Sum (Vp strainp) / Sum Vp
			// where Vp = J mp/rho0
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					Matrix3 F = mpm[p]->GetDeformationGradientMatrix();
					value += Vp*F.get(qid,1.);
					Vtot += Vp;;
				}
			}
			if(Vtot>0.) value /= Vtot;
			value *= UnitsController::Scaling(100.);
			break;
			
		// totalenergies (Volume*energy) (J in Legacy)
		case KINE_ENERGY:
		case WORK_ENERGY:
		case STRAIN_ENERGY:
		case HEAT_ENERGY:
        case ENTROPY_ENERGY:
        case INTERNAL_ENERGY:
        case HELMHOLZ_ENERGY:
		{	bool threeD = fmobj->IsThreeD();
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	switch(quantity)
                    {   case KINE_ENERGY:
						{	value += 0.5*mpm[p]->mp*(mpm[p]->vel.x*mpm[p]->vel.x
																+ mpm[p]->vel.y*mpm[p]->vel.y);
							if(threeD)
								value += 0.5*mpm[p]->mp*(mpm[p]->vel.z*mpm[p]->vel.z);
                            break;
						}
						case WORK_ENERGY:
							value += mpm[p]->mp*mpm[p]->GetWorkEnergy();
							break;
						case STRAIN_ENERGY:
							value += mpm[p]->mp*mpm[p]->GetStrainEnergy();
							break;
						case HEAT_ENERGY:
							value += mpm[p]->mp*mpm[p]->GetHeatEnergy();
							break;
						case ENTROPY_ENERGY:
							value += mpm[p]->mp*mpm[p]->GetEntropy();
                            break;
						case INTERNAL_ENERGY:
							value += mpm[p]->mp*(mpm[p]->GetWorkEnergy()+mpm[p]->GetHeatEnergy());
                            break;
						case HELMHOLZ_ENERGY:
							value += mpm[p]->mp*(mpm[p]->GetWorkEnergy()+mpm[p]->GetHeatEnergy()
                                                 - mpm[p]->pPreviousTemperature*mpm[p]->GetEntropy());
                            break;
						default:
							break;
					}
				}
			}
			value *= UnitsController::Scaling(1.e-9);
			break;
		}
		
		// interface energy (J in Legacy)
		case INTERFACE_ENERGY:
			value = NodalPoint::interfaceEnergy*UnitsController::Scaling(1.e-9);
			break;
			
		// fricitonal work (J in Legacy)
		case FRICTION_WORK:
			value = NodalPoint::frictionWork*UnitsController::Scaling(1.e-9);
			break;
			
		// energies (Volume*energy) (J in Legacy)
		case PLAS_ENERGY:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
                {   value+=mpm[p]->mp*mpm[p]->GetPlastEnergy();
                }
			}
			value *= UnitsController::Scaling(1.e-9);
			break;
		
		// velocity x
		case AVG_VELX:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					value += Vp*mpm[p]->vel.x;
					Vtot += Vp;
				}
			}
			if(Vtot>0.) value /= Vtot;
 			break;
			
		// velocity y
		case AVG_VELY:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					value += Vp*mpm[p]->vel.y;
					Vtot += Vp;
				}
			}
			if(Vtot>0.) value /= Vtot;
			break;
			
		// velocity z
		case AVG_VELZ:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					value += Vp*mpm[p]->vel.z;
					Vtot += Vp;
				}
			}
			if(Vtot>0.) value /= Vtot;
			break;
			
		// x displacement
		case AVG_DISPX:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not verage those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					value += Vp*(mpm[p]->pos.x-mpm[p]->origpos.x);
					Vtot += Vp;
				}
			}
			if(Vtot>0.) value /= Vtot;
			break;
			
		// y displacement
		case AVG_DISPY:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					value += Vp*(mpm[p]->pos.y-mpm[p]->origpos.y);
					Vtot += Vp;
				}
			}
			if(Vtot>0.) value /= Vtot;
			break;
			
		// z displacement
		case AVG_DISPZ:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					value += Vp*(mpm[p]->pos.z-mpm[p]->origpos.z);
					Vtot += Vp;
				}
			}
			if(Vtot>0.) value /= Vtot;
			break;
			
		// temperature
		case AVG_TEMP:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					value += Vp*mpm[p]->pPreviousTemperature;
					Vtot += Vp;
				}
			}
			if(Vtot>0.) value /= Vtot;
			break;
		
		case WTFRACT_CONC:
		{   if(subcode>0)
            {   // other diffusion tasks
                DiffusionTask *otherDiffusion = DiffusionTask::FindDiffusionTaskByOrder(subcode);
                if(otherDiffusion!=NULL)
                {   // Volume weighted average is Sum (Vp pDiff->prevConc) / Sum Vp
                    // where Vp = J mp/rho0
                    for(p=p0;p<pend;p++)
                    {   if(mpm[p]->InReservoir()) continue;        // do not average those in resevoir
                        matid=mpm[p]->MatID();
                        if(IncludeThisMaterial(matid,singleParticle))
                        {   rho0 = theMaterials[matid]->GetRho(mpm[p]);
                            Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
                            Vp = Jp*mpm[p]->mp/rho0;
                            value += Vp*otherDiffusion->GetPrevParticleValue(mpm[p]);
                            Vtot += Vp;
                        }
                    }
                    if(Vtot>0.) value /= Vtot;
                }
            }
            else if(fmobj->HasDiffusion())
			{	// get total solvent content divided by total mass for total weight fraction
                double totalWeight=0.;
                for(p=p0;p<pend;p++)
                {	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
                    matid = mpm[p]->MatID();
                    if(IncludeThisMaterial(matid,singleParticle))
                    {	double csat = mpm[p]->GetConcSaturation();
                        value += diffusion->GetPrevParticleValue(mpm[p])*csat*mpm[p]->mp;
                        totalWeight += mpm[p]->mp;
                    }
                }
                if(totalWeight>0.) value /= totalWeight;
			}
#ifdef POROELASTICITY
			else if(fmobj->HasPoroelasticity())
			{	// Volume weighted average is Sum (Vp pp) / Sum Vp
				// where Vp = J mp/rho0
				for(p=p0;p<pend;p++)
				{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
					matid=mpm[p]->MatID();
					if(IncludeThisMaterial(matid,singleParticle))
					{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
						Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
						Vp = Jp*mpm[p]->mp/rho0;
						value += Vp*diffusion->GetPrevParticleValue(mpm[p]);
						Vtot += Vp;
					}
				}
				if(Vtot>0.) value /= Vtot;
				value *= UnitsController::Scaling(1.e-6);
			}
#endif
			break;
		}
        		
		case STEP_NUMBER:
			value=(double)fmobj->mstep;
			break;
		
		case RESERVOIR_SIZE:
			value=(double)mpmReservoir->currentCount();;
			break;
		
		// always in seconts
		case CPU_TIME:
			value=fmobj->CPUTime();
			break;

		// always in seconds
		case ELAPSED_TIME:
			value=fmobj->ElapsedTime();
			break;

		case GRID_ALPHA:
			value=bodyFrc.GetGridDamping(mtime);
			break;
		
		case PARTICLE_ALPHA:
			value=bodyFrc.GetParticleDamping(mtime);
			break;
            
		case HISTORY_VARIABLE:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					value += Vp*theMaterials[mpm[p]->MatID()]->GetHistory(subcode,mpm[p]->GetHistoryPtr(0));
					Vtot += Vp;
				}
			}
			if(Vtot>0.) value /= Vtot;
			break;
        
        case CRACK_LENGTH:
        case DEBONDED_CRACK_LENGTH:
        {   // get desired crack (or 0 if  no such crack)
            int cnum = 0;
            value = 0.;
            CrackHeader *nextCrack=firstCrack;
            while(nextCrack!=NULL)
            {   cnum++;
                if(subcode==cnum)
                {   // found crack, check the length
                    value = quantity==CRACK_LENGTH ? nextCrack->Length() :
                                                    nextCrack->DebondedLength() ;
                    break;
                }
                nextCrack = (CrackHeader *)nextCrack->GetNextObject();
            }
            break;
        }
		
		case TOT_FCONX:
		case TOT_FCONY:
		case TOT_FCONZ:
		{	// this vector will be filled
			Vector ftotal;
			ZeroVector(&ftotal);
			// First step has zero force
			bool hasForce = fmobj->mstep==0 ? true : false;
			
			if(!hasForce)
			{	// Three options
				//   1. VTK active, but not doing contact
				//   2. VTK inactive
				//   3. VTK active and archiving contact
				// When 1 and 2, all contact stuff here, which means must
				//     a. update lastArchiveContactStep, which is done in GetArchiveContactStepInterval
				//     b. clear force after reading
				//     c. Store last contact force in case another component is called next
				//     d. totalSteps will never be zero, so no need to track last contact force
				// When 3
				//     a. VTK archiving tracks lastArchiveContactStep
				//     b. Do not clear force (it is cleared on each VTK archive)
				//     c. VTK archiving will also save contact forces in case get here just after VTK archiving
				//          (since global archiving done after step is done
				
				// time steps since last cleared
				int totalSteps = archiver->GetArchiveContactStepInterval();
				
				// true if VTK archve inactive or it is not archiving contact on the grid
				bool clearForces = !archiver->GetDoingArchiveContact();
				
				// get array of vectors to store contact force for each material
				Vector *forces = archiver->GetLastContactForcePtr();
				
				// if needed sum forces on all nodes
				if(totalSteps>0)
				{	// non-zero steps, may or may not be doing VTKArchive
					for(int im=0;im<maxMaterialFields;im++) ZeroVector(&forces[im]);
					double scale = -1./(double)totalSteps;
					for(p=1;p<=nnodes;p++)
					{	nd[p]->AddGetContactForce(clearForces,forces,scale,NULL);
					}
				}
				
				// extract proper force (sum or one value)
				for(int im=0;im<maxMaterialFields;im++)
				{	if(whichMat==0)
					{	AddVector(&ftotal,&forces[im]);
					}
					else if(whichMat==MaterialBase::GetFieldMatID(im)+1)
					{	AddVector(&ftotal,&forces[im]);
						break;
				   }
				}
			}
				
 			// pick the component
			if(quantity==TOT_FCONX)
				value=ftotal.x;
			else if(quantity==TOT_FCONY)
				value=ftotal.y;
			else
				value=ftotal.z;
			break;
		}
		
		case TOT_REACTX:
		case TOT_REACTY:
		case TOT_REACTZ:
		{	// find force for BCs with provided ID (N in Legacy)
			Vector freaction = NodalVelBC::TotalReactionForce(whichMat);
			
			// pick the component
			if(quantity==TOT_REACTX)
				value = freaction.x;
			else if(quantity==TOT_REACTY)
				value = freaction.y;
			else
				value = freaction.z;
			value *= UnitsController::Scaling(1.e-6);
			break;
		}
		
		case TOT_REACTQ:
		{	// find heat flow for BCs with provided ID (J/sec in Legacy)
			double qreaction = NodalTempBC::TotalHeatReaction(whichMat);
			value = qreaction*UnitsController::Scaling(1.e-9);
			break;
		}
		
        // diffusion flux
        case TOT_FLUX:
        {   // get flow for BCs with provided ID and diffusion style
            // convert to diffusion task, and exit with zero if none
            DiffusionTask *qdiff = DiffusionTask::FindDiffusionTaskByOrder(subcode);
            if(qdiff==NULL) break;
            
            // get reaction flow for diffusion task
            value = NodalConcBC::TotalConcReaction(whichMat,qdiff);
            
            // scale concentration and poroelasticity
            if(qdiff==diffusion)
            {
                double csat = 1.;
#ifdef POROELASTICITY
                if(fmobj->HasPoroelasticity())
                {    // use for units in poroelasticity (convert to MPa in Legacy)
                    csat = UnitsController::Scaling(1.e-6);
                }
#endif
                value*=csat;
            }
            break;
        }
            
		// linear momentum (Legacy N-sec)
		case LINMOMX:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
					value += mpm[p]->mp*mpm[p]->vel.x;
			}
			value *= UnitsController::Scaling(1.e-6);
			break;
			
		case LINMOMY:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
					value += mpm[p]->mp*mpm[p]->vel.y;
			}
			value *= UnitsController::Scaling(1.e-6);
			break;
			
		case LINMOMZ:
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
					value += mpm[p]->mp*mpm[p]->vel.z;
			}
			value *= UnitsController::Scaling(1.e-6);
			break;
			
		// total angular momentum (Legacy J-sec)
		case ANGMOMX:
		case ANGMOMY:
		case ANGMOMZ:
		{	Vector Ltot = MakeVector(0.,0.,0.);
			Vector cp;
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	// get Mp Xp X Vp
					CrossProduct(&cp,&mpm[p]->pos,&mpm[p]->vel);
					AddScaledVector(&Ltot,&cp,mpm[p]->mp);
                    
					// add particle spin angular momentum
					Vector Lp = mpm[p]->GetParticleAngMomentum();
					AddVector(&Ltot,&Lp);
				}
			}
			if(quantity==ANGMOMX)
				value = Ltot.x;
			else if(quantity==ANGMOMY)
				value = Ltot.y;
			else
				value = Ltot.z;
			value *= UnitsController::Scaling(1.e-9);
			break;
		}

		// particle spin angular momentum (Legacy J-sec)
		case LPMOMX:
		case LPMOMY:
		case LPMOMZ:
		{	Vector Ltot = MakeVector(0.,0.,0.);
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	Vector Lp = mpm[p]->GetParticleAngMomentum();
					AddVector(&Ltot,&Lp);
				}
			}
			if(quantity==LPMOMX)
				value = Ltot.x;
			else if(quantity==LPMOMY)
				value = Ltot.y;
			else
				value = Ltot.z;
			value *= UnitsController::Scaling(1.e-9);
			break;
		}
			
		// average particle angular velocity (radians per sec)
		case ANGVELX:
		case ANGVELY:
		case ANGVELZ:
		{	// sum over all particles
			Vector wptot = MakeVector(0.,0.,0.);
			for(p=p0;p<pend;p++)
			{	if(mpm[p]->InReservoir()) continue;		// do not average those in resevoir
				matid = mpm[p]->MatID();
				if(IncludeThisMaterial(matid,singleParticle))
				{	// get volume
					rho0 = theMaterials[matid]->GetRho(mpm[p]);
					Jp = theMaterials[matid]->GetCurrentRelativeVolume(mpm[p],0);
					Vp = Jp*mpm[p]->mp/rho0;
					Vtot += Vp;
					
					// angular spatial velocity gradient
					Matrix3 spatialGradVp = mpm[p]->GetParticleGradVp(true,false);
						
                    // Extract angular velocity antisymmetric grad Vp
					if(fmobj->IsThreeD())
					{	Vector wp = MakeVector(0.5*(spatialGradVp(2,1)-spatialGradVp(1,2)),
											   0.5*(spatialGradVp(0,2)-spatialGradVp(2,0)),
											   0.5*(spatialGradVp(1,0)-spatialGradVp(0,1)));
						AddScaledVector(&wptot,&wp,Vp);
					}
					else
					{	double wpz = 0.5*(spatialGradVp(1,0)-spatialGradVp(0,1));
						wptot.z += Vp*wpz;
					}
				}
			}
			if(quantity==ANGVELX)
				value = wptot.x/Vtot;
			else if(quantity==ANGVELY)
				value = wptot.y/Vtot;
			else
				value = wptot.z/Vtot;
			break;
		}

		// grid kinetic energy (J in Legacy)
		case GRID_KINE_ENERGY:
		{	double totalMass;
			for(p=1;p<=nnodes;p++)
				nd[p]->AddKineticEnergyAndMass(value,totalMass);
			value *= UnitsController::Scaling(1.e-9);
			break;
		}

		// skip decoheion and unknown
		case DECOHESION:
		case UNKNOWN_QUANTITY:
			return nextGlobal;
		
		// zero if not programmed yet
		default:
			break;
	}
	
	toArchive.push_back(value);
	
	// return next one
	return nextGlobal;
}

// append tab and color string
GlobalQuantity *GlobalQuantity::AppendColor(char *fline)
{
	if(quantity==UNKNOWN_QUANTITY || quantity==DECOHESION)
		return nextGlobal;
	
	strcat(fline,"\t");
	strcat(fline,PickColor(colorID));
	
	// return next one
	return nextGlobal;
}

// Class method to pick 10 colors using buNum mod 10
const char *GlobalQuantity::PickColor(int byNum)
{
	int numID = byNum % 10;
	switch(numID)
	{   case 0:
			return "black";
		case 1:
			return "blue";
		case 2:
			return "red";
		case 3:
			return "green";
		case 4:
			return "brown";
		case 5:
			return "cyan";
		case 6:
			return "magenta";
		case 7:
			return "orange";
		case 8:
			return "purple";
		case 9:
			return "yellow";
		default:
			return "black";
	}
	return "black";
}

#pragma mark GlobalQuantity::ACCESSORS

// decide if archiving this material
// return true is matches material (rigid or nonrigid) or is for a single partle.
// But ff whichMat is zero only return true if it is a non-rigid material.
//      Thus one material can average only any material, but average over all
//		particles is non rigid only
bool GlobalQuantity::IncludeThisMaterial(int matid,bool singleParticle)
{
	// accept any specified material
	if(matid+1==whichMat || singleParticle) return (bool)true;
	
	// otherwise only allow non-rigid materials
	if(whichMat==0 && !theMaterials[matid]->IsRigid()) return (bool)true;
	
	return (bool)false;
}

// set the next Global
GlobalQuantity *GlobalQuantity::GetNextGlobal(void) { return nextGlobal; }
void GlobalQuantity::SetNextGlobal(GlobalQuantity *newGlobal) { nextGlobal=newGlobal; }

// compare to settings
bool GlobalQuantity::IsSameQuantity(int qval,int qcode,int qmat)
{	if(quantity==qval && subcode==qcode && qmat==whichMat) return true;
	return false;
}

// return the global quantity by integer ID
int GlobalQuantity::GetQuantity(void) { return quantity; }

// true if is tracer particle
bool GlobalQuantity::IsTracerParticle(void) { return ptPos!=NULL; }

// Set name when known
GlobalQuantity *GlobalQuantity::SetTracerParticle(void)
{	if(ptPos!=NULL)
	{	char pnum[25];
		size_t psize=25;
		snprintf(pnum,psize,"%d",ptNum+1);
		strcat(name,pnum);
	}
	return nextGlobal;
}


