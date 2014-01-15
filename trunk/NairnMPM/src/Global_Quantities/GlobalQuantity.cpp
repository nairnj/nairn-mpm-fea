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

#include "Global_Quantities/GlobalQuantity.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "System/ArchiveData.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"

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

GlobalQuantity::GlobalQuantity(char *quant,int whichOne)
{
	char nameStr[200];
	whichMat=whichOne;
	
	// set quantity and subcode
	subcode=0;
	if(strcmp(quant,"sxx")==0 || strcmp(quant,"sRR")==0)
		quantity=AVG_SXX;
	else if(strcmp(quant,"syy")==0 || strcmp(quant,"sZZ")==0)
		quantity=AVG_SYY;
	else if(strcmp(quant,"sxy")==0 || strcmp(quant,"sRZ")==0)
		quantity=AVG_SXY;
	else if(strcmp(quant,"szz")==0 || strcmp(quant,"sTT")==0)
		quantity=AVG_SZZ;
	else if(strcmp(quant,"sxz")==0)
		quantity=AVG_SXZ;
	else if(strcmp(quant,"syz")==0)
		quantity=AVG_SYZ;
	else if(strcmp(quant,"exx")==0 || strcmp(quant,"eRR")==0)
		quantity=AVG_EXX;
	else if(strcmp(quant,"eyy")==0 || strcmp(quant,"eZZ")==0)
		quantity=AVG_EYY;
	else if(strcmp(quant,"exy")==0 || strcmp(quant,"eRZ")==0)
		quantity=AVG_EXY;
	else if(strcmp(quant,"ezz")==0 || strcmp(quant,"eTT")==0)
		quantity=AVG_EZZ;
	else if(strcmp(quant,"exz")==0)
		quantity=AVG_EXZ;
	else if(strcmp(quant,"eyz")==0)
		quantity=AVG_EYZ;
	else if(strcmp(quant,"exxe")==0 || strcmp(quant,"eRRe")==0)
		quantity=AVG_EXXE;
	else if(strcmp(quant,"eyye")==0 || strcmp(quant,"eZZe")==0)
		quantity=AVG_EYYE;
	else if(strcmp(quant,"exye")==0 || strcmp(quant,"eRZe")==0)
		quantity=AVG_EXYE;
	else if(strcmp(quant,"ezze")==0 || strcmp(quant,"eTTe")==0)
		quantity=AVG_EZZE;
	else if(strcmp(quant,"exze")==0)
		quantity=AVG_EXZE;
	else if(strcmp(quant,"eyze")==0)
		quantity=AVG_EYZE;
	else if(strcmp(quant,"exxp")==0 || strcmp(quant,"eRRp")==0)
		quantity=AVG_EXXP;
	else if(strcmp(quant,"eyyp")==0 || strcmp(quant,"eZZp")==0)
		quantity=AVG_EYYP;
	else if(strcmp(quant,"exyp")==0 || strcmp(quant,"eRZp")==0)
		quantity=AVG_EXYP;
	else if(strcmp(quant,"ezzp")==0 || strcmp(quant,"eTTp")==0)
		quantity=AVG_EZZP;
	else if(strcmp(quant,"exzp")==0)
		quantity=AVG_EXZP;
	else if(strcmp(quant,"eyzp")==0)
		quantity=AVG_EYZP;
	else if(strcmp(quant,"Kinetic Energy")==0)
		quantity=KINE_ENERGY;
	else if(strcmp(quant,"Grid Kinetic Energy")==0)
		quantity=GRID_KINE_ENERGY;
	else if(strcmp(quant,"Strain Energy")==0)
		quantity=STRAIN_ENERGY;
	else if(strcmp(quant,"Heat Energy")==0 || strcmp(quant,"Thermal Energy")==0)
		quantity=HEAT_ENERGY;
	else if(strcmp(quant,"Entropy")==0)
		quantity=ENTROPY_ENERGY;
	else if(strcmp(quant,"Internal Energy")==0)
		quantity=INTERNAL_ENERGY;
	else if(strcmp(quant,"Helmholz Energy")==0)
		quantity=HELMHOLZ_ENERGY;
	else if(strcmp(quant,"Interface Energy")==0)
		quantity=INTERFACE_ENERGY;
	else if(strcmp(quant,"Work Energy")==0)
		quantity=WORK_ENERGY;
	else if(strcmp(quant,"Plastic Energy")==0)
		quantity=PLAS_ENERGY;
	else if(strcmp(quant,"velx")==0 || strcmp(quant,"velR")==0)
		quantity=AVG_VELX;
	else if(strcmp(quant,"vely")==0 || strcmp(quant,"velZ")==0)
		quantity=AVG_VELY;
	else if(strcmp(quant,"velz")==0)
		quantity=AVG_VELZ;
	else if(strcmp(quant,"dispx")==0 || strcmp(quant,"dispR")==0)
		quantity=AVG_DISPX;
	else if(strcmp(quant,"dispy")==0 || strcmp(quant,"dispZ")==0)
		quantity=AVG_DISPY;
	else if(strcmp(quant,"dispz")==0)
		quantity=AVG_DISPZ;
	else if(strcmp(quant,"temp")==0)
		quantity=AVG_TEMP;
	else if(strcmp(quant,"concentration")==0)
		quantity=WTFRACT_CONC;
	else if(strcmp(quant,"Step number")==0)
		quantity=STEP_NUMBER;
	else if(strcmp(quant,"CPU time")==0)
		quantity=CPU_TIME;
	else if(strcmp(quant,"Elapsed time")==0)
		quantity=ELAPSED_TIME;
	else if(strcmp(quant,"alpha")==0)
		quantity=FEEDBACK_ALPHA;
	else if(strcmp(quant,"contactx")==0 || strcmp(quant,"contactR")==0)
		quantity=TOT_FCONX;
	else if(strcmp(quant,"contacty")==0 || strcmp(quant,"contactZ")==0)
		quantity=TOT_FCONY;
	else if(strcmp(quant,"contactz")==0)
		quantity=TOT_FCONZ;
	else if(strcmp(quant,"reactionx")==0 || strcmp(quant,"reactionR")==0)
		quantity=TOT_REACTX;
	else if(strcmp(quant,"reactiony")==0 || strcmp(quant,"reactionZ")==0)
		quantity=TOT_REACTY;
	else if(strcmp(quant,"reactionz")==0)
		quantity=TOT_REACTZ;
	else
	{	quantity=UNKNOWN_QUANTITY;
	
		// possible a history variable which must be "history n"
		if(strlen(quant)>7)
		{	strcpy(nameStr,quant);
			nameStr[7]=0;
			if(strcmp(nameStr,"history")==0)
			{	sscanf(quant,"%*s %d",&subcode);
				quantity=HISTORY_VARIABLE;
			}
		}
	}
	
	// set name
	if(whichMat!=0)
		sprintf(nameStr,"%s mat %d",quant,whichMat);
	else
		strcpy(nameStr,quant);
	name=new char[strlen(nameStr)+1];
	strcpy(name,nameStr);
	
	// set color ID
	colorID=numGlobal % 10;
	
	// this object is current the last one
	SetNextGlobal(NULL);
	
	// adjust previous global quantity or set firstGlobal if this is the first one
	if(lastGlobal!=NULL)
		lastGlobal->SetNextGlobal(this);
	else
		firstGlobal=this;
	lastGlobal=this;
	
	// count the number of objects
	numGlobal++;
}

/*******************************************************************
	GlobalQuantity: Methods
*******************************************************************/

// append quantity
GlobalQuantity *GlobalQuantity::AppendName(char *fline)
{
	char nameStr[100];
	
	if(quantity==UNKNOWN_QUANTITY)
		return nextGlobal;
	
	sprintf(nameStr,"\t%c%s%c",quote,name,quote);
	strcat(fline,nameStr);
	return nextGlobal;
}

// append quantity
GlobalQuantity *GlobalQuantity::AppendQuantity(char *fline)
{
	int p,numAvged=0;
	double value=0.,rho,rho0;
	char numStr[100];
	int matid,qid=0;
    bool threeD;
	
	switch(quantity)
	{   // stresses in MPa
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
            threeD = fmobj->IsThreeD();
			for(p=0;p<nmpms;p++)
			{   matid=mpm[p]->MatID();
				if(IncludeThisMaterial(matid))
				{	rho0=theMaterials[matid]->rho;
                    rho = rho0/theMaterials[matid]->GetCurrentRelativeVolume(mpm[p]);
                    Tensor sp = mpm[p]->ReadStressTensor();
 					value+=rho*Tensor_i(&sp,qid);
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			value*=1.e-6;
			break;
		
		// elastic strain in %
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
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	Tensor *ep=mpm[p]->GetStrainTensor();
				    value+=Tensor_i(ep,qid);
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			value*=100.;
			break;

		// plastic strain in %
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
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	Tensor *eplast=mpm[p]->GetPlasticStrainTensor();
					value+=Tensor_i(eplast,qid);
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			value*=100.;
			break;

		// total strain
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
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	Tensor *ep=mpm[p]->GetStrainTensor();
                    if(theMaterials[mpm[p]->MatID()]->PartitionsElasticAndPlasticStrain())
					{   Tensor *eplast=mpm[p]->GetPlasticStrainTensor();
                        value += Tensor_i(ep,qid)+Tensor_i(eplast,qid);
                    }
                    else
                        value += Tensor_i(ep,qid);
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			value*=100.;
			break;
		
		// energies (Volume*energy) in J
		case KINE_ENERGY:
		case WORK_ENERGY:
		case STRAIN_ENERGY:
		case HEAT_ENERGY:
        case ENTROPY_ENERGY:
        case INTERNAL_ENERGY:
        case HELMHOLZ_ENERGY:
            threeD = fmobj->IsThreeD();
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
                {   switch(quantity)
                    {   case KINE_ENERGY:
							value += 0.5e-3*mpm[p]->mp*(mpm[p]->vel.x*mpm[p]->vel.x
																+ mpm[p]->vel.y*mpm[p]->vel.y);
							if(threeD)
								value += 0.5e-3*mpm[p]->mp*(mpm[p]->vel.z*mpm[p]->vel.z);
                            break;
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
			value*=1.e-6;
			break;
		
		// interface energy in J
		case INTERFACE_ENERGY:
			value=1.e-9*NodalPoint::interfaceEnergy;
			break;
			
		// energies (Volume*energy) in J
		case PLAS_ENERGY:
            threeD = fmobj->IsThreeD();
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
                {   value+=mpm[p]->mp*mpm[p]->GetPlastEnergy();
                }
			}
			value*=1.e-6;
			break;
		
		// velocity x in mm/sec
		case AVG_VELX:
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	value+=mpm[p]->vel.x;
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			break;
			
		// velocity y in mm/sec
		case AVG_VELY:
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	value+=mpm[p]->vel.y;
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			break;
			
		// velocity z in mm/sec
		case AVG_VELZ:
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	value+=mpm[p]->vel.z;
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			break;
			
		// x displacement
		case AVG_DISPX:
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	value+=(mpm[p]->pos.x-mpm[p]->origpos.x);
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			break;
			
		// y displacement
		case AVG_DISPY:
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	value+=(mpm[p]->pos.y-mpm[p]->origpos.y);
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			break;
			
		// z displacement
		case AVG_DISPZ:
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	value+=(mpm[p]->pos.z-mpm[p]->origpos.z);
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			break;
			
		// temperature
		case AVG_TEMP:
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	value+=mpm[p]->pPreviousTemperature;
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			break;
		
		case WTFRACT_CONC:
		{	double totalWeight=0.;
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	double csat=theMaterials[mpm[p]->MatID()]->concSaturation;
					value+=mpm[p]->pPreviousConcentration*csat*mpm[p]->mp;
					totalWeight+=mpm[p]->mp;
					numAvged++;
				}
			}
			if(numAvged>0) value/=totalWeight;
			break;
		}
		
		case STEP_NUMBER:
			value=(double)fmobj->mstep;
			break;
			
		case CPU_TIME:
			value=fmobj->CPUTime();
			break;

		case ELAPSED_TIME:
			value=fmobj->ElapsedTime();
			break;

		case FEEDBACK_ALPHA:
			value=bodyFrc.GetAlpha();
			break;
		
		case HISTORY_VARIABLE:
			for(p=0;p<nmpms;p++)
			{	if(IncludeThisMaterial(mpm[p]->MatID()))
				{	value+=theMaterials[mpm[p]->MatID()]->GetHistory(subcode,mpm[p]->GetHistoryPtr());
					numAvged++;
				}
			}
			if(numAvged>0) value/=(double)numAvged;
			break;
		
		case TOT_FCONX:
		case TOT_FCONY:
		case TOT_FCONZ:
        {   // Three options
            //   1. VTK active, but not doing contact
            //   2. VTK inactive
            //   3. VTK active and archiving contact
            // When 1 and 2, all contact stuff here, which means must
            //     a. update lastArchiveContactStep, which is done in GetArchiveContactStepInterval
            //     b. clear force after reading
            //     c. totalSteps will never be zero, so no need to track last contact force
            // When 3
            //     a. VTK archiving tracks lastArchiveContactStep
            //     b. Do clear force (it is cleared on each VTK archive)
            //     c. VTK archiving will also set lastContactForce in case get here just after VTK archiving
            //          (since global archiving done after step is done
			int totalSteps=archiver->GetArchiveContactStepInterval();
			bool clearForces=!archiver->GetDoingArchiveContact();
			Vector ftotal;
			ZeroVector(&ftotal);
			
			if(totalSteps>0)
			{	// non-zero steps, may or may not be doing VTKArchive
				for(p=1;p<=nnodes;p++)
				{	Vector fcontact=nd[p]->GetTotalContactForce(clearForces);
					AddVector(&ftotal,&fcontact);
				}
				ScaleVector(&ftotal,-1./(double)totalSteps);           // force of rigid particles on the object (per step)
			}
			else
			{	// VTK task just found contact force, so use it here
				ftotal=archiver->GetLastContactForce();
			}
			// if totalSteps==0 and clearForces, then must be zero, as initialized above
				
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
		{	// find force for BCs with provided ID
			Vector freaction = NodalVelBC::TotalReactionForce(whichMat);
			
			// pick the component
			if(quantity==TOT_REACTX)
				value = 1.e-6*freaction.x;
			else if(quantity==TOT_REACTY)
				value = 1.e-6*freaction.y;
			else
				value = 1.e-6*freaction.z;
			break;
		}
		
		// grid kinetic energy (J)
		case GRID_KINE_ENERGY:
		{	double totalMass;
			for(p=1;p<=nnodes;p++)
				nd[p]->AddKineticEnergyAndMass(value,totalMass);
			value *= 1.e-9;
			break;
		}

		// skip unknown
		case UNKNOWN_QUANTITY:
			return nextGlobal;
		
		// zero if not programmed yet
		default:
			break;
	}
	
	sprintf(numStr,"\t%e",value);
	strcat(fline,numStr);
	
	// return next one
	return nextGlobal;
}

// append tab and color string
GlobalQuantity *GlobalQuantity::AppendColor(char *fline)
{
	if(quantity==UNKNOWN_QUANTITY) return nextGlobal;
	
	strcat(fline,"\t");
	switch(colorID)
	{   case 0:
			strcat(fline,"black");
			break;
		case 1:
			strcat(fline,"blue");
			break;
		case 2:
			strcat(fline,"red");
			break;
		case 3:
			strcat(fline,"green");
			break;
		case 4:
			strcat(fline,"brown");
			break;
		case 5:
			strcat(fline,"cyan");
			break;
		case 6:
			strcat(fline,"magenta");
			break;
		case 7:
			strcat(fline,"orange");
			break;
		case 8:
			strcat(fline,"purple");
			break;
		case 9:
			strcat(fline,"yellow");
			break;
		default:
			strcat(fline,"black");
			break;
	}
	
	// return next one
	return nextGlobal;
}

// decide if archiving this material
bool GlobalQuantity::IncludeThisMaterial(int matid)
{
	// accept any specified material
	if(matid+1==whichMat) return (bool)TRUE;
	
	// otherwise only allow non-rigid materials
	if(whichMat==0 && !theMaterials[matid]->Rigid()) return (bool)TRUE;
	
	return (bool)FALSE;
}

// set the next Global
void GlobalQuantity::SetNextGlobal(GlobalQuantity *newGlobal) { nextGlobal=newGlobal; }

