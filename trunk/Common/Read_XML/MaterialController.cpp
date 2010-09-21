/********************************************************************************
    MaterialController.cpp
    NairnFEA
    
    Created by John Nairn on 6/27/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_XML/MaterialController.hpp"

// Material types created here
#include "Materials/IsotropicMat.hpp"
#include "Materials/TransIsotropic.hpp"
#include "Materials/Orthotropic.hpp"
#ifdef MPM_CODE
	#include "NairnMPM_Class/NairnMPM.hpp"
	#include "Materials/Viscoelastic.hpp"
	#include "Materials/VonMisesHardening.hpp"
	#include "Materials/HillPlastic.hpp"
	#include "Materials/WoodMaterial.hpp"
	#include "Materials/JohnsonCook.hpp"
	#include "Materials/MGSCGLMaterial.hpp"
	#include "Materials/SLMaterial.hpp"
	#include "Materials/Mooney.hpp"
	#include "Materials/BistableIsotropic.hpp"
	#include "Materials/RigidMaterial.hpp"
	#include "Materials/CohesiveZone.hpp"
	#include "Materials/LinearTraction.hpp"
	#include "Materials/CubicTraction.hpp"
	#include "Materials/TrilinearTraction.hpp"
#else
	#include "Materials/ImperfectInterface.hpp"
#endif

MaterialController *matCtrl=NULL;

/********************************************************************************
	MaterialController: methods
********************************************************************************/

/* Create new material
	When a new material is added to NairnMPM or NairnFEA, you must include the
	header above and add a case below to create a new material when the matID
	calls for the new type. All other issues on implementing the material
	should be possible in the new material class alone, without modifying
	an other core code.
*/
int MaterialController::AddMaterial(int matID,char *matName)
{
	MaterialBase *newMaterial;
	
	switch(matID)
	{   case ISOTROPIC:
			newMaterial=new IsotropicMat(matName);
			break;
		case TRANSISO1:
		case TRANSISO2:
			newMaterial=new TransIsotropic(matName,matID);
			break;
		case ORTHO:
			newMaterial=new Orthotropic(matName);
			break;
#ifdef MPM_CODE
		case VISCOELASTIC:
			newMaterial=new Viscoelastic(matName);
			break;
		case MOONEYRIVLIN:
			newMaterial=new Mooney(matName);
			break;
		case VONMISESHARDENING:
			newMaterial=new VonMisesHardening(matName);
			break;
		case BISTABLEISO:
			newMaterial=new BistableIsotropic(matName);
			break;
		case RIGIDMATERIAL:
			newMaterial=new RigidMaterial(matName);
			break;
		case COHESIVEZONEMATERIAL:
			newMaterial=new CohesiveZone(matName);
			break;
		case LINEARTRACTIONMATERIAL:
			newMaterial=new LinearTraction(matName);
			break;
		case CUBICTRACTIONMATERIAL:
			newMaterial=new CubicTraction(matName);
			break;
		case HILLPLASTIC:
			newMaterial=new HillPlastic(matName);
			break;
		case WOODMATERIAL:
			newMaterial=new WoodMaterial(matName);
			break;
		case JOHNSONCOOK:
			newMaterial=new JohnsonCook(matName);
			break;
		case MGSCGLMATERIAL:
			newMaterial=new MGSCGLMaterial(matName);
			break;
		case SLMATERIAL:
			newMaterial=new SLMaterial(matName);
			break;
		case TRILINEARTRACTIONMATERIAL:
			newMaterial=new TrilinearTraction(matName);
			break;
#else
		case INTERFACEPARAMS:
			newMaterial=new ImperfectInterface(matName);
			break;
#endif
		default:
			return FALSE;
	}
	AddObject(newMaterial);
	return TRUE;
}

// assemble into array used in the code
int MaterialController::SetMaterialArray(void)
{
	theMaterials=(MaterialBase **)MakeObjectArray(0);
	if(theMaterials==NULL) return FALSE;
	
	// fill the array
	MaterialBase *obj=(MaterialBase *)firstObject;
	nmat=0;
	while(obj!=NULL)
	{	theMaterials[nmat]=obj;
		nmat++;
		obj=(MaterialBase *)obj->GetNextObject(); 
			
	}
	return TRUE;
}

// pointer to read a material property
char *MaterialController::InputPointer(char *property,int &input)
{	return ((MaterialBase *)lastObject)->InputMat(property,input);
}

// set material color (optional)
void MaterialController::SetMatColor(float red,float green,float blue)
{	((MaterialBase *)lastObject)->red=red;
	((MaterialBase *)lastObject)->green=green;
	((MaterialBase *)lastObject)->blue=blue;
}

#ifdef MPM_CODE
// return current material
void MaterialController::SetCriterion(int matCriterion,int setIndex)
{	((MaterialBase *)lastObject)->criterion[setIndex]=matCriterion;
}

// set propagation direction
void MaterialController::SetDirection(int matDirection,int setIndex)
{	((MaterialBase *)lastObject)->matPropagateDirection[setIndex]=matDirection;
}

// set propagation traction material
void MaterialController::SetTractionMat(int mat,int setIndex)
{	((MaterialBase *)lastObject)->tractionMat[setIndex]=mat;
}

// when done with Friction command, create new material friction object
void MaterialController::SetMaterialFriction(void)
{	((MaterialBase *)lastObject)->SetFriction(friction,otherMatID,Dn,Dnc,Dt);
}

#endif

