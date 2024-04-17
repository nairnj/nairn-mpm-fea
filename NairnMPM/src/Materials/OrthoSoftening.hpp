/********************************************************************************
	OrthoSoftening.hpp
	nairn-mpm-fea
 
	Created by John Nairn on 30 JAN 2018
	Copyright (c) 2018 John A. Nairn, All rights reserved.
 
	Dependencies
		TransIsoSoftening.hpp TransIsotropic.hpp Elastic.hpp MaterialBase.hpp
********************************************************************************/

#ifndef ORTHOSOFTENING

#define ORTHOSOFTENING 54

#include "Materials/TransIsoSoftening.hpp"

class OrthoSoftening : public TransIsoSoftening
{
	public:
	
		// constructors and destructors
		OrthoSoftening(char *,int);
	
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
		virtual bool AcceptInitiationLaw(FailureSurface *,int);
		virtual bool AcceptSofteningLaw(SofteningLaw *,int);
		virtual const char *VerifyAndLoadProperties(int);
		virtual const char *VerifyAndLoadOrthoProperties(int np);
		virtual void PrintMechanicalProperties(void) const;
		virtual void PrintTransportProperties(void) const;
	
		// history data
		virtual char *InitHistoryData(char *,MPMBase *);
		virtual void ResetHistoryData(char *,MPMBase *);

		// methods
		virtual int DecodeDamageInitiation(int,Vector *,int,double *) const;
		virtual void LoadCrackAxisProperties(int,CrackAxisProperties *,int ,ElasticProperties *) const;
		virtual bool GetRToCrack(Matrix3 *,double *, bool, int) const;

		// accessors
		virtual const char *MaterialType() const;
		virtual double WaveSpeed(bool,MPMBase *) const;
		virtual double GetDiffZ(void) const;
		virtual double GetKcondZ(void) const;
		virtual SofteningLaw *softeningXX(void) const;
		virtual SofteningLaw *softeningYY(void) const;
		virtual SofteningLaw *softeningXYX(void) const;
		virtual SofteningLaw *softeningXYY(void) const;
		virtual SofteningLaw *softeningYZY(void) const;
	
	protected:
		double Ex,Ey,Ez,Gxy,Gyz,Gxz,ax,ay,az,betax,betay,betaz;
		double nuxy,nuyx,nuxz,nuzx,nuzy,nuyz;
		double Dz,kCondz;
		SofteningLaw *softeningZZ;
		SofteningLaw *softeningXZX;
		SofteningLaw *softeningXZZ;
		SofteningLaw *softeningYZZ;
};

#endif
