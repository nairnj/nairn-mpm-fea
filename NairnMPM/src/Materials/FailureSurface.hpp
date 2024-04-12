/********************************************************************************
	FailureSurface.hpp
	nairn-mpm-fea

	Created by John Nairn, June 26, 2015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
 
	Base class for all failure surfaces used in softening materials

	Dependencies
		none
********************************************************************************/

#ifndef _FailureSurface_
#define _FailureSurface_

#define MAXPRINCIPALSTRESS 1

// Angle to pick tension or shear
#define TANPIOVER8 0.414213562373095

// TENSILE_FAILURE is max principle stress
// SHEAR_FAILURE is max shear stress
// combined are those in areas between two failures
// (note these must be 51 to 149)
#define NO_FAILURE 0
#define TENSILE_TENSILE_FAILURE 80
#define TENSILE_FAILURE 90
#define TENSILE_SHEAR_FAILURE 100
#define SHEAR_FAILURE 110
#define SHEAR_SHEAR_FAILURE 120

// pressure dependent models
enum { P_INDEPENDENT=0,P_LINEAR,P_STEPLINEAR,P_STEPEXP,P_SIGMOIDALEXP };

class MaterialBase;

class FailureSurface
{
	public:
		// constructors and destructors
		FailureSurface();
		FailureSurface(MaterialBase *);
 	
		// initialize
		virtual char *InputInitationProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
        virtual void RemapProperties(int);
		virtual void PrintInitiationProperties(void) const;
	
		// methods
		virtual int ShouldInitiateFailure(Tensor *,Vector *,int,double &,GenADaMVariables *) const;
		virtual bool CloserToShear(Vector *,double,double,double) const;

        // pressure dependent methods
        virtual void InitPDTerms(void);
        virtual char *InputPDProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadPDTerms(double);
        virtual void PrintPDTerms(void) const;
        virtual double GetPDStrength(double,double,double) const;
        virtual double sigmaIStabilityScale(void) const;
        virtual double sigmaIIStabilityScale(void) const;
        virtual bool IsPDModelConstant(double,double,double &) const;
        virtual bool IsPressureDependent(void) const;

		// accessors
		virtual double sigmaI(void) const;
		virtual double sigmaII(void) const;
		virtual double sigmaII(double,double &,Tensor &,GenADaMVariables *,double,double,double,int) const;
		virtual const char *GetInitiationLawName(void) const;
		virtual void SetFailureSurface(int);

	protected:
		double criticalNormal,sigmacRed;
		double criticalShear,taucRed;
		MaterialBase *parent;
		int failureSurface;
        int pdModel;

        // pressure dependence properties
        double tauh,tauhRed;
        double parg1,parg2;
        double shearShift;
        double maxPressure;
        double P1Red,P2Red;
};

#endif
