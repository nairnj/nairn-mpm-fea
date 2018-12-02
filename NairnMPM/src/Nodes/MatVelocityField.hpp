/********************************************************************************
    MatVelocityField.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 3 April 2009.
    Copyright (c) 2009 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _MATVELOCITYFIELD_

#define _MATVELOCITYFIELD_

#define MAX_FIELDS_FOR_CRACKS 4
// Would like to change this to 2, but better make sure code never tries to
// read fields [2] and [3]. Find SCWarning for where there are possible problems
#define MAX_FIELDS_FOR_ONE_CRACK 4

#define RIGID_FIELD_BIT 1
#define IGORE_CRACKS_BIT 2
#define RIGID_BLOCK_BIT 4

#define VSTAR_VEC 0
#define VSTARPREV_VEC 1
#define VSTARNEXT_VEC 2

class NodalPoint;

class MatVelocityField
{
	public:
		// variables (changed in MPM time step)
		int numberPoints;			// number of material points in this field
		double mass;				// total mass of this field
		Vector pk;					// momentum
		Vector disp;				// displacement for contact calculations
		Vector *volumeGrad;			// volume gradient allocated in multimaterial mode
		Vector *xpic;				// For XPIC (OSParticulas) or copy momenta (nfm)
		static int pkCopy;			// which xpic vector to store pk
	
		// constants (not changed in MPM time step)
				
        // constructors and destructors
        MatVelocityField(int);
		~MatVelocityField();
		void Zero(void);
		
		// methods
		void AddMomentumTask1(Vector *,Vector *,int);
		void CopyMassAndMomentum(NodalPoint *,int,int);
        void CopyMassAndMomentumLast(NodalPoint *,int,int);
		void CopyGridForces(NodalPoint *,int,int);
		void RestoreMomenta(void);
		void ChangeMatMomentum(Vector *,int,double);
		void AddContactForce(Vector *);
		void CalcVelocityForStrainUpdate(void);
		void AddGravityAndBodyForceTask3(Vector *);
        void AddFtot(Vector *);
        void AddFtotScaled(Vector *,double);
        void UpdateMomentum(double);
		virtual void IncrementNodalVelAcc(double,GridToParticleExtrap *gp) const;
	
        // only called if ADJUST_EXTRAPOLATED_PK_FOR_SYMMETRY is defined
        void AdjustForSymmetryBC(int);
	
		// accessors
		void Describe(int) const;
		void AddContactVolume(double);
		void SetContactVolume(double);
		double GetContactVolume(void) const;
        Vector GetVelocity(void);
        void SetMomentVelocityDirection(Vector *,int);
        void AddMomentVelocityDirection(Vector *,double,int);
        void SetFtotDirection(Vector *,double,Vector *);
        void AddFtotDirection(Vector *,double,double,Vector *);
        Vector GetFtot(void) const;
        Vector *GetFtotPtr(void);
		bool IsRigidField(void) const;
		void SetRigidField(bool);
		bool IsFixedRigidField(void) const;
		bool IsRigidBlockField(void) const;
		int GetFlags(void) const;
		bool IgnoresCracks(void) const;
	
		// class methods
		static bool ActiveField(MatVelocityField *);
		static bool ActiveRigidField(MatVelocityField *mvf);
		static bool ActiveNonrigidField(MatVelocityField *mvf);
		static bool ActiveNonrigidSourceField(MatVelocityField *,int);
		static bool ActiveNonrigidSeesCracksField(MatVelocityField *,bool);
	
	protected:
		int flags;					// bitwise flags for some field properties
		double volume;				// only for contact (cracks or multimaterial)

        Vector vk;					// velocity
        Vector ftot;				// total force or contact force for rigid material
};

#endif
