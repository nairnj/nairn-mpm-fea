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
#define DELTA_VSTARPREV_VEC 3
#define DELTA_VSTARNEXT_VEC 4
#define DELTA_VSTORE_VEC 5

// calculations options
enum { INITIALIZE_XPIC=0,UPDATE_VSTAR,COPY_VSTARNEXT };

class NodalPoint;

// contact extrapolations
typedef struct {
	double cvolume;				// volume for contact
	Vector *terms;				// array of contact vectors to extrapolate
} ContactTerms;

class MatVelocityField
{
	public:
		// variables (changed in MPM time step)
		int numberPoints;			// number of material points in this field
		double mass;				// total mass of this field
		Vector pk;					// momentum
		Vector *vk;					// vk in [0], extra entries for XPIC and to copy pk
		ContactTerms *contactInfo;	// extrapolations for contact

		static int pkCopy;			// which vk vector to store pk
	
		// constants (not changed in MPM time step)
				
        // constructors and destructors
        MatVelocityField(int);
		~MatVelocityField();
		void Zero(void);
		void ZeroContactTerms(void);
	
		// methods
		void AddMomentumTask1(Vector *,Vector *,int);
		void CopyMassAndMomentum(NodalPoint *,int,int);
        void CopyMassAndMomentumLast(NodalPoint *,int,int);
		double GetTotalMassAndCount(void);
		void CopyGridForces(NodalPoint *,int,int);
		void RestoreMomenta(void);
		void ChangeMatMomentum(Vector *,int,double);
		void AddContactForce(Vector *);
		void GridValueCalculation(int);
		void AddGravityAndBodyForceTask3(Vector *);
#ifdef RESTART_OPTION
        bool IsTravelTooMuch(double,double) const;
#endif
		void AddPk(Vector *);
		void AddPkScaled(Vector *,double);
        void AddFtot(Vector *);
        void AddFtotScaled(Vector *,double);
        void UpdateMomentum(double);
		virtual void IncrementNodalVelAcc(double,GridToParticleExtrap *gp) const;
		void RezeroNodeTask6(void);
		const Vector *GetContactDispPtr(bool) const;
		void AddContactVector(int,Vector *,double);
		void AddContactVector(int,Vector *);

#if ADJUST_COPIED_PK == 1
        void AdjustForSymmetryBC(int);
#endif
	
		// accessors
		void Describe(int) const;
		void AddContactVolume(double);
		void SetContactVolume(double);
		double GetContactVolume(void) const;
        Vector GetVelocity(void);
		Vector *GetVStarPrev(void) const;
	
		// boundary conditions
        void ZeroVelocityBC(Vector *,int,double,Vector *);
        void AddVelocityBC(Vector *,double,int,double,Vector *);
	
        Vector GetFtot(void) const;
        Vector *GetFtotPtr(void);
		bool IsRigidField(void) const;
		void SetRigidField(bool);
		bool IsFixedRigidField(void) const;
		int GetFlags(void) const;
		bool IgnoresCracks(void) const;
	
		// class methods
		static bool ActiveField(MatVelocityField *);
		static bool ActiveRigidField(MatVelocityField *mvf);
		static bool ActiveNonrigidField(MatVelocityField *mvf);
		static bool ActiveNonrigidSourceField(MatVelocityField *,int);
		static bool ActiveNonrigidSeesCracksField(MatVelocityField *,bool);
	
		// XPIC
		virtual void XPICSupport(int,int,NodalPoint *,double,int,int,double);
		virtual void AddVStarNext(Vector *,double,double);

	protected:
		int flags;					// bitwise flags for some field properties

        Vector ftot;				// total force or contact force for rigid material
};

#endif
