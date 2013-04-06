/********************************************************************************
    MatVelocityField.hpp
    NairnMPM
    
    Created by John Nairn on 3 April 2009.
    Copyright (c) 2009 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _MATVELOCITYFIELD_

#define _MATVELOCITYFIELD_

#define MAX_FIELDS_FOR_CRACKS 4

class MatVelocityField
{
	public:
		// variables (changed in MPM time step)
		int numberPoints;			// number of material points in this field
		double mass;				// total mass of this field
		Vector pk;					// momentum
		Vector disp;				// displacement for contact calculations
		Vector *volumeGrad;			// mass gradient allocated in multimaterial mode
		bool rigidField;			// TRUE or FALSE if for rigid contact particles
	
		// constants (not changed in MPM time step)
				
        // constructors and destructors
        MatVelocityField(short);
		~MatVelocityField();
		void Zero(void);
		
		// methods
		void AddMomentumTask1(Vector *,Vector *);
		void ChangeMatMomentum(Vector *,bool,double);
		void AddContactForce(Vector *);
		void CalcVelocityForStrainUpdate(void);
		void CalcFtotTask3(double);
        void AddFint(Vector *);
        void AddFint(Vector *,double);
        void AddFext(Vector *);
        void AddFext(Vector *,double);
        void UpdateMomentum(double);
        void IncrementNodalVelAcc(double,Vector *,Vector *) const;
	
		// accessors
		void Describe(void);
		void AddContactVolume(double);
		void SetContactVolume(double);
		double GetContactVolume(void);
        void SetVelocity(Vector *);
        Vector GetVelocity(void);
        Vector *GetVelocityPtr(void);
        void SetMomentVelocityDirection(Vector *);
        void AddMomentVelocityDirection(Vector *,double);
        void SetFtotDirection(Vector *,double);
        void AddFtotDirection(Vector *,double,double);
        Vector GetFtot(void);
        Vector *GetFtotPtr(void);
	
		// class methods
		static bool ActiveField(MatVelocityField *);
		static bool ActiveNonrigidField(MatVelocityField *mvf);
		static bool ActiveRigidField(MatVelocityField *mvf);
	
	private:
		double volume;				// only for contact (cracks or multimaterial)

        Vector vk;					// velocity
        Vector ftot;				// total force or contact force for rigid material
        Vector fint;				// internal force (fint.x is last time step archived for rigid material)
        Vector fext;				// external force
};

#endif
