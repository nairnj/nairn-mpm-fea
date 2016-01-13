/********************************************************************************
    Interface2D.hpp
    NairnFEA
    
    Created by John Nairn on 01/07/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
	
	Dependencies
		ElementBase.hpp
********************************************************************************/

#ifndef _INTERFACE_

#define _INTERFACE_

#include "Elements/ElementBase.hpp"

class Interface2D : public ElementBase
{
    public:
		double thickness;
		
        // constructors
		Interface2D(int,int *,int,double,double);
		
        // prototypes
		virtual bool BulkElement(void);
        virtual short PtInElement(Vector &) const;
        virtual void SetThickness(double);  
		virtual void Stiffness(int);
		virtual void IncrementStiffnessElements(double,double *,double,double,double,double,double);
		virtual void ForceStress(double *,int,int);
        virtual void FindExtent(void);
	
		// const methods
		virtual double GetThickness(void) const;
		virtual double GetArea(void) const;
		virtual double GetVolume(void) const;
		virtual int NumberSides(void) const;
	
	private:
		void InterfaceTraction(int,int,double,double,double,double,double);

 };

#endif

