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
		Interface2D(long,long *,int,double,double);
		
        // prototypes
        virtual int NumberSides(void);
        virtual double GetArea(void);
		virtual double GetVolume(void);
		virtual bool BulkElement(void);
        virtual short PtInElement(Vector &);
        virtual double GetThickness(void);
		virtual void Stiffness(int);
		virtual void IncrementStiffnessElements(double,double *,double,double,double,double,double);
		virtual void ForceStress(double *,int,int);
	
	private:
		void InterfaceTraction(int,int,double,double,double,double,double);

 };

#endif

