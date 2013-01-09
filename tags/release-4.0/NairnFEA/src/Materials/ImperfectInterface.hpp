/********************************************************************************
    ImperfectInterface.hpp
    NairnMPM
    
    Created by John Nairn on Jan 08 2006.
    Copyright (c) 2061 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef _INTERFACEMAT_

#define _INTERFACEMAT_

#define INTERFACEPARAMS 5

#include "Materials/MaterialBase.hpp"

class ImperfectInterface : public MaterialBase
{
    public:
        double Dn,Dt;
        
        // constructors and destructors
        ImperfectInterface(char *);
        
        // methods
        virtual char *InputMat(char *,int &);
        virtual void PrintMechanicalProperties(void);
        virtual void InitialLoadMechProps(int,int);
		
		// accessors
		virtual int MaterialTag(void);
		virtual const char *MaterialType(void);

    private:
};

#endif
