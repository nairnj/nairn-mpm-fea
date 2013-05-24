/********************************************************************************
    MPMReadHandler.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 23 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.

	Dependencies
		CommonReadHandler.hpp
********************************************************************************/

#ifndef _MPMREADHANDLER_

#define _MPMREADHANDLER_

#include "Read_XML/CommonReadHandler.hpp"

class MPMBase;
class ParseController;

class MPMReadHandler : public CommonReadHandler
{
    public:
        //  Constructors and Destructor
        MPMReadHandler();
        ~MPMReadHandler();
    
        // Handlers for processing FEA data
		virtual bool myStartElement(char *xName,const Attributes& attrs);
		virtual void myEndElement(char *);
		virtual void myCharacters(char *,const unsigned int);
		virtual void TranslateBMPFiles(void);
        
        // My methods
        short GenerateInput(char *,const Attributes&);
        short EndGenerator(char *xName);
		short BMPFileInput(char *,const Attributes&);
 		void SetLevelVelocity(double,double,double);
        void MPMPts();							// Generate mpms
		void SetGIMPBorderAsHoles(void);		// implicit block edge GIMP elements from getting particles
		void MPMCracks(int,int,double,double,int,int);			// Generate cracks segments, Liping Xue
        void DisplacementBCs(void);				// Generate Displacement BCs
        void LoadBCs(void);						// Set load BCs
		char *LastBC(char *);
		void CreateAxisymetricBCs();
    
    private:
        //  Private data members
		ParseController *velocityBCs,*concBCs,*tempBCs,*mpLoadCtrl,*mpTractionCtrl,*mpConcFluxCtrl;
        char *currentTask;
        
		void grid(void);
};

void SetMptAnglesFromFunctions(int,Vector *,MPMBase *);		// 3D transformation angles
void ConvertToZYX(MPMBase *newMpt,double,double,double,double,double,double,double,double);

#endif


