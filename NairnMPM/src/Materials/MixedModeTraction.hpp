/********************************************************************************
    MixedModeTraction.hpp
    nairn-mpm-fea

    Created by John Nairn on 6/20/2020.
    Copyright (c) 2020 John A. Nairn, All rights reserved.

    Dependencies
        TrilinearTraction.hpp,CohesiveZone.hpp, TractionLaw.hpp, MaterialBase.hpp
 ********************************************************************************/

#ifndef MIXEDMODETRACTIONMATERIAL

#define MIXEDMODETRACTIONMATERIAL 33

#include "Materials/TrilinearTraction.hpp"

enum {MM_D=0,MM_DELN,MM_DELT,MM_DW,MM_UN,MM_UT };

#define MMD_FOR_FAILURE 0.999

class MixedModeTraction : public TrilinearTraction
{
    public:

        // constructors and destructors
        MixedModeTraction(char *,int);
    
        // methods
        virtual char *InputTractionLawProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintMechanicalProperties(void) const;
        virtual char *InitHistoryData(char *);
    
        // traciton law methods
        virtual void CrackTractionLaw(CrackSegment *,double,double,Vector *,Vector *,double);
        virtual void GetStrengths(double *,double *,double *);
        virtual void GetRdeltaUc(double *,double *,double *);
        virtual void GetRPrimedelta(double *,double *,double *);
        virtual void MMGetDFromDelta(double *,bool);
        virtual void MMGetDeltaFromD(double *h,double *,double *);
        virtual void GetDeltaPrimeFromD(double *h,double *,double *);
        virtual void GetVarphis(double *,double *,double *);
        virtual double CrackWorkEnergy(CrackSegment *,double,double);
        virtual void CrackDissipatedEnergy(CrackSegment *,double &,double &);

        // accessors
        virtual const char *MaterialType(void) const;
    
    protected:
        int modelI,modelII;
        int useNewtonsMethod;
};

#endif

