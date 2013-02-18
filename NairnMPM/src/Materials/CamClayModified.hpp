//
//  CamClayModified.hpp
//  NairnMPM
//
//  Created by Raydel Lorenzo on 6/11/12.
//  Copyright (c) 2012 __Geotecnia-UNB__. All rights reserved.
//

/********************************************************************************
 Dependencies
 GeoMaterials.hpp (IsotropicMat.hpp, Elastic.hpp, MaterialBase.hpp)
 ********************************************************************************/

#ifndef CAMCLAYMODIFIED

#define CAMCLAYMODIFIED 101

#include "Materials/GeoMaterials.hpp"

class CamClayModified : public GeoMaterials
{
public:

    double popCC;
    double lamdaCC;                 //lamda Cam Clay
    double kapaCC;                  //kapa Cam Clay
    double Mcc;                     //M Cam Clay
    double e0;                      //Initial Void Ratio
    double phi;                      //internal friction angle of critical state 

    // constructors and destructors
    CamClayModified();
    CamClayModified(char *);
    
    // initialize
    virtual char *InputMat(char *,int &);
    virtual void InitialLoadMechProps(int,int);
    virtual void PrintMechanicalProperties(void);
    virtual void PrintYieldProperties(void);
    virtual const char *VerifyProperties(int);
    
    // override plastic potential functions

    
    virtual double GetFtrial(MPMBase *,int,Tensor *);                            //method return the valor of f for current conditions
    virtual double GetPcc(MPMBase *,int,Tensor *);                               //return mean stress I1/3   
    virtual double Getqcc(MPMBase *,int, Tensor *);                              //return shear stress q=sqrt(3*J2)
    virtual Tensor GetDfDsigmaGeo(MPMBase *,int, Tensor *);                      //return vector df/dsigma
    virtual double GetDfDdefvp(MPMBase *,int, Tensor *);                         //return value of df/ddeformaVolumPlas
    virtual Tensor GetDqdsigma(Tensor *);                                       //return dqdsigma
    virtual Tensor GetDpDsigma();                                               //return [1/3,1/3,1/3,0,0,0]    
    virtual double GetMultPlast(MPMBase *, int , Tensor *, double, double, double, double, double, double); //Return Multiplicador Plastico
    virtual double GetPoActual(MPMBase *, int);                                               //Return Po for each state for de volumetric plastic strain
    virtual Tensor GetElasticIncrement(MPMBase *, int, Tensor *, double, double, double, double, double, double);   //Elastic Stress Increment for trial update.

    
    // accessors
    virtual const char *MaterialType(void);
    virtual int MaterialTag();
                            
protected:
    // specific values of yield strength & modulus 
    double Mccred, lamdaCCred, kapaCCred,popCCred;
    
    bool readpopCC;
    bool readlamdaCC;
    bool readkapaCC;
    bool readMcc;
    bool reade0;
    bool readphi;

};

#endif


