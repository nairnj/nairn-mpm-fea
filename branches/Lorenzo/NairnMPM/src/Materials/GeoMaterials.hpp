//
//  GeoMaterials.hpp
//  NairnMPM
//
//  Created by Raydel Lorenzo on 6/11/12.
//  Copyright (c) 2012 __Geotecnia-UNB__. All rights reserved.
//

#ifndef NairnMPM_GeoMaterials_h

#define NairnMPM_GeoMaterials_h



/********************************************************************************
 Dependencies
 IsotropicMat.hpp (Elastic.hpp, MaterialBase.hpp)
 ********************************************************************************/


#define SQRT_TWOTHIRDS 0.8164965809277260
#define TWOTHIRDS 0.6666666666666667
#define ONETHIRD 0.3333333333333333
#define SQRT_EIGHT27THS 0.5443310539518174
#define EULER 2.71828182845905


#include "Materials/IsotropicMat.hpp"

class GeoMaterials : public IsotropicMat
{
public:
    // one yield stress for isotropic, plastic materials
    double yield;
    double yref;                            //ground surface level (mm)
    double OCR;                             // Over Consolidation Ratio
    double k0;                              // earth's pressure coefficient at rest

    // constructors and destructors
    GeoMaterials();
    GeoMaterials(char *matName);
    
    // initialize
    virtual char *InputMat(char *,int &);
    virtual const char *VerifyProperties(int);
    virtual void InitialLoadMechProps(int,int);
    virtual char *InitHistoryData(void);
    virtual void PrintYieldProperties(void) = 0;							// subclass must provide
    virtual void SetInitialParticleState(MPMBase *,int);
    
    // methods
    virtual void MPMConstLaw(MPMBase *,double,double,double,double,double, double,int);
    virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
    
    // custom methods: Find yield function and solve for lambda
 
 

    virtual void ElasticUpdateFinished(MPMBase *,int,double);
    
    virtual double GetFtrial(MPMBase *,int, Tensor *) = 0;          // Evaluate in the yield surface, subclass must provide (mat,anali)
    virtual Tensor GetDfDsigmaGeo(MPMBase *,int,Tensor *)=0;        // return un vector with tud six derivatives with respect to sigma, subclass must provide (mat,anal,&stk)
    virtual double GetDfDdefvp(MPMBase *,int,Tensor *)=0;           //Return valou of df/ddefvolplas, subclass must provide
    virtual double GetMultPlast(MPMBase *,int,Tensor *, double, double, double, double, double, double)=0;          //Return plastic multiplicator, subclass must provide (mat,tensoes, total strain)
    virtual Tensor GetElasticIncrement(MPMBase *, int, Tensor *, double, double, double, double, double, double)=0;        // Get elastic stress increment for trial update.
    //virtual Tensor CorrectDrift(MPMBase *, int, Tensor *, Tensor *);
    
    // Default internal variable as cumlative plastic strain
    virtual void UpdateDefVolPlas(MPMBase *);                   //Update plastic volumentric strain GEO
    virtual void StoragePlasticInternal(MPMBase *);             //storage plastic volumentric strain    GEo
    
    // accessors

    
protected:
    double yldred,Gred,Kred;
    double psRed,psLr2G,psKred;
    bool readyref, readOCR, readk0, readYield;
    double ddvplas, dvplas, pp;              // delta and accumulate volumetric plastic strain...
    double ftrial;                           //yield function
};

#endif

