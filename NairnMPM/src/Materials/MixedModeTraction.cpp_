/********************************************************************************
    MixedModeTraction.cpp
    nairn-mpm-fea

    Created by John Nairn on 60/20/2020.
    Copyright (c) 2020 John A. Nairn, All rights reserved.
 
    My method for coupled, mixed-mode cohesive zone with two
    separate strength models (which can each be different styles)
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Materials/MixedModeTraction.hpp"
#include "Materials/CubicTraction.hpp"
#include "Cracks/CrackSegment.hpp"

//#define DEBUG_NEWTON
#define NEWTON_CONVERGE 1.e-9
#define MAX_NEWTON_STEPS 10

extern double mtime;

#pragma mark MixedModeTraction::Constructors and Destructors

// Constructor
MixedModeTraction::MixedModeTraction(char *matName,int matID) : TrilinearTraction(matName,matID)
{
    // mode I cohesive law style (properties set in superclasses)
    modelI = TRIANGULARTRACTIONMATERIAL;
    
    // mode II cohesive law style (properties set in superclasses)
    modelII = TRIANGULARTRACTIONMATERIAL;
    
    // Activate Newton's method
    useNewtonsMethod = 0;
}

#pragma mark MixedModeTraction::Initialization

// no properties to read
char *MixedModeTraction::InputTractionLawProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"modelI")==0)
    {   input=INT_NUM;
        return (char *)&modelI;
    }
    
    else if(strcmp(xName,"modelII")==0)
    {   input=INT_NUM;
        return (char *)&modelII;
    }
    
    else if(strcmp(xName,"NewtonsMethod")==0)
    {   input=INT_NUM;
        return (char *)&useNewtonsMethod;
    }
    
    return TrilinearTraction::InputTractionLawProperty(xName,input,gScaling);
}

// Calculate properties used in analyses - here trilinear law
const char *MixedModeTraction::VerifyAndLoadProperties(int np)
{
    // mode I model
    const char *msg;
    switch(modelI)
    {   case TRIANGULARTRACTIONMATERIAL:
            msg = SetSawToothTractionLaw(stress1,kI1,delIc,JIc,umidI,phiI1,RI1);
            break;
        case TRILINEARTRACTIONMATERIAL:
        {   msg = SetTLTractionLaw(stress1,kI1,umidI,sI2,uI2,delIc,JIc,break1is2I,JI_1c,JI_2c);
            // See if break points are the same
            if(break1is2I)
            {   msg = "Mixed-mode traction material does not all delpkI=delpk2I (for stability reasons)";
                break;
            }
            // constant terms for phi(delta) and R(delta) and break point in terms of D
            double arg = (uI2-sI2*umidI/stress1)/(uI2-umidI);
            phiI1 = stress1*arg;
            RI1 = arg*umidI/delIc;
            phiI2 = sI2*delIc/(delIc-uI2);
            RI2 = phiI2*umidI/(stress1*delIc);
            DIbreak = 1. - sI2*umidI/(stress1*uI2);
            break;
        }
        case CUBICTRACTIONMATERIAL:
            msg = SetCubicTractionLaw(stress1,delIc,JIc,kI1);
            // change to initial stiffness (which different fhow cubic law finds it)
            kI1 = 27.*stress1/(4.*delIc);
            // for phi(delta)
            phiI1 = 2.*kI1*delIc;
            break;
        case EXPONENTIALTRACTIONMATERIAL:
            if(alphaI<=0.)
                msg = "Exponential cohesive laws must specify alphaI > 0";
            else
            {   expalphaI = exp(-alphaI);
                romexpalphaI = 1./(1.-expalphaI);
                double fealphaI = 2.*(1./alphaI - expalphaI*romexpalphaI);
                msg=SetExponentialTractionLaw(stress1,kI1,delIc,JIc,umidI,fealphaI);
            }
            break;
        default:
            msg = "Mode I strength model is not a valid type";
    }
    if(msg!=NULL) return msg;
    
    // mode II model
    switch(modelII)
    {   case TRIANGULARTRACTIONMATERIAL:
            msg = SetSawToothTractionLaw(stress2,kII1,delIIc,JIIc,umidII,phiII1,RII1);
            break;
        case TRILINEARTRACTIONMATERIAL:
        {   msg = SetTLTractionLaw(stress2,kII1,umidII,sII2,uII2,delIIc,JIIc,break1is2II,JII_1c,JII_2c);
            // See if break points are the same
            if(break1is2II)
            {   msg = "Mixed-mode traction material does not allow delpkII=delpk2II (for stability reasons)";
                break;
            }
            // constant terms for phi(delta) and R(delta) and break point in terms of D
            double arg = (uII2-sII2*umidII/stress2)/(uII2-umidII);
            phiII1 = stress2*arg;
            RII1 = arg*umidII/delIIc;
            phiII2 = sII2*delIIc/(delIIc-uII2);
            RII2 = phiII2*umidII/(stress2*delIIc);
            DIIbreak = 1. - sII2*umidII/(stress2*uII2);
            break;
        }
        case CUBICTRACTIONMATERIAL:
            msg = SetCubicTractionLaw(stress2,delIIc,JIIc,kII1);
            // change to initial stiffness (which different fhow cubic law finds it)
            kII1 = 27.*stress2/(4.*delIIc);
            // for phi(delta)
            phiII1 = 2.*kII1*delIIc;
            break;
        case EXPONENTIALTRACTIONMATERIAL:
            if(alphaII<=0.)
                msg = "Exponential cohesive laws must specify alphaII > 0";
            else
            {   expalphaII = exp(-alphaII);
                romexpalphaII = 1./(1.-expalphaII);
                double fealphaII = 2.*(1./alphaII - expalphaII*romexpalphaII);
                msg=SetExponentialTractionLaw(stress2,kII1,delIIc,JIIc,umidII,fealphaII);
            }
            break;
        default:
            msg = "Mode II strength model is not a valid type";
    }
    
    // check for only one direction cubic, if so, then force Newtons method
    if(modelI==CUBICTRACTIONMATERIAL || modelII==CUBICTRACTIONMATERIAL)
    {   if(modelI!=modelII) useNewtonsMethod = 1;
    }
    
    // exponential not set up for Newton's (could be if needed)
    if(useNewtonsMethod!=0 && msg!=NULL)
    {   if(modelI==EXPONENTIALTRACTIONMATERIAL || modelII==EXPONENTIALTRACTIONMATERIAL)
            msg = "Exponential cohesive laws not set up for use of Newton's method";
    }
    
    // exit no error
    if(msg!=NULL) return msg;
    
    // skip to parent class
    return TractionLaw::VerifyAndLoadProperties(np);
}

// print to output window
void MixedModeTraction::PrintMechanicalProperties(void) const
{
    switch(modelI)
    {   case TRIANGULARTRACTIONMATERIAL:
            cout << "Mode I: Trianglar" << endl;
            PrintSawToothModel("I",JIc,stress1,delIc,kI1,umidI,-1.);
            break;
        case TRILINEARTRACTIONMATERIAL:
            cout << "Mode I: Trilinear" << endl;
            PrintTriLinearModel("I",JIc,stress1,kI1,umidI,sI2,uI2,delIc);
            break;
        case EXPONENTIALTRACTIONMATERIAL:
            cout << "Mode I: Exponential" << endl;
            PrintSawToothModel("I",JIc,stress1,delIc,kI1,umidI,alphaI);
            break;
        default:
            cout << "Mode I: Cubic" << endl;
            PrintCubicModel("I",JIc,stress1,delIc,kI1);
            break;
    }
    switch(modelII)
    {   case TRIANGULARTRACTIONMATERIAL:
            cout << "Mode II: Triangular" << endl;
            PrintSawToothModel("II",JIIc,stress2,delIIc,kII1,umidII,-1.);
            break;
        case TRILINEARTRACTIONMATERIAL:
            cout << "Mode II: Trilinear" << endl;
            PrintTriLinearModel("II",JIIc,stress2,kII1,umidII,sII2,uII2,delIIc);
            break;
        case EXPONENTIALTRACTIONMATERIAL:
            cout << "Mode II: Exponential" << endl;
            PrintSawToothModel("II",JIIc,stress2,delIIc,kII1,umidII,alphaII);
            break;
        default:
            cout << "Mode II: Cubic" << endl;
            PrintCubicModel("II",JIIc,stress2,delIIc,kII1);
            break;
    }
    if(useNewtonsMethod!=0)
        cout << "Damage increments by Newton's method" << endl;

}

#pragma mark MixedModeTraction::History Data Methods

// history variables:
// h[0=MM_D] is D (<0 until initiation if both cubic)
// h[1=MM_DELN] is deltan
// h[2=MM_DELT] is deltat
// h[3=MM_DW] is dW or work increment
// h[4=MM_UN] is nCod stored to get dun on next step
// h[5=MM_UT] is tCod stored to get dut on next step
char *MixedModeTraction::InitHistoryData(char *pchr,MPMBase *mptr)
{
    numTractionHistory = 6;
    double *p = CreateAndZeroDoubles(pchr,numTractionHistory);
    if(modelI!=CUBICTRACTIONMATERIAL) p[MM_DELN]=umidI;
    if(modelII!=CUBICTRACTIONMATERIAL) p[MM_DELT]=umidII;
    return (char *)p;
}

#pragma mark MixedModeTraction::Traction Law

// Traction law - assume trianglar shape with unloading down slope back to the origin
void MixedModeTraction::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,Vector *n,Vector *t,double area)
{
    double Tn=0.,Tt=0.,dGI=0.,dGII=0.;
    double *h = (double *)cs->GetHistoryData();
    
    // get trial tractions (note that input nCod = un+dun, tCod=ut+dut)
    double D = h[MM_D];
    double Sn=0,St,testFailure;
    
    // transverse direction (note that dut and tCod may be negative)
    double tStiff = (1.-D)*kII1;
    double Tt0 = tStiff*h[MM_UT];
    Tt = tStiff*tCod;
    double dut = tCod-h[MM_UT];
    h[MM_UT] = tCod;

    // normal direction (will be zero if last step was compression and set h[MM_UN]=0)
    double nStiff = (1.-D)*kI1;
    double Tn0 = nStiff*h[MM_UN];
    double dun;
    if(nCod>0)
    {   // tensile half plane
        Tn = nStiff*nCod;
        dun = nCod-h[MM_UN];
        h[MM_UN] = nCod;
         
        // get strengths and test parameter
        GetStrengths(h,&Sn,&St);
        double tn = Tn/Sn,tt = Tt/St;
        testFailure = tn*tn + tt*tt;
    }
    else
    {   // compression half plane
        dun = 0.;
        h[MM_UN] = 0.;
        
        // get shear strength and test parameter
        GetStrengths(h,NULL,&St);
        testFailure = fabs(Tt)/St;
    }
    
    // is it elastic
    if(testFailure<=1.)
    {   // work increment (loading and unloading)
        h[MM_DW] += 0.5*((Tn0+Tn)*dun + (Tt0+Tt)*dut);
        
        // force is traction times area projected onto plane of unit vectors (units F)
        // tract = -area*(Tn*n + Tt*t)
        // In 2D, if t=(dx,dy), then n=(-dy,dx)
        cs->tract.x = -area*(Tn*n->x + Tt*t->x);
        cs->tract.y = -area*(Tn*n->y + Tt*t->y);
        cs->tract.z = -area*(Tn*n->z + Tt*t->z);

        return;
    }
    
    // partition into elastic increment and damage increment
    if(nCod>0 && Sn>0. && St>0.)
    {   double T0mag = Tn0*Tn0 + Tt0*Tt0;
        double Tn0hat = Tn0/T0mag;
        double Tt0hat = Tt0/T0mag;
        double TnSn = Tn0hat/Sn;
        double TtSt = Tt0hat/St;
        double that = TnSn*TnSn + TtSt*TtSt;
        double phiT0mag = 1./sqrt(that) - T0mag;
        double dune = phiT0mag*Tn0hat/nStiff;
        double dute = phiT0mag*Tt0hat/tStiff;
        
        // tractions on the failure surface T = (1+phi)T0
        Tn = Tn0hat*(T0mag + phiT0mag);
        Tt = Tt0hat*(T0mag + phiT0mag);
        
        // elastic work increment by midpoint rule
        h[MM_DW] += 0.5*((Tn0+Tn)*dune + (Tt0+Tt)*dute);
        
        // damage strain increment dutot = dud + due
        dun -= dune;
        dut -= dute;
        
        // New initial traction on the ellipse
        Tn0 = Tn;
        Tt0 = Tt;
    }
    else if(nCod<=0.)
    {   // elastic update when shear only (careful of signs)
        double absdute = (St - fabs(Tt0))/tStiff;
        
        // elastic work by midpoint rule
        h[MM_DW] += 0.5*(fabs(Tt0)+St)*absdute;
        
        // get increment causing damage
        if(dut>0.)
        {   dut -= absdute;
            Tt0 = St;
        }
        else
        {   dut += absdute;
            Tt0 = -St;
        }
    }
    
    // Now starting on failure surface with (Tn0,Tt0) and cod increment (dun,dut)
    // causing damage.
    double ddeltat,ddeltan;
    double varphin0,varphit0,varphin,varphit;

    // Bicubic is handled as special case
    if(modelI==CUBICTRACTIONMATERIAL && modelII==CUBICTRACTIONMATERIAL)
    {   if(nCod>0.)
        {   // tensile half plane - get initial dissipation rates for midpoint rule
            GetVarphis(h,&varphin0,&varphit0);
            
            if(h[MM_DELN]>0.)
            {   // Scale by (u/delta)^2
                double undn = nCod/h[MM_DELN];
                varphin0 *= undn*undn;
                double utdt = tCod/h[MM_DELT];
                varphit0 *= utdt*utdt;
            }
            
            // get dimensionless damage increment
            double unbar = nCod/delIc;
            double utbar = tCod/delIIc;
            double deltaBar = sqrt(unbar*unbar+utbar*utbar);
            if(deltaBar>=1.)
            {   deltaBar = 1.;
                cs->SetMatID(0);
            }

            // Set deltas and D
            double deltat = delIIc*deltaBar;
            ddeltat = deltat-h[MM_DELT];
            h[MM_DELT] = deltat;
            double deltan = delIc*deltaBar;
            ddeltan = deltan - h[MM_DELN];
            h[MM_DELN] = deltan;
            h[MM_D] = 1. - (1.-deltaBar)*(1.-deltaBar);
            
            // final traction
            Tn = (1-h[MM_D])*kI1*nCod;
            Tt = (1-h[MM_D])*kII1*tCod;
            
            // Work increment by midpoint rule (initial traction was zero)
            h[MM_DW] += 0.5*((Tn0+Tn)*dun + (Tt0+Tt)*dut);
            
            // energy release rate increments by midpoint rule
            GetVarphis(h,&varphin,&varphit);
            double arg = h[MM_UN]/h[MM_DELN];
            varphin *= arg*arg;
            arg = h[MM_UT]/h[MM_DELT];
            varphit *= arg*arg;
            dGI = 0.25*(varphin0+varphin)*ddeltan;
            dGII = 0.25*(varphit0+varphit)*ddeltat;
        }
        else
        {   // compression half plane - get initial varphi
            GetVarphis(h,NULL,&varphit0);
            
            // new deltaBar (watch for sign)
            double deltaBar = fabs(tCod)/delIIc;
            if(deltaBar>=1.)
            {   deltaBar = 1.;
                cs->SetMatID(0);
            }

            // Set deltas and D
            double deltat = delIIc*deltaBar;
            ddeltat = deltat-h[MM_DELT];
            h[MM_DELT] = deltat;
            double deltan = delIc*deltaBar;
            ddeltan = deltan - h[MM_DELN];
            h[MM_DELN] = deltan;
            h[MM_D] = 1. - (1.-deltaBar)*(1.-deltaBar);
            
            // final traction
            Tt = (1-h[MM_D])*kII1*tCod;
            
            // Work increment by midpoint rule (initial traction zero and dun=0)
            h[MM_DW] += 0.5*(Tt0+Tt)*dut;
            
            // energy release rate increments
            GetVarphis(h,NULL,&varphit);
            dGII = 0.25*(varphit0+varphit)*ddeltat;
        }
    }
    
    else
    {   // Damage update
        if(nCod>0)
        {   // tensile half plane
            
            if(useNewtonsMethod==0)
            {   // Solve with Taylor series expansion
                
                // Get initial dissipation rates for midpoint rule
                GetVarphis(h,&varphin0,&varphit0);
                
                // Scale by (u/delta)^2
                double undn = nCod/h[MM_DELN];
                varphin0 *= undn*undn;
                double utdt = tCod/h[MM_DELT];
                varphit0 *= utdt*utdt;

                // find most stable update
                double Rnunc,Rtutc;
                GetRdeltaUc(h,&Rnunc,&Rtutc);
                // Get u/delta^2 and update numerator for all
                double undn2 = undn/h[MM_DELN];
                double utdt2 = utdt/h[MM_DELT];
                double numer = undn2*dun + utdt2*dut;
                if(fmax(Rnunc,Rtutc)<1.)
                {   // update D
                    double denom = (delIc/Rnunc)*undn2*undn + (delIIc/Rtutc)*utdt2*utdt;
                    h[MM_D] += numer/denom;
                    if(h[MM_D]<MMD_FOR_FAILURE)
                    {   // calculate deltan and deltat from new D in h[0]
                        ddeltan = -h[MM_DELN];
                        ddeltat = -h[MM_DELT];
                        MMGetDeltaFromD(h,&h[MM_DELN],&h[MM_DELT]);
                        ddeltan += h[MM_DELN];
                        ddeltat += h[MM_DELT];
                    }
                    else
                    {   cs->SetMatID(0);
                        h[MM_D] = 1.;
                        ddeltan = delIc-h[MM_DELN];
                        h[MM_DELN] = delIc;
                        ddeltat = delIIc-h[MM_DELT];
                        h[MM_DELT] = delIIc;
                    }
                }
                else if(Rnunc>Rtutc)
                {   // update deltan
                    double denom = undn2*undn + (Rnunc*delIIc/(Rtutc*delIc))*utdt2*utdt;
                    ddeltan = numer/denom;
                    if(h[MM_DELN]+ddeltan<delIc)
                    {   h[MM_DELN] += ddeltan;
                        
                        // calculate D with updated h[MM_DELN]
                        MMGetDFromDelta(h,true);
                        
                        // calculate deltat from new D in h[MM_D]
                        ddeltat = -h[MM_DELT];
                        MMGetDeltaFromD(h,NULL,&h[MM_DELT]);
                        ddeltat += h[MM_DELT];
                    }
                    else
                    {   cs->SetMatID(0);
                        ddeltan = delIc-h[MM_DELN];
                        h[MM_DELN] = delIc;
                        h[MM_D] = 1.;
                        ddeltat = delIIc-h[MM_DELT];
                        h[MM_DELT] = delIIc;
                    }
                }
                else
                {   // update deltat
                    double denom = (Rtutc*delIc/(Rnunc*delIIc))*undn2*undn + utdt2*utdt;
                    ddeltat = numer/denom;
                    if(h[MM_DELT]+ddeltat<delIIc)
                    {   h[MM_DELT] += ddeltat;
                        
                        // calculate D with updated h[MM_DELT]
                        MMGetDFromDelta(h,false);
                        
                        // calculate deltan from new D in h[MM_D]
                        ddeltan = -h[MM_DELN];
                        MMGetDeltaFromD(h,&h[MM_DELN],NULL);
                        ddeltan += h[MM_DELN];
                    }
                    else
                    {   cs->SetMatID(0);
                        ddeltat = delIIc-h[MM_DELT];
                        h[MM_DELT] = delIIc;
                        h[MM_D] = 1.;
                        ddeltan = delIc-h[MM_DELN];
                        h[MM_DELN] = delIc;
                    }
                }
                
                // New tractions
                Tn = kI1*(1.-h[MM_D])*h[MM_UN];
                Tt = kII1*(1.-h[MM_D])*h[MM_UT];

                // work increment by midpoint rule
                h[MM_DW] += 0.5*((Tn0+Tn)*dun + (Tt0+Tt)*dut);
                
                // energy release rate increments by midpoint rule
                GetVarphis(h,&varphin,&varphit);
                double arg = h[MM_UN]/h[MM_DELN];
                varphin *= arg*arg;
                arg = h[MM_UT]/h[MM_DELT];
                varphit *= arg*arg;
                dGI = 0.25*(varphin0+varphin)*ddeltan;
                dGII = 0.25*(varphit0+varphit)*ddeltat;
            }
            
            else
            {   // Use Newton's Method

                // copy starting values
                double Dinit = h[MM_D];
                
                // solution bracketed by dD=0 and dD=1-D
                
                // take initial guess of dD based on dun or dut
                double deln,delt,Rn,Rt,delnp,deltp,Rnp,Rtp;
                GetRdeltaUc(h,&Rn,&Rt);
                double dD = fmax(Rn*dun/delIc,Rt*dut/delIIc);

#ifdef DEBUG_NEWTON
#pragma omp critical (output)
                {
                    cout << "Bracket 0 < " << dD << " < " <<
                        " based on dDn=" << Rn*dun/delIc << " and dDt=" << Rt*dut/delIIc << endl;
                    cout << "nCod=" << nCod << ", dun=" << dun << ", tCod=" << tCod << ", dut=" << dut << endl;
                    cout << "Rn=" << Rn << ", Rt=" << Rt << ", dun=" << dun << ", dut=" << dut <<
                        ", deln=" << h[MM_DELN] << ", delt=" << h[MM_DELT] << endl;
                }
#endif
                // Newton's method until done
                int step=1;
                while(true)
                {   // get deltas and R's
                    h[MM_D] = Dinit+dD;
                    MMGetDeltaFromD(h,&deln,&delt);
                    h[MM_DELN] = deln;
                    h[MM_DELT] = delt;
                    GetRdeltaUc(h,&Rn,&Rt);
                    Rn /= delIc;
                    Rt /= delIIc;
                    
                    // get function
                    double n1=nCod/deln,n2=dD/deln;
                    double narg = n1/(1. + n2/Rn);
                    double t1 = tCod/delt,t2=dD/delt;
                    double targ = t1/(1. + t2/Rt);
                    double gdD = narg*narg + targ*targ - 1.;
                    GetDeltaPrimeFromD(h,&delnp,&deltp);
                    GetRPrimedelta(h,&Rnp,&Rtp);
                    double nargp = -2.*narg*narg/(Rn*(Rn*deln+dD));
                    double targp = -2.*targ*targ/(Rt*(Rt*delt+dD));
                    double slope = nargp*(Rn + (Rn*Rn - dD*Rnp)*delnp)
                                    + targp*(Rt + (Rt*Rt - dD*Rtp)*deltp);
            
                    // new value (bisect if crosses zero) and check convergence
                    double dx = gdD/slope;
                    if(dx>dD) dx = 0.5*dD;
                    dD = dD - dx;
                
#ifdef DEBUG_NEWTON
#pragma omp critical (output)
                    {
                            cout << "...step " << step << " Newton to dD = " << dD << ", |dx| = " << dx << endl;
                    }
#endif
                    // watch for non-physical values and see if converged
                    if(dD>1.-Dinit) break;
                    if(fabs(dx)<NEWTON_CONVERGE || step>MAX_NEWTON_STEPS) break;
                    step++;
                }
            
#ifdef DEBUG_NEWTON
#pragma omp critical (output)
                {
                    cout << "...Converged in " << step << " steps to dD = " << dD << endl;
                }
#endif
                h[MM_D] = Dinit+dD;
                if(h[MM_D]<MMD_FOR_FAILURE)
                {   // calculate deltan and deltat from new D in h[0]
                    MMGetDeltaFromD(h,&h[MM_DELN],&h[MM_DELT]);
                }
                else
                {   cs->SetMatID(0);
                    dD = 1.-Dinit;
                    h[MM_D] = 1.;
                    h[MM_DELN] = delIc;
                    h[MM_DELT] = delIIc;
                }
                
                // New tractions
                Tn = kI1*(1.-h[MM_D])*h[MM_UN];
                Tt = kII1*(1.-h[MM_D])*h[MM_UT];

                // work increment by midpoint rule
                h[MM_DW] += 0.5*((Tn0+Tn)*dun + (Tt0+Tt)*dut);
                
                // energy release rate increments by midpoint rule
                dGI = 0.25*kI1*((h[MM_UN]-dun)*(h[MM_UN]-dun)+h[MM_UN]*h[MM_UN])*dD;
                dGII = 0.25*kII1*((h[MM_UT]-dut)*(h[MM_UT]-dut)+h[MM_UT]*h[MM_UT])*dD;
            }
        }
        else
        {   // compression half plane - get initial varphi
            GetVarphis(h,NULL,&varphit0);
            
            // update deltat (accounting for sign)
            ddeltat = fabs(dut);
            if(h[MM_DELT]+ddeltat<delIIc)
            {   h[MM_DELT] += ddeltat;
            
                // calculate D with update h[MM_DELT]
                MMGetDFromDelta(h,false);
            
                // calculate deltan from new D in h[MM_D]
                MMGetDeltaFromD(h,&h[MM_DELN],NULL);
            }
            else
            {   cs->SetMatID(0);
                ddeltat = delIIc-h[MM_DELT];
                h[MM_DELT] = delIIc;
                h[MM_D] = 1.;
                h[MM_DELN] = delIc;
            }
            
            // New traction
            Tt = kII1*(1.-h[MM_D])*h[MM_UT];
            
            // work increment (normal work zero because dun=0)
            h[MM_DW] += 0.5*(Tt0+Tt)*dut;
            
            // energy release rate increments
            GetVarphis(h,NULL,&varphit);
            dGII = 0.25*(varphit0+varphit)*ddeltat;
        }
    }
    
    // update dissipated energy
    cs->czmdG.x += dGI;
    cs->czmdG.y += dGII;
    cs->czmdG.z = 1.;
    
    // Report debond
    if(cs->MatID()<0)
    {   // it failed above in pure mode
        ReportDebond(mtime, cs, cs->czmdG.x/(cs->czmdG.x+cs->czmdG.y),cs->czmdG.x+cs->czmdG.y);
        Tn = 0.;                                       // turn off in tractions, if calculated
        Tt = 0.;
    }
    
    // force is traction times area projected onto plane of unit vectors (units F)
    // tract = -area*(Tn*n + Tt*t)
    // In 2D, if t=(dx,dy), then n=(-dy,dx)
    cs->tract.x = -area*(Tn*n->x + Tt*t->x);
    cs->tract.y = -area*(Tn*n->y + Tt*t->y);
    cs->tract.z = -area*(Tn*n->z + Tt*t->z);
}

// Return current traction law work energy (Int T.du).
//    This energy is needed for J integral (and only used in J Integral)
// units of F/L
double MixedModeTraction::CrackWorkEnergy(CrackSegment *cs,double nCod,double tCod)
{   double *h = (double *)cs->GetHistoryData();
    return h[MM_DW];
}

// Return mode I and II energy that has been released by current segment.
void MixedModeTraction::CrackDissipatedEnergy(CrackSegment *cs,double &GI,double &GII)
{   GI = cs->czmdG.x;
    GII = cs->czmdG.y;
}

#pragma mark MixedMode::Support Function Methods

// Traction law - get current normal and shear tractions
// Normal is skipped if Sn=NULL
void MixedModeTraction::GetStrengths(double *h,double *Sn,double *St)
{
    // mode I strength (optionally)
    if(Sn!=NULL)
    {   switch(modelI)
        {   case TRIANGULARTRACTIONMATERIAL:
                *Sn = CohesiveZone::Strength(1,h[MM_DELN]);
                break;
            case TRILINEARTRACTIONMATERIAL:
                *Sn = Strength(1,h[MM_DELN]);       // Trilinear is parent class
                break;
            case EXPONENTIALTRACTIONMATERIAL:
                *Sn = ExponentialTraction::Strength(1,h[MM_DELN]);
                break;
            default:
            {   // cubic
                double delbar = h[MM_DELN]/delIc;
                *Sn = kI1*h[MM_DELN]*(1.-delbar)*(1.-delbar);
                break;
            }
        }
    }
    
    // mode II strength
    if(St!=NULL)
    {   switch(modelII)
        {   case TRIANGULARTRACTIONMATERIAL:
                *St = CohesiveZone::Strength(2,h[MM_DELT]);
                break;
            case TRILINEARTRACTIONMATERIAL:
                *St = Strength(2,h[MM_DELT]);       // Trilinear is parent class
                break;
            case EXPONENTIALTRACTIONMATERIAL:
                *St = ExponentialTraction::Strength(2,h[MM_DELT]);
                break;
            default:
            {   // cubic
                double delbar = h[MM_DELT]/delIIc;
                *St = kII1*h[MM_DELT]*(1.-delbar)*(1.-delbar);
                break;
            }
        }
    }
}

// Traction law - Get D from one of the deltas and put it in h[MM_D]
void MixedModeTraction::MMGetDFromDelta(double *h,bool useModeI)
{
    double delbar;
    
    // Use mode I delta or mode II
    if(useModeI)
    {   switch(modelI)
        {   case TRIANGULARTRACTIONMATERIAL:
                h[MM_D] = CohesiveZone::GetDFromDelta(1,h[MM_DELN]);
                break;
            case TRILINEARTRACTIONMATERIAL:
                h[MM_D] = GetDFromDelta(1,h[MM_DELN]);              // Trilinear is parent class
                break;
            case EXPONENTIALTRACTIONMATERIAL:
                h[MM_D] = ExponentialTraction::GetDFromDelta(1,h[MM_DELN]);
                break;
            default:
                // cubic
                delbar = 1.-h[MM_DELN]/delIc;
                h[MM_D] = 1. - delbar*delbar;
                break;
        }
    }
    else
    {   switch(modelII)
        {   case TRIANGULARTRACTIONMATERIAL:
                h[MM_D] = CohesiveZone::GetDFromDelta(2,h[MM_DELT]);
                break;
            case TRILINEARTRACTIONMATERIAL:
                h[MM_D] = GetDFromDelta(2,h[MM_DELT]);              // Trilinear is parent class
                break;
            case EXPONENTIALTRACTIONMATERIAL:
                h[MM_D] = ExponentialTraction::GetDFromDelta(2,h[MM_DELT]);
                break;
            default:
                // cubic
                delbar = 1.-h[MM_DELT]/delIIc;
                h[MM_D] = 1. - delbar*delbar;
                break;
        }
    }
}

// Traction law - one or both deltas from current D
void MixedModeTraction::MMGetDeltaFromD(double *h,double *deln,double *delt)
{
    // mode I strength (optionally)
    if(deln!=NULL)
    {   switch(modelI)
        {   case TRIANGULARTRACTIONMATERIAL:
                *deln = CohesiveZone::GetDeltaFromD(1,h[MM_D]);
                break;
            case TRILINEARTRACTIONMATERIAL:
                *deln = GetDeltaFromD(1,h[MM_D]);       // Trilinear is parent class
                break;
            case EXPONENTIALTRACTIONMATERIAL:
                *deln = ExponentialTraction::GetDeltaFromD(1,h[MM_D]);
                break;
            default:
            {   // cubic
                double delProv = 1.-sqrt(1-h[MM_D]);
                if(delProv<1.e-5) delProv = 0.5*h[MM_D];
                *deln = delIc*delProv;
                break;
            }
        }
    }
    
    // mode II strength always
    if(delt!=NULL)
    {   switch(modelII)
        {   case TRIANGULARTRACTIONMATERIAL:
                *delt = CohesiveZone::GetDeltaFromD(2,h[MM_D]);
                break;
            case TRILINEARTRACTIONMATERIAL:
                *delt = GetDeltaFromD(2,h[MM_D]);       // Trilinear is parent class
                break;
            case EXPONENTIALTRACTIONMATERIAL:
                *delt = ExponentialTraction::GetDeltaFromD(2,h[MM_D]);
                break;
            default:
            {   // cubic
                double delProv = 1.-sqrt(1-h[MM_D]);
                if(delProv<1.e-5) delProv = 0.5*h[MM_D];
                *delt = delIIc*delProv;
                break;
            }
        }
    }
}

// Dissipation functions - one or both deltas from current delta's
void MixedModeTraction::GetVarphis(double *h,double *varphin,double *varphit)
{
    // mode I strength (optionally)
    if(varphin!=NULL)
    {   switch(modelI)
        {   case TRIANGULARTRACTIONMATERIAL:
                *varphin = phiI1;
                break;
            case TRILINEARTRACTIONMATERIAL:
                if(h[MM_DELN]<uI2)
                    *varphin = phiI1;
                else
                    *varphin = phiI2;
                break;
            case EXPONENTIALTRACTIONMATERIAL:
                *varphin = ExponentialTraction::DissipationRate(1,h[MM_DELN]);
                break;
            default:
            {   // cubic
                double xi = h[MM_DELN]/delIc;
                *varphin = phiI1*xi*xi*(1-xi);
                break;
            }
        }
    }
    
    // mode II strength always
    if(varphit!=NULL)
    {   switch(modelII)
        {   case TRIANGULARTRACTIONMATERIAL:
                *varphit = phiII1;
                break;
            case TRILINEARTRACTIONMATERIAL:
                if(h[MM_DELT]<uII2)
                    *varphit = phiII1;
                else
                    *varphit = phiII2;
                break;
            case EXPONENTIALTRACTIONMATERIAL:
                *varphit = ExponentialTraction::DissipationRate(2,h[MM_DELT]);
                break;
            default:
            {   // cubic
                double xi = h[MM_DELT]/delIIc;
                *varphit = phiII1*xi*xi*(1-xi);
                break;
            }
        }
    }
}

// Traction law - Get current R(delta)*uc (neither can be NULL)
// For linear R(delta)uc = R1/xi^2
// For trilinear R(delta)uc = R1/xi^2 for xi<u2/uc and R2/xi^2 otherwise
// For cubic R(delta)uc = 2(1-xi)
void MixedModeTraction::GetRdeltaUc(double *h,double *Rnunc,double *Rtutc)
{
    // mode I Rn*unc
    double xi = h[MM_DELN]/delIc;
    switch(modelI)
    {   case TRIANGULARTRACTIONMATERIAL:
            *Rnunc = RI1/(xi*xi);
            break;
        case TRILINEARTRACTIONMATERIAL:
            if(h[MM_DELN]<uI2)
                *Rnunc = RI1/(xi*xi);
            else
                *Rnunc = RI2/(xi*xi);
            break;
        case EXPONENTIALTRACTIONMATERIAL:
            *Rnunc = ExponentialTraction::RatioFunction(1,h[MM_DELN]);
            break;
        default:
            // cubic
            *Rnunc = 2.*(1-xi);
            break;
    }
    
    // mode II Rt*utc
    xi = h[MM_DELT]/delIIc;
    switch(modelII)
    {   case TRIANGULARTRACTIONMATERIAL:
            *Rtutc = RII1/(xi*xi);
            break;
        case TRILINEARTRACTIONMATERIAL:
            if(h[MM_DELT]<uII2)
                *Rtutc = RII1/(xi*xi);
            else
                *Rtutc = RII2/(xi*xi);
            break;
        case EXPONENTIALTRACTIONMATERIAL:
            *Rtutc = ExponentialTraction::RatioFunction(2,h[MM_DELT]);
            break;
        default:
            // cubic
            *Rtutc = 2.*(1-xi);
            break;
    }
}

#pragma mark MixedMode::Methods for Newtons Method

// Traction law - Get current R'(delta) (neither can be NULL)
// For linear R(delta) = R1 uc/del^2 or R'(delta) = -2*uc*R1/del^3
// For trilinear R(delta) = R1 uc/del^2 for del<u2 and R2 uc/del^2 otherwise
//         leading to R'(delta) = -2 R1 uc/del^3 and -2 R2 uc/del^3
// For cubic R(delta) = 2(uc-del)/uc^2 or R'(delta) = -2/uc^2
// Only used if try Newton's method for updates
void MixedModeTraction::GetRPrimedelta(double *h,double *Rnp,double *Rtp)
{
    // mode I Rn'
    double del = h[MM_DELN];
    switch(modelI)
    {   case TRIANGULARTRACTIONMATERIAL:
            *Rnp = -2.*delIc*RI1/(del*del*del);
            break;
        case TRILINEARTRACTIONMATERIAL:
            if(del<uI2)
                *Rnp = -2.*delIc*RI1/(del*del*del);
            else
                *Rnp = -2.*delIc*RI2/(del*del*del);
            break;
        default:
            // cubic
            *Rnp = -2./(delIc*delIc);
            break;
    }
    
    // mode II Rt'
    del = h[MM_DELT];
    switch(modelII)
    {   case TRIANGULARTRACTIONMATERIAL:
            *Rtp = -2.*delIIc*RII1/(del*del*del);
            break;
        case TRILINEARTRACTIONMATERIAL:
            if(del<uII2)
                *Rtp = -2.*delIIc*RII1/(del*del*del);
            else
                *Rtp = -2.*delIIc*RII2/(del*del*del);
            break;
        default:
            // cubic
            *Rtp = -2./(delIIc*delIIc);
            break;
    }
}

// Traction law - get both delta'(D) (both must be non-null)
// Only called when using newton's methog
void MixedModeTraction::GetDeltaPrimeFromD(double *h,double *delnp,double *deltp)
{
    // mode I
   switch(modelI)
    {   case TRIANGULARTRACTIONMATERIAL:
            *delnp = GetSawToothDeltaPrimeFromD(h[MM_D],delIc,umidI);
            break;
        case TRILINEARTRACTIONMATERIAL:
            *delnp = GetTLDeltaPrimeFromD(h[MM_D],DIbreak,stress1,sI2,umidI,uI2,delIc,break1is2I);
            break;
        default:
            // cubic
            *delnp = delIc/(2.*sqrt(1-h[MM_D]));
            break;
    }
    
    // mode II
    switch(modelII)
    {   case TRIANGULARTRACTIONMATERIAL:
            *deltp = GetSawToothDeltaPrimeFromD(h[MM_D],delIIc,umidII);
            break;
        case TRILINEARTRACTIONMATERIAL:
            *deltp = GetTLDeltaPrimeFromD(h[MM_D],DIIbreak,stress2,sII2,umidII,uII2,delIIc,break1is2II);
            break;
        default:
            // cubic
            *deltp = delIIc/(2.*sqrt(1-h[MM_D]));
            break;
    }
}

#pragma mark MixedModeTraction::Accessors

// return material type
const char *MixedModeTraction::MaterialType(void) const { return "Coupled Mixed Mode Cohesive Zone"; }

