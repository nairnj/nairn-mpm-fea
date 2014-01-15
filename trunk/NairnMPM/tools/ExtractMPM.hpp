/*********************************************************************
    ExtractMPM.hpp
    Nairn Research Group MPM and FEA Code
	MPM results extractor
    
    Created by John Nairn on Oct 23 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
*********************************************************************/

#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

using namespace std;

// quantity options
enum { SXX=0,SYY,SZZ,SXY,SXZ,SYZ,
		EXX,EYY,EZZ,EXY,EXZ,EYZ,
		PEXX,PEYY,PEZZ,PEXY,PEXZ,PEYZ,STRENG,PLENG,TEMP,CONC,
		DISPX,DISPY,DISPZ,VELX,VELY,VELZ,
		MAT,J1,J2,KI,KII,MASS,XPOS,YPOS,ZPOS,PRESSURE,EQUIVSTRESS,
        HIST1,HIST2,HIST3,HIST4,WORKENERGY,HEATENERGY };

// error codes
enum { noErr=0, NoInputFileErr, BadOptionErr, FileAccessErr, MemoryErr };

#define ARCH_ByteOrder 0
#define ARCH_Defaults 1

// Archiving options for material points
enum { ARCH_Velocity=2,ARCH_Stress,ARCH_Strain,ARCH_PlasticStrain,
            ARCH_OldOrigPosition,ARCH_WorkEnergy,ARCH_DeltaTemp,ARCH_PlasticEnergy,
            ARCH_ver2Empty, ARCH_ShearComponents, ARCH_StrainEnergy,
            ARCH_History, ARCH_Concentration,ARCH_HeatEnergy,ARCH_ElementCrossings,
            ARCH_RotStrain, ARCH_MAXMPMITEMS };

// Archiving options for crack segments
enum { ARCH_JIntegral=2,ARCH_StressIntensity,ARCH_BalanceResults,ARCH_MAXCRACKITEMS };

// prototypes
char *NextArgument(int,char * const [],int,char);
void Usage(const char *);
int ExtractMPMData(const char *,int,int);
void VTKLegacy(ostream &,unsigned char *,long,const char *);
void OutputQuantity(int,unsigned char *,ostream &,short,char);
short pointMatnum(unsigned char *);
bool skipThisPoint(short);
void OutputDouble(double *,int,char,bool,ostream &,int);
void OutputRecordEnd(ostream &,bool);
void BeginCrack(ostream &);
void EndCrack(ostream &);
void BeginMP(ostream &);
void EndMP(ostream &);
int CalcArchiveSize(int);
int Reverse(char *,int);

