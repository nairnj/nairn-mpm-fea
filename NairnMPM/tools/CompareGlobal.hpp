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

#define MAX_FILE_LINE 1000
#define MAX_COLUMNS 100
#define NUMBER_STATS 5
#define MAX_DIFF 1.0e-16
#define MAX_REL_DIFF 1.0e-7


using namespace std;

// error codes
enum { noErr=0, NoInputFileErr, BadOptionErr, FileAccessErr, MemoryErr, FileFormatErr };

// prototypes
char *NextArgument(int,char * const [],int,char);
void Usage(const char *);
void CompareGlobalFiles(const char *,int,int);
unsigned char *ReadLine(unsigned char *bptr,unsigned char *bend,unsigned char *rline,int *);
unsigned char *ReadGlobalFile(const char *,long &);
bool DbleEqual(double A, double B);

