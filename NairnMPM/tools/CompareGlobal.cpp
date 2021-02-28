/*********************************************************************
    CompareGlobal.cpp
    Nairn Research Group MPM and FEA Code
	Compare to global files numerically
    
    Created by John Nairn on Oct 13, 2018.
    Copyright (c) 2018 John A. Nairn, All rights reserved.
	
 	CompareGlobal -o filename (filename)
*********************************************************************/

#define _CRT_SECURE_NO_WARNINGS

#include "CompareGlobal.hpp"

// Global variable settings
char *origfile=NULL;
double minMean = 0.;

#pragma mark MAIN AND INPUT PARAMETERS

// main entry point
int main(int argc, char * const argv[])
{
    // Requires one argument after options for file name
    if(argc<2)
    {	Usage("Input file name is missing");
        return NoInputFileErr;
    }
	
	// initialize
	
    // Check for options
	int parmInd;
	unsigned long opt;
	char *parm;
    for(parmInd=1;parmInd<argc && argv[parmInd][0]=='-';parmInd++)
	{	// each character in the next argument (skipping the '-')
		unsigned long optNum = strlen(argv[parmInd]);
		for(opt=1;opt<optNum;opt++)
		{	// Help request
			if(argv[parmInd][opt]=='H' || argv[parmInd][opt]=='h')
			{	Usage(NULL);
				return noErr;
			}
			
			// original file name
			else if(argv[parmInd][opt]=='o')
			{	parm=NextArgument(++parmInd,argv,argc,'o');
				if(parm==NULL) return 1;
				if(origfile!=NULL) delete [] origfile;
				origfile=new char[strlen(parm)+1];
				strcpy(origfile,parm);
			}
			
			// ignore values with SM < 1e-(value)
			else if(argv[parmInd][opt]=='m')
			{	parm=NextArgument(++parmInd,argv,argc,'m');
				if(parm==NULL) return BadOptionErr;
				int mpow;
				sscanf(parm,"%d",&mpow);
				if(mpow<=0)
				{   cerr << "CompareGlobal option 'm' must be a positive integer" << endl;
					return BadOptionErr;
				}
				minMean = 2.*pow(10.,(double)(-mpow));
				break;
			}
			
			else
			{   cerr << "Unknown CompareGlobal option '" << argv[parmInd][opt] << "' was used\n";
				return BadOptionErr;
			}
		}
    }
	
	//  Original file is needed
	if(origfile == NULL)
	{	Usage("Original file name is missing");
		return NoInputFileErr;
	}
	
    //  Last parameter must be the file name
    if(parmInd>=argc)
    {	Usage("Comparison file name is missing");
        return NoInputFileErr;
    }
	
	// compare each files in the list
	int fileNum;
    for(fileNum=parmInd;fileNum<argc;fileNum++)
	{	CompareGlobalFiles(argv[fileNum],fileNum-parmInd,argc-parmInd-1);
	}
	
    // insert code here...
    return noErr;
}

// grab next argument or quit with error is not found
char *NextArgument(int parmInd,char * const argv[],int argc,char option)
{
	// error if not there
	if(parmInd>=argc)
	{   cerr << "CompareGlobal option '" << option << "' is missing its required argument.\n";
		return NULL;
	}
	return argv[parmInd];
}
	

// Explain usage of this program
void Usage(const char *msg)
{
	if(msg!=NULL)
		cout << "\nERROR: " << msg << endl;
		
	cout << "\nCompareGlobal\n    version 1.0.0" << endl;
    cout << "\nUsage:\n"
        "    CompareGlobal -o <OriginalFile> <CompareFile>\n\n"
        "This program compares MPM global results in <CompareFile> numerically\n"
        "to the MPM global results in <OriginalFile>\n\n"
        "  Options:\n"
		"    -o path            Original file name (required)\n"
		"		                       (quoted if name has spaces)\n"
		"    -m pow             Ignore values with symmetric mean < 1.e(-pow)\n"
        "    -H (or -?)         Show this help and exit\n"
		"\n"
          <<  endl;
}

#pragma mark EXTRACTION CODE

// main entry to extract results from one file
void CompareGlobalFiles(const char *mpmFile,int fileIndex,int lastIndex)
{
	cout << "Original File: " << origfile << endl;
	cout << "Comparison File: " << mpmFile << endl;
	
	// read the two files
	long fileLength;
	unsigned char *orig = ReadGlobalFile(origfile,fileLength);
	if(orig == NULL) return;
	unsigned char *origEnd = orig+fileLength;
	
	unsigned char *cmp = ReadGlobalFile(mpmFile,fileLength);
	if(cmp == NULL) return;
	unsigned char *cmpEnd = cmp+fileLength;
	
	// line buffer
	unsigned char origCols[MAX_FILE_LINE];
	int nCols = 0;
	int origColOffs[MAX_COLUMNS];
	
	// read to column labels in original file
	while(true)
	{	orig = ReadLine(orig,origEnd,origCols,origColOffs);
		if(orig==NULL)
		{	cout << "Could not find column labels in original file" << endl;
			return;
		}
		if(strcmp((const char *)origCols,"#setName")==0)
		{	// count columns
			// for i for 0 to < numCols, col[i(0 based)] at origCols+origColOffs[i]
			// time is at zero offset
			while(true)
			{	if(origColOffs[nCols]==0) break;
				nCols++;
			}
			break;
		}
	}
	
	// more line buffers
	unsigned char cmpLine[MAX_FILE_LINE],origLine[MAX_FILE_LINE];
	int cmpOffs[MAX_COLUMNS],origOffs[MAX_COLUMNS],map[MAX_COLUMNS],nCmpCols=0;
	bool hasMatch = false;
	
	// read target columns and map to original
	while(true)
	{	cmp = ReadLine(cmp,cmpEnd,cmpLine,cmpOffs);
		if(cmp==NULL)
		{	cout << "Could not find column labels in comparison file" << endl;
			return;
		}
		if(strcmp((const char *)cmpLine,"#setName")==0)
		{	// map columns
			while(true)
			{	if(cmpOffs[nCmpCols]==0)
				{	// non left in comparison file
					break;
				}
				map[nCmpCols] = -1;		// set not mapped
				for(int i=0;i<nCols;i++)
				{	if(strcmp((const char *)(cmpLine+cmpOffs[nCmpCols]),
							  	(const char *)(origCols+origColOffs[i]))==0)
					{	// found mapping
						map[nCmpCols] = i;
						hasMatch = true;
						break;
					}
				}
				if(map[nCmpCols]<0)
				{	cout << "MISMATCH: column labeled " << (cmpLine+cmpOffs[nCmpCols])
						<< " not found in the original file" << endl;
				}
				nCmpCols++;
			}
			break;
		}
	}
	
	// print error messages, exit if not match
	if(nCols<nCmpCols)
	{	cout << "MISMATH: comparison file has more columns than original file" << endl;
	}
	else if(nCols>nCmpCols)
	{	cout << "MISMATH: comparison file has fewer columns than original file" << endl;
	}
	if(!hasMatch)
	{	cout << "MISMATCH: no comparison file column matches any original file column" << endl;
		return;
	}
	
	// Evaulate matching columns
	double stats[MAX_COLUMNS+1][NUMBER_STATS+1];
	for(int i=0;i<=nCols;i++)
	{	for(int j=0;j<NUMBER_STATS+1;j++)
			stats[i][j]=0.;
	}
	
	// read each row
	int nrows = 0;
	int nmismatches = 0,nevaluated = 0;
	while(true)
	{	// read one row each file
		orig = ReadLine(orig,origEnd,origLine,origOffs);
		cmp = ReadLine(cmp,cmpEnd,cmpLine,cmpOffs);
		
		// if one is done then don
		if(orig==NULL || cmp==NULL)
		{	if(orig!=NULL)
				cout << "MISMATCH: comparison file has fewer rows than original file" << endl;
			else if(cmp!=NULL)
				cout << "MISMATCH: comparison file has more rows than original file" << endl;
			break;
		}
		
		// check each from
		double origVal,cmpVal;
		int i;
		int stringMatch;
		for(int j=-1;j<nCmpCols;j++)
		{	if(j==-1)
			{	// get time column
				sscanf((const char *)origLine,"%lf",&origVal);
				sscanf((const char *)cmpLine,"%lf",&cmpVal);
				
				// store in 0
				i=0;
				
				// are strings the same
				stringMatch = strcmp((const char *)origLine,(const char *)cmpLine);
			}
			else
			{	if(map[j]<0) continue;
			
				// compare column j of comparison to i fo original
				i = map[j];
			
				// problem if file lacks the columns
				sscanf((const char *)(origLine+origOffs[i]),"%lf",&origVal);
				sscanf((const char *)(cmpLine+cmpOffs[j]),"%lf",&cmpVal);
				
				// are strings the same
				stringMatch = strcmp((const char *)(origLine+origOffs[i]),
									 	(const char *)(cmpLine+cmpOffs[j]));
				if(stringMatch!=0) nmismatches++;
				nevaluated++;
				
				// store in i+1;
				i++;
			}
			
			// get terms
			double diff = fabs(origVal-cmpVal);
			double mean = fabs(origVal)+fabs(cmpVal);
			
			// [4] = max(mean)
			if(mean>stats[i][4])
				stats[i][4] = mean;
			
			// stats if different
			// [0] = sum fabs(diff)
			// [1] = max(diff)
			// [2] = sum fabs(diff/mean)
			// [3] = max(diff/mean)
			if(stringMatch!=0)
			{	stats[i][0] += diff;
				if(diff>stats[i][1])
					stats[i][1] = diff;
				// ignore points much smaller than mean (won't help initial points)
				if(mean>minMean && mean>1.e-12*stats[i][4])
				{	double relDiff = 2.*diff/mean;
					stats[i][2] += relDiff;
					if(relDiff>stats[i][3])
						stats[i][3] = relDiff;
				}
			}
		}
		
		nrows++;
	}
	
	if(nrows==0)
	{	cout << "No rows of data found for comparison" << endl;
		return;
	}
	
	cout << "Tabulated " << nrows << " rows of results" << endl;
	if(nevaluated>0)
	{	cout << "Found " << (100.*(double)nmismatches/(double)nevaluated)
					<< "% of entries mismatched" << endl;
	}
	if(minMean>0.)
		cout << "Ignore values with symmetric mean < " << 0.5*minMean << endl;
	cout << "-----------------------------------------------------------------------------------\n";
	cout << "Column Name          Max SM         MAE        Max AE      SMAPE (%)   Max APE (%)\n";
	cout << "-----------------------------------------------------------------------------------\n";

	// print summary
	char rep[200];
	for(int i=0;i<=nCols;i++)
	{	if(i==0)
			strcpy(rep,"Time ");
		else
			strcpy(rep,(const char *)(origCols+origColOffs[i-1]+1));
		long clen = strlen(rep)-1;
		for(int i=(int)clen;i<17;i++)
			rep[i] = ' ';
		rep[17]=0;
		
		// name
		cout << rep;
		
		// get results (max, MAE, MaxAE, SMAPE(%), MaxAPE(%)
		double MAE = stats[i][0]>0. ? stats[i][0]/((double)nrows) : 0.;
		double SMAPE = stats[i][2]>0. ? 100.*stats[i][2]/((double)nrows) : 0.;
		sprintf(rep," %12g %12g %12g %12g %12g",0.5*stats[i][4],MAE,stats[i][1],SMAPE,100.*stats[i][3]);
		cout << rep << endl;
	}
}

// Read next line (to next CR or LF)
// look for tabs and set offsets to column data
unsigned char *ReadLine(unsigned char *bptr,unsigned char *bend,unsigned char *rline,int *offset)
{
	// skip leading returns
	while(bptr<bend)
	{	if(*bptr!='\n' && *bptr!='\r') break;
		bptr++;
	}
	
	// return NULL if nothing left
	if(bptr>=bend) return NULL;
	
	// file the line
	int lineLength = 0,numTabs=0;
	while(bptr<bend && lineLength<MAX_FILE_LINE)
	{	// end of the line
		if(*bptr=='\n' || *bptr=='\r') break;
		*rline = *bptr++;
		lineLength++;
		
		// was it a tab
		if(*rline == '\t')
		{	offset[numTabs++] = lineLength;
			*rline = 0;
			if(numTabs==MAX_COLUMNS)
			{	cout << "Line in file has more than " << MAX_COLUMNS << " columns." << endl;
				return NULL;
			}
		}
		
		// next character
		rline++;
	}
	
	if(lineLength>=MAX_FILE_LINE)
	{	cout << "Line in file has more than " << MAX_FILE_LINE << " characters." << endl;
		return NULL;
	}
	
	// line is full, make it a string, and indicate last tab
	*rline = 0;
	offset[numTabs] = 0;
	
	// return new pointer
	return bptr;
}

// open file, read to  buffer, close file, return the buffer
unsigned char *ReadGlobalFile(const char *filename,long &fileLength)
{
	FILE *fo;
	
	// open the original file
	if((fo=fopen(filename,"r"))==NULL)
	{	cerr << "Target file '" << filename << "' could not be opened" << endl;
		return NULL;
	}
	
	// get file length and read it
	if(fseek(fo,0L,SEEK_END)!=0)
	{	cerr << "Target file '" << filename << "' access error (fseek())" << endl;
		fclose(fo);
		return NULL;
	}
	fileLength=ftell(fo);
	rewind(fo);
	
	// create buffer
	unsigned char *buffer=(unsigned char *)malloc(fileLength);
	if(buffer==NULL)
	{	cerr << "Out of memory creating file reading buffer" << endl;
		fclose(fo);
		return NULL;
	}
	
	// read the file
	if(fread(buffer,fileLength,1,fo)!=1)
	{	cerr << "Target file '" << filename << "' access error (not all read)" << endl;
		fclose(fo);
		return NULL;
	}
	
	return buffer;
	
}

bool DbleEqual(double A, double B)
{
	// Check if the numbers are really close -- needed
	// when comparing numbers near zero.
	double diff = fabs(A - B);
	if (diff <= MAX_DIFF)
		return true;
	
	// Check relative difference
	A = fabs(A);
	B = fabs(B);
	double largest = (B > A) ? B : A;
	if (diff <= largest * MAX_REL_DIFF)
		return true;
	
	// assuming unequal
	return false;
}
