/*********************************************************************
    ExtractMPM.cpp
    Nairn Research Group MPM and FEA Code
	MPM results extractor
    
    Created by John Nairn on Oct 23 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	To Do:
		stress (for all), strain (for all)
		more quantities
	
	Notes:
		Reads header of ver4 files, written since 10/24/07
		Requires -b for ver3 files
		Does not support ver2 or earlier files
*********************************************************************/

#define _CRT_SECURE_NO_WARNINGS

#include "ExtractMPM.hpp"

// Global variable settings
char fileFormat='T';
char *outfile=NULL;
bool stepExtension=false;
char *headerName=NULL;
char fileExtension[20];
char mpmOrder[50];
char crackOrder[50];
char stepNum[50];
vector< int > includeMaterial;
vector< int > excludeMaterial;
vector< int > quantity;
vector< char * > quantityName;
int singleMaterial=0;
char thisEndian;
bool threeD=false;
bool header=false;
bool crackDataOnly=false;
bool hasTempQ=false;
bool hasConcQ=false;
bool isStructured=false;
float archiveTimeMs=-1.;
double zeroDouble=0.;
int recSize;
int vectorSize,tensorSize;
bool reverseFromInput;
int angleOffset=-1,angleYOffset,angleXOffset,posOffset=-1,stressOffset=-1,strainOffset=-1;
int crackPosOffset=-1,jIntOffset=-1,kSifOffset=-1;
int velocityOffset=-1,origPosOffset=-1,plStrainOffset=-1;
int tempOffset=-1,concOffset=-1,strainEnergyOffset=-1,plasticEnergyOffset=-1;
int history1Offset=-1,history2Offset=-1,history3Offset=-1,history4Offset=-1;
int workEnergyOffset=-1,heatEnergyOffset=-1,rotStrainOffset=-1;

// global file variables
unsigned char *buffer;
long fileLength,seekOffset,blockSize;
FILE *fp;
unsigned char *ap,*apNextBlock;

// 50 MB chunks
#define MAX_READ_BLOCK 50000000

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
	mpmOrder[0]=0;
	crackOrder[0]=0;
	strcpy(fileExtension,"txt");
	
    // Check for options
	int parmInd;
	unsigned int optInd;
	char *parm;
    for(parmInd=1;parmInd<argc && argv[parmInd][0]=='-';parmInd++)
	{	// each character in the next argument (skipping the '-')
		for(optInd=1;optInd<strlen(argv[parmInd]);optInd++)
		{	// Help request
			if(argv[parmInd][optInd]=='H' || argv[parmInd][optInd]=='?')
			{	Usage(NULL);
				return noErr;
			}
			
			// quantity (must be q, space, text name of quantity)
			else if(argv[parmInd][optInd]=='q')
			{	parm=NextArgument(++parmInd,argv,argc,'q');
				if(parm==NULL) return BadOptionErr;
				int q;
				if(strcmp(parm,"sxx")==0)
					q=SXX;
				else if(strcmp(parm,"syy")==0)
					q=SYY;
				else if(strcmp(parm,"szz")==0)
					q=SZZ;
				else if(strcmp(parm,"sxy")==0 || strcmp(parm,"syx")==0)
					q=SXY;
				else if(strcmp(parm,"sxz")==0 || strcmp(parm,"szx")==0)
					q=SXZ;
				else if(strcmp(parm,"syz")==0 || strcmp(parm,"szy")==0)
					q=SYZ;
				else if(strcmp(parm,"stress")==0)
					q=STRESS;
				else if(strcmp(parm,"pressure")==0)
                    q=PRESSURE;
                else if(strcmp(parm,"vonmises")==0)
                    q=EQUIVSTRESS;
                else if(strcmp(parm,"equivstrain")==0)
                    q=EQUIVSTRAIN;
				else if(strcmp(parm,"exx")==0)
					q=EXX;
				else if(strcmp(parm,"eyy")==0)
					q=EYY;
				else if(strcmp(parm,"ezz")==0)
					q=EZZ;
				else if(strcmp(parm,"exy")==0 || strcmp(parm,"eyx")==0)
					q=EXY;
				else if(strcmp(parm,"exz")==0 || strcmp(parm,"ezx")==0)
					q=EXZ;
				else if(strcmp(parm,"eyz")==0 || strcmp(parm,"ezy")==0)
					q=EYZ;
				else if(strcmp(parm,"strain")==0)
					q=STRAIN;
				else if(strcmp(parm,"pexx")==0)
					q=PEXX;
				else if(strcmp(parm,"peyy")==0)
					q=PEYY;
				else if(strcmp(parm,"pezz")==0)
					q=PEZZ;
				else if(strcmp(parm,"pexy")==0 || strcmp(parm,"peyx")==0)
					q=PEXY;
				else if(strcmp(parm,"pexz")==0 || strcmp(parm,"pezx")==0)
					q=PEXZ;
				else if(strcmp(parm,"peyz")==0 || strcmp(parm,"pezy")==0)
					q=PEYZ;
				else if(strcmp(parm,"plasticstrain")==0)
					q=PLASTICSTRAIN;
				else if(strcmp(parm,"strerg")==0)
					q=STRENG;
				else if(strcmp(parm,"plerg")==0)
					q=PLENG;
				else if(strcmp(parm,"temp")==0)
				{	q=TEMP;
					hasTempQ=true;
				}
				else if(strcmp(parm,"conc")==0)
				{	q=CONC;
					hasConcQ=true;
				}
				else if(strcmp(parm,"velx")==0)
					q=VELX;
				else if(strcmp(parm,"vely")==0)
					q=VELY;
				else if(strcmp(parm,"velz")==0)
					q=VELZ;
				else if(strcmp(parm,"velocity")==0)
					q=VELOCITY;
				else if(strcmp(parm,"dispx")==0)
					q=DISPX;
				else if(strcmp(parm,"dispy")==0)
					q=DISPY;
				else if(strcmp(parm,"dispz")==0)
					q=DISPZ;
				else if(strcmp(parm,"displacement")==0)
					q=DISPLACEMENT;
				else if(strcmp(parm,"mat")==0)
					q=MAT;
				else if(strcmp(parm,"j1")==0)
					q=J1;
				else if(strcmp(parm,"j2")==0)
					q=J2;
				else if(strcmp(parm,"ki")==0)
					q=KI;
				else if(strcmp(parm,"kii")==0)
					q=KII;
				else if(strcmp(parm,"mass")==0)
					q=MASS;
				else if(strcmp(parm,"hist1")==0)
					q=HIST1;
				else if(strcmp(parm,"hist2")==0)
					q=HIST2;
				else if(strcmp(parm,"hist3")==0)
					q=HIST3;
				else if(strcmp(parm,"hist4")==0)
					q=HIST4;
				else if(strcmp(parm,"work")==0)
					q=WORKENERGY;
				else if(strcmp(parm,"heat")==0)
					q=HEATENERGY;
				else if(strcmp(parm,"angz")==0)
					q=MATANGLEZ;
				else if(strcmp(parm,"angy")==0)
					q=MATANGLEY;
				else if(strcmp(parm,"angx")==0)
					q=MATANGLEX;
				else
				{   cerr << "ExtractMPM option 'q' argument of '" << parm << "' is not recognized" << endl;
					return BadOptionErr;
				}
				
				quantity.push_back(q);
				quantityName.push_back(parm);
				break;
			}

			// exclude material (must be m, space, number)
			else if(argv[parmInd][optInd]=='m')
			{	parm=NextArgument(++parmInd,argv,argc,'m');
				if(parm==NULL) return BadOptionErr;
				int excluded;
				sscanf(parm,"%d",&excluded);
				if(excluded<=0)
				{   cerr << "ExtractMPM option 'm' must be a positive integer" << endl;
					return 1;
				}
				excludeMaterial.push_back(excluded);
				break;
			}
			
			// select one material (must be M, space, number)
			else if(argv[parmInd][optInd]=='M')
			{	parm=NextArgument(++parmInd,argv,argc,'M');
				if(parm==NULL) return BadOptionErr;
				int included;
				sscanf(parm,"%d",&included);
				if(included<=0)
				{   cerr << "ExtractMPM option 'M' must be a positive integer" << endl;
					return BadOptionErr;
				}
				includeMaterial.push_back(included);
				break;
			}
			
			// input file format for material points
			else if(argv[parmInd][optInd]=='b')
			{	parm=NextArgument(++parmInd,argv,argc,'b');
				if(parm==NULL) return 1;
				if(strlen(parm)>ARCH_MAXMPMITEMS)
				{   cerr << "ExtractMPM option 'b' is too long" << endl;
					return BadOptionErr;
				}
				strcpy(mpmOrder,parm);
				break;
			}
			
			// input file format for material points
			else if(argv[parmInd][optInd]=='c')
			{	parm=NextArgument(++parmInd,argv,argc,'c');
				if(parm==NULL) return 1;
				if(strlen(parm)>ARCH_MAXCRACKITEMS)
				{   cerr << "ExtractMPM option 'c' is too long" << endl;
					return BadOptionErr;
				}
				strcpy(crackOrder,parm);
				break;
			}
			
			// output file name
			else if(argv[parmInd][optInd]=='o')
			{	parm=NextArgument(++parmInd,argv,argc,'o');
				if(parm==NULL) return 1;
				if(outfile!=NULL) delete [] outfile;
				outfile=new char[strlen(parm)+1];
				strcpy(outfile,parm);
				break;
			}
			
			// text format (the default)
			else if(argv[parmInd][optInd]=='s' || argv[parmInd][optInd]=='S')
			{	stepExtension=true;
			}
			
			// title for the header
			else if(argv[parmInd][optInd]=='n')
			{	parm=NextArgument(++parmInd,argv,argc,'o');
				if(parm==NULL) return 1;
				if(headerName!=NULL) delete [] headerName;
				headerName=new char[strlen(parm)+1];
				strcpy(headerName,parm);
				break;
			}
			
			// binary file format
			else if(argv[parmInd][optInd]=='f' || argv[parmInd][optInd]=='F' ||
						argv[parmInd][optInd]=='d' || argv[parmInd][optInd]=='D')
			{	fileFormat=argv[parmInd][optInd];
				strcpy(fileExtension,"data");
			}
			
			// text format (the default)
			else if(argv[parmInd][optInd]=='t' || argv[parmInd][optInd]=='T')
			{	fileFormat='T';
				strcpy(fileExtension,"txt");
			}
			
			// text format (the default)
			else if(argv[parmInd][optInd]=='x' || argv[parmInd][optInd]=='X')
			{	fileFormat='X';
				strcpy(fileExtension,"xml");
			}
			
			// vtk
			else if(argv[parmInd][optInd]=='v' || argv[parmInd][optInd]=='V')
			{	fileFormat='V';
				strcpy(fileExtension,"vtk");
			}
			
			// xyz file
			else if(argv[parmInd][optInd]=='z' || argv[parmInd][optInd]=='Z')
			{	fileFormat='Z';
				strcpy(fileExtension,"xyz");
			}
			
			// 3D file
			else if(argv[parmInd][optInd]=='3')
			{	threeD=true;
			}
			
			// 2D file (the default)
			else if(argv[parmInd][optInd]=='2')
			{	threeD=false;
			}
			
			// exporting just crack data
			else if(argv[parmInd][optInd]=='C')
			{	crackDataOnly=true;
			}
			
			// exporting just particle data
			else if(argv[parmInd][optInd]=='P')
			{	crackDataOnly=false;
			}
			
			// include a header
			else if(argv[parmInd][optInd]=='h')
			{	header=true;
			}
			
			else
			{   cerr << "Unknown ExtractMPM option '" << argv[parmInd][optInd] << "' was used\n";
				return BadOptionErr;
			}
		}
    }
	
    //  Last parameter must be the file name
    if(parmInd>=argc)
    {	Usage("Input file name is missing");
        return NoInputFileErr;
    }
	
    /* byte order this machine
		Intel chips use little endian in which int 1 will have 0x01 in first byte.
		G5, etc, use big endian in which int 1 will have 0x01 in last byte.
		When done, byte order will be set to 'm' if output file is big endian
			or to 'i' if output file is little endian.
	*/
    int test=1;
    char *testPtr=(char *)&test;			// point to first byte
	thisEndian = *testPtr==1 ? 'i' : 'm' ;
	
	// extract the data
	int fileNum;
    for(fileNum=parmInd;fileNum<argc;fileNum++)
	{	int result=ExtractMPMData(argv[fileNum],fileNum-parmInd,argc-parmInd-1);
		if(result!=noErr) return result;
	}
	
    // insert code here...
    return noErr;
}

// grab next argument or quit with error is not found
char *NextArgument(int parmInd,char * const argv[],int argc,char option)
{
	// error if not there
	if(parmInd>=argc)
	{   cerr << "ExtractMPM option '" << option << "' is missing its required argument.\n";
		return NULL;
	}
	return argv[parmInd];
}
	

// Explain usage of this program
void Usage(const char *msg)
{
	if(msg!=NULL)
		cout << "\nERROR: " << msg << endl;
		
	cout << "\nExtractMPM\n    version 1.0.0" << endl;
    cout << "\nUsage:\n"
        "    ExtractMPM [-options] <InputFile>\n\n"
        "This program reads the <InputFile>, extracts the desired data as\n"
        "selected in the options and writes to standard output or to a file\n"
        "selected in the options.\n\n"
        "Options:\n"
		"  Input file information (only used for 'ver3' files)\n"
		"    -2                 Input file has 2D data (the default)\n"
		"    -3                 Input file has 3D data\n"
		"    -b format          MPM archive file format\n"
        "    -c format          MPM archive file crack format\n"
		"  Select output data to extract\n"
		"    -q name            Add column with named quantity\n"
		"    -P                 Extract particle data only (the default)\n"
		"    -C                 Extract crack data only\n"
		"    -m num             Exclude this material\n"
		"    -M num             Include only this material\n"
		"    -h                 Write header in file (default: false)\n"
		"    -n name            Title or name to be in the header\n"
		"    -s                 Include step number from input file numeric\n"
		"                              extension in the output file name\n"
		"  Select output file name and format\n"
		"    -o path            Output file name with no extension\n"
		"		                       (quoted if name has spaces)\n"
		"    -T                 Output as tab-delimited text file (the default)\n"
		"    -V                 Output as VTK Legacy file (particle date only)\n"
		"    -X                 Output as XML file\n"
		"    -d                 Output as little Endian binary doubles file\n"
        "    -D                 Output as big Endian binary doubles file\n"
        "    -f                 Output as little Endian binary floats file\n"
        "    -F                 Output as big Endian binary floats file\n"
		"  Other options\n"
        "    -H (of -?)         Show this help and exit\n"
		"\n"
        "See http://osupdocs.forestry.oregonstate.edu/index.php/ExtractMPM for documentation.\n\n"
          <<  endl;
}

#pragma mark EXTRACTION CODE

// main entry to extract results from one file
int ExtractMPMData(const char *mpmFile,int fileIndex,int lastIndex)
{
	int i;
	
	// skip non archive (all archives end in .# where # is integer)
	i = (int)strlen(mpmFile)-1;
	while(i>0 && mpmFile[i]!='.')
	{	if(mpmFile[i]<'0' || mpmFile[i]>'9')
		{	cout << "File '" << mpmFile <<"' is not an MPM archive file" << endl;
			return noErr;
		}
		i--;
	}
	if(i==0)
	{	cout << "File '" << mpmFile <<"' is not an MPM archive file" << endl;
		return noErr;
	}
	
	// open the file
	if((fp=fopen(mpmFile,"r"))==NULL)
	{	cerr << "Input file '" << mpmFile << "' could not be opened" << endl;
		return FileAccessErr;
	}
	
	// get file length
	if(fseek(fp,0L,SEEK_END)!=0)
	{	cerr << "Input file '" << mpmFile << "' access error (fseek())" << endl;
		fclose(fp);
		return FileAccessErr;
	}
	fileLength=ftell(fp);
	rewind(fp);
	
	// read first  block of file into buffer
    blockSize = fileLength>MAX_READ_BLOCK ? MAX_READ_BLOCK : fileLength;
	buffer=(unsigned char *)malloc(blockSize);
	if(buffer==NULL)
	{	cerr << "Out of memory creating file reading buffer" << endl;
		fclose(fp);
		return MemoryErr;
	}
    
    // read first block
	if(fread(buffer,blockSize,1,fp)!=1)
	{	cerr << "Input file '" << mpmFile << "' access error (not all read)" << endl;
		fclose(fp);
		return FileAccessErr;
	}
    
    // set up point
    ap=buffer;
    seekOffset = 0;
    
	// look for file header
    ap+=3;							// jump to version number
	int vernum = *ap-'0';
	if(vernum>=4)
	{	// read file format info in header and to replace set values
		ap++;
		unsigned char mpmLength=*ap++;
		for(i=0;i<(int)mpmLength;i++) mpmOrder[i]=*ap++;
		mpmOrder[(int)mpmLength]=0;

		unsigned char crackLength=*ap++;
		for(i=0;i<(int)crackLength;i++) crackOrder[i]=*ap++;
		crackOrder[(int)crackLength]=0;
		
		unsigned char fileDimension=*ap++;
		threeD = (fileDimension=='3') ? true : false ;
		
		if(vernum>=5)
		{	unsigned char structured=*ap++;
			isStructured = (structured=='1') ? true : false ;
			
			float *timePtr=(float *)ap;
			if(thisEndian!=mpmOrder[0]) Reverse((char *)ap, sizeof(float));
			archiveTimeMs = *timePtr;
			ap+=sizeof(float);
		}
		
		ap=buffer;
		ap+=64;			// jump past header
	}
	else if(vernum==3)
		ap++;			// jump past header
	else
	{   cerr << "ExtractMPM does not support archive file versions before 'ver3'\n";
		return FileAccessErr;
	}
	
	// check required parameters
	if(strlen(mpmOrder)<1)
	{   cerr << "Input file format for particle data (option 'b') is required when reading 'ver3' archive files\n";
		return BadOptionErr;
	}
	
	if(crackDataOnly && strlen(crackOrder)<1)
	{   cerr << "Input file format for crack data (option 'c') is required when reading 'ver3' archive files\n";
		return BadOptionErr;
	}
    
	// if this machine differs from input file, then must reverse bytes
	reverseFromInput = thisEndian==mpmOrder[0] ? false : true ;

	// get record size
	if(CalcArchiveSize(vernum)!=noErr)
	{	cerr << "Input file format too old for this tool (missing current defaults) or to new (unknown record size)" << endl;
		return FileAccessErr;
	}
    
    // adjust block end pointer
    apNextBlock = fileLength>MAX_READ_BLOCK ? ap+MAX_READ_BLOCK-2*recSize : ap+MAX_READ_BLOCK+1;

	// output file
	ofstream fout;
	if(outfile!=NULL)
	{	char fname[256];
		
		// get stepNum string
		int ext=0,dot=(int)strlen(mpmFile)-1;
		while(dot>=0 && mpmFile[dot]!='.') dot--;
		if(dot>=0)
		{	stepNum[ext++]='_';
			dot++;
			while(mpmFile[dot]!=0)
			{	if(mpmFile[dot]<'0' || mpmFile[dot]>'9')
				{	// skip if not a number
					ext=0;
					break;
				}
				stepNum[ext++]=mpmFile[dot++];
			}
		}
		stepNum[ext]=0;
		
		if(lastIndex==0)
		{	if(stepExtension)
				sprintf(fname,"%s%s.%s",outfile,stepNum,fileExtension);
			else
				sprintf(fname,"%s.%s",outfile,fileExtension);
		}
		else
		{	if(stepExtension && strlen(stepNum)>0)
				sprintf(fname,"%s%s.%s",outfile,stepNum,fileExtension);
			else
				sprintf(fname,"%s-%d.%s",outfile,fileIndex,fileExtension);
		}
		fout.open(fname);
		if(!outfile)
		{	cerr << "Output file '" << fname << "' could not be created" << endl;
			return FileAccessErr;
		}
		
		cout << "Writing file '" << fname << "'" << endl;
	}
	
	// get the stream
	ostream os((outfile!=NULL) ? fout.rdbuf() : cout.rdbuf());
	
	// special case for VTK legacy files
	if(fileFormat=='V')
	{	if(crackDataOnly)
		{	cerr << "Exports to VTK Legacy file format cannot be for crack data" << endl;
			return FileAccessErr;
		}
		int vtkResult = VTKLegacy(os,mpmFile);
        fclose(fp);
		free(buffer);
		return vtkResult;
	}
	
	// special case for xyz files
	else if(fileFormat=='Z')
	{	if(crackDataOnly)
		{	cerr << "Exports to XYZ file format cannot be for crack data" << endl;
			return FileAccessErr;
		}
		int xyzResult = XYZExport(os,mpmFile);
		fclose(fp);
		free(buffer);
		return xyzResult;
	}
	
	// optional header
	if(header)
	{	char *headBuffer=new char[2500];
		char headLine[50];
		headBuffer[0]=0;
		
		if(fileFormat=='X')
			strcat(headBuffer,"<!-- ");
		if(headerName!=NULL)
		{	strcat(headBuffer,"Name ");
			strcat(headBuffer,headerName);
			strcat(headBuffer,"\n");
		}
		strcat(headBuffer,"Source ");
		strcat(headBuffer,mpmFile);
		strcat(headBuffer,"\n");
		if(archiveTimeMs>-0.5)
		{	sprintf(headLine,"Time %g ms",archiveTimeMs);
			strcat(headBuffer,headLine);
			strcat(headBuffer,"\n");
		}
		
		if(crackDataOnly)
		{	strcat(headBuffer,"Export Crack_Data\n");
		}
		else
		{	strcat(headBuffer,"Export Particle_Data\n");
			if(includeMaterial.size()>0)
			{	strcat(headBuffer,"Included_Materials");
				for(i=0;i<(int)includeMaterial.size();i++)
				{	sprintf(headLine," %d",includeMaterial[i]);
					strcat(headBuffer,headLine);
				}
				strcat(headBuffer,"\n");
			}
			if(excludeMaterial.size()>0)
			{	strcat(headBuffer,"Excluded_Materials");
				for(i=0;i<(int)excludeMaterial.size();i++)
				{	sprintf(headLine," %d",excludeMaterial[i]);
					strcat(headBuffer,headLine);
				}
				strcat(headBuffer,"\n");
			}
		}
		
		strcat(headBuffer,"Data ");
		if(crackDataOnly) strcat(headBuffer,"# ");
		strcat(headBuffer,"x y");
		if(threeD) strcat(headBuffer," z");
		if(quantity.size()>0)
		{	for(i=0;i<(int)quantity.size();i++)
			{	strcat(headBuffer," ");
				strcat(headBuffer,quantityName[i]);
			}
		}
		strcat(headBuffer,"\n");
		
		strcat(headBuffer,"Format ");
		switch(fileFormat)
		{	case 'T':
				strcat(headBuffer,"text\n");
				break;
			case 'X':
				strcat(headBuffer,"xml\n");
				break;
			case 'd':
				strcat(headBuffer,"double\n");
				strcat(headBuffer,"Endian little\n");
				break;
			case 'D':
				strcat(headBuffer,"double\n");
				strcat(headBuffer,"Endian big\n");
				break;
			case 'f':
				strcat(headBuffer,"float\n");
				strcat(headBuffer,"Endian little\n");
				break;
			case 'F':
				strcat(headBuffer,"float\n");
				strcat(headBuffer,"Endian big\n");
				break;
			default:
				break;
		}
		strcat(headBuffer,"EndHeader");
		if(fileFormat=='X')
			strcat(headBuffer," -->");
		
		// pad for 8-byte alignment if binary
		if(fileFormat!='T' && fileFormat!='X')
		{	int pad= 7- (strlen(headBuffer) % 8);
			int i;
			for(i=0;i<pad;i++) strcat(headBuffer," ");
		}
		strcat(headBuffer,"\n");
		
		// write it
		os.write(headBuffer,strlen(headBuffer));
		delete [] headBuffer;
	}
	
	// the file
	int nummpms=(int)(fileLength/recSize);
	int p;
	short *mptr;
	BeginMP(os);
	for(p=0;p<nummpms;p++)
    {   // read next block when needed
        if(!GetNextFileBlock(mpmFile)) return FileAccessErr;
        
		short matnum=pointMatnum(ap);
		if(matnum<0) break;
		
		// skip if crack data export
		if(crackDataOnly)
		{	ap+=recSize;
			continue;
		}
		
		// if skip is true here, that skip this material
		if(skipThisPoint(matnum))
		{	ap+=recSize;
			continue;
		}
		
		// write position, special for XML which starts the <mp> element
		if(fileFormat!='X')
		{	OutputDouble((double *)(ap+posOffset),0,0,reverseFromInput,os,XPOS);
		    OutputDouble((double *)(ap+posOffset),1,'\t',reverseFromInput,os,YPOS);
		    if(threeD) OutputDouble((double *)(ap+posOffset),2,'\t',reverseFromInput,os,ZPOS);
		}
		else
		{	// begin mp element (mat, angle, thickness always output)
			os << "  <mp matl='" << matnum << "'";
			if(!threeD)
			{	double *angle=(double *)(ap+angleOffset);
				double *thickness=(double *)(ap+angleOffset+sizeof(double));
				if(reverseFromInput)
				{	Reverse((char *)angle, sizeof(double));
					Reverse((char *)thickness, sizeof(double));
				}
				os << " angle='" << *angle << "' thick='" << *thickness << "'";
			}
			if(tempOffset>0 && hasTempQ)
			{	double *temp=(double *)(ap+tempOffset);
				if(reverseFromInput) Reverse((char *)temp,sizeof(double));
				os << " temp='" << *temp << '"';
			}
			if(concOffset>0 && hasConcQ)
			{	double *conc=(double *)(ap+concOffset);
				if(reverseFromInput) Reverse((char *)conc,sizeof(double));
				os << " wtconc='" << *conc << '"';
			}
			// future: add conc and temp attributes if requested
			os << ">" << endl;
			
			// subordinate mp element for position
			double *xpos=(double *)(ap+posOffset);
			double *ypos=xpos+1;
			if(reverseFromInput)
			{	Reverse((char *)xpos, sizeof(double));
				Reverse((char *)ypos, sizeof(double));
			}
			os << "    <pt x='" << *xpos << "' y='" << *ypos << "'";
			if(threeD)
			{	double *zpos=xpos+2;
				if(reverseFromInput)  Reverse((char *)zpos, sizeof(double));
			    os << " z='" << *zpos << "'";
            }
			os << "/>" << endl;
		}
		
		// write quantities
		for(i=0;i<(int)quantity.size();i++)
			OutputQuantity(i,ap,os,matnum,'\t');
		
		OutputRecordEnd(os,true);
		ap+=recSize;
	}
	EndMP(os);
	
	// is this a crack data Export?
	if(crackDataOnly)
	{	int initCrack=p;
		int crackNumber=0,tipMatnum;
		for(p=initCrack;p<nummpms;p++)
        {   // read next block when needed
            if(!GetNextFileBlock(mpmFile)) return FileAccessErr;
            
			// read tipMatNum
			mptr=(short *)(ap+sizeof(int));
			if(reverseFromInput) Reverse((char *)mptr,sizeof(int));
			tipMatnum=*(int *)mptr;
			
			// read marker (-1 for new crack, -2 for same crack
			mptr=(short *)(ap+sizeof(int)+sizeof(double));
			if(reverseFromInput) Reverse((char *)mptr,sizeof(short));
			short matnum=*mptr;
			if(matnum==-1)
			{	if(p>initCrack) EndCrack(os);
				BeginCrack(os);
				crackNumber++;
			}
			
			// write number and position
			if(fileFormat!='X')
			{	double matDouble=(double)crackNumber;
				OutputDouble(&matDouble,0,0,false,os,-1);
				OutputDouble((double *)(ap+crackPosOffset),0,'\t',reverseFromInput,os,-1);
				OutputDouble((double *)(ap+crackPosOffset),1,'\t',reverseFromInput,os,-1);
				if(threeD) OutputDouble((double *)(ap+crackPosOffset),2,'\t',reverseFromInput,os,-1);
			
				// write quantities
				for(i=0;i<(int)quantity.size();i++)
				{	switch(quantity[i])
					{	case J1:
						case J2:
							if(jIntOffset>0)
								OutputDouble((double *)(ap+jIntOffset),quantity[i]-J1,'\t',reverseFromInput,os,quantity[i]);
							else
								OutputDouble(&zeroDouble,0,'\t',false,os,quantity[i]);
							break;
							
						case KI:
						case KII:
							if(kSifOffset>0)
								OutputDouble((double *)(ap+kSifOffset),quantity[i]-KI,'\t',reverseFromInput,os,quantity[i]);
							else
								OutputDouble(&zeroDouble,0,'\t',false,os,quantity[i]);
							break;
						
						default:
							OutputDouble(&zeroDouble,0,'\t',false,os,quantity[i]);
							break;
					}
				}
			}
			else
			{	double *xpos=(double *)(ap+crackPosOffset);
				double *ypos=xpos+1;
				if(reverseFromInput)
				{	Reverse((char *)xpos,sizeof(double));
					Reverse((char *)ypos,sizeof(double));
				}
				os << "  <pt x='" << *xpos << "' y='" << *ypos;
				if(tipMatnum==-2 || (tipMatnum>0 && tipMatnum<25))
					os << "' tip='" << tipMatnum;
				os << "'/>" << endl;
			}
			
			OutputRecordEnd(os,false);
			ap+=recSize;
		}
		EndCrack(os);
	}
	
	free(buffer);
    fclose(fp);
	return noErr;
}

// called when start MP output - only needed for XML output
int VTKLegacy(ostream &os,const char *mpmFile)
{
	os << "# vtk DataFile Version 4.0" << endl;
	if(headerName!=NULL)
		os << headerName;
	else
		os << "OSParticulas/NairnMPM particle data";
	os << " from file ";
	os << mpmFile ;
	if(archiveTimeMs>-0.5)
		os << " at time " << archiveTimeMs << " ms" ;
	os << endl;
	os << "ASCII" << endl;
	os << "DATASET POLYDATA" << endl;
	
	// time and step number
	if(strlen(stepNum)>1)
	{	os << "FIELD FieldData 2" << endl;
		os << "TIME 1 1 double" << endl;
		os << archiveTimeMs << endl;
		os << "STEP 1 1 int" << endl;
		os << &stepNum[1] << endl;
	}
	else
	{	os << "FIELD FieldData 1" << endl;
		os << "TIME 1 1 double" << endl;
		os << archiveTimeMs << endl;
	}

	// count points to extract
	int nummpms=(int)(fileLength/recSize);
	int p,numExtract=0;
    long origOffset = (long)(ap-buffer);
	for(p=0;p<nummpms;p++)
    {   // read next block when needed
        if(!GetNextFileBlock(mpmFile)) return FileAccessErr;
        
		short matnum=pointMatnum(ap);
		if(matnum<0) break;
		if(!skipThisPoint(matnum)) numExtract++;
		ap+=recSize;
	}
	
    // output the number of points
	os << "POINTS " << numExtract << " double" << endl;
	if(numExtract==0) return noErr;
    
    // back to start of the file
    if(!RestartFileBlocks(origOffset,mpmFile)) return FileAccessErr;
	
	// extract the point positions
	for(p=0;p<nummpms;p++)
    {   // read next block when needed
        if(!GetNextFileBlock(mpmFile)) return FileAccessErr;
        
		short matnum=pointMatnum(ap);
		if(matnum<0) break;
		if(!skipThisPoint(matnum))
		{	OutputDouble((double *)(ap+posOffset),0,0,reverseFromInput,os,XPOS);
		    OutputDouble((double *)(ap+posOffset),1,' ',reverseFromInput,os,YPOS);
		    if(threeD)
                OutputDouble((double *)(ap+posOffset),2,' ',reverseFromInput,os,ZPOS);
            else
                os << " 0";         // so ParaView can do 2D data
			os << endl;
		}
		ap+=recSize;
	}
	
	// extract quantities
	if(quantity.size()==0) return noErr;
	os << "POINT_DATA " << numExtract << endl;
	int i;
	for(i=0;i<(int)quantity.size();i++)
    {   // back to start of the file
        if(!RestartFileBlocks(origOffset,mpmFile)) return FileAccessErr;
		
		switch(quantity[i])
		{	case VELOCITY:
			case DISPLACEMENT:
				os << "VECTORS " << quantityName[i] << " double" << endl;
				break;
			case STRESS:
			case STRAIN:
			case PLASTICSTRAIN:
				os << "TENSORS " << quantityName[i] << " double" << endl;
				break;
			default:
				os << "SCALARS " << quantityName[i] << " double 1" << endl;
				os << "LOOKUP_TABLE default" << endl;
				break;
		}
		
		for(p=0;p<nummpms;p++)
        {   // read next block when needed
            if(!GetNextFileBlock(mpmFile)) return FileAccessErr;
            
			short matnum=pointMatnum(ap);
			if(matnum<0) break;
			if(!skipThisPoint(matnum))
			{	OutputQuantity(i,ap,os,matnum,0);
				os << endl;
			}
			ap+=recSize;
		}
	}
    
    // done
    return noErr;
}

// called when start MP output - only needed for XML output
int XYZExport(ostream &os,const char *mpmFile)
{
	// count points to extract
	int nummpms=(int)(fileLength/recSize);
	int p,numExtract=0;
	long origOffset = (long)(ap-buffer);
	for(p=0;p<nummpms;p++)
	{   // read next block when needed
		if(!GetNextFileBlock(mpmFile)) return FileAccessErr;
		
		short matnum=pointMatnum(ap);
		if(matnum<0) break;
		if(!skipThisPoint(matnum)) numExtract++;
		ap+=recSize;
	}
	
	// number of points
	//os << numExtract << endl;

	// one line comment
	//os << "Particle data from file " << mpmFile;
	//if(archiveTimeMs>-0.5)
	//	os << " at time " << archiveTimeMs << " ms" ;
	//os << endl;
	
	// back to start of the file
	if(!RestartFileBlocks(origOffset,mpmFile)) return FileAccessErr;
	
	// extract material and point positions
	for(p=0;p<nummpms;p++)
	{   // read next block when needed
		if(!GetNextFileBlock(mpmFile)) return FileAccessErr;
		
		short matnum=pointMatnum(ap);
		if(matnum<0) break;
		if(!skipThisPoint(matnum))
		{	// material number
			//os << matnum;
			
			// data
			if(!skipThisPoint(matnum))
			{	OutputDouble((double *)(ap+posOffset),0,0,reverseFromInput,os,XPOS);
				//OutputDouble((double *)(ap+posOffset),0,' ',reverseFromInput,os,XPOS);
				OutputDouble((double *)(ap+posOffset),1,' ',reverseFromInput,os,YPOS);
				if(threeD)
					OutputDouble((double *)(ap+posOffset),2,' ',reverseFromInput,os,ZPOS);
				else
					os << " 0.0";         // so ParaView can do 2D data
				os << endl;
			}
		}
		/*
		char pos[20];
		if(!skipThisPoint(matnum))
		{	// output material x y z
			if(matnum<10)
				os << "  ";
			else if(matnum<100)
				os << " ";
			os << matnum;
			
			// X pos
			double *data = (double *)(ap+posOffset);
			if(reverseFromInput) Reverse((char *)data,sizeof(double));
			sprintf(pos," %12.5f",*data);
			os << pos;
			
			// Y pos
			data+=1;
			if(reverseFromInput) Reverse((char *)data,sizeof(double));
			sprintf(pos," %12.5f",*data);
			os << pos;
			
			// Z pos
			if(threeD)
			{	data+=1;
				if(reverseFromInput) Reverse((char *)data,sizeof(double));
				sprintf(pos," %12.5f",*data);
			}
			else
				sprintf(pos," %12.5f",(double)0.);
			os << pos;
			
			// line odne
			os << endl;
		}
		*/
		ap+=recSize;
	}
	
	// done
	return noErr;
}

// If needed, read next block of data from the file
bool GetNextFileBlock(const char *mpmFile)
{
    if(ap<=apNextBlock) return true;
    
    // set new offset
    seekOffset += (long)(ap-buffer);
    if(fseek(fp,seekOffset,SEEK_SET)!=0)
    {	cerr << "Input file '" << mpmFile << "' access error (not able to reset file position)" << endl;
        fclose(fp);
        return false;
    }
        
    // read next block
    blockSize = seekOffset+MAX_READ_BLOCK > fileLength ? fileLength-seekOffset : MAX_READ_BLOCK;
    if(fread(buffer,blockSize,1,fp)!=1)
    {	cerr << "Input file '" << mpmFile << "' access error (next block not read)" << endl;
        fclose(fp);
        return false;
    }
        
    // reset pointers
    ap=buffer;
    apNextBlock = fileLength>seekOffset+MAX_READ_BLOCK ? ap+MAX_READ_BLOCK-2*recSize : ap+MAX_READ_BLOCK+1;
    return true;
}

// When reading VTK file, need several passes through file and restart before each one
bool RestartFileBlocks(long origOffset,const char *mpmFile)
{
    if(fileLength>MAX_READ_BLOCK)
    {   // go back and read first block
        rewind(fp);
        blockSize = MAX_READ_BLOCK;
        if(fread(buffer,blockSize,1,fp)!=1)
        {	cerr << "Input file '" << mpmFile << "' access error (not all read)" << endl;
            fclose(fp);
            return false;
        }
        ap=buffer;
        apNextBlock = ap+MAX_READ_BLOCK-2*recSize;
        seekOffset = 0;
    }
    
    // set to point offset
    ap = buffer+origOffset;
    return true;
}

// output one quantity to the selected file type
void OutputQuantity(int i,unsigned char *ap,ostream &os,short matnum,char delim)
{
	switch(quantity[i])
	{	case SXX:
		case SYY:
		case SZZ:
		case SXY:
		case SXZ:
		case SYZ:
            // Stress in Pa
			if(stressOffset>0 && (threeD || quantity[i]<SXZ))
				OutputDouble((double *)(ap+stressOffset),quantity[i]-SXX,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
		
		case STRESS:
			if(fileFormat=='V')
			{	// VTK legacy only
				if(stressOffset>0)
				{	// 9 elements of stress - SXX, SXY, SXZ
					OutputDouble((double *)(ap+stressOffset),0,0,reverseFromInput,os,quantity[i]);
					OutputDouble((double *)(ap+stressOffset),3,' ',reverseFromInput,os,quantity[i]);
					if(threeD)
						OutputDouble((double *)(ap+stressOffset),4,' ',reverseFromInput,os,quantity[i]);
					else
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
					
					// SXY, SYY, SYZ
					OutputDouble((double *)(ap+stressOffset),3,' ',reverseFromInput,os,quantity[i]);
					OutputDouble((double *)(ap+stressOffset),1,' ',reverseFromInput,os,quantity[i]);
					if(threeD)
						OutputDouble((double *)(ap+stressOffset),5,' ',reverseFromInput,os,quantity[i]);
					else
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
					
					// SXZ, SYZ, SZZ
					if(threeD)
					{	OutputDouble((double *)(ap+stressOffset),4,' ',reverseFromInput,os,quantity[i]);
						OutputDouble((double *)(ap+stressOffset),5,' ',reverseFromInput,os,quantity[i]);
					}
					else
					{	OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
					}
					OutputDouble((double *)(ap+stressOffset),2,' ',reverseFromInput,os,quantity[i]);
				}
				else
				{	// 9 zeros
					OutputDouble(&zeroDouble,0,0,false,os,quantity[i]);
					for(int z=0;z<8;z++)
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);

				}
			}
			else if(stressOffset>0)
			{	// each available element of stress SXX, SYY, SZZ, SXY, SXZ, SYZ
				OutputDouble((double *)(ap+stressOffset),0,delim,reverseFromInput,os,quantity[i]);
				OutputDouble((double *)(ap+stressOffset),1,delim,reverseFromInput,os,quantity[i]);
				OutputDouble((double *)(ap+stressOffset),2,delim,reverseFromInput,os,quantity[i]);
				OutputDouble((double *)(ap+stressOffset),3,delim,reverseFromInput,os,quantity[i]);
				if(threeD)
				{	OutputDouble((double *)(ap+stressOffset),4,delim,reverseFromInput,os,quantity[i]);
					OutputDouble((double *)(ap+stressOffset),5,delim,reverseFromInput,os,quantity[i]);
				}
				else
				{	OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
					OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
				}
			}
			else
			{	// 6 zeros
				for(int z=0;z<6;z++)
					OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			}
			break;

        case PRESSURE:
            // Pressure in Pa
			if(stressOffset>0)
            {   double *sxx=(double *)(ap+stressOffset);
                double *syy=(sxx+SYY-SXX);
                double *szz=(sxx+SZZ-SXX);
				if(reverseFromInput)
				{	Reverse((char *)sxx,sizeof(double));
					Reverse((char *)syy,sizeof(double));
					Reverse((char *)szz,sizeof(double));
				}
                double pressure = -(*sxx+*syy+*szz)/3.;
				OutputDouble(&pressure,0,delim,false,os,quantity[i]);
            }
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
            
        case EQUIVSTRESS:
            // Equivalent stress = sqrt(3 J2) in Pa (aka vonmises stress)
            if(stressOffset>0)
            {   double *sxx=(double *)(ap+stressOffset);
                double *syy=(sxx+SYY-SXX);
                double *szz=(sxx+SZZ-SXX);
                double *sxy=(sxx+SXY-SXX);
				double *sxz = (sxx + SXZ - SXX);
				double *syz = (sxx + SYZ - SXX);
                if(reverseFromInput)
                {	Reverse((char *)sxx,sizeof(double));
                    Reverse((char *)syy,sizeof(double));
                    Reverse((char *)szz,sizeof(double));
                    Reverse((char *)sxy,sizeof(double));
					if(threeD)
					{	Reverse((char *)sxz, sizeof(double));
						Reverse((char *)syz, sizeof(double));
					}
				}
                double se = pow(*sxx-*syy,2.) + pow(*syy-*szz,2.) + pow(*sxx-*szz,2.);
                se += 6.*(*sxy)*(*sxy);
                if(threeD) se += 6.*((*sxz)*(*sxz) + (*syz)*(*syz));
                se = sqrt(0.5*se);
                OutputDouble(&se,0,delim,false,os,quantity[i]);
            }
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;

        case EQUIVSTRAIN:
            // Equivalent strain = sqrt(2 s.s / 3) where s is  deviatoric strain (absolute)
            if(strainOffset>0)
            {   double *exx=(double *)(ap+strainOffset);
                double *eyy=(exx+EYY-EXX);
                double *ezz=(exx+EZZ-EXX);
                double *exy=(exx+EXY-EXX);
                if(reverseFromInput)
                {	Reverse((char *)exx,sizeof(double));
                    Reverse((char *)eyy,sizeof(double));
                    Reverse((char *)ezz,sizeof(double));
                    Reverse((char *)exy,sizeof(double));
                }
                double dxz=0.,dyz=0.;
                if(threeD)
                {   double *exz=(exx+EXZ-EXX);
                    double *eyz=(exx+EYZ-EXX);
                    if(reverseFromInput)
                    {	Reverse((char *)exz,sizeof(double));
                        Reverse((char *)eyz,sizeof(double));
                    }
                    dxz = *exz;
                    dyz = *eyz;
                }
                double tre = (*exx+*eyy+*ezz)/3.;
                double dxx = *exx - tre;
                double dyy = *eyy - tre;
                double dzz = *ezz - tre;
                double dxy = *exy;
                double se = dxx*dxx + dyy*dyy + dzz*dzz + 0.5*(dxy*dxy + dxz*dxz + dyz*dyz);
                se = sqrt(2.*se/3.);
                OutputDouble(&se,0,delim,false,os,quantity[i]);
            }
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
		
		case EXX:
		case EYY:
		case EZZ:
		case EXY:
		case EXZ:
		case EYZ:
			if(strainOffset>0 && (threeD || quantity[i]<EXZ))
				OutputDouble((double *)(ap+strainOffset),quantity[i]-EXX,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
		case STRAIN:
			if(fileFormat=='V')
			{	// VTK legacy only
				if(strainOffset>0)
				{	// 9 elements of stress - EXX, EXY, EXZ
					OutputDouble((double *)(ap+strainOffset),0,0,reverseFromInput,os,quantity[i]);
					OutputDouble((double *)(ap+strainOffset),3,' ',reverseFromInput,os,quantity[i]);
					if(threeD)
						OutputDouble((double *)(ap+strainOffset),4,' ',reverseFromInput,os,quantity[i]);
					else
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
					
					// EXY, EYY, EXZ
					OutputDouble((double *)(ap+strainOffset),3,' ',reverseFromInput,os,quantity[i]);
					OutputDouble((double *)(ap+strainOffset),1,' ',reverseFromInput,os,quantity[i]);
					if(threeD)
						OutputDouble((double *)(ap+strainOffset),5,' ',reverseFromInput,os,quantity[i]);
					else
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
					
					// EXZ, EYZ, EZZ
					if(threeD)
					{	OutputDouble((double *)(ap+strainOffset),4,' ',reverseFromInput,os,quantity[i]);
						OutputDouble((double *)(ap+strainOffset),5,' ',reverseFromInput,os,quantity[i]);
					}
					else
					{	OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
					}
					OutputDouble((double *)(ap+strainOffset),2,' ',reverseFromInput,os,quantity[i]);
					break;
				}
				else
				{	// 9 zeros
					OutputDouble(&zeroDouble,0,0,false,os,quantity[i]);
					for(int z=0;z<8;z++)
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
					
				}
			}
			else if(strainOffset>0)
			{	// each available element of stress EXX, EYY, EZZ, EXY, EXZ, EYZ
				OutputDouble((double *)(ap+strainOffset),0,delim,reverseFromInput,os,quantity[i]);
				OutputDouble((double *)(ap+strainOffset),1,delim,reverseFromInput,os,quantity[i]);
				OutputDouble((double *)(ap+strainOffset),2,delim,reverseFromInput,os,quantity[i]);
				OutputDouble((double *)(ap+strainOffset),3,delim,reverseFromInput,os,quantity[i]);
				if(threeD)
				{	OutputDouble((double *)(ap+strainOffset),4,delim,reverseFromInput,os,quantity[i]);
					OutputDouble((double *)(ap+strainOffset),5,delim,reverseFromInput,os,quantity[i]);
				}
				else
				{	OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
					OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
				}
			}
			else
			{	// 6 zeros
				for(int z=0;z<6;z++)
					OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			}
			break;

		case PEXX:
		case PEYY:
		case PEZZ:
		case PEXY:
		case PEXZ:
		case PEYZ:
			if(plStrainOffset>0 && (threeD || quantity[i]<PEXZ))
				OutputDouble((double *)(ap+plStrainOffset),quantity[i]-PEXX,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
		
		case PLASTICSTRAIN:
			if(fileFormat=='V')
			{	// VTK legacy only
				if(plStrainOffset>0)
				{	// 9 elements of stress - EPXX, EPXY, EPXZ
					OutputDouble((double *)(ap+plStrainOffset),0,0,reverseFromInput,os,quantity[i]);
					OutputDouble((double *)(ap+plStrainOffset),3,' ',reverseFromInput,os,quantity[i]);
					if(threeD)
						OutputDouble((double *)(ap+plStrainOffset),4,' ',reverseFromInput,os,quantity[i]);
					else
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
					
					// EXY, EYY, EXZ
					OutputDouble((double *)(ap+plStrainOffset),3,' ',reverseFromInput,os,quantity[i]);
					OutputDouble((double *)(ap+plStrainOffset),1,' ',reverseFromInput,os,quantity[i]);
					if(threeD)
						OutputDouble((double *)(ap+plStrainOffset),5,' ',reverseFromInput,os,quantity[i]);
					else
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
					
					// EPXZ, EPYZ, EPXZ
					if(threeD)
					{	OutputDouble((double *)(ap+plStrainOffset),4,' ',reverseFromInput,os,quantity[i]);
						OutputDouble((double *)(ap+plStrainOffset),5,' ',reverseFromInput,os,quantity[i]);
					}
					else
					{	OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
					}
					OutputDouble((double *)(ap+plStrainOffset),2,' ',reverseFromInput,os,quantity[i]);
				}
				else
				{	// 9 zeros
					OutputDouble(&zeroDouble,0,0,false,os,quantity[i]);
					for(int z=0;z<8;z++)
						OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
					
				}
			}
			else if(plStrainOffset>0)
			{	// each available element of stress EPXX, EPYY, EPZZ, EPXY, EPXZ, EPYZ
				OutputDouble((double *)(ap+plStrainOffset),0,delim,reverseFromInput,os,quantity[i]);
				OutputDouble((double *)(ap+plStrainOffset),1,delim,reverseFromInput,os,quantity[i]);
				OutputDouble((double *)(ap+plStrainOffset),2,delim,reverseFromInput,os,quantity[i]);
				OutputDouble((double *)(ap+plStrainOffset),3,delim,reverseFromInput,os,quantity[i]);
				if(threeD)
				{	OutputDouble((double *)(ap+plStrainOffset),4,delim,reverseFromInput,os,quantity[i]);
					OutputDouble((double *)(ap+plStrainOffset),5,delim,reverseFromInput,os,quantity[i]);
				}
				else
				{	OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
					OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
				}
			}
			else
			{	// 6 zeros
				for(int z=0;z<6;z++)
					OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			}
			break;

		case VELX:
		case VELY:
		case VELZ:
			if(velocityOffset>0 && (threeD || quantity[i]<VELZ))
				OutputDouble((double *)(ap+velocityOffset),quantity[i]-VELX,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
		
		case VELOCITY:
			if(fileFormat=='V')
			{	// VTK legacy only
				OutputDouble((double *)(ap+velocityOffset),0,0,reverseFromInput,os,quantity[i]);
				OutputDouble((double *)(ap+velocityOffset),1,' ',reverseFromInput,os,quantity[i]);
				if(threeD)
					OutputDouble((double *)(ap+velocityOffset),2,' ',reverseFromInput,os,quantity[i]);
				else
					OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
			}
			else
			{	OutputDouble((double *)(ap+velocityOffset),0,delim,reverseFromInput,os,quantity[i]);
				OutputDouble((double *)(ap+velocityOffset),1,delim,reverseFromInput,os,quantity[i]);
				if(threeD)
					OutputDouble((double *)(ap+velocityOffset),2,delim,reverseFromInput,os,quantity[i]);
				else
					OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			}
			break;
			
		case DISPX:
		case DISPY:
		case DISPZ:
			if(threeD || quantity[i]<DISPZ)
			{	double *pPtr=(double *)(ap+posOffset);
				pPtr+=(quantity[i]-DISPX);
				double *opPtr=(double *)(ap+origPosOffset);
				opPtr+=(quantity[i]-DISPX);
				if(reverseFromInput)
				{	Reverse((char *)pPtr,sizeof(double));
					Reverse((char *)opPtr,sizeof(double));
				}
				double displacement=*pPtr-*opPtr;
				OutputDouble(&displacement,0,delim,false,os,quantity[i]);
			}
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
		case DISPLACEMENT:
		{	double *pPtr=(double *)(ap+posOffset);
			double *opPtr=(double *)(ap+origPosOffset);
			if(reverseFromInput)
			{	Reverse((char *)pPtr,sizeof(double));
				Reverse((char *)opPtr,sizeof(double));
			}
			double displacement=*pPtr-*opPtr;
			if(fileFormat=='V')
				OutputDouble(&displacement,0,0,false,os,quantity[i]);
			else
				OutputDouble(&displacement,0,delim,false,os,quantity[i]);
			
			pPtr++;
			opPtr++;
			if(reverseFromInput)
			{	Reverse((char *)pPtr,sizeof(double));
				Reverse((char *)opPtr,sizeof(double));
			}
			displacement=*pPtr-*opPtr;
			if(fileFormat=='V')
				OutputDouble(&displacement,0,' ',false,os,quantity[i]);
			else
				OutputDouble(&displacement,0,delim,false,os,quantity[i]);

			if(threeD)
			{	pPtr++;
				opPtr++;
				if(reverseFromInput)
				{	Reverse((char *)pPtr,sizeof(double));
					Reverse((char *)opPtr,sizeof(double));
				}
				displacement=*pPtr-*opPtr;
				if(fileFormat=='V')
					OutputDouble(&displacement,0,' ',false,os,quantity[i]);
				else
					OutputDouble(&displacement,0,delim,false,os,quantity[i]);
			}
			else if(fileFormat=='V')
				OutputDouble(&zeroDouble,0,' ',false,os,quantity[i]);
			else
				OutputDouble(&displacement,0,delim,false,os,quantity[i]);
			break;
		}
			
		case STRENG:
			if(strainEnergyOffset>0)
				OutputDouble((double *)(ap+strainEnergyOffset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
		case PLENG:
			if(plasticEnergyOffset>0)
				OutputDouble((double *)(ap+plasticEnergyOffset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
		case TEMP:
			if(tempOffset>0)
				OutputDouble((double *)(ap+tempOffset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
		case CONC:
			if(concOffset>0)
				OutputDouble((double *)(ap+concOffset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
		
		case MAT:
		{	double matDouble=(double)matnum;
			OutputDouble(&matDouble,0,delim,false,os,quantity[i]);
			break;
		}
		
		case MASS:
			OutputDouble((double *)(ap+sizeof(int)),0,delim,reverseFromInput,os,quantity[i]);
			break;
		
		case MATANGLEZ:
			OutputDouble((double *)(ap+angleOffset),0,delim,reverseFromInput,os,quantity[i]);
			break;
			
		case MATANGLEY:
			if(angleYOffset>0)
				OutputDouble((double *)(ap+angleYOffset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
		case MATANGLEX:
			if(angleXOffset>0)
				OutputDouble((double *)(ap+angleXOffset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
        case HIST1:
            if(history1Offset>0)
                OutputDouble((double *)(ap+history1Offset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
        case HIST2:
            if(history2Offset>0)
                OutputDouble((double *)(ap+history2Offset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
        case HIST3:
            if(history3Offset>0)
                OutputDouble((double *)(ap+history3Offset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
        case HIST4:
            if(history4Offset>0)
                OutputDouble((double *)(ap+history4Offset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
        
        case WORKENERGY:
            if(workEnergyOffset>0)
                OutputDouble((double *)(ap+workEnergyOffset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
        case HEATENERGY:
            if(heatEnergyOffset>0)
                OutputDouble((double *)(ap+heatEnergyOffset),0,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
			
		default:
			OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
			break;
	}
}

// get material number for next material point or crack point
short pointMatnum(unsigned char *ap)
{	short *mptr=(short *)(ap+sizeof(int)+sizeof(double));
	if(reverseFromInput) Reverse((char *)mptr,sizeof(short));
	short matnum=*mptr;
	if(matnum<0)
	{	if(reverseFromInput) Reverse((char *)mptr,sizeof(short));
	}
	return matnum;
}

// return true to skip this point or false to include it
bool skipThisPoint(short matnum)
{
	// is it included?
	int i;
	bool skip=false;
	if(includeMaterial.size()>0)
	{	// if have list of included, make sure on the list
		skip=true;
		for(i=0;i<(int)includeMaterial.size();i++)
		{	if(includeMaterial[i]==matnum)
			{	skip=false;
				break;
			}
		}
	}
	// skip is false here if no list of included materials or if is in that list
		
	// can override with exclusion
	if(!skip)
	{	for(i=0;i<(int)excludeMaterial.size();i++)
		{	if(excludeMaterial[i]==matnum)
			{	skip=true;
				break;
			}
		}
	}
	
	return skip;
}

/* output double pointed to by *(data+offset)
	delim 0 for first column or 1 char delimiter for other columns
	mustReverse==reverseFromInput if pointer is to raw data
		or false if pointer is a local double
*/
void OutputDouble(double *data,int offset,char delim,bool mustReverse,ostream &os,int q)
{	// offset to component in vector or tensor
	data+=offset;
	
	// output in desired format
	float fdata;
	switch(fileFormat)
	{	case 'T':
		case 'V':
		case 'Z':
			if(mustReverse) Reverse((char *)data,sizeof(double));
			if(delim!=0) os << delim;
			os << *data;
			break;
		
		case 'd':
			if((thisEndian=='m' && !mustReverse) || (thisEndian=='i' && mustReverse))
				Reverse((char *)data,sizeof(double));
			os.write((char *)data,sizeof(double));
			break;
			
		case 'D':
			if((thisEndian=='i' && !mustReverse) || (thisEndian=='m' && mustReverse))
				Reverse((char *)data,sizeof(double));
			os.write((char *)data,sizeof(double));
			break;
		
		case 'f':
			if(mustReverse) Reverse((char *)data,sizeof(double));
			fdata=(float)(*data);
			if(thisEndian=='m') Reverse((char *)&fdata,sizeof(float));
			os.write((char *)&fdata,sizeof(float));
			break;
			
		case 'F':
			if(mustReverse) Reverse((char *)data,sizeof(double));
			fdata=(float)(*data);
			if(thisEndian=='i' ) Reverse((char *)&fdata,sizeof(float));
			os.write((char *)&fdata,sizeof(float));
			break;
			
		case 'X':
			if(mustReverse) Reverse((char *)data,sizeof(double));
			switch(q)
			{	case MASS:
					os << "    <mass m='" << *data << "'/>" << endl;
					break;
				
				case VELX:
					os << "    <vel x='" << *data << "'";
					data++;
					if(mustReverse) Reverse((char *)data,sizeof(double));
					os << "    y='" << *data << "'";
					if(threeD)
					{	data++;
						if(mustReverse) Reverse((char *)data,sizeof(double));
						os << "    z='" << *data << "'";
					}
					os << "/>" << endl;
					break;
				
				default:
					break;
			}
			break;
		
		default:
			break;
	}
}


// called when record is done. mpData true for end of material point record of false for crack record
void OutputRecordEnd(ostream &os, bool mpData)
{
	switch(fileFormat)
	{	case 'T':
			os << endl;
			break;
		
		case 'X':
			if(mpData)
				os << "  </mp>" << endl;
			break;
			
		default:
			break;
	}
}

// called when start MP output - only needed for XML output
void BeginMP(ostream &os)
{
	switch(fileFormat)
	{	case 'X':
	        if(!crackDataOnly)
			{	os << "<PointList>" << endl;
			}
			break;
		
		default:
			break;
	}
}

// called when end MP output - only needed for XML output
void EndMP(ostream &os)
{
	switch(fileFormat)
	{	case 'X':
	        if(!crackDataOnly)
			{	os << "</PointList>" << endl;
            }
			break;
		
		default:
			break;
	}
}


// called when start new crack - only needed for XML output
void BeginCrack(ostream &os)
{
	switch(fileFormat)
	{	case 'X':
			os << "<CrackList>" << endl;
			break;
		
		default:
			break;
	}
}

// called when start new crack - only needed for XML output
void EndCrack(ostream &os)
{
	switch(fileFormat)
	{	case 'X':
			os << "</CrackList>" << endl;
			break;
		
		default:
			break;
	}
}

#pragma mark UTILITIES

/**********************************************************
    get archive record size using current save order
	and save offsets to possible data
*/
int CalcArchiveSize(int vernum)
{
    int i;
    
    // pad if needed, and truncate if unknowns are all 'N'
    if(strlen(mpmOrder)<2) strcpy(mpmOrder,"mY");
    if(strlen(mpmOrder)<ARCH_MAXMPMITEMS)
    {	for(i=(int)strlen(mpmOrder);i<ARCH_MAXMPMITEMS;i++)
            mpmOrder[i]='N';
    }
	else if(strlen(mpmOrder)>ARCH_MAXMPMITEMS)
    {	// File might have more data, but OK if all unknowns are 'N'
		for(i=ARCH_MAXMPMITEMS;i<(int)strlen(mpmOrder);i++)
		{	if(mpmOrder[i]!='N')
				return FileAccessErr;
		}
    }
	// cut of at known items
    mpmOrder[ARCH_MAXMPMITEMS]=0;
	
    // pad if needed, and truncate if unknowns are all 'N'
    if(strlen(crackOrder)<2) strcpy(crackOrder,"mY");
    if(strlen(crackOrder)<ARCH_MAXCRACKITEMS)
    {	for(i=(int)strlen(crackOrder);i<ARCH_MAXCRACKITEMS;i++)
            crackOrder[i]='N';
    }
	else if(strlen(crackOrder)>ARCH_MAXCRACKITEMS)
    {	// File might have more data, but OK if all unknowns are 'N'
		for(i=ARCH_MAXCRACKITEMS;i<(int)strlen(crackOrder);i++)
		{	if(crackOrder[i]!='N')
				return FileAccessErr;
		}
    }
	// cut of at known items
    crackOrder[ARCH_MAXCRACKITEMS]=0;
	
	// check it has defaults
	if(mpmOrder[ARCH_Defaults]=='N' || mpmOrder[ARCH_OldOrigPosition]=='Y' || mpmOrder[ARCH_ver2Empty]=='Y')
		return FileAccessErr;
	if(crackOrder[ARCH_Defaults]=='N')
		return FileAccessErr;
    
	// sizes
	int thirdAngle=0;
	if(threeD)
	{	vectorSize=3*sizeof(double);
		tensorSize=6*sizeof(double);
		if(vernum>=6) thirdAngle=sizeof(double);
	}
	else
	{	vectorSize=2*sizeof(double);
		tensorSize=4*sizeof(double);
	}
    
    // check what is there for material points
	
	/* ARCH_Defaults are
		2D: elemID (int), mass (double), matId (short) angle (double), thickness (double),
					pos (Vector), origPos (Vector)
		3D: thickness replaced by second rotation and Vectors longer and add extra double 
					in vernum==6 or higher for third rotation angle
	*/
	int mpmRecSize=sizeof(int)+3*sizeof(double)+thirdAngle+2*vectorSize+sizeof(short)+2;
	angleOffset=sizeof(int)+sizeof(double)+sizeof(short)+2;
	if(threeD)
	{	angleYOffset=angleOffset+sizeof(double);
		angleXOffset=angleYOffset+sizeof(double);
	}
	posOffset=angleOffset+2*sizeof(double)+thirdAngle;
	origPosOffset=posOffset+vectorSize;
	
    if(mpmOrder[ARCH_Velocity]=='Y')
	{	velocityOffset=mpmRecSize;
        mpmRecSize+=vectorSize;
	}
    if(mpmOrder[ARCH_Stress]=='Y')
	{	stressOffset=mpmRecSize;
        mpmRecSize+=tensorSize;
	}
    if(mpmOrder[ARCH_Strain]=='Y')
	{	strainOffset=mpmRecSize;
        mpmRecSize+=tensorSize;
	}
    if(mpmOrder[ARCH_PlasticStrain]=='Y')
	{	plStrainOffset=mpmRecSize;
        mpmRecSize+=tensorSize;
	}
    if(mpmOrder[ARCH_WorkEnergy]=='Y')
    {   workEnergyOffset=mpmRecSize;
        mpmRecSize+=sizeof(double);
    }
    if(mpmOrder[ARCH_DeltaTemp]=='Y')
	{	tempOffset=mpmRecSize;
        mpmRecSize+=sizeof(double);
	}
    if(mpmOrder[ARCH_PlasticEnergy]=='Y')
	{	plasticEnergyOffset=mpmRecSize;
        mpmRecSize+=sizeof(double);
	}
    if(mpmOrder[ARCH_ShearComponents]=='Y')
        mpmRecSize+=2*sizeof(double);
    if(mpmOrder[ARCH_StrainEnergy]=='Y')
	{	strainEnergyOffset=mpmRecSize;
        mpmRecSize+=sizeof(double);
	}
    if(mpmOrder[ARCH_History]=='Y')
    {   history1Offset=mpmRecSize;
        mpmRecSize+=sizeof(double);
    }
	else if(mpmOrder[ARCH_History]!='N')
	{	if(mpmOrder[ARCH_History]&0x01)
        {   history1Offset=mpmRecSize;
            mpmRecSize+=sizeof(double);
        }
		if(mpmOrder[ARCH_History]&0x02)
        {   history2Offset=mpmRecSize;
            mpmRecSize+=sizeof(double);
        }
		if(mpmOrder[ARCH_History]&0x04)
        {   history3Offset=mpmRecSize;
            mpmRecSize+=sizeof(double);
        }
		if(mpmOrder[ARCH_History]&0x08)
        {   history4Offset=mpmRecSize;
            mpmRecSize+=sizeof(double);
        }
	}
    if(mpmOrder[ARCH_Concentration]=='Y')
	{	concOffset=mpmRecSize;
        mpmRecSize+=sizeof(double)+vectorSize;
	}
    if(mpmOrder[ARCH_HeatEnergy]=='Y')
    {   heatEnergyOffset=mpmRecSize;
        mpmRecSize+=sizeof(double);
    }
    if(mpmOrder[ARCH_ElementCrossings]=='Y')
        mpmRecSize+=sizeof(int);
    if(mpmOrder[ARCH_RotStrain]=='Y')
	{	rotStrainOffset=mpmRecSize;
		if(threeD)
			mpmRecSize+=3*sizeof(double);
		else
			mpmRecSize+=sizeof(double);
	}
    if(mpmOrder[ARCH_DamageNormal]=='Y')
	{	if(threeD)
			mpmRecSize+=3*sizeof(double);
		else
			mpmRecSize+=2*sizeof(double);
	}
    if(mpmOrder[ARCH_SpinMomentum]=='Y')
	{	if(threeD)
			mpmRecSize+=3*sizeof(double);
		else
			mpmRecSize+=sizeof(double);
	}
    if(mpmOrder[ARCH_SpinVelocity]=='Y')
	{	if(threeD)
			mpmRecSize+=3*sizeof(double);
		else
			mpmRecSize+=sizeof(double);
	}
		
    // check what will be there for crack segments
 	crackPosOffset=sizeof(int)+sizeof(double)+sizeof(short)+2;
	
	/* ARCH_Defaults are
		plane element (int), empty (double), newCrack (short+2), pos (Vector), origPos (Vector),
			above element (int), above (Vector) below element (int), below (Vector)
	*/
    int crackRecSize=sizeof(int)+sizeof(double)+sizeof(short)+2;
	//crackPosOffset=crackRecSize;
	crackRecSize+=2*sizeof(int)+8*sizeof(double);
    if(crackOrder[ARCH_JIntegral]=='Y')
	{	jIntOffset=crackRecSize;
        crackRecSize+=2*sizeof(double);
	}
    if(crackOrder[ARCH_StressIntensity]=='Y')
	{	kSifOffset=crackRecSize;
        crackRecSize+=2*sizeof(double);
	}
    if(crackOrder[ARCH_BalanceResults]=='Y')
        crackRecSize+=sizeof(int)+2*sizeof(double);
    
    // record is max of these two sizes
    recSize=mpmRecSize>crackRecSize ? mpmRecSize : crackRecSize;
	
	return noErr;
}

// reverse bytes
int Reverse(char *nptr,int nlen)
{
    int i,imax;
    char *rptr,temp;
    
    rptr=nptr+nlen-1;
    imax=nlen>>1;
    for(i=0;i<imax;i++)
    {	temp=*nptr;
        *nptr++=*rptr;
        *rptr--=temp;
    }
    return nlen;
}
