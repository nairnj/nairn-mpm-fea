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

#include "ExtractMPM.hpp"

// Global variable settings
char fileFormat='T';
char *outfile=NULL;
bool stepExtension=false;
char *headerName=NULL;
char fileExtension[20];
char mpmOrder[50];
char crackOrder[50];
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
int angleOffset=-1,posOffset=-1,stressOffset=-1,strainOffset=-1;
int crackPosOffset=-1,jIntOffset=-1,kSifOffset=-1;
int velocityOffset=-1,origPosOffset=-1,plStrainOffset=-1;
int tempOffset=-1,concOffset=-1,strainEnergyOffset=-1,plasticEnergyOffset=-1;

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
			if(argv[parmInd][optInd]=='H')
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
                else if(strcmp(parm,"pressure")==0)
                    q=PRESSURE;
                else if(strcmp(parm,"vonmises")==0)
                    q=EQUIVSTRESS;
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
				else if(strcmp(parm,"dispx")==0)
					q=DISPX;
				else if(strcmp(parm,"dispy")==0)
					q=DISPY;
				else if(strcmp(parm,"dispz")==0)
					q=DISPZ;
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
			{	parm=NextArgument(++parmInd,argv,argc,'m');
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
            "    -b format          MPM archive file format (only used for 'ver3' files)\n"
            "    -c format          MPM archive file crack format (only used for 'ver3' files)\n"
            "    -C                 Extract crack data only\n"
            "    -d                 Output as little Endian binary doubles file\n"
            "    -D                 Output as big Endian binary doubles file\n"
            "    -f                 Output as little Endian binary floats file\n"
            "    -F                 Output as big Endian binary floats file\n"
            "    -h                 Write header in file (default: false)\n"
            "    -H                 Show this help and exit\n"
            "    -m num             Exclude this material\n"
            "    -M num             Include only this material\n"
            "    -o path            Output file name with no extension\n"
			"		                       (quoted if name has spaces)\n"
			"    -s                 Include step number from input file numeric\n"
			"                              extension in the output file name\n"
            "    -P                 Extract particle data only (the default)\n"
            "    -q name            Add column with named quantity\n"
            "    -T                 Output as tab-delimited text file (the default)\n"
            "    -V                 Output as VTK Legacy file (particle date only)\n"
            "    -X                 Output as XML file\n"
            "    -2                 Input file has 2D data (the default)\n"
            "    -3                 Input file has 3D data (only used for 'ver3' files)\n"
			"\n"
            "See http://oregonstate.edu/nairnj for documentation.\n\n"
          <<  endl;
}

#pragma mark EXTRACTION CODE

// main entry to extract results from one file
int ExtractMPMData(const char *mpmFile,int fileIndex,int lastIndex)
{
	int i;
	
	FILE *fp;
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
	long fileLength=ftell(fp);
	rewind(fp);
	
	// read entire file into a buffer
	unsigned char *buffer=(unsigned char *)malloc(fileLength);
	if(buffer==NULL)
	{	cerr << "Out of memory creating file reading buffer" << endl;
		fclose(fp);
		return MemoryErr;
	}
	if(fread(buffer,fileLength,1,fp)!=1)
	{	cerr << "Input file '" << mpmFile << "' access error (not all read)" << endl;
		fclose(fp);
		return FileAccessErr;
	}
	fclose(fp);
	
	// look for file header
    unsigned char *ap=buffer;
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
	{	cerr << "Input file format too old for this tool (missing current defaults)" << endl;
		return FileAccessErr;
	}

	// output file
	ofstream fout;
	if(outfile!=NULL)
	{	char fname[256],step[256];
		if(stepExtension)
		{	int ext=0,dot=strlen(mpmFile)-1;
			while(dot>=0 && mpmFile[dot]!='.') dot--;
			if(dot>=0)
			{	step[ext++]='_';
				dot++;
				while(mpmFile[dot]!=0)
				{	if(mpmFile[dot]<'0' || mpmFile[dot]>'9')
					{	// skip if not a number
						ext=0;
						break;
					}
					step[ext++]=mpmFile[dot++];
				}
			}
			step[ext]=0;
		}
		if(lastIndex==0)
		{	if(stepExtension)
				sprintf(fname,"%s%s.%s",outfile,step,fileExtension);
			else
				sprintf(fname,"%s.%s",outfile,fileExtension);
		}
		else
		{	if(stepExtension && strlen(step)>0)
				sprintf(fname,"%s%s.%s",outfile,step,fileExtension);
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
		{	cerr << "Exports ot VTK Legacy file format cannot be for crack data" << endl;
			return FileAccessErr;
		}
		VTKLegacy(os,ap,fileLength,mpmFile);
		free(buffer);
		return noErr;
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
				for(i=0;i<includeMaterial.size();i++)
				{	sprintf(headLine," %d",includeMaterial[i]);
					strcat(headBuffer,headLine);
				}
				strcat(headBuffer,"\n");
			}
			if(excludeMaterial.size()>0)
			{	strcat(headBuffer,"Excluded_Materials");
				for(i=0;i<excludeMaterial.size();i++)
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
		{	for(i=0;i<quantity.size();i++)
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
	int nummpms=fileLength/recSize;
	int p;
	short *mptr;
	BeginMP(os);
	for(p=0;p<nummpms;p++)
	{	short matnum=pointMatnum(ap);
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
		for(i=0;i<quantity.size();i++)
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
		{	// read tipMatNum
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
				for(i=0;i<quantity.size();i++)
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
	return noErr;
}

// called when start MP output - only needed for XML output
void VTKLegacy(ostream &os,unsigned char *ap,long fileLength,const char *mpmFile)
{
	os << "# vtk DataFile Version 4.2" << endl;
	if(headerName!=NULL)
		os << headerName;
	else
		os << "NairnMPM particle data";
	os << " from file ";
	os << mpmFile ;
	if(archiveTimeMs>-0.5)
		os << " at time " << archiveTimeMs << " ms" ;
	os << endl;
	os << "ASCII" << endl;
	os << "DATASET POLYDATA" << endl;
	
	// count points to extract
	int nummpms=fileLength/recSize;
	int p,numExtract=0;
	short *mptr;
	unsigned char *origap=ap;
	for(p=0;p<nummpms;p++)
	{	short matnum=pointMatnum(ap);
		if(matnum<0) break;
		if(!skipThisPoint(matnum)) numExtract++;
		ap+=recSize;
	}
	
	os << "POINTS " << numExtract << " double" << endl;
	if(numExtract==0) return;
	
	// extract the point positions
	ap=origap;
	for(p=0;p<nummpms;p++)
	{	short matnum=pointMatnum(ap);
		if(matnum<0) break;
		if(!skipThisPoint(matnum))
		{	OutputDouble((double *)(ap+posOffset),0,0,reverseFromInput,os,XPOS);
		    OutputDouble((double *)(ap+posOffset),1,' ',reverseFromInput,os,YPOS);
		    if(threeD) OutputDouble((double *)(ap+posOffset),2,' ',reverseFromInput,os,ZPOS);
			os << endl;
		}
		ap+=recSize;
	}
	
	// extract quantities
	if(quantity.size()==0) return;
	os << "POINT_DATA " << numExtract << endl;
	int i;
	double matDouble,zeroDouble=0.;
	for(i=0;i<quantity.size();i++)
	{	os << "SCALARS " << quantityName[i] << " double 1" << endl;
		os << "LOOKUP_TABLE default" << endl;
		
		ap=origap;
		for(p=0;p<nummpms;p++)
		{	short matnum=pointMatnum(ap);
			if(matnum<0) break;
			if(!skipThisPoint(matnum))
			{	OutputQuantity(i,ap,os,matnum,0);
				os << endl;
			}
			ap+=recSize;
		}
	}
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
                double pressure = (*sxx+*syy+*szz)/3.;
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
                if(reverseFromInput)
                {	Reverse((char *)sxx,sizeof(double));
                    Reverse((char *)syy,sizeof(double));
                    Reverse((char *)szz,sizeof(double));
                    Reverse((char *)sxy,sizeof(double));
                }
                double *sxz,*syz;
                if(threeD)
                {   double *sxz=(sxx+SXZ-SXX);
                    double *syz=(sxx+SYZ-SXX);
                    if(reverseFromInput)
                    {	Reverse((char *)sxz,sizeof(double));
                        Reverse((char *)syz,sizeof(double));
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
		
		case VELX:
		case VELY:
		case VELZ:
			if(velocityOffset>0 && (threeD || quantity[i]<VELZ))
				OutputDouble((double *)(ap+velocityOffset),quantity[i]-VELX,delim,reverseFromInput,os,quantity[i]);
			else
				OutputDouble(&zeroDouble,0,delim,false,os,quantity[i]);
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
		for(i=0;i<includeMaterial.size();i++)
		{	if(includeMaterial[i]==matnum)
			{	skip=false;
				break;
			}
		}
	}
	// skip is false here if no list of included materials or if is in that list
		
	// can override with exclusion
	if(!skip)
	{	for(i=0;i<excludeMaterial.size();i++)
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
			fdata=*data;
			if(thisEndian=='m') Reverse((char *)&fdata,sizeof(float));
			os.write((char *)&fdata,sizeof(float));
			break;
			
		case 'F':
			if(mustReverse) Reverse((char *)data,sizeof(double));
			fdata=*data;
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


// called when record is done. mpTrue for end of material point record of false for crack record
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
    
    // pad if needed, and truncate if not recognized
    if(strlen(mpmOrder)<2) strcpy(mpmOrder,"mY");
    if(strlen(mpmOrder)<ARCH_MAXMPMITEMS)
    {	for(i=strlen(mpmOrder);i<ARCH_MAXMPMITEMS;i++)
            mpmOrder[i]='N';
    }
    mpmOrder[ARCH_MAXMPMITEMS]=0;
	
    if(strlen(crackOrder)<2) strcpy(crackOrder,"mY");
    if(strlen(crackOrder)<ARCH_MAXCRACKITEMS)
    {	for(i=strlen(crackOrder);i<ARCH_MAXCRACKITEMS;i++)
            crackOrder[i]='N';
    }
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
    if(mpmOrder[ARCH_ExtWork]=='Y')
        mpmRecSize+=sizeof(double);
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
        mpmRecSize+=sizeof(double);
	else if(mpmOrder[ARCH_History]!='N')
	{	if(mpmOrder[ARCH_History]&0x01) mpmRecSize+=sizeof(double);
		if(mpmOrder[ARCH_History]&0x02) mpmRecSize+=sizeof(double);
		if(mpmOrder[ARCH_History]&0x04) mpmRecSize+=sizeof(double);
		if(mpmOrder[ARCH_History]&0x08) mpmRecSize+=sizeof(double);
	}
    if(mpmOrder[ARCH_Concentration]=='Y')
	{	concOffset=mpmRecSize;
        mpmRecSize+=sizeof(double)+vectorSize;
	}
    if(mpmOrder[ARCH_ThermalEnergy]=='Y')
        mpmRecSize+=sizeof(double);
    if(mpmOrder[ARCH_ElementCrossings]=='Y')
        mpmRecSize+=sizeof(int);
    if(mpmOrder[ARCH_RotStrain]=='Y')
	{	if(threeD)
			mpmRecSize+=3*sizeof(double);
		else
			mpmRecSize+=sizeof(double);
	}
		
    // check what will be there for crack segments
 	crackPosOffset=sizeof(int)+sizeof(double)+sizeof(short)+2;
	
	/* ARCH_Defaults are
		planeInElem (int), empty (double), newCrack (short+2), pos (Vector), origPos (Vector),
			aboveInElem (int), above (Vector) belowInElem (int), below (Vector)
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
