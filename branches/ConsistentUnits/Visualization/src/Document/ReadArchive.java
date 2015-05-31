/*******************************************************************
	ReadArchive.java
	NairnFEAMPMViz

	Created by John Nairn on 2/28/07.
	Copyright 2007 RSAC Software. All rights reserved.
	
	No instance, only static variables and methods for code to
	read MPM archive files
*******************************************************************/

import java.io.*;
import java.nio.*;

public class ReadArchive
{
	// archving material points
	static final int ARCH_Defaults=1;
	static final int ARCH_Velocity=2;
	static final int ARCH_Stress=3;
	static final int ARCH_Strain=4;
	static final int ARCH_PlasticStrain=5;
	static final int ARCH_OldOrigPosition=6;
	static final int ARCH_WorkEnergy=7;
	static final int ARCH_DeltaTemp=8;
	static final int ARCH_PlasticEnergy=9;
	static final int ARCH_ver2Empty=10;
	static final int ARCH_ShearComponents=11;
	static final int ARCH_StrainEnergy=12;
	static final int ARCH_History=13;
	static final int ARCH_Concentration=14;
	static final int ARCH_HeatEnergy=15;
	static final int ARCH_ElementCrossings=16;
	static final int ARCH_RotStrain=17;
	static final int ARCH_MAXMPMITEMS=18;

	// Archiving options for crack segments
	static final int ARCH_JIntegral=2;
	static final int ARCH_StressIntensity=3;
	static final int ARCH_BalanceResults=4;
	static final int ARCH_MAXCRACKITEMS=5;
	
	// FEA archiving options
	static final int ARCH_FEADisplacements=0;
	static final int ARCH_FEAAvgStress=1;
	static final int ARCH_FEAElemStress=2;
	static final int ARCH_FEAElemForce=3;
	static final int ARCH_FEAElemEnergy=4;
	static final int ARCH_Interfaces=5;
	
	// assumed binary sizes
	public static final int sizeofInt=4;
	public static final int sizeofLong=4;
	public static final int sizeofDouble=8;
	public static final int sizeofShort=2;
	
	// read buffer for a file
	public static ByteBuffer openBuffer(ResultsDocument doc,File archivePath) throws Exception
	{	// determine if need to reverse the bytes from this file
		ByteOrder fileSystem = doc.archFormat.charAt(0)=='m' ? ByteOrder.BIG_ENDIAN : ByteOrder.LITTLE_ENDIAN;
		
		// read the file
		ByteBuffer bb;
		try
		{	// load ByteBuffer and set byte order
			bb=ByteBuffer.allocate((int)archivePath.length());
			FileInputStream fs=new FileInputStream(archivePath);
			fs.read(bb.array(),0,(int)archivePath.length());
			fs.close();
			bb.order(fileSystem);
		}
		catch (Exception fe)
		{	throw new Exception("Error reading archive file: "+fe.getMessage());
		}
		
		return bb;
	}
	
	// load data into an archive files into the ResultsDocument
	public static void load(ResultsDocument doc,File archivePath) throws Exception
	{
		// read the file
		try
		{	// determine if need to reverse the bytes from this file
			ByteOrder fileSystem = doc.archFormat.charAt(0)=='m' ? ByteOrder.BIG_ENDIAN : ByteOrder.LITTLE_ENDIAN;
		
			// formats
			char[] mpmOrder=new char[ARCH_MAXMPMITEMS];
			char[] crackOrder=new char[ARCH_MAXCRACKITEMS];
			doc.archFormat.getChars(0,ARCH_MAXMPMITEMS,mpmOrder,0);
			doc.crackFormat.getChars(0,ARCH_MAXCRACKITEMS,crackOrder,0);
		
			// load ByteBuffer and set byte order
			ByteBuffer bb=ByteBuffer.allocate((int)archivePath.length());
			FileInputStream fs=new FileInputStream(archivePath);
			fs.read(bb.array(),0,(int)archivePath.length());
			fs.close();
			bb.order(fileSystem);
			
			// version
			byte[] version=new byte[4];
			bb.get(version);
			int headerLength=4;
			int vernum=0;
			if(version[2]>='1' && version[2]<='9')
				vernum+=10*((int)version[2]-'0');
			if(version[3]>='0' && version[3]<='9')
				vernum+=(int)(version[3]-'0');
			if(vernum>=4)
			{	headerLength=64;
				bb.position(headerLength);
			}
			else if(vernum!=3)
				throw new Exception("Archive file is too old for this tool");
			
			// number of records (some may be cracks)
			int newMpms=(int)((archivePath.length()-headerLength)/doc.recSize);
			int newsegs=0;
			doc.mpmPoints.clear();
			doc.mpmPoints.ensureCapacity(newMpms);
			
			// loop over material points
			int p;
			for(p=0;p<newMpms;p++)
			{	// create material point
				int pos=bb.position();
				MaterialPoint mpm=new MaterialPoint(p+1);
				mpm.readRecord(bb,mpmOrder,doc.units,doc.is3D());
				
				// A negative material number means start of crack particles
				if(mpm.material<0)
				{	newsegs=newMpms-p;
					bb.position(pos);
					break;
				}
					
				// add to ResultsDocument
				doc.mpmPoints.add(mpm);
				
				// next record
				bb.position(pos+doc.recSize);
			}
			
			doc.mpmCracks.clear();
			doc.mpmCracks.ensureCapacity(newsegs);
			
			// read crack segments
			if(newsegs>0)
			{	if(crackOrder[ReadArchive.ARCH_Defaults]!='Y')
					throw new Exception("Old crack archive format not supported by this tool");
				CrackHeader ch=null;
				for(p=0;p<newsegs;p++)
				{	// read matnum
					int pos=bb.position();
					bb.position(pos+12);
					int matnum=bb.getShort();
					bb.position(pos);
					
					// start new crack
					if(matnum==-1)
					{	ch=new CrackHeader();
						doc.mpmCracks.add(ch);
					}

					// add the segment
					CrackSegment cs=new CrackSegment();
					cs.readRecord(bb,crackOrder,doc.units);
					ch.add(cs);
					
					// next record
					bb.position(pos+doc.recSize);
				}
			}
				
		}
		catch(Exception e)
		{	throw new Exception("Archive File Error:  " + e.getMessage());
		}
	}
	
	// calculate size of each binary record in the archive files
	public static int getRecordSize(ResultsDocument doc)
	{
		char[] mpmOrder,crackOrder;
		mpmOrder=new char[ARCH_MAXMPMITEMS];
		crackOrder=new char[ARCH_MAXCRACKITEMS];
		
		// these were padded when read
		doc.archFormat.getChars(0,ARCH_MAXMPMITEMS,mpmOrder,0);
		doc.crackFormat.getChars(0,ARCH_MAXCRACKITEMS,crackOrder,0);
		
		// initialize
		int mpmRecSize=0;
		int crackRecSize=0;
		
		int vectorSize,tensorSize;
		if(doc.is3D())
		{	vectorSize=3*sizeofDouble;
			tensorSize=6*sizeofDouble;
		}
		else
		{	vectorSize=2*sizeofDouble;
			tensorSize=4*sizeofDouble;
		}
		// check what will be there for material points
		if(mpmOrder[ARCH_Defaults]=='Y')
		{	mpmRecSize+=sizeofInt+3*sizeofDouble+2*vectorSize+sizeofShort+2;
			if(doc.is3D()) mpmRecSize+=sizeofDouble;		// extra double for 3rd angle
		}
		else
			mpmRecSize+=sizeofInt+5*sizeofDouble+sizeofShort+2;
		
		if(mpmOrder[ARCH_Velocity]=='Y')
			mpmRecSize+=vectorSize;
		if(mpmOrder[ARCH_Stress]=='Y')
			mpmRecSize+=tensorSize;
		if(mpmOrder[ARCH_Strain]=='Y')
			mpmRecSize+=tensorSize;
		if(mpmOrder[ARCH_PlasticStrain]=='Y')
			mpmRecSize+=tensorSize;
		if(mpmOrder[ARCH_OldOrigPosition]=='Y')
			mpmRecSize+=vectorSize;
		if(mpmOrder[ARCH_WorkEnergy]=='Y')
			mpmRecSize+=sizeofDouble;
		if(mpmOrder[ARCH_DeltaTemp]=='Y')
			mpmRecSize+=sizeofDouble;
		if(mpmOrder[ARCH_PlasticEnergy]=='Y')
			mpmRecSize+=sizeofDouble;
		if(mpmOrder[ARCH_ver2Empty]=='Y')
			mpmRecSize+=sizeofDouble;
		if(mpmOrder[ARCH_ShearComponents]=='Y')
			mpmRecSize+=2*sizeofDouble;
		if(mpmOrder[ARCH_StrainEnergy]=='Y')
			mpmRecSize+=sizeofDouble;
		if(mpmOrder[ARCH_History]=='Y')
			mpmRecSize+=sizeofDouble;
		else if(mpmOrder[ARCH_History]!='N')
		{	int history=(int)mpmOrder[ReadArchive.ARCH_History];
			if((history & 0x01) !=0) mpmRecSize+=sizeofDouble;
			if((history & 0x02) !=0) mpmRecSize+=sizeofDouble;
			if((history & 0x04) !=0) mpmRecSize+=sizeofDouble;
			if((history & 0x08) !=0) mpmRecSize+=sizeofDouble;
		}
		if(mpmOrder[ARCH_Concentration]=='Y')
			mpmRecSize+=vectorSize+sizeofDouble;
		if(mpmOrder[ARCH_HeatEnergy]=='Y')
			mpmRecSize+=sizeofDouble;
		if(mpmOrder[ARCH_ElementCrossings]=='Y')
			mpmRecSize+=sizeofInt;
		if(mpmOrder[ARCH_RotStrain]=='Y')
		{	mpmRecSize+=sizeofDouble;
			if(doc.is3D()) mpmRecSize+=2*sizeofDouble;
		}
			   
		// check what will be there for crack segments
		crackRecSize+=sizeofInt+sizeofDouble+sizeofShort+2;
		if(crackOrder[ARCH_Defaults]=='Y')
			crackRecSize+=2*sizeofInt+8*sizeofDouble;
		else
			crackRecSize+=20*sizeofDouble;
		if(crackOrder[ARCH_JIntegral]=='Y')
			crackRecSize+=2*sizeofDouble;
		if(crackOrder[ARCH_StressIntensity]=='Y')
			crackRecSize+=2*sizeofDouble;
		if(crackOrder[ARCH_BalanceResults]=='Y')
			crackRecSize+=sizeofInt+2*sizeofDouble;
		
		// record is max of these two sizes
		int recSize=mpmRecSize>crackRecSize ? mpmRecSize : crackRecSize;
		return recSize;
	}
}
