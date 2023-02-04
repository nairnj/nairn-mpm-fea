/*******************************************************************
	CrackSegment.java
	NairnFEAMPMViz

	Created by John Nairn on 3/11/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.nio.*;
import java.awt.geom.*;

public class CrackSegment
{
	// variables and constants
	public static int PLANE_POS=0;
	public static int ORIG_POS=1;
	public static int ABOVE_POS=2;
	public static int BELOW_POS=3;
	
	public int[] inElem;
	public short startFlag;
	public short tractionMaterial;
	public double[] xpos;
	public double[] ypos;
	public double[] zpos;
	public double J1,J2;
	public double KI,KII;
	public int czmHasDisp;
	public double czmGI,czmGII;
	public double[] tractionData;

	// initialize
	CrackSegment()
	{	inElem=new int[4];
		xpos=new double[4];
		ypos=new double[4];
		zpos=new double[4];
		tractionData = null;
	}
	
	//---------------------------------------------------------------------
	// read record from and archive file into this material point
	//---------------------------------------------------------------------
	
	public void readRecord(ByteBuffer bb,char[] crackOrder,JNUnits units,boolean has3D)
	{	// required elements
		inElem[PLANE_POS]=bb.getInt()-1;			// in element number (zero based)
		bb.position(bb.position()+8);				// skip empty double
		startFlag=bb.getShort();					// -1 to start a new crack
		tractionMaterial=bb.getShort();				// >0 if traction law
		
		// position in length units
		xpos[PLANE_POS]=bb.getDouble()*units.lengthScale();
		ypos[PLANE_POS]=bb.getDouble()*units.lengthScale();
		if(has3D) zpos[PLANE_POS]=bb.getDouble()*units.lengthScale();
		
		// original position
		xpos[ORIG_POS]=bb.getDouble()*units.lengthScale();
		ypos[ORIG_POS]=bb.getDouble()*units.lengthScale();
		if(has3D) zpos[ORIG_POS]=bb.getDouble()*units.lengthScale();
		
		// above the crack
		inElem[ABOVE_POS]=bb.getInt();
		xpos[ABOVE_POS]=bb.getDouble()*units.lengthScale();
		ypos[ABOVE_POS]=bb.getDouble()*units.lengthScale();
		if(has3D) zpos[ABOVE_POS]=bb.getDouble()*units.lengthScale();
		
		// above the crack
		inElem[BELOW_POS]=bb.getInt();
		xpos[BELOW_POS]=bb.getDouble()*units.lengthScale();
		ypos[BELOW_POS]=bb.getDouble()*units.lengthScale();
		if(has3D) zpos[BELOW_POS]=bb.getDouble()*units.lengthScale();
		
		// J1 and J2 in J/m^2
		if(crackOrder[ReadArchive.ARCH_JIntegral]=='Y')
		{	J1=bb.getDouble()*units.outputERRScale();
			J2=bb.getDouble()*units.outputERRScale();
		}
		
		// KI ad KII in MPa m^1/2
		if(crackOrder[ReadArchive.ARCH_StressIntensity]=='Y')
		{	KI=bb.getDouble()*units.outputSIScale();
			KII=bb.getDouble()*units.outputSIScale();
		}
		
		// Energy flow in J/m^2
		if(crackOrder[ReadArchive.ARCH_CZMDeltaG]=='Y')
		{	czmHasDisp=bb.getInt();
			czmGI=bb.getDouble()*units.outputERRScale();
			czmGII=bb.getDouble()*units.outputERRScale();
		}
		
		// traction history
		char history=crackOrder[ReadArchive.ARCH_Traction15];
		if(history=='Y')
			setTractionHistory(1,bb.getDouble());
		else if(history!='N')
		{	if((history&0x01)!=0) setTractionHistory(1,bb.getDouble());
			if((history&0x02)!=0) setTractionHistory(2,bb.getDouble());
			if((history&0x04)!=0) setTractionHistory(3,bb.getDouble());
			if((history&0x08)!=0) setTractionHistory(4,bb.getDouble());
			if((history&0x10)!=0) setTractionHistory(5,bb.getDouble());
		}
		history=crackOrder[ReadArchive.ARCH_Traction610];
		if(history=='Y')
			setTractionHistory(6,bb.getDouble());
		else if(history!='N')
		{	if((history&0x01)!=0) setTractionHistory(6,bb.getDouble());
			if((history&0x02)!=0) setTractionHistory(7,bb.getDouble());
			if((history&0x04)!=0) setTractionHistory(8,bb.getDouble());
			if((history&0x08)!=0) setTractionHistory(9,bb.getDouble());
			if((history&0x10)!=0) setTractionHistory(10,bb.getDouble());
		}
	}
	
	// set history variable (create array if needed)
	private void setTractionHistory(int num,double value)
	{	if(tractionData==null) tractionData=new double[10];
		tractionData[num-1]=value;
	}
	
	// set history variable (create array if needed)
	public double getTractionHistory(int num)
	{	if(tractionData==null) return 0.;
		if(num>tractionData.length) return 0;
		return tractionData[num-1];
	}
	
	//---------------------------------------------------------------------
	// accessors
	//---------------------------------------------------------------------
	
	public Point2D.Double getPt()
	{	return new Point2D.Double(xpos[PLANE_POS],ypos[PLANE_POS]);
	}
	
	public Point2D.Double getMedianPosition()
	{	return new Point2D.Double((xpos[ABOVE_POS]+xpos[BELOW_POS])/2.,(ypos[ABOVE_POS]+ypos[BELOW_POS])/2.);
	}

	public Point2D.Double getCOD()
	{	return new Point2D.Double(xpos[BELOW_POS]-xpos[ABOVE_POS],ypos[BELOW_POS]-ypos[ABOVE_POS]);
	}
}
