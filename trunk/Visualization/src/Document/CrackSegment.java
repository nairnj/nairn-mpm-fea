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
	public double J1,J2;
	public double KI,KII;
	public int crackIncrements;
	public double release,absorb;

	// initialize
	CrackSegment()
	{	inElem=new int[4];
		xpos=new double[4];
		ypos=new double[4];
	}
	
	//---------------------------------------------------------------------
	// read record from and archive file into this material point
	//---------------------------------------------------------------------
	
	public void readRecord(ByteBuffer bb,char[] crackOrder,double lengthScale,double timeScale)
	{	// required elements
		inElem[PLANE_POS]=bb.getInt()-1;			// in element number (zero based)
		bb.position(bb.position()+8);				// skip empty double
		startFlag=bb.getShort();					// -1 to start a new crack
		tractionMaterial=bb.getShort();				// >0 if traction law
		
		// position in length units
		xpos[PLANE_POS]=bb.getDouble()*lengthScale;
		ypos[PLANE_POS]=bb.getDouble()*lengthScale;
		
		// original position
		xpos[ORIG_POS]=bb.getDouble()*lengthScale;
		ypos[ORIG_POS]=bb.getDouble()*lengthScale;
		
		// above the crack
		inElem[ABOVE_POS]=bb.getInt();
		xpos[ABOVE_POS]=bb.getDouble()*lengthScale;
		ypos[ABOVE_POS]=bb.getDouble()*lengthScale;
		
		// above the crack
		inElem[BELOW_POS]=bb.getInt();
		xpos[BELOW_POS]=bb.getDouble()*lengthScale;
		ypos[BELOW_POS]=bb.getDouble()*lengthScale;
		
		// J1 and J2 in J/m^2
		if(crackOrder[ReadArchive.ARCH_JIntegral]=='Y')
		{	J1=bb.getDouble();
			J2=bb.getDouble();
		}
		
		// KI ad KII in MPa mm^1/2
		if(crackOrder[ReadArchive.ARCH_StressIntensity]=='Y')
		{	KI=bb.getDouble();
			KII=bb.getDouble();
		}
		
		// Energy flow in J/m^2
		if(crackOrder[ReadArchive.ARCH_BalanceResults]=='Y')
		{	crackIncrements=bb.getInt();
			release=bb.getDouble();
			absorb=bb.getDouble();
		}
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
