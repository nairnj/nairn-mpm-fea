/*******************************************************************
	BoundaryCondition.java
	NairnFEAMPMViz

	Created by John Nairn on 2/27/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.geom.*;

public class BoundaryCondition
{
	// variables and constants
	static final double BC_SIZE=0.035;
	static final int CONSTANT_VALUE=1;
	static final int LINEAR_VALUE=2;
	static final int SINE_VALUE=3;
	static final int COSINE_VALUE=4;
	static final int SILENT=5;
	static final int FUNCTION_VALUE=6;
	static final int FEA_DISPBC=7;
	static final int FEA_ROTATEBC=8;
	static final int FEA_LOADBC=9;
	static final int FEA_ELEMSTRESSBC=10;
	static final int SKEWXY_DIRECTION=0;	
	static final int X_DIRECTION=1;
	static final int Y_DIRECTION=2;
		
	protected int node,dof,style;
	protected double value,argument,angle;
	
	// initialize
	BoundaryCondition(int nodeNum,int bcDof,int theID,double theVal,double theArg,double theAngle)
	{
		node=nodeNum-1;
		dof=bcDof;
		style=theID;
		value=theVal;
		argument=theArg;
		angle=theAngle;
	}
	
	// utility to make an arrow return length in relLength
	// pointing in the x direction
	public GeneralPath makeArrow(Rectangle2D.Double bounds,float relLength)
	{
		// get size
		float length=(float)(Math.max(bounds.getWidth(),bounds.getHeight())*BC_SIZE*relLength);
		float awidth=length/6.f;
		float adepth=length/4.f;
		float ainset=length/5.f;
		
		// Path and generic arrow (horizontal, starting at (0,0)  ---->
		GeneralPath arrow=new GeneralPath();
		arrow.moveTo(0.f,0.f);
		arrow.lineTo(length,0.f);
		arrow.lineTo(length-adepth,awidth);
		arrow.lineTo(length-ainset,0.f);
		arrow.lineTo(length-adepth,-awidth);
		arrow.lineTo(length,(float)0.);
		
		return arrow;
	}
	
	// accessors
	public double getValue() { return value; }
	
}
