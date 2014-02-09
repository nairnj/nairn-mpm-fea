/*
 * MPMGrid.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 27 Dec 2012.
 * Copyright (c) 2012 RSAC Software. All rights reserved.
 */

import java.util.*;

public class MPMGrid
{
	private CmdViewer doc;
	private int[] ncells;
	private double[] ratios;
	private double[] symmin;
	private double[] symmax;
	private boolean[] hasmin;
	private boolean[] hasmax;
	private double xmin,xmax,ymin,ymax,zmin,zmax;
	private boolean hasGrid;
	double thickness;

	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public MPMGrid(CmdViewer cmdDoc)
	{	// save parent CmdViewer
		doc = cmdDoc;
		ncells = new int[3];
		ratios = new double[3];
		symmin = new double[3];
		symmax = new double[3];
		hasmin = new boolean[3];
		hasmax = new boolean[3];
	}
	
	public void initRunSettings()
	{	ncells[0] = 0;
		ncells[1] = 0;
		ncells[2] = 0;
		ratios[0] = 1.;
		ratios[1] = 1.;
		ratios[2] = 1.;
		thickness = -1.;
		hasmin[0] = false;
		hasmin[1] = false;
		hasmin[2] = false;
		hasmax[0] = false;
		hasmax[1] = false;
		hasmax[2] = false;
		hasGrid = false;
	}
	
	//----------------------------------------------------------------------------
	// Methods
	//----------------------------------------------------------------------------
	
	// GridHoriz, GridVert, or GridDepth #1 (number of cells) #2 (ratio, not used)
	public void doGridAxis(ArrayList<String> args,int axis) throws Exception
	{
	    // MPM Only
	    doc.requiresMPM(args);
	    
	    // needs one
	    if(args.size()<2)
	    	throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
	    
	    // number of cells
	    ncells[axis] = doc.readIntArg(args.get(1));
	    ratios[axis] = 1.;
	    
	    // symmetry planes
	    if(args.size()>2)
	    {	double sym = doc.readDoubleArg(args.get(2));
	    	int symdir = -1;
	    	if(args.size()>3) symdir = doc.readIntArg(args.get(3));
	    	if(symdir==-1)
	    	{	hasmin[axis]=true;
	    		symmin[axis]=sym;
	    	}
	    	else if(symdir==1)
	    	{	hasmax[axis]=true;
	    		symmax[axis]=sym;
	    	}
	    	else
	    		throw new Exception("'"+args.get(0)+"' has invalid symmetry direction:\n"+args);
	    	if(args.size()>4)
	    	{	sym = doc.readDoubleArg(args.get(4));
	    		if(symdir==-1)
	    		{	hasmax[axis]=true;
	    			symmax[axis]=sym;
	    		}
	    		else if(symdir==1)
	    		{	hasmin[axis]=true;
	    			symmin[axis]=sym;
	    		}
	    	}
	    }
	}
	
	// GridRect xmin,xmax,ymin,ymax (zmin,zmax if 3D)
	public void doGridRect(ArrayList<String> args) throws Exception
	{
	    // FEA Only
		doc.requiresMPM(args);

	    // needs at least 5 arguments
	    if(args.size()<5)
	    	throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
	    
	    // limits
	    hasGrid = true;
	    double temp;
	    xmin = doc.readDoubleArg(args.get(1));
	    xmax = doc.readDoubleArg(args.get(2));
	    if(xmax < xmin)
	    {	temp = xmin;
	    	xmin = xmax;
	    	xmax = temp;
	    }
	    ymin = doc.readDoubleArg(args.get(3));
	    ymax = doc.readDoubleArg(args.get(4));
	    if(ymax < ymin)
	    {	temp = ymin;
	    	ymin = ymax;
	    	ymax = temp;
	    }
	    
	    // z axis
	    if(doc.isMPM3D())
		{	if(args.size()<7)
		    	throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
	    	zmin = doc.readDoubleArg(args.get(5));
	    	zmax = doc.readDoubleArg(args.get(6));
		    if(zmax < zmin)
		    {	temp = zmin;
		    	zmin = zmax;
		    	zmax = temp;
		    }
	    }
	}
	
	// GridThickness #1
	public void doGridThickness(ArrayList<String> args) throws Exception
	{
	    // MPM Only
	    doc.requiresMPM(args);
	    
	    // needs one
	    if(args.size()<2)
	    	throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
	    
	    thickness = doc.readDoubleArg(args.get(1));
	    if(thickness<=0.)
	    	throw new Exception("The grid thickness must be positive:\n"+args);
	}
	    
	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	public String toXMLString()
	{
		if(!hasGrid) return "";
		
		// Grid element
		StringBuffer xml = new StringBuffer("    <Grid");
		xml.append(" xmin='"+xmin+"' xmax='"+xmax+"'");
		xml.append(" ymin='"+ymin+"' ymax='"+ymax+"'");
		if(doc.isMPM3D())
			xml.append(" zmin='"+zmin+"' zmax='"+zmax+"'");
		else if(thickness>0.)
			xml.append(" thickness='"+thickness+"'");
		xml.append(">\n");
		
		// Horiz, Vert, Depth elements
		xml.append("      <Horiz nx='"+ncells[0]+"'");
		if(hasmin[0]) xml.append(" symmin='"+symmin[0]+"'");
		if(hasmax[0]) xml.append(" symmax='"+symmax[0]+"'");
		xml.append("/>\n");
		xml.append("      <Vert ny='"+ncells[1]+"'");
		if(hasmin[1]) xml.append(" symmin='"+symmin[1]+"'");
		if(hasmax[1]) xml.append(" symmax='"+symmax[1]+"'");
		xml.append("/>\n");
		if(doc.isMPM3D())
		{	xml.append("      <Depth nz='"+ncells[2]+"'");
			if(hasmin[2]) xml.append(" symmin='"+symmin[2]+"'");
			if(hasmax[2]) xml.append(" symmax='"+symmax[2]+"'");
			xml.append("/>\n");
		}
		
		// End Grid
		xml.append("    </Grid>\n");
		
		return xml.toString();
	}
}
