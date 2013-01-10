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
	}
	
	public void initRunSettings()
	{	ncells[0] = 0;
		ncells[1] = 0;
		ncells[2] = 0;
		ratios[0] = 1.;
		ratios[1] = 1.;
		ratios[2] = 1.;
		thickness = -1.;
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
	    	throw new Exception("'"+args.get(0)+"' has too few parameters: "+args);
	    
	    // number of cells
	    ncells[axis] = doc.readIntArg(args.get(1));
	    
	    // ratios
	    if(args.size()>2)
	    	ratios[axis] = doc.readDoubleArg(args.get(2));
	    
	    // validity?
	    if(ncells[axis]<1 || ratios[axis]==0.)
	    	throw new Exception("The grid parameters are invalid: "+args);
	}
	
	// GridRect xmin,xmax,ymin,ymax (zmin,zmax if 3D)
	public void doGridRect(ArrayList<String> args) throws Exception
	{
	    // FEA Only
		doc.requiresMPM(args);

	    // needs at least 5 arguments
	    if(args.size()<5)
	    	throw new Exception("'"+args.get(0)+"' has too few parameters: "+args);
	    
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
		    	throw new Exception("'"+args.get(0)+"' has too few parameters: "+args);
	    	zmin = doc.readDoubleArg(args.get(3));
	    	zmax = doc.readDoubleArg(args.get(4));
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
	    	throw new Exception("'"+args.get(0)+"' has too few parameters: "+args);
	    
	    thickness = doc.readDoubleArg(args.get(1));
	    if(thickness<=0.)
	    	throw new Exception("The grid thickness must be positive: "+args);
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
			xml.append(" zmin='"+xmin+"' xmax='"+zmax+"'");
		else if(thickness>0.)
			xml.append(" thickness='"+thickness+"'");
		xml.append(">\n");
		
		// Horiz, Vert, Depth elements
		xml.append("      <Horiz nx='"+ncells[0]+"' rx='1'/>\n");
		xml.append("      <Vert ny='"+ncells[1]+"' ry='1'/>\n");
		if(doc.isMPM3D())
			xml.append("      <Depth nz='"+ncells[2]+"' rz='1'/>\n");
		
		// End Grid
		xml.append("    </Grid>\n");
		
		return xml.toString();
	}
}
