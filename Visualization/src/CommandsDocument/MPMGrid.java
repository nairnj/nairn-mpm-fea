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
	private int[] rstyle;
	private double[] symmin;
	private double[] symmax;
	private boolean[] hasmin;
	private boolean[] hasmax;
	private double xmin,xmax,ymin,ymax,zmin,zmax;
	
	private String tartanBorder;
	private String tartanAOIs;
	private boolean isTartan;
	
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
		rstyle = new int[3];
		symmin = new double[3];
		symmax = new double[3];
		hasmin = new boolean[3];
		hasmax = new boolean[3];
	}
	
	public void initRunSettings()
	{	ncells[0] = 0;
		ncells[1] = 0;
		ncells[2] = 0;
		ratios[0] = -1.;
		ratios[1] = -1.;
		ratios[2] = -1.;
		rstyle[0] = -1;
		rstyle[1] = -1;
		rstyle[2] = -1;
		thickness = -1.;
		hasmin[0] = false;
		hasmin[1] = false;
		hasmin[2] = false;
		hasmax[0] = false;
		hasmax[1] = false;
		hasmax[2] = false;
		hasGrid = false;
		isTartan = false;
		tartanBorder = null;
		tartanAOIs = "";
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
	    
		// if #2 is "R-Style" then get tartan grid parameters
		hasmin[axis] = false;
		hasmax[axis] = false;;
	    ratios[axis] = -1.;
	    rstyle[axis] = -1;
	    int argnum = 2;
	    if(args.size()>2)
	    {	Object rarg = doc.readStringOrDoubleArg(args.get(2));
			if(rarg.getClass().equals(String.class))
			{	String theText = (String)rarg;
				isTartan=true;
				if(!theText.toLowerCase().equals("r-style"))
			    	throw new Exception("The second parameter is neither a number nor 'R-Style':\n"+args);
					
				// read next two if there, for Tartan R and style
				if(args.size()<4)
					throw new Exception("Entering 'R-Style' options but no R provided in next parameter:\n"+args);
				ratios[axis] = doc.readDoubleArg(args.get(3));
				if(ratios[axis]<1.)
					throw new Exception("'R-Style' R value must be >=1:\n"+args);
				
				// optional style
				if(args.size()>4)
				{	HashMap<String,Integer> options = new HashMap<String,Integer>(10);
					options.put("geometric", new Integer(0));
					options.put("linear", new Integer(1));
					rstyle[axis] = doc.readIntOption(args.get(4),options,"MPM update method");
				}
				
				// symmetry block start
				argnum = 5;
			}
		}

		// symmetry planes
	    if(args.size()>argnum)
	    {	double sym = doc.readDoubleArg(args.get(argnum));
	    	int symdir = -1;
	    	if(args.size()>argnum+1) symdir = doc.readIntArg(args.get(argnum+1));
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
	    	if(args.size()>argnum+2)
	    	{	sym = doc.readDoubleArg(args.get(argnum+2));
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
	
	// Tartan Grid border
	public void doTartanBorder(ArrayList<String> args) throws Exception
	{
	    // MPM Only
		doc.requiresMPM(args);
		
		int x1=1,x2=1,y1=1,y2=1,z1=1,z2=1;
		
		if(args.size()>2)
	    {	x1 = doc.readIntArg(args.get(1));
	    	x2 = doc.readIntArg(args.get(2));
	    }
		if(args.size()>4)
	    {	y1 = doc.readIntArg(args.get(3));
	    	y2 = doc.readIntArg(args.get(4));
	    }
		if(args.size()>6)
	    {	z1 = doc.readIntArg(args.get(5));
	    	z2 = doc.readIntArg(args.get(6));
	    }
		
		if(doc.isMPM3D())
		{	tartanBorder = "      <Border xmin='"+x1+"' xmax='"+x2
										+"' ymin='"+y1+"' ymax='"+y2
										+"' zmin='"+z1+"' ymax='"+z2+"'/>\n";
		}
		else
		{	tartanBorder = "      <Border xmin='"+x1+"' xmax='"+x2
										+"' ymin='"+y1+"' ymax='"+y2+"'/>\n";
		}
	}
	
	// GridRect xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz
	public void doTartanAOI(ArrayList<String> args) throws Exception
	{
	    // FEA Only
		doc.requiresMPM(args);

	    // needs at least 7 arguments
	    if(args.size()<7 || (doc.isMPM3D() && args.size()<10))
	    	throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
	    
	    // x limits
	    double temp;
	    double x1 = doc.readDoubleArg(args.get(1));
	    double x2 = doc.readDoubleArg(args.get(2));
	    if(x2 < x1)
	    {	temp = x1;
	    	x1 = x2;
	    	x2 = temp;
	    }
	    int nx = doc.readIntArg(args.get(3));
	    
	    // y limits
	    double y1 = doc.readDoubleArg(args.get(4));
	    double y2 = doc.readDoubleArg(args.get(5));
	    if(y2 < y1)
	    {	temp = y1;
	    	y1 = y2;
	    	y2 = temp;
	    }
	    int ny = doc.readIntArg(args.get(6));
	    
	    if(doc.isMPM3D())
	    {	double z1 = doc.readDoubleArg(args.get(7));
	    	double z2 = doc.readDoubleArg(args.get(8));
	    	if(z2 < z1)
	    	{	temp = z1;
	    		z1 = z2;
	    		z2 = temp;
	    	}
	    	int nz = doc.readIntArg(args.get(9));
	    	
	    	tartanAOIs = tartanAOIs+"      <AreaOfInterest x1='"+x1+"' x2='"+x2+"' nx='"+nx
									+"' y1='"+y1+"' y2='"+y2+"' ny='"+ny
									+"' z1='"+z1+"' z2='"+z2+"' nz='"+nz+"'/>\n";

	    }
	    else
	    {	tartanAOIs = tartanAOIs+"      <AreaOfInterest x1='"+x1+"' x2='"+x2+"' nx='"+nx
	    							+"' y1='"+y1+"' y2='"+y2+"' ny='"+ny+"'/>\n";
	    }
	}
	    
	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	public String toXMLString()
	{
		if(!hasGrid) return "";
		
		// Grid element
		StringBuffer xml = new StringBuffer("    <Grid");
		xml.append(" xmin='"+doc.formatDble(xmin)+"' xmax='"+doc.formatDble(xmax)+"'");
		xml.append(" ymin='"+doc.formatDble(ymin)+"' ymax='"+doc.formatDble(ymax)+"'");
		if(doc.isMPM3D())
			xml.append(" zmin='"+doc.formatDble(zmin)+"' zmax='"+doc.formatDble(zmax)+"'");
		else if(thickness>0.)
			xml.append(" thickness='"+doc.formatDble(thickness)+"'");
		xml.append(">\n");
		
		// Horiz, Vert, Depth elements
		xml.append("      <Horiz nx='"+ncells[0]+"'");
		if(ratios[0]>0.) xml.append(" rx='"+ratios[0]+"'");
		if(rstyle[0]>=0) xml.append(" style='"+rstyle[0]+"'");
		if(hasmin[0]) xml.append(" symmin='"+doc.formatDble(symmin[0])+"'");
		if(hasmax[0]) xml.append(" symmax='"+doc.formatDble(symmax[0])+"'");
		xml.append("/>\n");
		xml.append("      <Vert ny='"+ncells[1]+"'");
		if(ratios[1]>0.) xml.append(" ry='"+ratios[1]+"'");
		if(rstyle[1]>=0) xml.append(" style='"+rstyle[1]+"'");
		if(hasmin[1]) xml.append(" symmin='"+doc.formatDble(symmin[1])+"'");
		if(hasmax[1]) xml.append(" symmax='"+doc.formatDble(symmax[1])+"'");
		xml.append("/>\n");
		if(doc.isMPM3D())
		{	xml.append("      <Depth nz='"+ncells[2]+"'");
			if(ratios[2]>0.) xml.append(" rz='"+ratios[2]+"'");
			if(rstyle[2]>=0) xml.append(" style='"+rstyle[2]+"'");
			if(hasmin[2]) xml.append(" symmin='"+doc.formatDble(symmin[2])+"'");
			if(hasmax[2]) xml.append(" symmax='"+doc.formatDble(symmax[2])+"'");
			xml.append("/>\n");
		}
		if(isTartan)
		{	if(tartanBorder!=null) xml.append(tartanBorder);
			if(tartanAOIs!=null) xml.append(tartanAOIs);
		}
		
		// check for more xml
		String more = doc.getXMLData("grid");
		if(more != null) xml.append(more);
		
		// End Grid
		xml.append("    </Grid>\n");
		
		return xml.toString();
	}
}
