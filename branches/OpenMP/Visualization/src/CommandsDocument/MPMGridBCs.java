/*
 * MPMGridBCs.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 7 Sep 2013.
 * Copyright (c) 2013 RSAC Software. All rights reserved.
 * 
 * Start BC: MoveLine, MoveArv, MoveBox
 * Subordinate: Velocity, Temperature, Concentration
 */

import java.util.*;

public class MPMGridBCs
{
	private CmdViewer doc;
	private String bcAttrs;
	private String bcCmd;
	private StringBuffer bcSettings;
	private StringBuffer xmlbcs = null;
	private int boundaryID;

	private int inBC;

	public static final int MOVELINE_BC=1;
	public static final int MOVEARC_BC=2;
	public static final int MOVEBOX_BC=3;
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public MPMGridBCs(CmdViewer cmdDoc)
	{	// save parent CmdViewer
		doc = cmdDoc;
	}
	
	public void initRunSettings()
	{	inBC = 0;
		xmlbcs = new StringBuffer("");
		boundaryID = 0;
	}

	//----------------------------------------------------------------------------
	// Methods
	//----------------------------------------------------------------------------
	
	// start grid BC line
	// MoveLine x1,y1,x2,y2,(tolerance)
	public void StartMoveLine(ArrayList<String> args) throws Exception
	{
	    // MPM Only
		doc.requiresMPM(args);

	    // verify not nested
	    if(inBC != 0)
	    	throw new Exception("MoveLine, MoveArc, and MoveBox cannot be nested:\n"+args);
	    
	    // needs at least 4 arguments
	    if(args.size()<5)
	    	throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
	    
	    // get x1,y1,x2,y2
	    double x1 = doc.readDoubleArg(args.get(1));
	    double y1 = doc.readDoubleArg(args.get(2));
	    double x2 = doc.readDoubleArg(args.get(3));
	    double y2 = doc.readDoubleArg(args.get(4));
	    	
	    // get optional tolerance
	    double tolerance = -1.;
	    if(args.size()>5)
	    	tolerance = doc.readDoubleArg(args.get(5));
	    
    	bcAttrs = "<BCLine x1='"+x1+"' y1='"+y1+"' x2='"+x2+"' y2='"+y2+"'";
    	if(tolerance > 0.)
    		bcAttrs = bcAttrs + " tolerance='" + tolerance + "'>\n";
    	else
    		bcAttrs = bcAttrs + ">\n";
        bcSettings = new StringBuffer("");
    	
	    inBC = MOVELINE_BC;
	    bcCmd = "BCLine";
 	}
	
	// start grid BC line
	public void StartMoveArc(ArrayList<String> args) throws Exception
	{	throw new Exception("MoveArc command not implemented yet.");
	}
	
	// start grid BC line
	public void StartMoveBox(ArrayList<String> args) throws Exception
	{	
		// MPM Only
		doc.requiresMPM(args);
		if(!doc.isMPM3D())
			throw new Exception("MoveBoxonly allowedin 3D MPM:\n"+args);

		// verify not nested
		if(inBC != 0)
			throw new Exception("MoveLine, MoveArc, and MoveBox cannot be nested:\n"+args);
    
		// needs at least 6 arguments
		if(args.size()<7)
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
    
		// get x1,y1,x2,y2,z1,z2
		double x1 = doc.readDoubleArg(args.get(1));
		double y1 = doc.readDoubleArg(args.get(2));
		double z1 = doc.readDoubleArg(args.get(3));
		double x2 = doc.readDoubleArg(args.get(4));
		double y2 = doc.readDoubleArg(args.get(5));
		double z2 = doc.readDoubleArg(args.get(6));
    	
		// get optional axis
		int axis = -1;
		if(args.size()>7)
		{	HashMap<String,Integer> options = new HashMap<String,Integer>(3);
			options.put("x", new Integer(1));
			options.put("y", new Integer(2));
			options.put("z", new Integer(3));
			axis = doc.readIntOption(args.get(7),options,"Cylinder axis");
			if(axis<1 || axis>3)
				throw new Exception("'MoveBox' cylinder axis is not valis:\n"+args);
		}
    
		bcAttrs = "<BCBox xmin='"+x1+"' ymin='"+y1+"' zmin='"+z1+"' xmax='"+x2+"' ymax='"+y2+"' zmax='"+z2+"'";
		if(axis > 0)
			bcAttrs = bcAttrs + " axis='" + axis + "'>\n";
		else
			bcAttrs = bcAttrs + ">\n";
		bcSettings = new StringBuffer("");
	
		inBC = MOVEBOX_BC;
		bcCmd = "BCBox";
	}

	// MoveLine, MoveArc, or MoveBox is done
	public void EndMoveBlock(ArrayList<String> args,int endType) throws Exception
	{
		if(inBC != endType)
			throw new Exception("'"+args.get(0)+"' does not match current boundary condition block:\n"+args);
		
		// append block
		xmlbcs.append("    "+bcAttrs+bcSettings+"    </"+bcCmd+">\n");
		
		inBC = 0;
	}
	
	// add velocity condition
	// Velocity (x or y or z),type,<arg1>,<arg2>
	// Velocity (skewed),type,arg1,arg2,angle1,<angle2>
	public void AddVelocity(ArrayList<String> args) throws Exception
	{
		// FEA only
		doc.requiresMPM(args);
		
		if(inBC == 0)
			throw new Exception("'Velocity' command must by in 'MoveLine', 'MoveArc', or 'MoveBox' block:\n"+args);
		
		// always needs #1 and #2
		if(args.size()<3)
	    	throw new Exception("'Velocity' has too few parameters:\n"+args);
		
		// read direction
		HashMap<String,Integer> options = new HashMap<String,Integer>(8);
		options.put("x", new Integer(1));
		options.put("y", new Integer(2));
		options.put("z", new Integer(3));
		options.put("skewxy", new Integer(12));
		options.put("skewrz", new Integer(12));
		options.put("skewxz", new Integer(13));
		options.put("skewyz", new Integer(23));
		options.put("skewxyz", new Integer(123));
		int dof = doc.readIntOption(args.get(1),options,"Velocity direction");
		
		if(dof>10 && args.size()<6)
	    	throw new Exception("Skewed 'Velocity' has too few parameters:\n"+args);
		
		// read style
		options = new HashMap<String,Integer>(5);
		options.put("constant", new Integer(1));
		options.put("linear", new Integer(2));
		options.put("sine", new Integer(3));
		options.put("cosine", new Integer(4));
		options.put("function", new Integer(6));
		int style = doc.readIntOption(args.get(2),options,"Velocity style");
		
		// all need a arg1 except constant
		if(style!=1 && args.size()<4)
	    	throw new Exception("'Velocity' has too few parameters:\n"+args);
		
		// read arg1 and arg2
		double arg1=0.,arg2=0.;
		String function = null;
		boolean hasArg2 = false;
		
		// arg1
		if(args.size()>3)
		{	if(style==6)
				function = doc.readStringArg(args.get(3));
			else
				arg1 = doc.readDoubleArg(args.get(3));
		}
		
		// arg2
		if(args.size()>4)
		{	hasArg2 = true;
			arg2 = doc.readDoubleArg(args.get(4));
		}
		
		// angles
		double angle1=0.,angle2=0.;
		if(dof>10)
		{	angle1 = doc.readDoubleArg(args.get(5));
			if(args.size()>6)
				angle2 = doc.readDoubleArg(args.get(6));
		}
		
		// add to xml
		bcSettings.append("      <DisBC dir='"+dof+"' style='"+style+"'");
		if(style==6)
			bcSettings.append(" function='"+function+"'");
		else
			bcSettings.append(" vel='"+arg1+"'");
		if(hasArg2) bcSettings.append(" time='"+arg2+"'");
		
		if(dof>10) bcSettings.append(" angle='"+angle1+"'");
		if(dof>100) bcSettings.append(" angle2='"+angle2+"'");
		
		if(boundaryID!=0) bcSettings.append(" id='"+boundaryID+"'");
		bcSettings.append("/>\n");
	}

	// add velocity condition
	public void AddTemperature(ArrayList<String> args) throws Exception
	{	throw new Exception("Temperature command not implemented yet");		
	}

	// add velocity condition
	public void AddConcentration(ArrayList<String> args) throws Exception
	{	throw new Exception("Concentration command not implemented yet");
	}
	
	// set boundary ID
	public void SetBoundaryID(ArrayList<String> args) throws Exception
	{
		// MPM only
		doc.requiresMPM(args);
		
		// no arg reverts to 0
		if(args.size()<2)
		{	boundaryID = 0;
			return;
		}
		
		// set
		boundaryID = doc.readIntArg(args.get(1));
	}

	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	// return xml data
	public String toXMLString()
	{	return xmlbcs.toString();
	}

}
