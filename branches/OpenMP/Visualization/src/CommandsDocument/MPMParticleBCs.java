/*
 * MPMParticleBCs.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 26 Sep 2013.
 * Copyright (c) 2013 RSAC Software. All rights reserved.
 * 
 * Start BC: LoadLine, LoadArc, LoadRect, MoveBox
 * Subordinate: Load, Traction ConcentrationFlux, HeatFlux, LoadType
 */

import java.util.*;

public class MPMParticleBCs
{
	private CmdViewer doc;
	private StringBuffer xmlbcs = null;
	private String bcAttrs;
	private String bcCmd;
	private StringBuffer bcSettings;

	private int inBC;

	public static final int LOADLINE_BC=1;
	public static final int LOADARC_BC=2;
	public static final int LOADRECT_BC=3;
	public static final int LOADBOX_BC=4;
	
	public static final int ADD_LOAD=1;
	public static final int ADD_TRACTION=2;
	public static final int ADD_HEATFLUX=3;
	public static final int ADD_CONCENTRATIONFLUX=4;
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public MPMParticleBCs(CmdViewer cmdDoc)
	{	// save parent CmdViewer
		doc = cmdDoc;
	}
	
	public void initRunSettings()
	{	inBC = 0;
		xmlbcs = new StringBuffer("");
	}
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	// start particle BC line
	// LoadLine x1,y1,x2,y2,(tolerance)
	public void StartLoadLine(ArrayList<String> args) throws Exception
	{
	    // MPM Only
		doc.requiresMPM(args);

	    // verify not nested
	    if(inBC != 0)
	    	throw new Exception("LoadLine, LoadArc, LoadRect, and LoadBox cannot be nested:\n"+args);
	    
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
    	
	    inBC = LOADLINE_BC;
	    bcCmd = "BCLine";
 	}

	// MoveLine, MoveArc, or MoveBox is done
	public void EndLoadBlock(ArrayList<String> args,int endType) throws Exception
	{
		if(inBC != endType)
			throw new Exception("'"+args.get(0)+"' does not match current boundary condition block:\n"+args);
		
		// append block
		xmlbcs.append("    "+bcAttrs+bcSettings+"    </"+bcCmd+">\n");
		
		inBC = 0;
	}
	
	// Add on of the following boundary conditions
	// Load dir,style,arg1,arg2 
	// Traction dir,face,style,arg1,arg2
	// HeatFlux "external",face,style,arg1,arg2
	// ConcentrationFlux "external",face,style,arg1,arg2
	public void AddCondition(ArrayList<String> args,int theType) throws Exception
	{
		// MPM only
		doc.requiresMPM(args);
		
		if(inBC == 0)
			throw new Exception("'"+args.get(0)+"' command must by in 'LoadLine',\n'LoadArc', 'LoadRect', or 'LoadBox' block:\n"+args);
		
		// always needs #1, #2, and #3 (those with face need #4 to)
		if(args.size()<5 || (theType==ADD_LOAD && args.size()<4))
	    	throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		// read direction
		HashMap<String,Integer> options = new HashMap<String,Integer>(7);
		if(theType==ADD_CONCENTRATIONFLUX || theType==ADD_HEATFLUX)
		{	options.put("external", new Integer(1));
		}
		else
		{	options.put("x", new Integer(1));
			options.put("y", new Integer(2));
			options.put("z", new Integer(3));
			options.put("R", new Integer(1));
			options.put("Z", new Integer(2));
			if(theType==ADD_TRACTION)
			{	options.put("normal", new Integer(11));
				options.put("shear", new Integer(12));
			}
		}
		int dof = doc.readIntOption(args.get(1),options,"Load, traction, or flux style direction");
		
		// face if needed
		int face=0;
		int arg=2;
		if(theType!=ADD_LOAD)
		{	face = doc.readIntArg(args.get(2));
			if(face<1 || face>6)
		    	throw new Exception("'"+args.get(0)+"' has ionvalid face:\n"+args);
			arg++;
		}
		
		// read style
		options = new HashMap<String,Integer>(5);
		options.put("constant", new Integer(1));
		options.put("linear", new Integer(2));
		options.put("sine", new Integer(3));
		options.put("cosine", new Integer(4));
		options.put("function", new Integer(6));
		int style = doc.readIntOption(args.get(arg),options,"Load, traction, or flux style");
		arg++;
		
		// read arg1 and arg2
		double arg1=0.,arg2=0.;
		String function = null;
		boolean hasArg2 = false;
		
		// arg1
		if(args.size()>arg)
		{	if(style==6)
				function = doc.readStringArg(args.get(arg));
			else
				arg1 = doc.readDoubleArg(args.get(arg));
			arg++;
		}
		
		// arg2
		if(args.size()>arg)
		{	hasArg2 = true;
			arg2 = doc.readDoubleArg(args.get(arg));
		}
		
		if(theType==ADD_LOAD)
		{	// <LoadBC dir='1' style='1' load='400' time='0.0' function='x*t'/>
			bcSettings.append("      <LoadBC dir='"+dof+"' style='"+style+"'");
			if(style==6)
				bcSettings.append(" function='"+function+"'");
			else
				bcSettings.append(" load='"+arg1+"'");
		}
		else if(theType==ADD_TRACTION)
		{	// <TractionBC dir="2" face="3" style="1" stress="1" time='0.0' function='x*t'/>
			bcSettings.append("      <TractionBC dir='"+dof+"' face='"+face+"' style='"+style+"'");
			if(style==6)
				bcSettings.append(" function='"+function+"'");
			else
				bcSettings.append(" stress='"+arg1+"'");
		}
		else if(theType==ADD_HEATFLUX)
		{	// <HeatFluxBC dir='1' face='1' style='1' value='0' time='0.0' function='sinh(t)'/>
			bcSettings.append("      <HeatFluxBC dir='"+dof+"' face='"+face+"' style='"+style+"'");
			if(style==6)
				bcSettings.append(" function='"+function+"'");
			else
				bcSettings.append(" value='"+arg1+"'");
			
		}
		else
		{	// <ConcFluxBC dir='1' face='1' style='1' value='0' time='0.0' function='sinh(t)'/>
			bcSettings.append("      <ConcFluxBC dir='"+dof+"' face='"+face+"' style='"+style+"'");
			if(style==6)
				bcSettings.append(" function='"+function+"'");
			else
				bcSettings.append(" value='"+arg1+"'");
		}
		
		// time arg
		if(hasArg2) bcSettings.append(" time='"+arg2+"'");
		
		// end command
		bcSettings.append("/>\n");
	}
	
	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	// return xml data
	public String toXMLString()
	{	return xmlbcs.toString();
	}

}
