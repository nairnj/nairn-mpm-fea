/*
 * FEABCs.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 17 Dec 2012.
 * Copyright (c) 2012 RSAC Software. All rights reserved.
 */

import java.util.*;

public class FEABCs
{
	private CmdViewer doc;
	
	private int inBC;
	private int select;
	private String bcAttrs;
	private StringBuffer xmlbcs = null;
	private StringBuffer bcSettings;
	private String rsKeyID;
	private double rsx;
	private double rsy;
	
	private static final int FIXLINE_BC=1;
	private static final int FIXPOINT_BC=2;
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public FEABCs(CmdViewer cmdDoc)
	{	// save parent CmdViewer
		doc = cmdDoc;
	}
	
	public void initRunSettings()
	{	inBC = 0;
		xmlbcs = new StringBuffer("");
		rsKeyID = null;
	}
	
	//----------------------------------------------------------------------------
	// Methods
	//----------------------------------------------------------------------------
	
	// start FEA mesh area with Area #1,#2,<#3>
	// #1 is material, #2 is thickness, #3 is material angle
	public void StartFixLine(ArrayList<String> args) throws Exception
	{
	    // FEA Only
	    if(!doc.isFEA())
	    	throw new Exception("The 'FixLine' command can only be used in FEA commands.");

	    // verify not nested
	    if(inBC != 0)
	    	throw new Exception("FixLine and FixPt cannot be nested.");
	    
	    // one means a path
	    if(args.size()==2)
	    {	String pathID = doc.readStringArg(args.get(1));
	    	if(!doc.areas.hasPath(pathID))
	    		throw new Exception("'FixLine' uses an undefined path: "+pathID);
	    	
	    	// start the command (leave room for option select)
	    	bcAttrs = "path='"+pathID+"'";
	    }
	    
	    // needs at least 5 arguments
	    else if(args.size()<5)
	    	throw new Exception("'FixLine' has too few parameters.");
	    
	    else
	    {	// get x1,y1,x2,y2
	    	double x1 = doc.readDoubleArg(args.get(1));
	    	double y1 = doc.readDoubleArg(args.get(2));
	    	double x2 = doc.readDoubleArg(args.get(3));
	    	double y2 = doc.readDoubleArg(args.get(4));
	    	
	    	// get optional tolerance
	    	double tolerance = -1.;
	    	if(args.size()>5)
	    		tolerance = doc.readDoubleArg(args.get(5));
	    	
	    	bcAttrs = "x1='"+x1+"' y1='"+y1+"' x2='"+x2+"' y2='"+y2+"'";
	    	if(tolerance > 0.)
	    		bcAttrs = bcAttrs + " tolerance='" + tolerance + "'";
	    }
	    
	    // now in a BC
	    select = 0;
	    bcSettings = new StringBuffer("");
	    inBC = FIXLINE_BC;
	}
	
	// FixLine done, add to xml
	public void EndFixLine(ArrayList<String> args) throws Exception
	{
		if(inBC != FIXLINE_BC)
			throw new Exception("'EndFixLine' not matched by 'FixLine' command.");
		
		// append block 
		xmlbcs.append("    <BCLine "+bcAttrs);
		if(select==1) xmlbcs.append(" select='1'");
		xmlbcs.append(">\n"+bcSettings+"    </BCLine>\n");
		
		inBC = 0;
	}
	    
	// start FEA mesh area with Area #1,#2,<#3>
	// #1 is material, #2 is thickness, #3 is material angle
	public void StartFixPoint(ArrayList<String> args) throws Exception
	{
	    // FEA Only
	    if(!doc.isFEA())
	    	throw new Exception("The 'FixPoint' command can only be used in FEA commands.");

	    // verify not nested
	    if(inBC != 0)
	    	throw new Exception("FixLine and FixPt cannot be nested.");
	    
	    // one means a keypoint
	    if(args.size()<3)
	    {	String keyID = doc.readStringArg(args.get(1));
	    	if(!doc.areas.hasKeypoint(keyID))
	    		throw new Exception("'FixPoint' uses an undefined keypoint: "+keyID);
	    	
	    	// start the command (leave room for option select)
	    	bcAttrs = "keypt='"+keyID+"'";
	    }
	    
	    // must have at least two arguments
	    else
	    {	// get x,y
	    	double x = doc.readDoubleArg(args.get(1));
	    	double y = doc.readDoubleArg(args.get(2));
	    	
	    	bcAttrs = "x='"+x+"' y='"+y+"'";
	    }
	    
	    // now in a BC
	    select = 0;
	    bcSettings = new StringBuffer("");
	    inBC = FIXPOINT_BC;
	}

	// FixLine done, add to xml
	public void EndFixPoint(ArrayList<String> args) throws Exception
	{
		if(inBC != FIXPOINT_BC)
			throw new Exception("'EndFixPoint' not matched by 'FixPoint' command.");
		
		// append block 
		xmlbcs.append("    <BCPt "+bcAttrs);
		if(select==1) xmlbcs.append(" select='1'");
		xmlbcs.append(">\n"+bcSettings+"    </BCPt>\n");
		
		inBC = 0;
	}
	    
	// Displace nodes (Displacement #1,#2)
	//	#1 is direction, #2 is value
	public void AddDisplacement(ArrayList<String> args) throws Exception
	{
		if(inBC == 0)
			throw new Exception("'Displacement' command must by in 'FixLine' or 'FixPoint' block.");
		
		// read direction
		if(args.size()<2)
	    	throw new Exception("'Displacement' has too few parameters.");
		int dof = readDirection(args.get(1));
		bcSettings.append("      <DisBC dof='"+dof+"'");
		
		// read value
		String value="0";
		if(args.size()>2)
			value = doc.readStringArg(args.get(2));
		bcSettings.append(" disp='"+value+"'");
		
		// end it
		bcSettings.append("/>\n");
	}
	
	// Load nodes (Load #1,#2)
	//	#1 is direction, #2 is value
	public void AddLoad(ArrayList<String> args) throws Exception
	{
		if(inBC == 0)
			throw new Exception("'Load' command must by in 'FixLine' or 'FixPoint' block.");
		
		// read direction
		if(args.size()<3)
	    	throw new Exception("'Load' has too few parameters.");
		int dof = readDirection(args.get(1));
		bcSettings.append("      <LoadBC dof='"+dof+"'");
		
		// read value
		try
		{	double value = doc.readDoubleArg(args.get(2));
			bcSettings.append(" load='"+value+"'");
		}
		catch(Exception e)
		{	String func = doc.readStringArg(args.get(2));
			bcSettings.append(" function='"+func+"'");
		}
		
		// end it
		bcSettings.append("/>\n");
	}
	
	// Resequence Command x and y or keypoint
	public void Resequence(ArrayList<String> args) throws Exception
	{
	    // FEA Only
	    if(!doc.isFEA())
	    	throw new Exception("The 'Resequence' command can only be used in FEA commands.");

	    // one means a keypoint
	    if(args.size()<3)
	    {	rsKeyID = doc.readStringArg(args.get(1));
	    	if(!doc.areas.hasKeypoint(rsKeyID))
	    		throw new Exception("'Resequence' uses an undefined keypoint: "+rsKeyID);
	    }
	    
	    // must have at least two arguments
	    else
	    {	// get x,y
	    	rsx = doc.readDoubleArg(args.get(1));
	    	rsy = doc.readDoubleArg(args.get(2));
	    	rsKeyID = "";
	    }
	}
	
	// read argument and convert to FEA direction
	public int readDirection(String arg) throws Exception
	{	// options
		HashMap<String,Integer> options = new HashMap<String,Integer>(5);
		options.put("x", new Integer(1));
		options.put("r", new Integer(1));
		options.put("y", new Integer(2));
		options.put("z", new Integer(2));
		
		// read it
		int ndir = doc.readIntOption(arg,options,"FEA BC direction");
		return ndir;
	}

	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	// return xml data
	public String toXMLString()
	{	if(rsKeyID == null)
			return xmlbcs.toString();
	
		String reseq;
		if(rsKeyID.length()==0)
			reseq = "    <Resequence x='"+rsx+"' y='"+rsy+"'/>\n";
		else
			reseq = "    <Resequence keypt='"+rsKeyID+"'/>\n";
		return xmlbcs.toString()+reseq;
	}

}
