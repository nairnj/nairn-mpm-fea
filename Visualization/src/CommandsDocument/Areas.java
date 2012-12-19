/*
 * Areas.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 17 Dec 2012.
 * Copyright (c) 2012 RSAC Software. All rights reserved.
 */

import java.util.*;

public class Areas
{
	private CmdViewer doc;
	private boolean inArea;
	private boolean inPath;
	private double originX;
	private double originY;
	
	// current area
	private int matnum;
	private double thickness;
	private String angle;
	private int elType;
	private int flipTriangles;
	private ArrayList<String> paths;			// in current area
	private StringBuffer xmlareas;
	
	// path
	private double ratio;
	private int intervals;
	private String pathID;
	private HashMap<String,Integer> pathIDs;
	private ArrayList<String> keys;				// in current path
	private StringBuffer xmlpaths;

	// keypoints
	private HashMap<String,Integer> keyIDs;
	private StringBuffer xmlkeys;
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public Areas(CmdViewer cmdDoc)
	{	// save parent CmdViewer
		doc = cmdDoc;
	}
	
	public void initRunSettings()
	{	inArea = false;
		inPath = false;
		flipTriangles = 0;
		elType = ElementBase.EIGHT_NODE_ISO;
		pathIDs = new HashMap<String,Integer>(20);
		keyIDs = new HashMap<String,Integer>(50);
		originX = 0.;
		originY = 0.;
		xmlkeys = new StringBuffer("");
		xmlpaths = new StringBuffer("");
		xmlareas = new StringBuffer("");
	}
	
	//----------------------------------------------------------------------------
	// Methods
	//----------------------------------------------------------------------------
	
	// start FEA mesh area with Area #1,#2,<#3>
	// #1 is material, #2 is thickness, #3 is material angle
	public void StartArea(ArrayList<String> args) throws Exception
	{
	    // FEA Only
	    if(!doc.isFEA())
	    	throw new Exception("The 'Area' command can only be used in FEA commands.");

	    // verify not nested
	    if(inArea)
	    	throw new Exception("FEA mesh Areas cannot be nested.");
	    if(inPath)
	    	throw new Exception("FEA mesh Areas cannot be inside Paths.");
	    
	    // activate
	    inArea = true;
	    paths = new ArrayList<String>(4);
	    
	    // read material by ID
	    if(args.size()<2)
		    throw new Exception("'Area' command missing material ID.");
	    String matID = doc.readStringArg(args.get(1));
		matnum=0;
	    if(!matID.equals("_NONE_"))
	    {	matnum = doc.mats.getMatID(matID);
			if(matnum<=0)
				throw new Exception("'Area' command has unknown material ID.");
		}
		
		// read thickness
	    if(args.size()<3)
		    throw new Exception("'Area' command missing thickness.");
	    thickness = doc.readDoubleArg(args.get(2));
		
		// read angle
	    if(args.size()>3)
	    	angle = doc.readStringArg(args.get(3));
	    else
	    	angle = null;
	    
	}
	
	// EndPath command
	public void EndArea(ArrayList<String> args)  throws Exception
	{
	    // check in a path (which must be in an area)
		if(!inArea)
	    	throw new Exception("'EndArea' not matched by 'Area' command.");
		
		// is it acceptable?
	    
	    // create XML for the current path
		xmlareas.append("    <Area mat='"+matnum+"' thick='"+thickness+"'");
		xmlareas.append(" type='"+elType+"' flip='"+flipTriangles+"'");
		if(angle!=null)
			xmlareas.append(" angle='"+angle+"'");
		xmlareas.append(">\n");
		
		int i;
		for(i=0;i<paths.size();i++)
			xmlareas.append("      <path id='"+paths.get(i)+"'/>\n");
		
		xmlareas.append("    </Area>\n");
	    
	    // now done
	    inArea = false;
	}
	
	// Start Path for FEA Mesh - Path #1,#2,<#3>
	// #1 is name or ID, #2 is intervals, #3 is ratio
	public void StartPath(ArrayList<String> args)  throws Exception
	{
	    // FEA Only
	    if(!doc.isFEA())
	    	throw new Exception("The 'Path' command can only be used in FEA commands.");
	    
	    // verify not nested
	    if(inPath)
	    	throw new Exception("FEA mesh Paths cannot be nested.");
	    
	    // activate
	    ratio = 1.;
	    keys = new ArrayList<String>(4);
		
	    // read path name
	    if(args.size()<2)
		    throw new Exception("'Path' command missing path ID.");
	    pathID = doc.readStringArg(args.get(1));
	    Integer existingPath = pathIDs.get(pathID);
		
		// pre-existing path
		if(existingPath != null)
		{	// error if data there to define the path
			if(args.size()>2)
				throw new Exception("Duplicate Path name: "+pathID);
			
			// ... and must be in an Area
			if(!inArea)
				throw new Exception("FEA Path reference (to "+pathID+") must be within an Area.");
			
			// add
			paths.add(pathID);
		}
		
		// new path
		else
		{	// read intervals
			if(args.size()<3)
				throw new Exception("Path ("+pathID+") missing the number of intervals.");
			intervals = doc.readIntArg(args.get(2));
		
			// read ratio
			if(args.size()>3)
				ratio = doc.readDoubleArg(args.get(3));
			
			// add to list
			pathIDs.put(pathID, new Integer(pathIDs.size()));
			
			// if area add there to
			if(inArea) paths.add(pathID);
			
			// now in a path
		    inPath = true;
		}

	}

	// EndPath command
	public void EndPath(ArrayList<String> args)  throws Exception
	{
	    // check in a path (which must be in an area)
		if(!inPath)
	    	throw new Exception("'EndPath' not matched by 'Path' command.");
		
		// is it acceptable?
	    
	    // create XML for the current path
		xmlpaths.append("    <Path id='"+pathID+"' intervals='"+intervals+"'");
		if(ratio!=1.) xmlpaths.append(" ratio='"+ratio+"'");
		xmlpaths.append(">\n");
		
		int i;
		for(i=0;i<keys.size();i++)
			xmlpaths.append("      <keypt id='"+keys.get(i)+"'/>\n");
		
		xmlpaths.append("    </Path>\n");
	    
	    // now done
	    inPath = false;
		
	}
	
	// Add Paths (any number) to current area
	public void AddPaths(ArrayList<String> args)  throws Exception
	{
	    // FEA Only
	    if(!inArea)
	    	throw new Exception("The 'Paths' command can only be used in Area commands.");
	    if(inPath)
	    	throw new Exception("The 'Paths' command cannot be in another 'Path' command.");
	    
	    // add each one
	    int i;
	    for(i=1;i<args.size();i++)
	    {	String nextPath = doc.readStringArg(args.get(i));
	    	Integer existingPath = pathIDs.get(nextPath);
	    	if(existingPath == null)
	    		throw new Exception("Undefined path ("+nextPath+") referenced in a 'Paths' command.");
	    	paths.add(nextPath);
	    }
	}
	    
	// Define keypoint for FEA Mesh - Keypoint #1,#2,#3,<#4>
	// #1 is name or ID, #2,#3 are x and y
	// #4 can be polar for polar coordinate entry
	public void AddKeypoint(ArrayList<String> args)  throws Exception
	{
	    // FEA Only
	    if(!doc.isFEA())
	    	throw new Exception("The 'Keypoint' command can only be used in FEA commands.");

	    // read keypoint name
	    if(args.size()<2)
		    throw new Exception("'Keypoint' command missing keypoint ID.");
	    String keyID = doc.readStringArg(args.get(1));
	    Integer existingKey = keyIDs.get(keyID);
		
		// pre-existing keypoint
		if(existingKey != null)
		{	// cannot define keypoint here using more arguments
			if(args.size()>2)
				throw new Exception("Duplicate keypoint name: "+keyID);
			
			// ... and must be in a Path
			if(!inPath)
				throw new Exception("FEA Keypoint reference (to "+keyID+") must be within a Path.");
			
			// add
			keys.add(keyID);
		}
		
		// new keypoint
		else
		{	// x value
			if(args.size()<3)
				throw new Exception("Keypoint ("+keyID+") missing x value.");
			double x = doc.readDoubleArg(args.get(2));
					
			// y value
			if(args.size()<4)
				throw new Exception("Keypoint ("+keyID+") missing y value.");
			double y = doc.readDoubleArg(args.get(3));
			
			// look for polar
			if(args.size()>4)
			{	String polar = doc.readStringArg(args.get(4));
				if(!polar.toLowerCase().equals("polar"))
					throw new Exception("Keypoint ("+keyID+") has invalid parameter #4.");
					
				double theta=Math.PI*y/180.;
				y=originY+x*Math.sin(theta);
				x=originX+x*Math.cos(theta);
			}
			
			// add to list
			keyIDs.put(keyID, new Integer(keyIDs.size()));

			// add to path
			if(inPath) keys.add(keyID);
		
			// add to XML
			xmlkeys.append("      <pt x='"+x+"' y = '"+y+"' id='"+keyID+"'/>\n");
		}
	}

	// Add Paths (any number) to current area
	public void AddKeypoints(ArrayList<String> args)  throws Exception
	{
	    // FEA Only
	    if(!inPath)
	    	throw new Exception("The 'Keypoints' command can only be used in Path commands.");
	    
	    // add each one
	    int i;
	    for(i=1;i<args.size();i++)
	    {	String nextKey = doc.readStringArg(args.get(i));
	    	Integer existingKey = keyIDs.get(nextKey);
	    	if(existingKey == null)
	    		throw new Exception("Undefined keypoint ("+existingKey+") referenced in a 'Keypoints' command.");
	    	keys.add(nextKey);
	    }
	}

	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------

	// set whenever changed in main commands
	public void setElementType(int lnameEl) { elType = lnameEl; }
	
	// combine to xml data
	public String toXMLString()
	{
		StringBuffer xml = new StringBuffer("");
		xml.append("    <Keypoints>\n"+xmlkeys+"    </Keypoints>\n\n");
		xml.append(xmlpaths+"\n");
		xml.append(xmlareas);
		
		return xml.toString();
	}
	
	// FlipTriangles <#1>
	// If #1=yes flip, otherwise do not flip, 
	public void setFlipTriangles(ArrayList<String> args) throws Exception
	{	// if omitted, toggle option
		if(args.size()<2)
		{	flipTriangles = 1-flipTriangles;
			return;
		}
		
		String setting = doc.readStringArg(args.get(1));
		if(setting.toLowerCase().equals("yes"))
			flipTriangles = 1;
		else
			flipTriangles = 0;
	}
	

	// see if has defined path
	public boolean hasPath(String anID)
	{	Integer exists = pathIDs.get(anID);
		return exists == null ? false : true ;
	}
	
	// see if has defined keypoint
	public boolean hasKeypoint(String anID)
	{	Integer exists = keyIDs.get(anID);
		return exists == null ? false : true ;
	}
	
	// Origin #1, #2
	public void setOrigin(ArrayList<String> args) throws Exception
	{
		if(args.size()<3)
			throw new Exception("'Origin' command requires two coorinates.");
		
		originX = doc.readDoubleArg(args.get(1));
		originY = doc.readDoubleArg(args.get(2));
	}

}
	
