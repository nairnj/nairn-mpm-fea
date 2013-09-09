/*
 * Regions.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 21 Dec 2012.
 * Copyright (c) 2012 RSAC Software. All rights reserved.
 */

import java.util.*;

public class Regions
{
	private CmdViewer doc;
	private int inRegion;
	private StringBuffer xmlRegions;
	private String indent;
	private boolean inPoly;
	
	private static int REGION_BLOCK=1;
	private static int HOLE_BLOCK=1;
	//private static int BMPREGION_BLOCK=1;
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public Regions(CmdViewer cmdDoc)
	{	// save parent CmdViewer
		doc = cmdDoc;
	}
	
	public void initRunSettings()
	{	inRegion = 0;
		xmlRegions = new StringBuffer("");
		indent = "";
		inPoly = false;
	}
	
	//----------------------------------------------------------------------------
	// Methods
	//----------------------------------------------------------------------------
	
	// start FEA Region #1,#2,<#3>
	// 		#1 is material, #2 is thickness, #3 is material angle function
	// or MPM Region #1,#2,#3,#4,(#5,#6 pairs)
	//   #1 is material, (#2,#3)=(vx,vy), #5=thickness or vz
	public void StartRegion(ArrayList<String> args) throws Exception
	{
		// verify not nested
		if(inRegion != 0)
		    throw new Exception("Regions, Holes, and BMPRegions cannot be nested:\n"+args);
		
	    // activate
	    inRegion = REGION_BLOCK;
	    
	    // read material by ID
	    if(args.size()<2)
		    throw new Exception("'Region' command missing material ID:\n"+args);
	    String matID = doc.readStringArg(args.get(1));
	    int matnum = doc.mats.getMatID(matID);
		if(matnum<=0)
			throw new Exception("'Region' command has unknown material ID:\n"+args);
		
	    // MPM or FEA
		if(doc.isMPM())
		{	indent = "    ";
		
			// read two velocities
			if(args.size()<4)
		    	throw new Exception("'Region' command missing two few arguments:\n"+args);
			double velx = doc.readDoubleArg(args.get(2));
			double vely = doc.readDoubleArg(args.get(3));
			
		    // start tag
		    xmlRegions.append("    <Body mat='"+matnum+"' vx='"+velx+"' vy='"+vely+"'");
		    
			// thickness or velocity z
		    double thick = 0.;
			if(args.size()>4)
			{	thick = doc.readDoubleArg(args.get(4));
			}
			if(doc.isMPM3D())
				xmlRegions.append(" vz='"+thick+"'");
			else if(args.size()>4)
				xmlRegions.append(" thick='"+thick+"'");
			
			// extra pairs (angle, temp, conc)
			int nextSize = 5;
			while(args.size()>nextSize)
			{	// get property
				String prop = doc.readStringArg(args.get(nextSize)).toLowerCase();
				if(!prop.equals("angle") && !prop.equals("temp") && !prop.equals("conc"))
					throw new Exception("Region '"+prop+"' property not recogonized:\n"+args);
				
				// need next one
				if(args.size()<nextSize+2)
					throw new Exception("Region '"+prop+"' property missing a value:\n"+args);
				double dvalue = doc.readDoubleArg(args.get(nextSize+1));
				
				// add it and increment
				xmlRegions.append(" "+prop+"='"+dvalue+"'");
				nextSize += 2;
			}
			
			// end region tag
			xmlRegions.append(">\n");
		}
		
		else if(doc.isFEA())
		{	indent = "  ";
		    
			// read thickness
		    if(args.size()<3)
			    throw new Exception("'Region' command missing thickness:\n"+args);
		    double thickness = doc.readDoubleArg(args.get(2));
		    
		    // start tag
		    xmlRegions.append("  <Body mat='"+matnum+"' thick='"+thickness+"'");
			
			// read angle
		    if(args.size()>3)
		    {	String angle = doc.readStringArg(args.get(3));
		    	xmlRegions.append(" angle='"+angle+"'");
		    }
		    
		    // end the line
		    xmlRegions.append(">\n");
		}
		
		else
			throw new Exception("'Region' command not allowed before analysis type is set:\n"+args);
	}

	// end current region
	public void EndRegion(ArrayList<String> args) throws Exception
	{
		// must be in region
		if(inRegion != REGION_BLOCK)
			throw new Exception("'EndRegion' command when not in a region:\n"+args);
		
		// active polygon
		if(inPoly == true)
		{	xmlRegions.append(indent+"  </Polygon>\n");
			inPoly = false;
		}
		
		// end the body
		xmlRegions.append(indent+"</Body>\n\n");
		inRegion = 0;
	}
	
	// start Hole (FEA or MPM)
	public void StartHole(ArrayList<String> args) throws Exception
	{
		// verify not nested
		if(inRegion != 0)
			throw new Exception("Regions, Holes, and BMPRegions cannot be nested:\n"+args);
		    
		// activate
		inRegion = HOLE_BLOCK;
		indent = doc.isMPM() ? "    " : "  ";
		    
		// start tag
		xmlRegions.append("  <Hole>\n");
	}

	// end current region
	public void EndHole(ArrayList<String> args) throws Exception
	{
		// must be in region
		if(inRegion != HOLE_BLOCK)
			throw new Exception("'EndHole' command when not in a hole:\n"+args);
		
		// active polygon
		if(inPoly == true)
		{	xmlRegions.append(indent+"  </Polygon>\n");
			inPoly = false;
		}
		
		// end the body
		xmlRegions.append(indent+"</Hole>\n\n");
		inRegion = 0;
	}
	
	// add shape for Rect #1,#2,#3,#4
	public void AddRectOrOval(ArrayList<String> args,String shape) throws Exception
	{	// times not allowed
		if(inRegion == 0)
			throw new Exception("'"+shape+"' command is only allowed within a region block:\n"+args);
		if(inPoly == true)
			throw new Exception("'"+shape+"' command is not allowed in a polygon block:\n"+args);
		if(doc.isMPM3D())
			throw new Exception("'PolyPt' command is only allowed within 2D MPM:\n"+args);
		
		// four numbers
		if(args.size()<5)
			throw new Exception("'"+shape+"' command has too few parameters:\n"+args);
		double xmin = doc.readDoubleArg(args.get(1));
		double xmax = doc.readDoubleArg(args.get(2));
		double ymin = doc.readDoubleArg(args.get(3));
		double ymax = doc.readDoubleArg(args.get(4));
		
		// add it
		xmlRegions.append(indent+"  <"+shape+" units='mm' xmin='"+xmin+"' xmax='"+xmax+"'");
		xmlRegions.append(" ymin='"+ymin+"' ymax='"+ymax+"'/>\n");
	}

	// add point to polygon
	public void AddPolypoint(ArrayList<String> args) throws Exception
	{	// times not allowed
		if(inRegion == 0)
			throw new Exception("'PolyPt' command is only allowed within a polygon sequence:\n"+args);
		if(inPoly==false && args.size()<2)
			throw new Exception("Empty 'PolyPt' command only allowed in a polygon sequence:\n"+args);
		if(doc.isMPM3D())
			throw new Exception("'PolyPt' command is only allowed within 2D MPM:\n"+args);
		
		// start a polygon
		if(inPoly == false)
		{	xmlRegions.append(indent+"  <Polygon>\n");
			inPoly = true;
		}
		
		// end the polygon
		if(args.size()==1)
		{	xmlRegions.append(indent+"  </Polygon>\n");
			inPoly = false;
			return;
		}
		
		// needs two arguments
		if(args.size()<3)
			throw new Exception("'PolyPt' command has too few parameters:\n"+args);
		double x = doc.readDoubleArg(args.get(1));
		double y = doc.readDoubleArg(args.get(2));
		
		// add it
		xmlRegions.append(indent+"    <pt units='mm' x='"+x+"' y='"+y+"'/>\n");
	}

	// add shape for Rect #1,#2,#3,#4
	public void AddBox(ArrayList<String> args) throws Exception
	{	// times not allowed
		if(inRegion == 0)
			throw new Exception("'Box' command is only allowed within a region block:\n"+args);
		if(!doc.isMPM3D())
			throw new Exception("'Box' command is only allowed within 3D MPM:\n"+args);
		
		// four numbers
		if(args.size()<7)
			throw new Exception("'Box' command has too few parameters: "+args);
		double xmin = doc.readDoubleArg(args.get(1));
		double xmax = doc.readDoubleArg(args.get(2));
		double ymin = doc.readDoubleArg(args.get(3));
		double ymax = doc.readDoubleArg(args.get(4));
		double zmin = doc.readDoubleArg(args.get(5));
		double zmax = doc.readDoubleArg(args.get(6));
		
		// add it
		xmlRegions.append(indent+"  <Box units='mm' xmin='"+xmin+"' xmax='"+xmax+"'");
		xmlRegions.append(" ymin='"+ymin+"' ymax='"+ymax+"'");
		xmlRegions.append(" zmin='"+zmin+"' zmax='"+zmax+"'/>\n");
	}
	
	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------

	// combine to xml data
	public String toXMLString() { return xmlRegions.toString(); }

}
