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
	private static int HOLE_BLOCK=2;
	private static int BMPREGION_BLOCK=3;
	
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

	// start FEA Region #1,#2,<#3>
	// 		#1 is material, #2 is thickness, #3 is material angle function
	// or MPM Region #1,#2,#3,#4,(#5,#6 pairs)
	//   #1 is material, (#2,#3)=(vx,vy), #5=thickness or vz
	public void StartBMPRegion(ArrayList<String> args) throws Exception
	{
		// verify not nested
		if(inRegion != 0)
		    throw new Exception("Regions, Holes, and BMPRegions cannot be nested:\n"+args);
		
	    // activate
	    inRegion = BMPREGION_BLOCK;
	    
	    // read path and width
	    if(args.size()<3)
		    throw new Exception("'BMP' has too few parameters:\n"+args);
	    
	    String filePath = doc.readStringArg(args.get(1));
	    double width = doc.readDoubleArg(args.get(2));
	    
	    // optional height
	    double height = -1.e8;
	    if(args.size()>3) height = doc.readDoubleArg(args.get(3));
	    
	    // optional angles path
	    String anglesPath = null;
	    if(args.size()>4) anglesPath = doc.readStringArg(args.get(4));
	    
	    // set indent
	    if(doc.isMPM())
			indent = "    ";
	    else if(doc.isFEA())
	    	indent = "  ";
		else
			throw new Exception("'BMPRegion' command not allowed before analysis type is set:\n"+args);

	    // create the command
	    xmlRegions.append(indent+"<BMP name='"+filePath+"' width='"+width+"'");
	    
	    // optional attrributes
	    if(height>-0.99e8) xmlRegions.append(" height='"+height+"'");
	    if(anglesPath!=null)
	    {	if(anglesPath.length()>0)
	    		xmlRegions.append(" angles='"+anglesPath+"'");
	    }
	    
	    // end it
	    xmlRegions.append(">\n");
	    
}
	
	// end current region
	public void EndRegion(ArrayList<String> args) throws Exception
	{
		// must be in region
		if(inRegion != REGION_BLOCK && inRegion!=BMPREGION_BLOCK)
			throw new Exception("'EndRegion' command when not in a Region or BMPRegion:\n"+args);
		
		if(inRegion==BMPREGION_BLOCK)
		{	xmlRegions.append(indent+"</BMP>\n\n");
		}
		else
		{	// active polygon
			if(inPoly == true)
			{	xmlRegions.append(indent+"  </Polygon>\n");
				inPoly = false;
			}
		
			// end the body
			xmlRegions.append(indent+"</Body>\n\n");
		}
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
		if(inRegion == 0 || inRegion==BMPREGION_BLOCK)
			throw new Exception("'"+shape+"' command is only allowed within a Region or Hole block:\n"+args);
		if(doc.isMPM3D())
			throw new Exception("'"+shape+"' command is only allowed within 2D MPM:\n"+args);
		if(inPoly == true)
		{	xmlRegions.append(indent+"  </Polygon>\n");
			inPoly = false;
		}
		
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
		if(inRegion == 0 || inRegion==BMPREGION_BLOCK)
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
	
	// Origin #1,#2,<#3>,<#4>
	public void setOrigin(ArrayList<String> args) throws Exception
	{
	    // must have at least two arguments
	    if(args.size()<3)
			throw new Exception("'Origin' command requires at leat two coordinates: "+args);
	    
	    double xo = doc.readDoubleArg(args.get(1));
	    double yo = doc.readDoubleArg(args.get(2));
	    
	    double zo = 0;
	    if(args.size()>3) zo = doc.readDoubleArg(args.get(3));
	    
	    String flip = null;
	    if(args.size()>4) flip = doc.readStringArg(args.get(4)).toLowerCase();
	    
	    // add the command
	    xmlRegions.append(indent+"  <Origin units='mm' x='"+xo+"' y='"+yo+"'");
	    if(doc.isMPM3D()) xmlRegions.append(" z='"+zo+"'");
	    if(flip!=null) xmlRegions.append(" flipped='"+flip+"'");
	    xmlRegions.append("/>\n");
	}

	// Intensity #1,#2,#3,<#4,#5>,...
	// Intensity "angles",#2,#3,#4,#5
	public void AddIntensity(ArrayList<String> args) throws Exception
	{
		// verify not nested
		if(inRegion != BMPREGION_BLOCK)
		    throw new Exception("'Intensity' command only allowed in a BMPRegion block:\n"+args);
		
	    // read material by ID
	    if(args.size()<2)
		    throw new Exception("'Intensity' command has no arguments:\n"+args);
	    String matID = doc.readStringArg(args.get(1));
	    int matnum = doc.mats.getMatID(matID);
	    
		// get gray range
		if(args.size()<4)
			throw new Exception("'Intensity' command is missing gray scale range parameters:\n"+args);
		int gray1 = doc.readIntArg(args.get(2));
		int gray2 = doc.readIntArg(args.get(3));
		
	    // map materials or angles
		if(matnum>0)
		{	// get gray range
			if(args.size()<4)
				 throw new Exception("'Intensity' command has no arguments:\n"+args);
			
			// the command
			xmlRegions.append(indent+"  <Intensity mat='"+matnum+"' imin='"+gray1+"' imax='"+gray2+"'>\n");
			
			// process each pair
			int propNum=5;
			double vx=0.,vy=0.,vz=0.;
			boolean hasVx=false,hasVy=false,hasVz=false;
			while(args.size()>propNum)
			{	String prop = doc.readStringArg(args.get(propNum-1)).toLowerCase();
				double value = doc.readDoubleArg(args.get(propNum));
				propNum += 2;
				
				if(prop.equals("thick"))
					xmlRegions.append(indent+"    <Thickness>"+value+"</Thickness>\n");
				else if(prop.equals("temp"))
					xmlRegions.append(indent+"    <Temperature>"+value+"</Temperature>\n");
				else if(prop.equals("angle"))
					xmlRegions.append(indent+"    <Angle>"+value+"</Angle>\n");
				else if(prop.equals("conc") && doc.isMPM())
					xmlRegions.append(indent+"    <Concentration>"+value+"</Concentration>\n");
				else if(prop.equals("vx") && doc.isMPM())
				{	vx = value;
					hasVx = true;
				}
				else if(prop.equals("vy") && doc.isMPM())
				{	vy = value;
					hasVy = true;
				}
				else if(prop.equals("vz") && doc.isMPM())
				{	vz = value;
					hasVz = true;
				}
				else
					throw new Exception("'Intensity' command invalid property options:\n"+args);
			}
			
			// velocity
			if(hasVx || hasVy || hasVz)
			{	xmlRegions.append(indent+"    <vel");
				if(hasVx) xmlRegions.append(" x='"+vx+"'");
				if(hasVy) xmlRegions.append(" y='"+vy+"'");
				if(hasVz && doc.isMPM3D()) xmlRegions.append(" z='"+vz+"'");
				xmlRegions.append("/>\n");
			}
			
			// finish up
			xmlRegions.append(indent+"  </Intensity>\n");
		}
		else
		{	// get angle range
			if(args.size()<6)
				throw new Exception("'Intensity' command is missing angle range parameters:\n"+args);
			double angle1 = doc.readDoubleArg(args.get(4));
			double angle2 = doc.readDoubleArg(args.get(5));
			
			// the commmad
			xmlRegions.append(indent+"  <Intensity imin='"+gray1+"' imax='"+gray2+"'");
			xmlRegions.append(" minAngle='"+angle1+"' maxAngle='"+angle2+"'/>\n");
		}
	}
	
	// Rotate #1,#2,<#3,#4>,<#5,$6>
	public void AddRotate(ArrayList<String> args) throws Exception
	{	// times not allowed
		if(inRegion==0 || inRegion==HOLE_BLOCK)
			throw new Exception("'Rotate' command is only allowed within a Region or BMPRegion block:\n"+args);
		
		// check for reset
		if(args.size()<2)
			throw new Exception("'Rotate' command has too few parameters:\n"+args);
		
		// Is is "reset"
		String reset = doc.readStringArg(args.get(1)).toLowerCase();
		if(reset.equals("reset"))
		{	xmlRegions.append(indent+"  <Unrotate/>\n");
			return;
		}
		
		// need at least one axis
		if(args.size()<3)
			throw new Exception("'Rotate' command has too few parameters:\n"+args);
		
		// up to three pairs (3D only)
		int axisNum=2;
		while(args.size()>axisNum && axisNum<8)
		{	// get axis
			int axis=0;
			Object axisArg = doc.readStringOrDoubleArg(args.get(axisNum-1));
			if(axisArg.getClass().equals(Double.class))
				axis = ((Double)axisArg).intValue();
			else if(((String)axisArg).toLowerCase().equals("x"))
				axis = 1;
			else if(((String)axisArg).toLowerCase().equals("y"))
				axis = 2;
			else if(((String)axisArg).toLowerCase().equals("z"))
				axis = 3;
			if(axis<1 || axis>3)
				throw new Exception("'Rotate' command has invalid rotation axis:\n"+args);
			if(axis!=3 && !doc.isMPM3D())
				throw new Exception("'Rotate' command axis must be z axis in 2D simulations:\n"+args);
			
			// get angle (can be function)
			String angle = doc.readStringArg(args.get(axisNum));
			
			if(axis==1)
				xmlRegions.append(indent+"  <RotateX>"+angle+"</RotateX>\n");
			else if(axis==2)
				xmlRegions.append(indent+"  <RotateY>"+angle+"</RotateY>\n");
			else
				xmlRegions.append(indent+"  <RotateZ>"+angle+"</RotateZ>\n");
			
			// next pair
			axisNum += 2;
			if(!doc.isMPM3D()) break;
		}
	}
	
	// add shape for Box #1,#2,#3,#4,#5,#6 or Sphere
	// Cylinder #1,#2,#3,#4,#5,#6,#7,#8
	public void AddBox(ArrayList<String> args,String shape) throws Exception
	{	// times not allowed
		if(inRegion == 0 || inRegion==BMPREGION_BLOCK)
			throw new Exception("'"+shape+"' command is only allowed within a Region or Hole block:\n"+args);
		if(!doc.isMPM3D())
			throw new Exception("'"+shape+"' command is only allowed within 3D MPM:\n"+args);
		
		// four numbers
		if(args.size()<7)
			throw new Exception("'"+shape+"' command has too few parameters: "+args);
		double xmin = doc.readDoubleArg(args.get(1));
		double xmax = doc.readDoubleArg(args.get(2));
		double ymin = doc.readDoubleArg(args.get(3));
		double ymax = doc.readDoubleArg(args.get(4));
		double zmin = doc.readDoubleArg(args.get(5));
		double zmax = doc.readDoubleArg(args.get(6));
		
		// add it
		xmlRegions.append(indent+"  <"+shape+" units='mm' xmin='"+xmin+"' xmax='"+xmax+"'");
		xmlRegions.append(" ymin='"+ymin+"' ymax='"+ymax+"'");
		xmlRegions.append(" zmin='"+zmin+"' zmax='"+zmax+"'");
		if(shape.equals("Cylinder"))
		{	if(args.size()<8)
				throw new Exception("'"+shape+"' command has too few parameters: "+args);
			xmlRegions.append(" axis='"+doc.readIntArg(args.get(7))+"'");
			if(args.size()>8)
				xmlRegions.append(" radius='"+doc.readDoubleArg(args.get(8))+"'");
		}
		xmlRegions.append("/>\n");
	}
	
	// insert XML
	public void AddXML(String rawXML)
	{	xmlRegions.append(rawXML);
	}

	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------

	// combine to xml data
	public String toXMLString() { return xmlRegions.toString(); }
	
	public boolean isInBMPRegion() { return inRegion==BMPREGION_BLOCK; }

}
