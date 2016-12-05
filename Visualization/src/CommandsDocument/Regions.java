/*
 * Regions.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 21 Dec 2012.
 * Copyright (c) 2012 RSAC Software. All rights reserved.
 */

import geditcom.JNFramework.JNEvaluator;

import java.util.*;

public class Regions
{
	private CmdViewer doc;
	private int inRegion;
	private StringBuffer xmlRegions;
	private String indent;
	private ArrayList<RegionPiece> pieces;
	private RegionPiece currentPiece;
	
	private static int REGION_BLOCK=1;
	private static int HOLE_BLOCK=2;
	private static int BMPREGION_BLOCK=3;
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public Regions(CmdViewer cmdDoc)
	{	// save parent CmdViewer
		doc = cmdDoc;
		pieces = null;
	}
	
	public void initRunSettings()
	{	inRegion = 0;
		xmlRegions = new StringBuffer("");
		indent = "";
		pieces=new ArrayList<RegionPiece>(20);
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
	    pieces.clear();
	    
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
			double velx=0.,vely=0.,velz=0.,thick=0.;
			String vel0X=null,vel0Y=null,vel0Z=null;
			if(args.size()<4)
		    	throw new Exception("'Region' command missing two few arguments:\n"+args);
			Object velFxn = doc.readStringOrDoubleArg(args.get(2));
			if(velFxn.getClass().equals(Double.class))
				velx = ((Double)velFxn).doubleValue();
			else
				vel0X = "      <vel0X>"+(String)velFxn+"</vel0X>\n";
			velFxn = doc.readStringOrDoubleArg(args.get(3));
			if(velFxn.getClass().equals(Double.class))
				vely = ((Double)velFxn).doubleValue();
			else
				vel0Y = "      <vel0Y>"+(String)velFxn+"</vel0Y>\n";
			
		    // start tag
		    xmlRegions.append("    <Body mat='"+matnum+"' vx='"+doc.formatDble(velx)+"' vy='"+doc.formatDble(vely)+"'");
		    
			// thickness or velocity z
			if(args.size()>4)
			{	if(doc.isMPM3D())
				{	velFxn = doc.readStringOrDoubleArg(args.get(4));
					if(velFxn.getClass().equals(Double.class))
						velz = ((Double)velFxn).doubleValue();
					else
						vel0Z = "      <vel0Z>"+(String)velFxn+"</vel0Z>\n";
				}
				else
					thick = doc.readDoubleArg(args.get(4));
			}
			if(doc.isMPM3D())
				xmlRegions.append(" vz='"+doc.formatDble(velz)+"'");
			else if(args.size()>4)
				xmlRegions.append(" thick='"+doc.formatDble(thick)+"'");
			
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
				xmlRegions.append(" "+prop+"='"+doc.formatDble(dvalue)+"'");
				nextSize += 2;
			}
			
			// end region tag
			xmlRegions.append(">\n");
			
			// initial velocities
			if(vel0X!=null) xmlRegions.append(vel0X);
			if(vel0Y!=null) xmlRegions.append(vel0Y);
			if(vel0Z!=null) xmlRegions.append(vel0Z);
		}
		
		else if(doc.isFEA())
		{	indent = "  ";
		    
			// read thickness
		    if(args.size()<3)
			    throw new Exception("'Region' command missing thickness:\n"+args);
		    double thickness = doc.readDoubleArg(args.get(2));
		    
		    // start tag
		    xmlRegions.append("  <Body mat='"+matnum+"' thick='"+doc.formatDble(thickness)+"'");
			
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

	// end current region, bmpregion, or hole
	public void EndRegion(ArrayList<String> args,String theCmd) throws Exception
	{
		// must be in region
		if(theCmd.equals("endhole"))
		{	if(inRegion!=HOLE_BLOCK)
			throw new Exception("'EndHole' command when not in a Hole:\n"+args);
		}
		else
		{	if(inRegion != REGION_BLOCK && inRegion!=BMPREGION_BLOCK)
			throw new Exception("'EndRegion' command when not in a Region or BMPRegion:\n"+args);
		}
		
		// add XML data
		if(inRegion==BMPREGION_BLOCK)
		{	xmlRegions.append(indent+"</BMP>\n");
		}
		else
		{	// go through pieces
			int numPieces = pieces.size();
			if(numPieces>0)
			{	// initialize
				currentPiece = null;
				int currentLevel = 0;
				
				int i=0;
				while(i<numPieces)
				{	RegionPiece obj = pieces.get(i);
					int level = obj.getLevel();
					
					// check if this level is less than parent level
					if(level<=currentLevel && currentPiece!=null)
					{	// climb back up the tree
						currentLevel = insertPriorElements(level);
					}
					
					switch(obj.getType())
					{	// Standard shapes
						case RegionPiece.RECT_OR_OVAL:
						case RegionPiece.SHAPE_3D:
							obj.setParent(currentPiece);
							currentLevel = level;
							currentPiece = obj;
							break;
							
						// 2D Polygons
						case RegionPiece.POLY_PT:
							obj.setParent(currentPiece);
							currentLevel = level;
							currentPiece = obj;
							String ptIndent = indent;
							for(int ii=0;ii<level;ii++) ptIndent = ptIndent+"  ";
							
							// loop until done
							while(i<numPieces-1)
							{	obj = pieces.get(i+1);
								
								// exit if new level or not a polygon
								if(obj.getLevel()!=level || obj.getType()!=RegionPiece.POLY_PT) break;
							
								// add a polypt
								currentPiece.appendXmlStart(ptIndent+obj.getXmlStart());
								i++;
							}
							
							// skip a break piece
							if(i<numPieces-1 && obj.getType()==RegionPiece.END_POLYGON) i++;
							break;
						
						case RegionPiece.END_POLYGON:
							// POLYPT_PIECE should always handle this
							break;
						
						// non shape options (must be at level 0)
						case RegionPiece.COMMAND_PIECE:
							xmlRegions.append(obj.getXmlStart());
							break;
							
						default:
							break;
					}
					
					// next object
					i++;
				}
			}
		
			// finish current elements
			if(currentPiece!=null)
				insertPriorElements(0);
			
			// end the body or holr
			if(inRegion==REGION_BLOCK)
				xmlRegions.append(indent+"</Body>\n");
			else
				xmlRegions.append(indent+"</Hole>\n");
		}
		
		// region is done
		inRegion = 0;
	}
	
	// climb back tree
	public int insertPriorElements(int level)
	{
		int newLevel = 0;
		
		while(true)
		{	RegionPiece parent = currentPiece.getParent();
			if(parent==null)
			{	xmlRegions.append(currentPiece.xmlString(indent));
			}
			else
			{	// add child
				parent.addChild(currentPiece);
			}
			
			// up to parent
			currentPiece = parent;
			
			// exit if done
			if(currentPiece==null) break;
			newLevel = currentPiece.getLevel();
			if(level>newLevel) break;
		}
		
		return newLevel;
	}
	
	// start Hole (FEA or MPM)
	public void StartHole(ArrayList<String> args) throws Exception
	{
		// verify not nested
		if(inRegion != 0)
			throw new Exception("Regions, Holes, and BMPRegions cannot be nested:\n"+args);
		    
		// activate
		inRegion = HOLE_BLOCK;
		pieces.clear();
		indent = doc.isMPM() ? "    " : "  ";
		    
		// start tag
		xmlRegions.append("    <Hole>\n");
	}
	
	// Cut Cut ... shape args command
	public void AddCutShape(ArrayList<String> args) throws Exception
	{
		// in a region
		if(inRegion==0)
			throw new Exception("Cut shape commands must be within a region");
		if(args.size()<2)
			throw new Exception("Cut copmmand with no cut shape");
		
		// assume next argument is delimited with spaces
		String [] cutArgs = args.get(1).split(" ");
		
		// find the level
		int level = 1;
		String shape = "cut";
		while(cutArgs.length>level-1)
		{	shape = cutArgs[level-1].toLowerCase();
			if(!shape.equals("cut")) break;
			level++;
		}
		
		// set command to shape name and first paremeter to last cutArg
		args.set(0,shape);
		args.set(1,cutArgs[cutArgs.length-1]);
		
		// each type
		if(shape.equals("rect"))
			AddRectOrOval(args,"Rect");
		else if(shape.equals("oval"))
			AddRectOrOval(args,"Oval");
		else if(shape.equals("polypt"))
			AddPolypoint(args);
		else if(shape.equals("box"))
			AddBox(args,"Box");
		else if(shape.equals("sphere"))
			AddBox(args,"Sphere");
		else if(shape.equals("cylinder"))
			AddBox(args,"Cylinder");
		else if(shape.equals("torus"))
			AddBox(args,"Torus");
		else
			throw new Exception("An invalid shape ('"+shape+"') in a 'Cut' command");
		
		// set last piece level
		if(!setLastPieceLevel(level))
			throw new Exception("Incorrectly nested shape commands in 'Region' or 'Hole'");

	}
	
	// bool set last piece level
	public boolean setLastPieceLevel(int level)
	{
		int numPieces = pieces.size();
		if(numPieces<1) return false;
		RegionPiece cutPiece = pieces.get(numPieces-1);
		RegionPiece prevPiece = pieces.get(numPieces-2);
		
		// set cut piece level
		if(cutPiece.getType()==RegionPiece.COMMAND_PIECE) return false;
		if(level<1) return false;
		cutPiece.setLevel(level);
		
		// check previous piece
		// Must be shape with level one less or greater
		if(prevPiece.getType()==RegionPiece.COMMAND_PIECE) return false;
		if(prevPiece.getLevel()<level-1) return false;
		
		return true;
	}
	
	// add shape for Rect #1,#2,#3,#4
	public void AddRectOrOval(ArrayList<String> args,String shape) throws Exception
	{	// times not allowed
		if(inRegion == 0 || inRegion==BMPREGION_BLOCK)
			throw new Exception("'"+shape+"' command is only allowed within a Region or Hole block:\n"+args);
		if(doc.isMPM3D())
			throw new Exception("'"+shape+"' command is only allowed within 2D MPM:\n"+args);
		
		// four numbers
		if(args.size()<5)
			throw new Exception("'"+shape+"' command has too few parameters:\n"+args);
		double xmin = doc.readDoubleArg(args.get(1));
		double xmax = doc.readDoubleArg(args.get(2));
		double ymin = doc.readDoubleArg(args.get(3));
		double ymax = doc.readDoubleArg(args.get(4));
		
		// arc angles
		double arcStart=-1.,arcEnd=0.;
		if(args.size()>5)
		{	arcStart = doc.readDoubleArg(args.get(5));
			if(args.size()>6)
			{	arcEnd = doc.readDoubleArg(args.get(6));
				if(arcStart<0. || arcStart>360. || arcEnd<arcStart)
					throw new Exception("Invalid arc angles (need 0<=start<=360 and end>=start)");
			}
			else
				throw new Exception("Has arc start angle but no end angle");
		}
		
		// start string
		StringBuffer newShape = new StringBuffer("");
		newShape.append(indent+"  <"+shape+" xmin='"+doc.formatDble(xmin)+"' xmax='"+doc.formatDble(xmax)+"'");
		newShape.append(" ymin='"+doc.formatDble(ymin)+"' ymax='"+doc.formatDble(ymax)+"'");
		
		// add piece
		RegionPiece newPiece = new RegionPiece(RegionPiece.RECT_OR_OVAL,newShape.toString(),shape);
		if(arcStart>=0.) newPiece.setArcAngles(arcStart,arcEnd);
		pieces.add(newPiece);
	}

	// add point to polygon
	public void AddPolypoint(ArrayList<String> args) throws Exception
	{	// times not allowed
		if(inRegion == 0 || inRegion==BMPREGION_BLOCK)
			throw new Exception("'PolyPt' command is only allowed within a polygon sequence:\n"+args);
		if(doc.isMPM3D())
			throw new Exception("'PolyPt' command is only allowed within 2D MPM:\n"+args);
		
		// end the polygon
		if(args.size()==1)
		{	RegionPiece newPiece = new RegionPiece(RegionPiece.END_POLYGON,"","");
			pieces.add(newPiece);
			return;
		}
		
		// needs two arguments
		if(args.size()<3)
			throw new Exception("'PolyPt' command has too few parameters:\n"+args);
		double x = doc.readDoubleArg(args.get(1));
		double y = doc.readDoubleArg(args.get(2));
		
		// add it
		String ptStr = "    <pt x='"+doc.formatDble(x)+"' y='"+doc.formatDble(y)+"'/>\n";
		RegionPiece newPiece = new RegionPiece(RegionPiece.POLY_PT,ptStr,"Polygon");
		pieces.add(newPiece);
	}
	
	// add shape for Box #1,#2,#3,#4,#5,#6 or Sphere
	// Cylinder #1,#2,#3,#4,#5,#6,#7,<#8>  or Torus
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
		StringBuffer newShape = new StringBuffer("");
		newShape.append(indent+"  <"+shape+" xmin='"+doc.formatDble(xmin)+"' xmax='"+doc.formatDble(xmax)+"'");
		newShape.append(" ymin='"+doc.formatDble(ymin)+"' ymax='"+doc.formatDble(ymax)+"'");
		newShape.append(" zmin='"+doc.formatDble(zmin)+"' zmax='"+doc.formatDble(zmax)+"'");
		if(shape.equals("Cylinder") || shape.equals("Torus"))
		{	if(args.size()<8)
				throw new Exception("'"+shape+"' command has too few parameters: "+args);
			// axis
			HashMap<String,Integer> options = new HashMap<String,Integer>(10);
			options.put("x", new Integer(1));
			options.put("y", new Integer(2));
			options.put("z", new Integer(3));
			int axis = doc.readIntOption(args.get(7),options,"shape axis");
			newShape.append(" axis='"+axis+"'");
			if(args.size()>8)
				newShape.append(" radius='"+doc.formatDble(doc.readDoubleArg(args.get(8)))+"'");
		}
		
		// add piece
		RegionPiece newPiece = new RegionPiece(RegionPiece.SHAPE_3D,newShape.toString(),shape);
		pieces.add(newPiece);
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
		{	String rotStr = indent+"  <Unrotate/>\n";
			if(inRegion == REGION_BLOCK)
			{	RegionPiece newPiece = new RegionPiece(RegionPiece.COMMAND_PIECE,rotStr,"");
				pieces.add(newPiece);
			}
			else
				xmlRegions.append(rotStr);
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
			
			// add piece
			String newShape;
			if(axis==1)
				newShape = indent+"  <RotateX>"+angle+"</RotateX>\n";
			else if(axis==2)
				newShape = indent+"  <RotateY>"+angle+"</RotateY>\n";
			else
				newShape = indent+"  <RotateZ>"+angle+"</RotateZ>\n";
			if(inRegion == REGION_BLOCK)
			{	RegionPiece newPiece = new RegionPiece(RegionPiece.COMMAND_PIECE,newShape,"");
				pieces.add(newPiece);
			}
			else
				xmlRegions.append(newShape);
			
			// next pair
			axisNum += 2;
			if(!doc.isMPM3D()) break;
		}
	}
	
	// AngularMom0 Lpz (if 2D) or AngularMom0 Lpx,Lpy,Lpz (if 3D)
	public void AddAngularMom0(ArrayList<String> args) throws Exception
	{	// check allowed
		if(inRegion == 0 || inRegion!=REGION_BLOCK)
			throw new Exception("'AngularMom0' command is only allowed within a Region block:\n"+args);
		
		// 2D or 3D
		if(args.size()<2)
			throw new Exception("'AngulaMom0' command has too few parameters: "+args);
		
		String Lp;
		Object LpFxn = doc.readStringOrDoubleArg(args.get(1));
		if(doc.isMPM3D())
		{	Lp = "      <Lp0X>"+LpFxn+"</Lp0X>\n";
			if(args.size()>2)
				Lp = Lp + "      <Lp0Y>"+doc.readStringOrDoubleArg(args.get(2))+"</Lp0Y>\n";
			if(args.size()>3)
				Lp = Lp + "      <Lp0Z>"+doc.readStringOrDoubleArg(args.get(3))+"</Lp0Z>\n";
		}
		else
			Lp = "      <Lp0Z>"+LpFxn+"</Lp0Z>\n";
		xmlRegions.append(Lp);
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
	    pieces.clear();
	    
	    // read path and width
	    if(args.size()<2)
		    throw new Exception("'BMP' has too few parameters:\n"+args);
	    
	    String filePath = doc.readStringArg(args.get(1));
	    // optional width
	    double width = -1.e9;
	    if(args.size()>2) width = doc.readDoubleArg(args.get(2));
	    
	    // optional height
	    double height = -1.e9;
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
	    xmlRegions.append(indent+"<BMP name='"+filePath+"'");
	    
	    // optional attrributes
	    if(width>-1.e8) xmlRegions.append(" width='"+doc.formatDble(width)+"'");
	    if(height>-1.e8) xmlRegions.append(" height='"+doc.formatDble(height)+"'");
	    if(anglesPath!=null)
	    {	if(anglesPath.length()>0)
	    		xmlRegions.append(" angles='"+anglesPath+"'");
	    }
	    
	    // end it
	    xmlRegions.append(">\n");	    
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
	    xmlRegions.append(indent+"  <Origin x='"+doc.formatDble(xo)+"' y='"+doc.formatDble(yo)+"'");
	    if(doc.isMPM3D()) xmlRegions.append(" z='"+doc.formatDble(zo)+"'");
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
					xmlRegions.append(indent+"    <Thickness>"+doc.formatDble(value)+"</Thickness>\n");
				else if(prop.equals("temp"))
					xmlRegions.append(indent+"    <Temperature>"+doc.formatDble(value)+"</Temperature>\n");
				else if(prop.equals("angle"))
					xmlRegions.append(indent+"    <Angle>"+doc.formatDble(value)+"</Angle>\n");
				else if(prop.equals("conc") && doc.isMPM())
					xmlRegions.append(indent+"    <Concentration>"+doc.formatDble(value)+"</Concentration>\n");
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
				if(hasVx) xmlRegions.append(" x='"+doc.formatDble(vx)+"'");
				if(hasVy) xmlRegions.append(" y='"+doc.formatDble(vy)+"'");
				if(hasVz && doc.isMPM3D()) xmlRegions.append(" z='"+doc.formatDble(vz)+"'");
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
			xmlRegions.append(" minAngle='"+doc.formatDble(angle1)+"' maxAngle='"+doc.formatDble(angle2)+"'/>\n");
		}
	}
	
	// insert XML data in pieces or xml data
	public void AddXML(String rawXML)
	{	xmlRegions.append(rawXML);
	}

	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------

	// combine to xml data
	public String toXMLString() { return xmlRegions.toString(); }
	
	public boolean isInBMPRegion() { return inRegion==BMPREGION_BLOCK; }

	public boolean isInRegion() { return inRegion!=0; }
}
