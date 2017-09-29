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
	public ArrayList<RegionPiece> pieces;
	
	public static final int GRID_BC = 1;
	public static final int MOVELINE_BC = 2;
	public static final int MOVEARC_BC = 3;
	public static final int MOVEBOX_BC = 4;

	public static final int PARTICLE_BC = 10;
	public static final int LOADLINE_BC = 11;
	public static final int LOADARC_BC = 12;
	public static final int LOADRECT_BC = 14;
	public static final int LOADBOX_BC = 13;

	public static final int ADD_TEMPERATURE = 1;
	public static final int ADD_CONCENTRATION = 2;

	// ----------------------------------------------------------------------------
	// Initialize
	// ----------------------------------------------------------------------------

	public MPMGridBCs(CmdViewer cmdDoc)
	{	// save parent CmdViewer
		doc = cmdDoc;
	}

	public void initRunSettings()
	{	inBC = 0;
		xmlbcs = new StringBuffer("");
		boundaryID = 0;
		pieces=new ArrayList<RegionPiece>(20);
	}

	// ----------------------------------------------------------------------------
	// Methods
	// ----------------------------------------------------------------------------

	// start grid BC line or arc
	// MoveLine x1,y1,x2,y2,(tolerance)
	// MoveArc x1,y1,x2,y2,start,end,(tolerance)
	public void StartMoveLine(ArrayList<String> args, int theType) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);

		// verify not nested
		if (inBC != 0 && theType < PARTICLE_BC)
		{	if(theType==GRID_BC)
				throw new Exception("GridBC blocks cannot be nested:\n" + args);
			else
				throw new Exception("MoveLine, MoveArc, and MoveBox cannot be nested:\n" + args);
		}
		else if (doc.mpmParticleBCs.getInBC() != 0 && theType >= PARTICLE_BC)
		{	if(theType==PARTICLE_BC)
				throw new Exception("Particle blocks cannot be nested:\n" + args);
			else
			{	throw new Exception("LoadLine, LoadArc, LoadRect, and LoadBox cannot be nested:\n"
							+ args);
			}
		}
		
		if(theType==GRID_BC)
		{	// no arguments
			if(args.size()>2)
				throw new Exception("GridBC command should have not parameters:\n" + args);
			bcAttrs = "<BCShape>\n";
			bcSettings = new StringBuffer("");
			inBC = GRID_BC;
			bcCmd = "BCShape";
			pieces.clear();
		}
		
		else if(theType==PARTICLE_BC)
		{	// no arguments
			if(args.size()>2)
				throw new Exception("ParticleBC command should have not parameters:\n" + args);
			String theAttrs = "<BCShape>\n";
			doc.mpmParticleBCs.SetLoadLine(theAttrs,PARTICLE_BC,"BCShape");
		}

		else
		{	// needs at least 4 arguments
			if (((theType == MOVELINE_BC || theType == LOADLINE_BC) && args.size() < 5)
					|| ((theType == MOVEARC_BC || theType == LOADARC_BC) && args.size() < 7))
			{	throw new Exception("'" + args.get(0)
						+ "' has too few parameters:\n" + args);
			}
	
			// get x1,y1,x2,y2
			double x1 = doc.readDoubleArg(args.get(1));
			double y1 = doc.readDoubleArg(args.get(2));
			double x2 = doc.readDoubleArg(args.get(3));
			double y2 = doc.readDoubleArg(args.get(4));
	
			double tolerance = -1., startAng = 0., endAng = 0.;
			if (theType == MOVELINE_BC || theType == LOADLINE_BC)
			{ 	// get optional tolerance
				if (args.size() > 5)
					tolerance = doc.readDoubleArg(args.get(5));
				if (theType == MOVELINE_BC)
				{	bcAttrs = "<BCLine x1='" + doc.formatDble(x1) + "' y1='" + doc.formatDble(y1) +
								"' x2='" + doc.formatDble(x2) + "' y2='" + doc.formatDble(y2) + "'";
					if (tolerance > 0.)
						bcAttrs = bcAttrs + " tolerance='" + doc.formatDble(tolerance) + "'>\n";
					else
						bcAttrs = bcAttrs + ">\n";
					bcSettings = new StringBuffer("");
					inBC = MOVELINE_BC;
					bcCmd = "BCLine";
				}
				else
				{	String theAttrs = "<BCLine x1='" + doc.formatDble(x1) + "' y1='" + doc.formatDble(y1)
							+ "' x2='" + doc.formatDble(x2) + "' y2='" + doc.formatDble(y2) + "'";
					if (tolerance > 0.)
						theAttrs = theAttrs + " tolerance='" + tolerance + "'>\n";
					else
						theAttrs = theAttrs + ">\n";
					doc.mpmParticleBCs.SetLoadLine(theAttrs,LOADLINE_BC,"BCLine");
				}
			}
			else
			{	// angles
				startAng = doc.readDoubleArg(args.get(5));
				endAng = doc.readDoubleArg(args.get(6));
				if (args.size() > 7)
					tolerance = doc.readDoubleArg(args.get(7));
				if (theType == MOVEARC_BC)
				{	bcAttrs = "<BCArc x1='" + doc.formatDble(x1) + "' y1='" + doc.formatDble(y1)
							+ "' x2='" + doc.formatDble(x2) + "' y2='" + doc.formatDble(y2)
							+ "' start='" + doc.formatDble(startAng) + "' end='" + doc.formatDble(endAng) + "'";
					if (tolerance > 0.)
						bcAttrs = bcAttrs + " tolerance='" + doc.formatDble(tolerance) + "'>\n";
					else
						bcAttrs = bcAttrs + ">\n";
					bcSettings = new StringBuffer("");
					inBC = MOVEARC_BC;
					bcCmd = "BCArc";
				}
				else
				{	String theAttrs = "<BCArc x1='" + doc.formatDble(x1) + "' y1='" + doc.formatDble(y1)
							+ "' x2='" + doc.formatDble(x2) + "' y2='" + doc.formatDble(y2)
							+ "' start='" + doc.formatDble(startAng) + "' end='" + doc.formatDble(endAng) + "'";
					if (tolerance > 0.)
						theAttrs = theAttrs + " tolerance='" + doc.formatDble(tolerance) + "'>\n";
					else
						theAttrs = theAttrs + ">\n";
					doc.mpmParticleBCs.SetLoadLine(theAttrs,LOADARC_BC, "BCArc");
				}
			}
		}
	}

	// start grid BC line
	public void StartMoveBox(ArrayList<String> args, int theType) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);
		if (!doc.isMPM3D())
			throw new Exception("MoveBox and LoadBox only allowed in 3D MPM:\n"
					+ args);

		// verify not nested
		if (inBC != 0 && theType == MOVEBOX_BC)
			throw new Exception(
					"MoveLine, MoveArc, and MoveBox cannot be nested:\n" + args);
		else if (doc.mpmParticleBCs.getInBC() != 0 && theType == LOADBOX_BC)
			throw new Exception(
					"LoadLine, LoadArc, LoadRect, and LoadBox cannot be nested:\n"
							+ args);

		// needs at least 6 arguments
		if (args.size() < 7)
			throw new Exception("'" + args.get(0)
					+ "' has too few parameters:\n" + args);

		// get x1,y1,x2,y2,z1,z2
		double x1 = doc.readDoubleArg(args.get(1));
		double y1 = doc.readDoubleArg(args.get(2));
		double z1 = doc.readDoubleArg(args.get(3));
		double x2 = doc.readDoubleArg(args.get(4));
		double y2 = doc.readDoubleArg(args.get(5));
		double z2 = doc.readDoubleArg(args.get(6));

		// get optional axis
		int axis = -1;
		if (args.size() > 7)
		{	HashMap<String, Integer> options = new HashMap<String, Integer>(3);
			options.put("x", new Integer(1));
			options.put("y", new Integer(2));
			options.put("z", new Integer(3));
			axis = doc.readIntOption(args.get(7), options, "Cylinder axis");
			if (axis < 1 || axis > 3)
				throw new Exception("'" + args.get(0)
						+ "' cylinder axis is not valid:\n" + args);
		}

		if (theType == MOVEBOX_BC)
		{	bcAttrs = "<BCBox xmin='" + doc.formatDble(x1) + "' ymin='" + doc.formatDble(y1)
					+ "' zmin='" + doc.formatDble(z1) + "' xmax='" + doc.formatDble(x2)
					+ "' ymax='" + doc.formatDble(y2) + "' zmax='" + doc.formatDble(z2) + "'";
			if (axis > 0)
				bcAttrs = bcAttrs + " axis='" + axis + "'>\n";
			else
				bcAttrs = bcAttrs + ">\n";
			bcSettings = new StringBuffer("");
			inBC = MOVEBOX_BC;
			bcCmd = "BCBox";
		}
		else
		{	String theAttrs = "<BCBox xmin='" + doc.formatDble(x1) + "' ymin='" + doc.formatDble(y1)
					+ "' zmin='" + doc.formatDble(z1) + "' xmax='" + doc.formatDble(x2)
					+ "' ymax='" + doc.formatDble(y2) + "' zmax='" + doc.formatDble(z2) + "'";
			if (axis > 0)
				theAttrs = theAttrs + " axis='" + axis + "'>\n";
			else
				theAttrs = theAttrs + ">\n";
			doc.mpmParticleBCs.SetLoadLine(theAttrs,LOADBOX_BC,"BCBox");
		}
	}

	// MoveLine, MoveArc, or MoveBox is done
	public void EndMoveBlock(ArrayList<String> args, int endType) throws Exception
	{	if (inBC != endType)
		{	throw new Exception("'" + args.get(0)
					+ "' does not match current boundary condition block:\n"
					+ args);
		}
	
		// get shape
		StringBuffer shapeXML = new StringBuffer("");
		if(inBC == GRID_BC)
		{	doc.regions.compilePieces(pieces, shapeXML);
		}

		// append block
		xmlbcs.append("    " + bcAttrs + shapeXML + bcSettings + "    </" + bcCmd + ">\n");

		inBC = 0;
	}

	// add velocity condition
	// Velocity (x or y or z),type,<arg1>,<arg2>
	// Velocity (skewed),type,arg1,arg2,angle1,<angle2>
	public void AddVelocity(ArrayList<String> args) throws Exception
	{	// MPM only
		doc.requiresMPM(args);

		if (inBC == 0)
		{	throw new Exception(
					"'Velocity' command must by in 'MoveLine', 'MoveArc', or 'MoveBox' block:\n"
							+ args);
		}

		// always needs #1 and #2
		if (args.size() < 3)
			throw new Exception("'Velocity' has too few parameters:\n" + args);

		// read direction
		HashMap<String, Integer> options = new HashMap<String, Integer>(8);
		options.put("x", new Integer(1));
		options.put("y", new Integer(2));
		options.put("z", new Integer(3));
		options.put("skewxy", new Integer(12));
		options.put("skewrz", new Integer(12));
		options.put("skewxz", new Integer(13));
		options.put("skewyz", new Integer(23));
		options.put("skewxyz", new Integer(123));
		int dof = doc.readIntOption(args.get(1), options, "Velocity direction");

		if (dof > 10 && args.size() < 6)
		{	throw new Exception("Skewed 'Velocity' has too few parameters:\n"
					+ args);
		}

		// read style
		options = new HashMap<String, Integer>(5);
		options.put("constant", new Integer(1));
		options.put("linear", new Integer(2));
		options.put("sine", new Integer(3));
		options.put("cosine", new Integer(4));
		options.put("function", new Integer(6));
		int style = doc.readIntOption(args.get(2), options, "Velocity style");

		// all need a arg1 except constant
		if (style != 1 && args.size() < 4)
			throw new Exception("'Velocity' has too few parameters:\n" + args);

		// read arg1 and arg2
		double arg1 = 0., arg2 = 0.;
		String function = null;
		boolean hasArg2 = false;

		// arg1
		if (args.size() > 3)
		{	if (style == 6)
				function = doc.readStringArg(args.get(3));
			else
				arg1 = doc.readDoubleArg(args.get(3));
		}

		// arg2
		if (args.size() > 4)
		{	hasArg2 = true;
			arg2 = doc.readDoubleArg(args.get(4));
		}

		// angles
		double angle1 = 0., angle2 = 0.;
		if (dof > 10)
		{	angle1 = doc.readDoubleArg(args.get(5));
			if (args.size() > 6)
				angle2 = doc.readDoubleArg(args.get(6));
		}

		// add to xml
		bcSettings.append("      <DisBC dir='" + dof + "' style='" + style + "'");
		if (style == 6)
			bcSettings.append(" function='" + function + "'");
		else
			bcSettings.append(" vel='" + doc.formatDble(arg1) + "'");
		if (hasArg2)
			bcSettings.append(" time='" + doc.formatDble(arg2) + "'");

		if (dof > 10)
			bcSettings.append(" angle='" + doc.formatDble(angle1) + "'");
		if (dof > 100)
			bcSettings.append(" angle2='" + doc.formatDble(angle2) + "'");

		if (boundaryID != 0)
			bcSettings.append(" id='" + boundaryID + "'");
		bcSettings.append("/>\n");
	}

	// add velocity condition
	public void AddTempConc(ArrayList<String> args, int theType) throws Exception
	{	// MPM only
		doc.requiresMPM(args);

		if (inBC == 0)
		{	throw new Exception("'"+ args.get(0)
							+ "' command must by in 'MoveLine', 'MoveArc', or 'MoveBox' block:\n"
							+ args);
		}

		// always needs #1 and #2
		if (args.size() < 3)
		{	throw new Exception("'" + args.get(0)
					+ "' has too few parameters:\n" + args);
		}

		// read style
		HashMap<String, Integer> options = new HashMap<String, Integer>(5);
		options.put("constant", new Integer(1));
		options.put("linear", new Integer(2));
		options.put("sine", new Integer(3));
		options.put("cosine", new Integer(4));
		options.put("function", new Integer(6));
		int style = doc.readIntOption(args.get(1), options,
				"Temperature or concentration style");

		// read arg1 and arg2
		double arg1 = 0., arg2 = 0.;
		String function = null;
		boolean hasArg2 = false;

		// arg1
		if (args.size() > 2)
		{	if (style == 6)
				function = doc.readStringArg(args.get(2));
			else
				arg1 = doc.readDoubleArg(args.get(2));
		}

		// arg2
		if (args.size() > 3)
		{	hasArg2 = true;
			arg2 = doc.readDoubleArg(args.get(3));
		}

		// add to xml
		if (theType == ADD_TEMPERATURE)
			bcSettings.append("      <TempBC");
		else
			bcSettings.append("      <ConcBC");
		bcSettings.append(" style='" + style + "'");
		if (style == 6)
			bcSettings.append(" function='" + function + "'");
		else
			bcSettings.append(" value='" + doc.formatDble(arg1) + "'");
		if (hasArg2)
			bcSettings.append(" time='" + doc.formatDble(arg2) + "'");

		bcSettings.append("/>\n");
	}

	// add velocity condition
	// set boundary ID
	public void SetBoundaryID(ArrayList<String> args) throws Exception
	{	// MPM only
		doc.requiresMPM(args);

		// no arg reverts to 0
		if (args.size() < 2)
		{	boundaryID = 0;
			return;
		}

		// set
		boundaryID = doc.readIntArg(args.get(1));
	}
	
	// add region pience
	public void addPiece(RegionPiece newPiece)
	{	pieces.add(newPiece);
	}

	// insert XML
	public void AddXML(String rawXML)
	{	xmlbcs.append(rawXML);
	}

	// ----------------------------------------------------------------------------
	// Accessors
	// ----------------------------------------------------------------------------

	// return xml data
	public String toXMLString()
	{	return xmlbcs.toString();
	}

	public int getInBC()
	{	return inBC;
	}
	
	public boolean allowsShape(int level)
	{	if(inBC!=GRID_BC) return false;
		if(level>0) return true;  // it is checked later
		if(pieces.size()==0) return true;
		return false;
	}

}
