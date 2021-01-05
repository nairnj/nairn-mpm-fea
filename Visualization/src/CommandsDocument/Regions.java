
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
	private String deform;
	private String ptsPerElement;
	private ArrayList<RegionPiece> pieces;
	private RegionPiece currentPiece;

	private static int REGION_BLOCK = 1;
	private static int HOLE_BLOCK = 2;
	private static int BMPREGION_BLOCK = 3;

	// ----------------------------------------------------------------------------
	// Initialize
	// ----------------------------------------------------------------------------

	public Regions(CmdViewer cmdDoc)
	{ // save parent CmdViewer
		doc = cmdDoc;
		pieces = null;
	}

	public void initRunSettings()
	{
		inRegion = 0;
		xmlRegions = new StringBuffer("");
		indent = "    ";
		pieces = new ArrayList<RegionPiece>(20);
	}

	// ----------------------------------------------------------------------------
	// Methods
	// ----------------------------------------------------------------------------

	// start FEA Region #1,#2,<#3>
	// #1 is material, #2 is thickness, #3 is material angle function
	// or MPM Region #1,#2,#3,#4,(#5,#6 pairs)
	// #1 is material, (#2,#3)=(vx,vy), #4=thickness or vz
	public void StartRegion(ArrayList<String> args) throws Exception
	{
		// verify not nested
		if(inRegion != 0)
			throw new Exception("Regions, Holes, and BMPRegions cannot be nested:\n" + args);

		// activate
		inRegion = REGION_BLOCK;
		pieces.clear();
		deform = null;
		ptsPerElement = null;

		// read material by ID
		if(args.size() < 2)
			throw new Exception("'Region' command missing material ID:\n" + args);
		String matID = doc.readStringArg(args.get(1));
		int matnum = doc.mats.getMatID(matID);
		if(matnum <= 0)
			throw new Exception("'Region' command has unknown material ID:\n" + args);

		// MPM or FEA
		if(doc.isMPM())
		{
			indent = "    ";

			// read two velocities
			double velx = 0., vely = 0., velz = 0., thick = 0.;
			String vel0X = null, vel0Y = null, vel0Z = null;
			if(args.size() < 4)
				throw new Exception("'Region' command missing two few arguments:\n" + args);
			Object velFxn = doc.readStringOrDoubleArg(args.get(2));
			if(velFxn.getClass().equals(Double.class))
				velx = ((Double) velFxn).doubleValue();
			else
				vel0X = "      <vel0X>" + (String) velFxn + "</vel0X>\n";
			velFxn = doc.readStringOrDoubleArg(args.get(3));
			if(velFxn.getClass().equals(Double.class))
				vely = ((Double) velFxn).doubleValue();
			else
				vel0Y = "      <vel0Y>" + (String) velFxn + "</vel0Y>\n";

			// start tag
			xmlRegions.append("    <Body mat='" + matnum + "' vx='" + doc.formatDble(velx) + "' vy='"
					+ doc.formatDble(vely) + "'");

			// thickness or velocity z
			if(args.size() > 4)
			{
				if(doc.isMPM3D())
				{
					velFxn = doc.readStringOrDoubleArg(args.get(4));
					if(velFxn.getClass().equals(Double.class))
						velz = ((Double) velFxn).doubleValue();
					else
						vel0Z = "      <vel0Z>" + (String) velFxn + "</vel0Z>\n";
				}
				else
					thick = doc.readDoubleArg(args.get(4));
			}
			if(doc.isMPM3D())
				xmlRegions.append(" vz='" + doc.formatDble(velz) + "'");
			else if(args.size() > 4)
				xmlRegions.append(" thick='" + doc.formatDble(thick) + "'");

			// extra pairs (angle, temp, conc, res)
			int nextSize = 5;
			while (args.size() > nextSize)
			{ // get property
				String prop = doc.readStringArg(args.get(nextSize)).toLowerCase();
				if(!prop.equals("angle") && !prop.equals("temp") && !prop.equals("conc")
						&& !prop.equals("pp") && !prop.equals("res"))
				{	throw new Exception("Region '" + prop + "' property not recogonized:\n" + args);
				}

				// need next one
				if(args.size() < nextSize + 2)
					throw new Exception("Region '" + prop + "' property missing a value:\n" + args);
				double dvalue = doc.readDoubleArg(args.get(nextSize + 1));

				// add it and increment
				if(prop.equals("res"))
				{	int totalPoints;
					if(doc.isMPM3D())
						totalPoints = (int)(dvalue*dvalue*dvalue+0.1);
					else
						totalPoints = (int)(dvalue*dvalue+0.1);	
					ptsPerElement = "<MatlPtsPerElement>" + totalPoints + "</MatlPtsPerElement>";
				}
				else
					xmlRegions.append(" " + prop + "='" + doc.formatDble(dvalue) + "'");
				nextSize += 2;
			}

			// end region tag
			xmlRegions.append(">\n");

			// initial velocities
			if(vel0X != null)
				xmlRegions.append(vel0X);
			if(vel0Y != null)
				xmlRegions.append(vel0Y);
			if(vel0Z != null)
				xmlRegions.append(vel0Z);
		}

		else if(doc.isFEA())
		{
			indent = "  ";

			// read thickness
			if(args.size() < 3)
				throw new Exception("'Region' command missing thickness:\n" + args);
			double thickness = doc.readDoubleArg(args.get(2));

			// start tag
			xmlRegions.append("  <Body mat='" + matnum + "' thick='" + doc.formatDble(thickness) + "'");

			// read angle
			if(args.size() > 3)
			{
				String angle = doc.readStringArg(args.get(3));
				xmlRegions.append(" angle='" + angle + "'");
			}

			// end the line
			xmlRegions.append(">\n");
		}

		else
			throw new Exception("'Region' command not allowed before analysis type is set:\n" + args);
	}

	// end current region, bmpregion, or hole
	public void EndRegion(ArrayList<String> args,String theCmd) throws Exception
	{
		// must be in region
		if(theCmd.equals("endhole"))
		{
			if(inRegion != HOLE_BLOCK)
				throw new Exception("'EndHole' command when not in a Hole:\n" + args);
		}
		else
		{
			if(inRegion != REGION_BLOCK && inRegion != BMPREGION_BLOCK)
				throw new Exception("'EndRegion' command when not in a Region or BMPRegion:\n" + args);
		}

		// add XML data
		if(inRegion == BMPREGION_BLOCK)
		{	// pts per element
			if(ptsPerElement!=null)
				xmlRegions.append(indent + "  " + ptsPerElement + "\n");
			
			xmlRegions.append(indent + "</BMP>\n");
		}
		else
		{	// pts per element
			if(ptsPerElement!=null)
				xmlRegions.append(indent + "  " + ptsPerElement + "\n");
			// deform option
			if(deform != null)
				xmlRegions.append(indent + "  " + deform + "\n");
			// go through pieces
			compilePieces(pieces, xmlRegions);

			// end the body or holr
			if(inRegion == REGION_BLOCK)
				xmlRegions.append(indent + "</Body>\n");
			else
				xmlRegions.append(indent + "</Hole>\n");
		}

		// region is done
		inRegion = 0;
	}

	// compile pieces
	public void compilePieces(ArrayList<RegionPiece> myPieces,StringBuffer myXML)
	{
		// go through pieces
		currentPiece = null;
		int numPieces = myPieces.size();
		if(numPieces > 0)
		{ // initialize
			int currentLevel = 0;

			int i = 0;
			while (i < numPieces)
			{
				RegionPiece obj = myPieces.get(i);
				int level = obj.getLevel();

				// check if this level is less than parent level
				if(level <= currentLevel && currentPiece != null)
				{	// climb back up the tree
					currentLevel = insertPriorElements(level, myXML);
				}

				switch(obj.getType())
				{	// Standard shapes
					case RegionPiece.RECT_OR_OVAL:
					case RegionPiece.SHAPE_3D:
						obj.setParent(currentPiece);
						currentLevel = level;
						currentPiece = obj;
						break;

					// 2D Polygons and 3D Polyhedrons
					case RegionPiece.POLY_PT:
						obj.setParent(currentPiece);
						currentLevel = level;
						currentPiece = obj;
						String ptIndent = indent;
						for(int ii = 0; ii < level; ii++)
							ptIndent = ptIndent + "  ";

						// loop until done
						while (i < numPieces - 1)
						{	obj = myPieces.get(i + 1);

							// exit if new level or not a polygon
							if(obj.getLevel() != level || obj.getType() != RegionPiece.POLY_PT)
								break;

							// add a polypt
							currentPiece.appendXmlStart(ptIndent + obj.getXmlStart());
							i++;
						}

						// skip a break piece
						if(i < numPieces - 1 && obj.getType() == RegionPiece.END_POLYGON)
						{	// for 3D surround currentPiece with faces
							if(doc.isMPM3D())
							{	currentPiece.setXmlStart(ptIndent+obj.getXmlStart()+ptIndent
									+currentPiece.getXmlStart()+ptIndent+"    </faces>\n");
							}
							i++;
						}
						break;

					case RegionPiece.END_POLYGON:
						// POLYPT_PIECE should always handle this
						break;

						// non shape options (must be at level 0)
					case RegionPiece.COMMAND_PIECE:
						myXML.append(obj.getXmlStart());
						break;

					default:
						break;
				}

				// next object
				i++;
			}
		}

		// finish current elements
		if(currentPiece != null)
			insertPriorElements(0, myXML);

		// done with pieces
		myPieces.clear();
	}

	// climb back tree
	public int insertPriorElements(int level,StringBuffer myXML)
	{
		int newLevel = 0;

		while (true)
		{
			RegionPiece parent = currentPiece.getParent();
			if(parent == null)
			{
				myXML.append(currentPiece.xmlString(indent));
			}
			else
			{ // add child
				parent.addChild(currentPiece);
			}

			// up to parent
			currentPiece = parent;

			// exit if done
			if(currentPiece == null)
				break;
			newLevel = currentPiece.getLevel();
			if(level > newLevel)
				break;
		}

		return newLevel;
	}

	// start Hole (FEA or MPM)
	public void StartHole(ArrayList<String> args) throws Exception
	{
		// verify not nested
		if(inRegion != 0)
			throw new Exception("Regions, Holes, and BMPRegions cannot be nested:\n" + args);

		// activate
		inRegion = HOLE_BLOCK;
		pieces.clear();
		indent = doc.isMPM() ? "    " : "  ";

		// start tag
		xmlRegions.append("    <Hole>\n");
	}

	// Cut Cut ... shape args command
	public void AddCutShape(ArrayList<String> args) throws Exception
	{ // Must be in region or a BC
			// note that level is verified later
		if(inRegion == 0 && doc.mpmGridBCs.getInBC() == 0 && doc.mpmParticleBCs.getInBC() == 0)
			throw new Exception("Cut shape commands must be in a 'Region', 'Hole' or BC block");

		// needs more parameters
		if(args.size() < 2)
			throw new Exception("Cut copmmand with no cut shape");

		// assume next argument is delimited with spaces
		// It should be cut cut ... shape first_arg where initial cuts are
		// higher levels
		String[] cutArgs = args.get(1).split(" ");

		// find the level
		int level = 1;
		String shape = "cut";
		while (cutArgs.length > level - 1)
		{
			shape = cutArgs[level - 1].toLowerCase();
			if(!shape.equals("cut"))
				break;
			level++;
		}

		// set command to shape name and first parameter to last cutArg
		args.set(0, shape); // shape name
		// first argument to shape (if it has one, e.g. empty PolyPt command)
		if(cutArgs.length > 1)
			args.set(1, cutArgs[cutArgs.length - 1]);
		else
			args.remove(1);

		// each type
		if(shape.equals("rect"))
			AddRectOrOval(args, "Rect", level);
		else if(shape.equals("oval"))
			AddRectOrOval(args, "Oval", level);
		else if(shape.equals("polypt"))
			AddPolypoint(args, level);
		else if(shape.equals("arc"))
			AddRectOrOval(args, "Arc", level);
		else if(shape.equals("line"))
		{
			if(doc.isMPM3D())
				AddBox(args, "Line", level);
			else
				AddRectOrOval(args, "Line", level);
		}
		else if(shape.equals("box"))
			AddBox(args, "Box", level);
		else if(shape.equals("sphere"))
			AddBox(args, "Sphere", level);
		else if(shape.equals("cylinder"))
			AddBox(args, "Cylinder", level);
		else if(shape.equals("torus"))
			AddBox(args, "Torus", level);
		else
			throw new Exception("An invalid shape ('" + shape + "') in a 'Cut' command");

		// set last piece level
		boolean validLevel;
		if(doc.mpmGridBCs.getInBC() != 0)
			validLevel = setLastPieceLevel(level, doc.mpmGridBCs.pieces);
		else if(doc.mpmParticleBCs.getInBC() != 0)
			validLevel = setLastPieceLevel(level, doc.mpmParticleBCs.pieces);
		else
			validLevel = setLastPieceLevel(level, pieces);
		if(!validLevel)
			throw new Exception("Incorrectly nested shape commands in 'Region', 'Hole' or BC block");
	}

	// bool set last piece level
	public boolean setLastPieceLevel(int level,ArrayList<RegionPiece> myPieces)
	{
		// cut piece already added, need at least one more for potential parent
		int numPieces = myPieces.size();
		if(numPieces < 2)
			return false;
		RegionPiece cutPiece = myPieces.get(numPieces - 1);
		RegionPiece prevPiece = myPieces.get(numPieces - 2);

		// set cut piece level
		if(cutPiece.getType() == RegionPiece.COMMAND_PIECE)
			return false;
		if(level < 1)
			return false;
		cutPiece.setLevel(level);

		// check previous piece
		// Must be shape with level one less or greater
		if(prevPiece.getType() == RegionPiece.COMMAND_PIECE)
			return false;
		if(prevPiece.getLevel() < level - 1)
			return false;

		return true;
	}

	// add shape for Rect (or Oval) #1,#2,#3,#4,<#5>,<#6> (last two for arc
	// range)
	// or Line #1,#2,#3,#4,<#5> (last is tolerance)
	// or Arc #1,#2,#3,#4,#5,#6,<#7> (last is tolerance)
	public void AddRectOrOval(ArrayList<String> args,String shape,int level) throws Exception
	{
		// Allowed in GridBC and ParticleBC (maybe)
		if(doc.mpmGridBCs.getInBC() != 0)
		{
			if(!doc.mpmGridBCs.allowsShape(level))
				throw new Exception("'" + shape + "' command not allowed in current BC block:\n" + args);
		}
		else if(doc.mpmParticleBCs.getInBC() != 0)
		{
			if(!doc.mpmParticleBCs.allowsShape(level))
				throw new Exception("'" + shape + "' command not allowed in current BC block:\n" + args);
		}

		// otherwise must be in region
		else if(inRegion == 0 || inRegion == BMPREGION_BLOCK)
			throw new Exception("'" + shape + "' command is only allowed within a Region or Hole block:\n" + args);

		// here 2D shapes only
		if(doc.isMPM3D())
			throw new Exception("'" + shape + "' command is only allowed within 2D MPM:\n" + args);

		// four numbers
		if(args.size() < 5)
			throw new Exception("'" + shape + "' command has too few parameters:\n" + args);
		double xmin = doc.readDoubleArg(args.get(1));
		double xmax = doc.readDoubleArg(args.get(2));
		double ymin = doc.readDoubleArg(args.get(3));
		double ymax = doc.readDoubleArg(args.get(4));

		// arc angles
		double arcStart = -1., arcEnd = 0., tolerance = -1.;
		if(args.size() > 5)
		{
			arcStart = doc.readDoubleArg(args.get(5));

			if(!shape.equals("Line"))
			{
				if(args.size() > 6)
				{
					arcEnd = doc.readDoubleArg(args.get(6));
					if(arcStart < 0. || arcStart > 360. || arcEnd < arcStart)
						throw new Exception("Invalid arc angles (need 0<=start<=360 and end>=start)");
				}
				else
					throw new Exception("Has arc start angle but no end angle");

				// arc may also have tolerances
				if(shape.equals("Arc"))
				{
					if(args.size() > 7)
						tolerance = doc.readDoubleArg(args.get(7));
				}
			}
		}

		// start string
		StringBuffer newShape = new StringBuffer("");
		newShape.append(
				indent + "  <" + shape + " xmin='" + doc.formatDble(xmin) + "' xmax='" + doc.formatDble(xmax) + "'");
		newShape.append(" ymin='" + doc.formatDble(ymin) + "' ymax='" + doc.formatDble(ymax) + "'");

		// Finish up line or arc
		if(shape.equals("Arc"))
		{
			if(arcStart < 0. || arcStart > 360. || arcEnd < arcStart)
				throw new Exception("Invalid arc angles (0<=start<=360 and end>=start)");

			newShape.append(" start='" + doc.formatDble(arcStart) + "' end='" + doc.formatDble(arcEnd) + "'");
			if(tolerance > 0.)
				newShape.append(" tolerance='" + doc.formatDble(tolerance) + "'");
			arcStart = -1.;
		}
		else if(shape.equals("Line"))
		{
			if(arcStart > 0.)
				newShape.append(" tolerance='" + doc.formatDble(arcStart) + "'");
			arcStart = -1.;
		}

		// add piece
		RegionPiece newPiece = new RegionPiece(RegionPiece.RECT_OR_OVAL, newShape.toString(), shape, doc);
		if(arcStart >= 0.)
			newPiece.setArcAngles(arcStart, arcEnd);

		// add to BC block or current region
		if(doc.mpmGridBCs.getInBC() != 0)
			doc.mpmGridBCs.addPiece(newPiece);
		else if(doc.mpmParticleBCs.getInBC() != 0)
			doc.mpmParticleBCs.addPiece(newPiece);
		else
			pieces.add(newPiece);
	}

	// add point to polygon
	public void AddPolypoint(ArrayList<String> args,int level) throws Exception
	{
		// Allowed in GridBC and ParticleBC (maybe)
		if(doc.mpmGridBCs.getInBC() != 0)
		{
			if(!doc.mpmGridBCs.allowsShape(level))
				throw new Exception("'PolyPt' command not allowed in current BC block:\n" + args);
		}
		else if(doc.mpmParticleBCs.getInBC() != 0)
		{
			if(!doc.mpmParticleBCs.allowsShape(level))
				throw new Exception("'PolyPt' command not allowed in current BC block:\n" + args);
		}

		// otherwise must be in region
		else if(inRegion == 0 || inRegion == BMPREGION_BLOCK)
			throw new Exception("'PolyPt' command is only allowed within a polygon sequence:\n" + args);

		// this is a 2D shape
		//if(doc.isMPM3D())
		//	throw new Exception("'PolyPt' command is only allowed within 2D MPM:\n" + args);

		// create a piece
		RegionPiece newPiece;

		// end the polygon
		if(args.size() == 1)
		{	newPiece = new RegionPiece(RegionPiece.END_POLYGON, "", "", doc);
		}

		else if(!doc.isMPM3D())
		{ 	// 2D needs two arguments
			if(args.size() < 3)
				throw new Exception("'PolyPt' command in 2D has too few parameters:\n" + args);
			double x = doc.readDoubleArg(args.get(1));
			double y = doc.readDoubleArg(args.get(2));

			// add it
			String ptStr = "    <pt x='" + doc.formatDble(x) + "' y='" + doc.formatDble(y) + "'/>\n";
			newPiece = new RegionPiece(RegionPiece.POLY_PT, ptStr, "Polygon", doc);
		}
		
		else
		{	// 3D needs three points or (style),(details)
			Object xval = doc.readStringOrDoubleArg(args.get(1));
			if(xval.getClass().equals(Double.class))
			{	if(args.size() < 4)
					throw new Exception("'PolyPt' command in 3D has too few parameters:\n" + args);
				double x = ((Double) xval).doubleValue();
				double y = doc.readDoubleArg(args.get(2));
				double z = doc.readDoubleArg(args.get(3));
				// add it
				String ptStr = "  "+doc.formatDble(x)+" "+doc.formatDble(y)
								+" "+doc.formatDble(z)+"\n";
				newPiece = new RegionPiece(RegionPiece.POLY_PT, ptStr, "Polyhedron", doc);
			}
			else
			{	// box, pyramid, tripts, or trivectors
				String style = (String)xval;
				if(style.equals("box"))
				{	// needs details
					if(args.size() < 3)
						throw new Exception("Polyhedron box needs list of points:\n" + args);
					String details = doc.readStringArg(args.get(2));
					if(details.length()!=8)
						throw new Exception("Polyhedron box must specify order of corners in 8 integers:\n" + args);
					for(int ii=0;ii<details.length();ii++)
					{	char anum = details.charAt(ii);
						if(anum<'1' || anum>'8')
							throw new Exception("Polyhedron box must specify order of corners in 8 integers:\n" + args);
					}
					style = details;
				}
				style = "<faces style='"+style+"'>\n";
				newPiece = new RegionPiece(RegionPiece.END_POLYGON, style, "", doc);
			}
		}

		// add to BC block or current region
		if(doc.mpmGridBCs.getInBC() != 0)
			doc.mpmGridBCs.addPiece(newPiece);
		else if(doc.mpmParticleBCs.getInBC() != 0)
			doc.mpmParticleBCs.addPiece(newPiece);
		else
			pieces.add(newPiece);
	}

	// add shape for Box #1,#2,#3,#4,#5,#6 or Sphere
	// Cylinder #1,#2,#3,#4,#5,#6,#7,<#8> or Torus
	public void AddBox(ArrayList<String> args,String shape,int level) throws Exception
	{
		// Allowed in GridBC and ParticleBC (maybe)
		if(doc.mpmGridBCs.getInBC() != 0)
		{
			if(!doc.mpmGridBCs.allowsShape(level))
				throw new Exception("'" + shape + "' command not allowed in current BC block:\n" + args);
		}
		else if(doc.mpmParticleBCs.getInBC() != 0)
		{
			if(!doc.mpmParticleBCs.allowsShape(level))
				throw new Exception("'" + shape + "' command not allowed in current BC block:\n" + args);
		}

		// otherwise must be in region
		else if(inRegion == 0 || inRegion == BMPREGION_BLOCK)
			throw new Exception("'" + shape + "' command is only allowed within a Region or Hole block:\n" + args);

		// this is a 3D shape
		if(!doc.isMPM3D())
			throw new Exception("'" + shape + "' command is only allowed within 3D MPM:\n" + args);

		// four numbers
		if(args.size() < 7)
			throw new Exception("'" + shape + "' command has too few parameters: " + args);
		double xmin = doc.readDoubleArg(args.get(1));
		double xmax = doc.readDoubleArg(args.get(2));
		double ymin = doc.readDoubleArg(args.get(3));
		double ymax = doc.readDoubleArg(args.get(4));
		double zmin = doc.readDoubleArg(args.get(5));
		double zmax = doc.readDoubleArg(args.get(6));

		// add it
		StringBuffer newShape = new StringBuffer("");
		newShape.append(
				indent + "  <" + shape + " xmin='" + doc.formatDble(xmin) + "' xmax='" + doc.formatDble(xmax) + "'");
		newShape.append(" ymin='" + doc.formatDble(ymin) + "' ymax='" + doc.formatDble(ymax) + "'");
		newShape.append(" zmin='" + doc.formatDble(zmin) + "' zmax='" + doc.formatDble(zmax) + "'");
		if(shape.equals("Cylinder") || shape.equals("Torus"))
		{
			if(args.size() < 8)
				throw new Exception("'" + shape + "' command has too few parameters: " + args);
			// axis
			HashMap<String, Integer> options = new HashMap<String, Integer>(10);
			options.put("x", new Integer(1));
			options.put("y", new Integer(2));
			options.put("z", new Integer(3));
			int axis = doc.readIntOption(args.get(7), options, "shape axis");
			newShape.append(" axis='" + axis + "'");
			if(args.size() > 8)
				newShape.append(" radius='" + doc.formatDble(doc.readDoubleArg(args.get(8))) + "'");
		}
		else if(shape.equals("Line"))
		{
			if(args.size() > 7)
			{
				double tolerance = doc.readDoubleArg(args.get(7));
				if(tolerance > 0.)
					newShape.append(" tolerance='" + doc.formatDble(tolerance) + "'");
			}
		}

		// add piece
		RegionPiece newPiece = new RegionPiece(RegionPiece.SHAPE_3D, newShape.toString(), shape, doc);

		// add to BC block or current region
		if(doc.mpmGridBCs.getInBC() != 0)
			doc.mpmGridBCs.addPiece(newPiece);
		else if(doc.mpmParticleBCs.getInBC() != 0)
			doc.mpmParticleBCs.addPiece(newPiece);
		else
			pieces.add(newPiece);
	}

	public void AddTransform(ArrayList<String> args) throws Exception
	{
		// times not allowed
		if(inRegion == 0 || inRegion == HOLE_BLOCK || inRegion == BMPREGION_BLOCK)
			throw new Exception("'Transform' command is only allowed within a Region block:\n" + args);

		// only one deform in region
		if(deform != null)
			throw new Exception("Only one 'Transform' command within a Region block:\n" + args);

		// get properties (note that matrix is R-I)
		Matrix3 rotate = new Matrix3(1., 0., 0., 1., 1.);
		Vector3 Trans = new Vector3(0., 0., 0.);
		Vector3 Orig = new Vector3(0., 0., 0.);
		if(args.size() > 1)
		{
			double angle = doc.readDoubleArg(args.get(1));
			double rad = angle * Math.PI / 180.;
			rotate.set(Math.cos(rad), -Math.sin(rad), Math.sin(rad), Math.cos(rad), 1.);
		}
		if(args.size() > 2)
		{
			Trans.x = doc.readDoubleArg(args.get(2));
		}
		if(args.size() > 3)
		{
			Trans.y = doc.readDoubleArg(args.get(3));
		}
		if(args.size() > 4)
		{
			Orig.x = doc.readDoubleArg(args.get(4));
		}
		if(args.size() > 5)
		{
			Orig.y = doc.readDoubleArg(args.get(5));
		}

		// 3D variables
		if(doc.isMPM3D())
		{
			if(args.size() > 6)
			{
				double angle2 = doc.readDoubleArg(args.get(6));
				Matrix3 rotY = new Matrix3();
				double rad = angle2 * Math.PI / 180.;
				rotY.set(0, 0, Math.cos(rad) - 1.);
				rotY.set(0, 2, Math.sin(rad));
				rotY.set(1, 1, 1.);
				rotY.set(2, 0, -Math.sin(rad));
				rotY.set(2, 2, Math.cos(rad) - 1.);
				rotY.setIs2D(false);
				rotate.times(rotY);
			}
			if(args.size() > 7)
			{
				double angle3 = doc.readDoubleArg(args.get(7));
				double rad = angle3 * Math.PI / 180.;
				Matrix3 rotZ = new Matrix3(Math.cos(rad), -Math.sin(rad), Math.sin(rad), Math.cos(rad), 1.);
				rotate.times(rotZ);
			}
			if(args.size() > 8)
			{
				Trans.z = doc.readDoubleArg(args.get(8));
			}
			if(args.size() > 8)
			{
				Orig.z = doc.readDoubleArg(args.get(9));
			}
		}

		// subtract identity
		rotate.set(0, 0, rotate.get(0, 0) - 1.);
		rotate.set(1, 1, rotate.get(1, 1) - 1.);
		rotate.set(2, 2, rotate.get(2, 2) - 1.);

		// build deformation gradient from non-zero elements of rotate
		String dF = "";

		// Get (R-I)(x-Orig) for component x
		String dX = "";
		double m11 = rotate.get(0, 0);
		double m12 = rotate.get(0, 1);
		if(m11 != 0.)
		{
			if(Orig.x > 0.)
				dX = dX + doc.formatDble(m11) + "*(x-" + doc.formatDble(Orig.x) + ")";
			else if(Orig.x < 0.)
				dX = dX + doc.formatDble(m11) + "*(x+" + doc.formatDble(-Orig.x) + ")";
			else
				dX = dX + doc.formatDble(m11) + "*x";
			dF = dF + " F11='" + doc.formatDble(m11 + 1) + "'";
		}
		if(m12 != 0.)
		{
			if(Orig.y > 0.)
				dX = dX + "+(" + doc.formatDble(m12) + "*(y-" + doc.formatDble(Orig.y) + "))";
			else if(Orig.y < 0.)
				dX = dX + "+(" + doc.formatDble(m12) + "*(y+" + doc.formatDble(-Orig.y) + "))";
			else
				dX = dX + "+(" + doc.formatDble(m12) + "*y)";
			dF = dF + " F12='" + doc.formatDble(m12) + "'";
		}
		if(!rotate.getIs2D())
		{
			double m13 = rotate.get(0, 1);
			if(m13 != 0.)
			{
				if(Orig.z > 0.)
					dX = dX + "+(" + doc.formatDble(m13) + "*(z-" + doc.formatDble(Orig.z) + "))";
				else if(Orig.z < 0.)
					dX = dX + "+(" + doc.formatDble(m13) + "*(z+" + doc.formatDble(-Orig.z) + "))";
				else
					dX = dX + "+(" + doc.formatDble(m13) + "*z)";
				dF = dF + " F13='" + doc.formatDble(m13) + "'";
			}
		}

		// add X translation
		if(dX.length() == 0 && Trans.x != 0.)
			dX = dX + doc.formatDble(Trans.x);
		else if(Trans.x > 0.)
			dX = dX + "+" + doc.formatDble(Trans.x);
		else if(Trans.x < 0.)
			dX = dX + "-" + doc.formatDble(-Trans.x);
		String dXattr = "";
		if(dX.length() > 0)
			dXattr = " dX='" + dX + "'";

		// Get (R-I)(x-Orig) for component y
		String dY = "";
		double m21 = rotate.get(1, 0);
		double m22 = rotate.get(1, 1);
		if(m21 != 0.)
		{
			if(Orig.x > 0.)
				dY = dY + doc.formatDble(m21) + "*(x-" + doc.formatDble(Orig.x) + ")";
			else if(Orig.x < 0.)
				dY = dY + doc.formatDble(m21) + "*(x+" + doc.formatDble(-Orig.x) + ")";
			else
				dY = dY + doc.formatDble(m21) + "*x";
			dF = dF + " F21='" + doc.formatDble(m21) + "'";
		}
		if(m22 != 0.)
		{
			if(Orig.y > 0.)
				dY = dY + "+(" + doc.formatDble(m22) + "*(y-" + doc.formatDble(Orig.y) + "))";
			else if(Orig.y < 0.)
				dY = dY + "+(" + doc.formatDble(m22) + "*(y+" + doc.formatDble(-Orig.y) + "))";
			else
				dY = dY + "+(" + doc.formatDble(m22) + "*y)";
			dF = dF + " F22='" + doc.formatDble(m22 + 1.) + "'";
		}
		if(!rotate.getIs2D())
		{
			double m23 = rotate.get(1, 2);
			if(m23 != 0.)
			{
				if(Orig.z > 0.)
					dY = dY + "+(" + doc.formatDble(m23) + "*(z-" + doc.formatDble(Orig.z) + "))";
				else if(Orig.z < 0.)
					dY = dY + "+(" + doc.formatDble(m23) + "*(z+" + doc.formatDble(-Orig.z) + "))";
				else
					dY = dY + "+(" + doc.formatDble(m23) + "*z)";
				dF = dF + " F23='" + doc.formatDble(m23) + "'";
			}
		}

		// add Y translation
		if(dY.length() == 0 && Trans.y != 0.)
			dY = dY + doc.formatDble(Trans.y);
		else if(Trans.y > 0.)
			dY = dY + "+" + doc.formatDble(Trans.y);
		else if(Trans.y < 0.)
			dY = dY + "-" + doc.formatDble(-Trans.y);
		String dYattr = "";
		if(dY.length() > 0)
			dYattr = " dY='" + dY + "'";

		// Get (R-I)(x-Orig) for component z (only if 3D)
		String dZattr = "";
		if(!rotate.getIs2D())
		{
			String dZ = "";
			double m31 = rotate.get(2, 0);
			double m32 = rotate.get(2, 1);
			double m33 = rotate.get(2, 2);
			if(m31 != 0.)
			{
				if(Orig.x > 0.)
					dZ = dZ + doc.formatDble(m31) + "*(x-" + doc.formatDble(Orig.x) + ")";
				else if(Orig.x < 0.)
					dZ = dZ + doc.formatDble(m31) + "*(x+" + doc.formatDble(-Orig.x) + ")";
				else
					dZ = dZ + doc.formatDble(m31) + "*x";
				dF = dF + " F31='" + doc.formatDble(m31) + "'";
			}
			if(m32 != 0.)
			{
				if(Orig.y > 0.)
					dZ = dZ + "+(" + doc.formatDble(m32) + "*(y-" + doc.formatDble(Orig.y) + "))";
				else if(Orig.y < 0.)
					dZ = dZ + "+(" + doc.formatDble(m32) + "*(y+" + doc.formatDble(-Orig.y) + "))";
				else
					dZ = dZ + "+(" + doc.formatDble(m32) + "*y)";
				dF = dF + " F32='" + doc.formatDble(m32) + "'";
			}
			if(m33 != 0.)
			{
				if(Orig.z > 0.)
					dZ = dZ + "+(" + doc.formatDble(m33) + "*(z-" + doc.formatDble(Orig.z) + "))";
				else if(Orig.z < 0.)
					dZ = dZ + "+(" + doc.formatDble(m33) + "*(z+" + doc.formatDble(-Orig.z) + "))";
				else
					dZ = dZ + "+(" + doc.formatDble(m33) + "*z)";
				dF = dF + " F33='" + doc.formatDble(m33 + 1.) + "'";
			}

			// add Z translation
			if(dZ.length() == 0 && Trans.z != 0.)
				dZ = dZ + doc.formatDble(Trans.z);
			else if(Trans.z > 0.)
				dZ = dZ + "+" + doc.formatDble(Trans.z);
			else if(Trans.z < 0.)
				dZ = dZ + "-" + doc.formatDble(-Trans.z);
			if(dZ.length() > 0)
				dZattr = " dZ='" + dZ + "'";
		}

		// set final command
		deform = "<Deform" + dXattr + dYattr + dZattr + dF + "/>";

		// remove empty transform
		if(deform.length() <= 0)
			deform = null;
	}

	// Rotate #1,#2,<#3,#4>,<#5,$6>
	public void AddRotate(ArrayList<String> args) throws Exception
	{ // times
			// not
		// allowed
		if(inRegion == 0 || inRegion == HOLE_BLOCK)
			throw new Exception("'Rotate' command is only allowed within a Region or BMPRegion block:\n" + args);

		// check for reset
		if(args.size() < 2)
			throw new Exception("'Rotate' command has too few parameters:\n" + args);

		// Is is "reset"
		String reset = doc.readStringArg(args.get(1)).toLowerCase();
		if(reset.equals("reset"))
		{
			String rotStr = indent + "  <Unrotate/>\n";
			if(inRegion == REGION_BLOCK)
			{
				RegionPiece newPiece = new RegionPiece(RegionPiece.COMMAND_PIECE, rotStr, "", doc);
				pieces.add(newPiece);
			}
			else
				xmlRegions.append(rotStr);
			return;
		}

		// need at least one axis
		if(args.size() < 3)
			throw new Exception("'Rotate' command has too few parameters:\n" + args);

		// up to three pairs (3D only)
		int axisNum = 2;
		while (args.size() > axisNum && axisNum < 8)
		{ // get axis
			int axis = 0;
			Object axisArg = doc.readStringOrDoubleArg(args.get(axisNum - 1));
			if(axisArg.getClass().equals(Double.class))
				axis = ((Double) axisArg).intValue();
			else if(((String) axisArg).toLowerCase().equals("x"))
				axis = 1;
			else if(((String) axisArg).toLowerCase().equals("y"))
				axis = 2;
			else if(((String) axisArg).toLowerCase().equals("z"))
				axis = 3;
			if(axis < 1 || axis > 3)
				throw new Exception("'Rotate' command has invalid rotation axis:\n" + args);
			if(axis != 3 && !doc.isMPM3D())
				throw new Exception("'Rotate' command axis must be z axis in 2D simulations:\n" + args);

			// get angle (can be function)
			String angle = doc.readStringArg(args.get(axisNum));

			// add piece
			String newShape;
			if(axis == 1)
				newShape = indent + "  <RotateX>" + angle + "</RotateX>\n";
			else if(axis == 2)
				newShape = indent + "  <RotateY>" + angle + "</RotateY>\n";
			else
				newShape = indent + "  <RotateZ>" + angle + "</RotateZ>\n";
			if(inRegion == REGION_BLOCK)
			{
				RegionPiece newPiece = new RegionPiece(RegionPiece.COMMAND_PIECE, newShape, "", doc);
				pieces.add(newPiece);
			}
			else
				xmlRegions.append(newShape);

			// next pair
			axisNum += 2;
			if(!doc.isMPM3D())
				break;
		}
	}

	// AngularMom0 Lpz (if 2D) or AngularMom0 Lpx,Lpy,Lpz (if 3D)
	public void AddAngularMom0(ArrayList<String> args) throws Exception
	{ // check
			// allowed
		if(inRegion == 0 || inRegion != REGION_BLOCK)
			throw new Exception("'AngularMom0' command is only allowed within a Region block:\n" + args);

		// 2D or 3D
		if(args.size() < 2)
			throw new Exception("'AngulaMom0' command has too few parameters: " + args);

		String Lp;
		Object LpFxn = doc.readStringOrDoubleArg(args.get(1));
		if(doc.isMPM3D())
		{
			Lp = "      <Lp0X>" + LpFxn + "</Lp0X>\n";
			if(args.size() > 2)
				Lp = Lp + "      <Lp0Y>" + doc.readStringOrDoubleArg(args.get(2)) + "</Lp0Y>\n";
			if(args.size() > 3)
				Lp = Lp + "      <Lp0Z>" + doc.readStringOrDoubleArg(args.get(3)) + "</Lp0Z>\n";
		}
		else
			Lp = "      <Lp0Z>" + LpFxn + "</Lp0Z>\n";
		xmlRegions.append(Lp);
	}

	// start FEA Region #1,#2,<#3>
	// #1 is material, #2 is thickness, #3 is material angle function
	// or MPM BMPRegion #1,#2,#3,<#4>,(#5,#6 pairs)
	// #1 is file name, (#2,#3)=(width,height), #4=angles file
	public void StartBMPRegion(ArrayList<String> args) throws Exception
	{
		// verify not nested
		if(inRegion != 0)
			throw new Exception("Regions, Holes, and BMPRegions cannot be nested:\n" + args);

		// activate
		inRegion = BMPREGION_BLOCK;
		pieces.clear();
		ptsPerElement = null;

		// read path and width
		if(args.size() < 2)
			throw new Exception("'BMP' has too few parameters:\n" + args);
		String filePath = doc.readStringArg(args.get(1));
		
		// optional width
		double width = -1.e9;
		if(args.size() > 2)
			width = doc.readDoubleArg(args.get(2));

		// optional height
		double height = -1.e9;
		if(args.size() > 3)
			height = doc.readDoubleArg(args.get(3));

		// optional angles path
		String[] anglesPath = new String[3];
		anglesPath[0] = null;
		anglesPath[1] = null;
		anglesPath[2] = null;
		String scheme = "";
		int nextSize=4;
		if(args.size() > nextSize)
		{	// get parameter, which will be scheme or a file path (or empty to skip)
			anglesPath[0] = doc.readStringArg(args.get(nextSize));
			nextSize++;

			// MPM file can have anglesPath (old method) or
			// scheme,angle1,angle2,angle3
			// FEA file can have anglesPath (old method) or "Z",angle1
			// assume scheme if 3 or less characters; all files name have more than 3 characters
			int schemeLength = anglesPath[0].length();
			if((schemeLength<4) && (schemeLength>0))
			{	scheme = anglesPath[0].toUpperCase();
				if(!doc.isMPM3D() && !scheme.equals("Z"))
					throw new Exception("2D Simulations can only have a 'Z' rotation scheme:\n" + args);
				
				// read more files
				for(int i = 0; i < scheme.length(); i++)
				{	if(args.size() > nextSize)
					{	anglesPath[i] = doc.readStringArg(args.get(nextSize));
						nextSize++;
					}
					else
					{	throw new Exception(
								"'BMPRegion' command missing angle file name needed for rotation scheme:\n" + args);
					}
				}
			}
			else if(schemeLength>0)
			{	// accept one angle file to rotate about Z
				scheme = "Z";
			}
		}

		// extra pairs (res) (with future room for more)
		while (args.size() > nextSize)
		{ // get property
			String prop = doc.readStringArg(args.get(nextSize)).toLowerCase();
			if(!prop.equals("res"))
			{	throw new Exception("BMPRegion '" + prop + "' property not recogonized:\n" + args);
			}

			// need next one
			if(args.size() < nextSize + 2)
				throw new Exception("BMPRegion '" + prop + "' property missing a value:\n" + args);
			double dvalue = doc.readDoubleArg(args.get(nextSize + 1));

			// add it and increment
			if(prop.equals("res"))
			{	int totalPoints;
				if(doc.isMPM3D())
					totalPoints = (int)(dvalue*dvalue*dvalue+0.1);
				else
					totalPoints = (int)(dvalue*dvalue+0.1);	
				ptsPerElement = "<MatlPtsPerElement>" + totalPoints + "</MatlPtsPerElement>";
			}
			else
				xmlRegions.append(" " + prop + "='" + doc.formatDble(dvalue) + "'");
			nextSize += 2;
		}
		
		// set indent
		if(doc.isMPM())
			indent = "    ";
		else if(doc.isFEA())
			indent = "  ";
		else
			throw new Exception("'BMPRegion' command not allowed before analysis type is set:\n" + args);

		// create the command
		xmlRegions.append(indent + "<BMP name='" + filePath + "'");

		// optional attributes
		if(width > -1.e8)
			xmlRegions.append(" width='" + doc.formatDble(width) + "'");
		if(height > -1.e8)
			xmlRegions.append(" height='" + doc.formatDble(height) + "'");

		// 0 to 3 angle file path names
		for(int i = 0; i < scheme.length(); i++)
		{
			xmlRegions.append(" angles" + scheme.charAt(i) + "='" + anglesPath[i] + "'");
		}

		// end it
		xmlRegions.append(">\n");
	}

	// Origin #1,#2,<#3>,<#4>
	public void setOrigin(ArrayList<String> args) throws Exception
	{
		// must have at least two arguments
		if(args.size() < 3)
			throw new Exception("'Origin' command requires at leat two coordinates: " + args);

		double xo = doc.readDoubleArg(args.get(1));
		double yo = doc.readDoubleArg(args.get(2));

		double zo = 0;
		if(args.size() > 3)
			zo = doc.readDoubleArg(args.get(3));

		String flip = null;
		if(args.size() > 4)
			flip = doc.readStringArg(args.get(4)).toLowerCase();

		// add the command
		xmlRegions.append(indent + "  <Origin x='" + doc.formatDble(xo) + "' y='" + doc.formatDble(yo) + "'");
		if(doc.isMPM3D())
			xmlRegions.append(" z='" + doc.formatDble(zo) + "'");
		if(flip != null)
			xmlRegions.append(" flipped='" + flip + "'");
		xmlRegions.append("/>\n");
	}

	// Intensity #1,#2,#3,<#4,#5>,...
	// Intensity "angles",#2,#3,#4,#5
	public void AddIntensity(ArrayList<String> args) throws Exception
	{
		// verify not nested
		if(inRegion != BMPREGION_BLOCK)
			throw new Exception("'Intensity' command only allowed in a BMPRegion block:\n" + args);

		// read material by ID
		if(args.size() < 2)
			throw new Exception("'Intensity' command has no arguments:\n" + args);
		String matID = doc.readStringArg(args.get(1));
		int matnum = doc.mats.getMatID(matID);

		// get gray range
		if(args.size() < 4)
			throw new Exception("'Intensity' command is missing gray scale range parameters:\n" + args);
		int gray1 = doc.readIntArg(args.get(2));
		int gray2 = doc.readIntArg(args.get(3));

		// map materials or angles
		if(matnum > 0)
		{ // get gray range
			if(args.size() < 4)
				throw new Exception("'Intensity' command has no arguments:\n" + args);

			// the command
			xmlRegions
					.append(indent + "  <Intensity mat='" + matnum + "' imin='" + gray1 + "' imax='" + gray2 + "'>\n");

			// process each pair
			int propNum = 5;
			double vx = 0., vy = 0., vz = 0.;
			boolean hasVx = false, hasVy = false, hasVz = false;
			while (args.size() > propNum)
			{
				String prop = doc.readStringArg(args.get(propNum - 1)).toLowerCase();
				double value = doc.readDoubleArg(args.get(propNum));
				propNum += 2;

				if(prop.equals("thick"))
					xmlRegions.append(indent + "    <Thickness>" + doc.formatDble(value) + "</Thickness>\n");
				else if(prop.equals("temp"))
					xmlRegions.append(indent + "    <Temperature>" + doc.formatDble(value) + "</Temperature>\n");
				else if(prop.equals("angle"))
					xmlRegions.append(indent + "    <Angle>" + doc.formatDble(value) + "</Angle>\n");
				else if(prop.equals("conc") && doc.isMPM())
					xmlRegions.append(indent + "    <Concentration>" + doc.formatDble(value) + "</Concentration>\n");
				else if(prop.equals("pp") && doc.isMPM())
					xmlRegions.append(indent + "    <Concentration>" + doc.formatDble(value) + "</Concentration>\n");
				else if(prop.equals("vx") && doc.isMPM())
				{
					vx = value;
					hasVx = true;
				}
				else if(prop.equals("vy") && doc.isMPM())
				{
					vy = value;
					hasVy = true;
				}
				else if(prop.equals("vz") && doc.isMPM())
				{
					vz = value;
					hasVz = true;
				}
				else
					throw new Exception("'Intensity' command invalid property options:\n" + args);
			}

			// velocity
			if(hasVx || hasVy || hasVz)
			{
				xmlRegions.append(indent + "    <vel");
				if(hasVx)
					xmlRegions.append(" x='" + doc.formatDble(vx) + "'");
				if(hasVy)
					xmlRegions.append(" y='" + doc.formatDble(vy) + "'");
				if(hasVz && doc.isMPM3D())
					xmlRegions.append(" z='" + doc.formatDble(vz) + "'");
				xmlRegions.append("/>\n");
			}

			// finish up
			xmlRegions.append(indent + "  </Intensity>\n");
		}
		else
		{ // get angle range
			if(args.size() < 6)
				throw new Exception("'Intensity' command is missing angle range parameters:\n" + args);
			double angle1 = doc.readDoubleArg(args.get(4));
			double angle2 = doc.readDoubleArg(args.get(5));

			// the commmad
			xmlRegions.append(indent + "  <Intensity imin='" + gray1 + "' imax='" + gray2 + "'");
			xmlRegions
					.append(" minAngle='" + doc.formatDble(angle1) + "' maxAngle='" + doc.formatDble(angle2) + "'/>\n");
		}
	}

	// insert XML data in pieces or xml data
	public void AddXML(String rawXML)
	{
		xmlRegions.append(rawXML);
	}

	// ----------------------------------------------------------------------------
	// Accessors
	// ----------------------------------------------------------------------------

	// combine to xml data
	public String toXMLString()
	{
		return xmlRegions.toString();
	}

	public boolean isInBMPRegion()
	{
		return inRegion == BMPREGION_BLOCK;
	}

	public boolean isInRegion()
	{
		return inRegion != 0;
	}
}
