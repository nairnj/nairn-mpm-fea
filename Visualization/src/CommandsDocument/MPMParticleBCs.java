/*
 * MPMParticleBCs.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 26 Sep 2013.
 * Copyright (c) 2013 RSAC Software. All rights reserved.
 * 
 * Start BC: LoadLine, LoadArc, LoadRect, MoveBox
 * Subordinate: Load, Traction, ConcentrationFlux, PorePressureFlux, HeatFlux, LoadType
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
	public ArrayList<RegionPiece> pieces;

	public static final int ADD_LOAD=1;
	public static final int ADD_TRACTION=2;
	public static final int ADD_HEATFLUX=3;
	public static final int ADD_CONCENTRATIONFLUX=4;
	public static final int ADD_DAMAGE=5;
	public static final int ADD_PHASEFIELD=6;
	
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
		pieces=new ArrayList<RegionPiece>(20);
	}
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	// finish Load line and load arc decoded in grid BCs
	public void SetLoadLine(String theAttrs,int theType,String theCmd)
	{	bcAttrs = theAttrs;
		inBC = theType;
		bcCmd = theCmd;
		bcSettings = new StringBuffer("");
		pieces.clear();
	}

	// start grid BC line
	public void StartLoadRect(ArrayList<String> args) throws Exception
	{	
		// MPM Only
		doc.requiresMPM(args);

		// verify not nested
		if(inBC != 0)
			throw new Exception("LoadLine, LoadArc, LoadRect, and LoadBox cannot be nested:\n"+args);
    
		// needs at least 5 arguments
		if(args.size()<5)
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
    
		// get xmin,xmax,ymin,ymax
		double xmin = doc.readDoubleArg(args.get(1));
		double xmax = doc.readDoubleArg(args.get(2));
		double ymin = doc.readDoubleArg(args.get(3));
		double ymax = doc.readDoubleArg(args.get(4));
    	
		bcAttrs = "<LdRect xmin='"+doc.formatDble(xmin)+"' xmax='"+doc.formatDble(xmax)
					+"' ymin='"+doc.formatDble(ymin)+"' ymax='"+doc.formatDble(ymax)+"'>\n";
		bcSettings = new StringBuffer("");
		inBC = MPMGridBCs.LOADRECT_BC;
		bcCmd = "LdRect";
		pieces.clear();				// but never used
	}

	// LoadLine, LoadArc, LoadRect or LoadBox is done
	public void EndLoadBlock(ArrayList<String> args,int endType) throws Exception
	{
		if(inBC != endType)
			throw new Exception("'"+args.get(0)+"' does not match current boundary condition block:\n"+args);
		
		// get shape
		StringBuffer shapeXML = new StringBuffer("");
		if(inBC == MPMGridBCs.PARTICLE_BC)
		{	doc.regions.compilePieces(pieces, shapeXML);
		}

		// append block
		xmlbcs.append("    "+bcAttrs+ shapeXML + bcSettings+"    </"+bcCmd+">\n");
		
		inBC = 0;
	}
	
	// Add one of the following boundary conditions
	// ADD_LOAD: Load dir,style,arg1,arg2
	// ADD_TRACTION: Traction dir,face,style,arg1,arg2
	// ADD_HEATFLUX: HeatFlux "external",face,style,arg1,arg2
	// ADD_CONCENTRATIONFLUX: ConcentrationFlux "external",face,style,arg1,arg2,phaseStyle
	// ADD_DAMAGE: Damage (nx),(ny),<(nz)>,<(dn)>,<(dxy)>,<(dxz)>,<(mode)>
	// ADD_PHASEFIELD: PhaseField (phi)
	public void AddCondition(ArrayList<String> args,int theType) throws Exception
	{
		// MPM only
		doc.requiresMPM(args);
		
		if(inBC == 0)
			throw new Exception("'"+args.get(0)+"' command must by in a 'ParticleBC' block:\n"+args);
		
		// always needs #1, #2, and #3 (those with face need #4 to), add damage 2D only needs 2, but 3 in 3D
		if(theType==ADD_DAMAGE)
		{	if(args.size()<3 || (doc.isMPM3D() && args.size()<4))
	    		throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
			
		}
		else if(theType==ADD_PHASEFIELD)
		{	if(args.size()<2)
    			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
			String phi = doc.readStringArg(args.get(1));
			bcSettings.append("      <Damage phi='"+phi+"'/>\n");
			return;
		}
		else if((theType!=ADD_LOAD && args.size()<5) || (theType==ADD_LOAD && args.size()<4))
	    	throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		// read direction
		HashMap<String,Integer> options = new HashMap<String,Integer>(7);
		if(theType==ADD_CONCENTRATIONFLUX || theType==ADD_HEATFLUX)
		{	// external or coupled
			options.put("external", new Integer(1));
			options.put("coupled", new Integer(2));
		}
		else if(theType==ADD_DAMAGE)
		{	// (#1,#2,#3) is normal vector, (#4,#5,#6) are dn, dxy, and dxz
			// #7 is mode
			// process and then done with this BC type
			double nx = doc.readDoubleArg(args.get(1)); 
			double ny = doc.readDoubleArg(args.get(2));
			double nz = args.size()>=4 ? doc.readDoubleArg(args.get(3)) : 0.; 
			double dn = args.size()>=5 ? doc.readDoubleArg(args.get(4)) : 1.; 
			double dxy = args.size()>=6 ? doc.readDoubleArg(args.get(5)) : 1.; 
			double dxz = args.size()>=7 ? doc.readDoubleArg(args.get(6)) : 1.;
			double mode = args.size()>=8 ? doc.readDoubleArg(args.get(7)) : 1.;
			if(doc.isMPM3D())
			{	bcSettings.append("      <Damage nx='"+doc.formatDble(nx)+"' ny='"+doc.formatDble(ny)+
					"' nz='"+doc.formatDble(nz)+"' dn='"+doc.formatDble(dn)+
					"' dxy='"+doc.formatDble(dxy)+"' dxz='"+doc.formatDble(dxz)+
					"' mode='"+doc.formatDble(mode)+"'");
			}
			else
			{	bcSettings.append("      <Damage nx='"+doc.formatDble(nx)+"' ny='"+doc.formatDble(ny)+
					"' dn='"+doc.formatDble(dn)+"' dxy='"+doc.formatDble(dxy)+
					"' mode='"+doc.formatDble(mode)+"'");
			}
			bcSettings.append("/>\n");
			return;
		}
		else
		{	// directions for load and traction
			options.put("x", new Integer(1));
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
		
		// face if needed (all except load which is at particle center)
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
		// silent not allowed for traction of for concentrationflux with phaseStyle not solvent
		if(theType!=ADD_TRACTION) options.put("silent", new Integer(5));
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
			arg++;
		}
		
		// get optional phaseStyle
		int phaseStyle = 1;
		if(theType==ADD_CONCENTRATIONFLUX && args.size()>arg)
		{	options = new HashMap<String,Integer>(5);
			options.put("solvent", new Integer(1));
			options.put("moisture", new Integer(1));
			options.put("fracture", new Integer(3));
			options.put("battery", new Integer(4));
			options.put("conduction", new Integer(5));
			phaseStyle = doc.readIntOption(args.get(arg),options,"flux phaseStyle");
			arg++;
		}
		
		if(theType==ADD_LOAD)
		{	// <LoadBC dir='1' style='1' load='400' time='0.0' function='x*t'/>
			bcSettings.append("      <LoadBC dir='"+dof+"' style='"+style+"'");
			if(style==6)
				bcSettings.append(" function='"+function+"'");
			else
				bcSettings.append(" load='"+doc.formatDble(arg1)+"'");
		}
		else if(theType==ADD_TRACTION)
		{	// <TractionBC dir="2" face="3" style="1" stress="1" time='0.0' function='x*t'/>
			bcSettings.append("      <TractionBC dir='"+dof+"' face='"+face+"' style='"+style+"'");
			if(style==6)
				bcSettings.append(" function='"+function+"'");
			else
				bcSettings.append(" stress='"+doc.formatDble(arg1)+"'");
		}
		else if(theType==ADD_HEATFLUX)
		{	// <HeatFluxBC dir='1' face='1' style='1' value='0' time='0.0' function='sinh(t)'/>
			bcSettings.append("      <HeatFluxBC dir='"+dof+"' face='"+face+"' style='"+style+"'");
			if(dof==2 && style!=6)
				throw new Exception("Coupled HeatFlux must use a function");
			if(style==6)
				bcSettings.append(" function='"+function+"'");
			else
				bcSettings.append(" value='"+doc.formatDble(arg1)+"'");
			
		}
		else
		{	// <ConcFluxBC dir='1' face='1' style='1' value='0' time='0.0' function='sinh(t)'/>
			bcSettings.append("      <ConcFluxBC dir='"+dof+"' face='"+face+"' style='"+style+"'");
			if(dof==2 && style!=6)
				throw new Exception("Coupled ConcentrationFlux or PorePressureFlux must use a function");
			if(style==6)
				bcSettings.append(" function='"+function+"'");
			else
				bcSettings.append(" value='"+doc.formatDble(arg1)+"'");
			if(phaseStyle!=1)
				bcSettings.append(" phase='"+phaseStyle+"'");
		}	
		
		// time arg
		if(hasArg2) bcSettings.append(" time='"+doc.formatDble(arg2)+"'");
		
		// end command
		bcSettings.append("/>\n");
	}
	
	public void doLoadType(ArrayList<String> args) throws Exception
	{
		// MPM only
		doc.requiresMPM(args);
	
		if(inBC == 0)
			throw new Exception("'"+args.get(0)+"' command must by in 'LoadLine',\n'LoadArc', 'LoadRect', or 'LoadBox' block:\n"+args);
		
		// always needs #1, #2, and #3 (those with face need #4 to)
		if(args.size()<2)
	    	throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		String netType = doc.readStringArg(args.get(1)).toLowerCase();
		
		if(netType.equals("net"))
			bcSettings.append("      <net/>\n");
		else if(netType.equals("perparticle"))
			bcSettings.append("      <perParticle/>\n");
		else
	    	throw new Exception("'"+args.get(0)+"' has an invalid options:\n"+args);

	}
	
	// add region pience
	public void addPiece(RegionPiece newPiece)
	{	pieces.add(newPiece);
	}
	
	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	// return xml data
	public String toXMLString()
	{	return xmlbcs.toString();
	}
	
	public int getInBC()
	{	return inBC;
	}

	public boolean allowsShape(int level)
	{	if(inBC!=MPMGridBCs.PARTICLE_BC) return false;
		if(level>0) return true;  // it is checked later
		if(pieces.size()==0) return true;
		return false;
	}

}
