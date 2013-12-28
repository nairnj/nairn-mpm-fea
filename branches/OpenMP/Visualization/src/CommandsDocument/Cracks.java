/*
 * Cracks.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 23 SEP 2013.
 * Copyright (c) 2013 RSAC Software. All rights reserved.
 */

import java.util.*;

public class Cracks
{
	private CmdViewer doc;
	
	private StringBuffer settings = null;
	
	private StringBuffer crackList = null;
	private StringBuffer currentCrack = null;
	private boolean crackFixed;
	private String crackFriction;
	private String crackThickness;
	private String movePlane;
	private double cx,cy;
	private String Friction;
	private String altPropagate;
	private String propagate;
	private String propLength;
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public Cracks(CmdViewer cmdDoc)
	{	// save parent CmdViewer
		doc = cmdDoc;
	}
	
	public void initRunSettings()
	{	
		settings = new StringBuffer("");
		crackList = new StringBuffer("");
		currentCrack = null;
		Friction = null;
		altPropagate = null;
		propagate = null;
		movePlane = null;
		propLength = null;
	}
	
	//----------------------------------------------------------------------------
	// Commands
	//----------------------------------------------------------------------------
	
	// JContour #1,#2 (size and number of terms)
	public void doJContour(ArrayList<String> args) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);
		
		// read analysis type
		if(args.size()<2)
			throw new Exception("'JContour' has too few parameters:\n"+args);
		
		// contour size
		int cells = doc.readIntArg(args.get(1));
		settings.append("      <JContour type='1' size='"+cells+"'");
		
		// number of terms
		int terms=0;
		if(args.size()>2)
		{	terms = doc.readIntArg(args.get(2));
			if(terms<1 || terms>2)
				throw new Exception("'JContour' second parameter not valid:\n"+args);
			
			settings.append(" terms='"+terms+"'/>\n");
		}
		else
			settings.append("/>\n");
	}
	
	// NewCrack x,y,<#3>,<#4>,<#5>
	// #3 = material (tip=), exterior (tip=-2), or fixed (type='fixed'), or free (no tip)
	// #4 = friction setting (friction=) or "traction"
	// #5 - traction material (mat=)
	public void StartCrack(ArrayList<String> args) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);
		
		// save current crack
		appendCurrentCrack();
		
		// reset global crack properties
		currentCrack = new StringBuffer("");
		crackFixed = false;
		crackFriction = null;
		crackThickness = null;
		
		// read analysis type
		if(args.size()<3)
			throw new Exception("'NewCrack' has too few parameters:\n"+args);
		double ptx = doc.readDoubleArg(args.get(1));
		double pty = doc.readDoubleArg(args.get(2));
		
		// material ID or alternates
		int tip=-1;
		if(args.size()>3)
		{	String matID = doc.readStringArg(args.get(3));
			if(matID.toLowerCase().equals("fixed"))
			{	crackFixed = true;
			}
			else if(matID.toLowerCase().equals("exterior"))
			{	tip = -2;
			}
			else if(!matID.toLowerCase().equals("free"))
			{	// look for material ID (but can't check valid material)
				tip = doc.mats.getMatID(matID);
				if(tip==-1)
					throw new Exception("'NewCrack' has unknown crack tip material:\n"+args);
			}	
		
		}
		
		// friction or traction
		int mat=-1;
		if(args.size()>4)
		{	String frict = doc.readStringArg(args.get(4));
			if(frict.toLowerCase().equals("traction"))
			{	// get traction material (but can't check if traction law)
				if(args.size()>5)
					mat = doc.mats.getMatID(doc.readStringArg(args.get(5)));
				if(mat==-1)
					throw new Exception("'NewCrack' has unknown traction law material:\n"+args);
			}
			else
			{	// look for custom friction setting
				ArrayList<String> fargs = new ArrayList<String>();
				fargs.add("Friction");
				fargs.add(frict);
				crackFriction=doc.doFriction(fargs,3);
			}
		}
		
		// start current crack XML
		currentCrack.append("    <pt x='"+ptx+"' y='"+pty+"'");
		if(tip!=-1) currentCrack.append(" tip='"+tip+"'");
		if(mat!=-1) currentCrack.append(" mat='"+mat+"'");
		currentCrack.append("/>\n");
		
		cx = ptx;
		cy = pty;
	}
	
	// 0: GrowCrack x,y,<tip>,<mat>
	// 1: GrowCrackLine x,y,resolution,<tip>,<mat
	// <tip> = material ID (invalid OK if #5 is there) (tip= or end_tip=)
	// <mat> - traction material (mat=)
	public void GrowCrack(ArrayList<String> args,int option) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);
		
		// must be in a crack
		if(currentCrack==null)
			throw new Exception("'"+args.get(0)+"' command must be in an active crack:\n"+args);
		
		// read analysis type
		if(args.size()<3)
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		double ptx = doc.readDoubleArg(args.get(1));
		double pty = doc.readDoubleArg(args.get(2));
		
		// align params
		int param=3;
		
		// segments (GrowCrackLine)
		int segs=1;
		if(option>0)
		{	if(args.size()>param)
			{	segs = doc.readIntArg(args.get(param));
				param++;
			}
		}
		
		// material ID (error only if #3 is absent)
		int tip=-1;
		if(args.size()>param)
		{	String matID = doc.readStringArg(args.get(param));
			if(matID.toLowerCase().equals("exterior"))
			{	tip = -2;
			}
			else
			{	tip = doc.mats.getMatID(matID);
				if(tip==-1 && args.size()<param+2)
					throw new Exception("'"+args.get(0)+"' has unknown crack tip material:\n"+args);
			}
			param++;
		}
		
		// traction law
		int mat=-1;
		if(args.size()>param)
		{	String tract = doc.readStringArg(args.get(param));
			mat = doc.mats.getMatID(tract);
			if(mat==-1)
				throw new Exception("'"+args.get(0)+"' has unknown traction law material:\n"+args);
			param++;
		}
		
		// grow current crack XML
		if(option==0)
		{	currentCrack.append("    <Line x='"+ptx+"' y='"+pty+"'");
			if(tip!=-1) currentCrack.append(" tip='"+tip+"'");
		}
		else
		{	currentCrack.append("    <Line xmin='"+cx+"' ymin='"+cy+"'");
			currentCrack.append(" xmax='"+ptx+"' ymax='"+pty+"' resolution='"+segs+"'");
			if(tip!=-1) currentCrack.append(" end_tip='"+tip+"'");
		}
		if(mat!=-1) currentCrack.append(" mat='"+mat+"'");
		currentCrack.append("/>\n");
		
		// save location
		cx = ptx;
		cy = pty;
	}
	
	// set global friction command
	public void setFriction(String cmd) { Friction = cmd; }

	// set global friction command
	public void setCrackFriction(ArrayList<String> args,String cmd) throws Exception
	{	// must be in a crack
		if(currentCrack==null)
			throw new Exception("'CrackInteface' command must be in an active crack:\n"+args);
		crackFriction = cmd;
	}

	// CrackThickness #1
	public void doProagateLength(ArrayList<String> args) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);
		
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' command with too few arguments:\n"+args);
		
		propLength = "      <PropagateLength>"+doc.readDoubleArg(args.get(1))+"</PropagateLength>\n";
	}
	
	// CrackThickness #1
	public void doCrackThickness(ArrayList<String> args) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);
		
		// must be in a crack
		if(currentCrack==null)
			throw new Exception("'"+args.get(0)+"' command must be in an active crack:\n"+args);
		
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' command with too few arguments:\n"+args);
		
		crackThickness = "    <Thickness>"+doc.readDoubleArg(args.get(1))+"</Thickness>\n";
	}
	
	// Propagaate (crit),<(dir)>,<(traction)>
	// AltProagate (crit),<(dir)>,<(traction)>
	public void doPropagate(ArrayList<String> args,boolean isAlt) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);
		
		// must have criterion at least
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' command with too few arguments:\n"+args);
		
		// get criterion
		int critNum = decodeCriterion(args.get(1));
		
		// crack growth direction
		int dirNum = -1;
		if(args.size()>2)
			dirNum = decodeDirection(args.get(2));
		
		// traction material
		int mat = -1;
		if(args.size()>3)
		{	String tract = doc.readStringArg(args.get(3));
			mat = doc.mats.getMatID(tract);
			if(mat==-1)
				throw new Exception("'"+args.get(0)+"' has unknown traction law material:\n"+args);
		}
		
		// ouput
		StringBuffer cmd = new StringBuffer("");
		if(isAlt)
			cmd.append("      <AltPropagate");
		else
			cmd.append("      <Propagate");
		cmd.append(" criterion='"+critNum+"'");
		if(dirNum>=0) cmd.append(" direction='"+dirNum+"'");
		if(mat>=0) cmd.append(" traction='"+mat+"'");
		cmd.append("/>\n");
		
		if(isAlt)
			altPropagate = cmd.toString();
		else
			propagate = cmd.toString();
		
	}

	//  type,<prevent>
	public void doMovePlane(ArrayList<String> args) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);
		
		// must have criterion at least
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' command with too few arguments:\n"+args);
		
		// type
		String mtype = doc.readStringArg(args.get(1));
		if(!mtype.equals("avg") && !mtype.equals("cm"))
			throw new Exception("'"+args.get(0)+"' type argument is not valid:\n"+args);
		
		// prevent
		String prevent = null;
		if(args.size()>2)
		{	prevent = doc.readStringArg(args.get(2));
			if(prevent.equals("0")) prevent = "no";
			if(prevent.equals("1")) prevent = "yes";
			if(prevent.equals("true")) prevent = "yed";
			if(prevent.equals("false")) prevent = "no";
			if(!prevent.equals("yes") && !prevent.equals("no"))
				throw new Exception("'"+args.get(0)+"' prevent argument is not valid:\n"+args);
			
			movePlane = "      <MovePlane type='"+mtype+"' prevent='"+prevent+"'/>\n";
		}
		else
			movePlane = "      <MovePlane type='"+mtype+"'/>\n";
	}
	
	// when done or start new crack, append current one
	public void appendCurrentCrack()
	{
		if(currentCrack==null) return;
		
		crackList.append("  <CrackList");
		if(crackFixed) crackList.append(" type='fixed'");
		if(crackFriction!=null) crackList.append(crackFriction);
		
		// finish up
		crackList.append(">\n"+currentCrack.toString());
		if(crackThickness!=null) crackList.append(crackThickness);
		crackList.append("  </CrackList>\n\n");
	}

	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	// <Cracks> for MPM Header
	public String getSettings(int MMNormals,String ContactPosition)
	{	if(settings.length()==0 && Friction==null && (MMNormals>=0 || ContactPosition==null)
			&& propagate==null && altPropagate==null && movePlane==null && propLength==null) return null;
	
		StringBuffer cracks = new StringBuffer("    <Cracks>\n");
		if(settings.length()>0) cracks.append(settings);
		if(Friction!=null) cracks.append(Friction);
		if(MMNormals<0 && ContactPosition!=null) cracks.append(ContactPosition);
		if(propagate!=null) cracks.append(propagate);
		if(altPropagate!=null) cracks.append(altPropagate);
		if(movePlane!=null) cracks.append(movePlane);
		if(propLength!=null) cracks.append(propLength);
		cracks.append("    </Cracks>\n");
		return cracks.toString();
	}
	
	// return list or null if no cracks
	public String getCrackList()
	{	appendCurrentCrack();
		if(crackList.length()==0) return null;
		return crackList.toString();
	}
	
	// read propagation criterion
	public int decodeCriterion(String code) throws Exception
	{	HashMap<String,Integer> options = new HashMap<String,Integer>(8);
		options.put("none", new Integer(0));
		options.put("max energy release", new Integer(1));
		options.put("steady state", new Integer(2));
		options.put("energy balance", new Integer(3));
		options.put("energy density", new Integer(4));
		options.put("elliptical", new Integer(5));
		options.put("max ctod", new Integer(6));
		options.put("critical err", new Integer(7));
		return doc.readIntOption(code,options,"Propagation criterion");
	}
	
	// read propagation criterion
	public int decodeDirection(String code) throws Exception
	{	HashMap<String,Integer> options = new HashMap<String,Integer>(5);
		options.put("default", new Integer(0));
		options.put("self similar", new Integer(1));
		options.put("cod normal", new Integer(2));
		options.put("cod hoop", new Integer(3));
		options.put("initial", new Integer(4));
		return doc.readIntOption(code,options,"Propagation direction");
	}


}
