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
	private boolean crackFriction;
	private double cx,cy;
	
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
	// #4 = friction setting (friction=) of "traction"
	// #5 - traction material (mat=)
	public void StartCrack(ArrayList<String> args) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);
		
		// save current crack
		appendCurrentCrack();
		
		// reset global crack properties
		currentCrack = new StringBuffer("");
		crackFixed = false;
		crackFriction = false;
		
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
			{	// get traction material (but can't check is traction law)
				if(args.size()>5)
					mat = doc.mats.getMatID(doc.readStringArg(args.get(5)));
				if(mat==-1)
					throw new Exception("'NewCrack' has unknown traction law material:\n"+args);
			}
			else
			{	// look for custom friction setting
				throw new Exception("'NewCrack' does not yet read friction settings:\n"+args);
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

	// when done or start new crack, append current one
	public void appendCurrentCrack()
	{
		if(currentCrack==null) return;
		
		crackList.append("  <CrackList");
		if(crackFixed) crackList.append(" type='fixed'");
		if(crackFriction)
		{	// add friction and/or imperfect interface settings
			
		}
		
		// finish up
		crackList.append(">\n"+currentCrack.toString()+"  </CrackList>\n\n");
	}

	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	// <Cracks> for MPM Header
	public String getSettings()
	{	if(settings.length()==0) return null;
		return "    <Cracks>\n"+settings.toString()+"    </Cracks>\n";
	}
	
	// return list or null if no cracks
	public String getCrackList()
	{	appendCurrentCrack();
		if(crackList.length()==0) return null;
		return crackList.toString();
	}

}
