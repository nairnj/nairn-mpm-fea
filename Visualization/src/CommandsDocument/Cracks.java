/*
 * Cracks.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 23 SEP 2013.
 * Copyright (c) 2013 RSAC Software. All rights reserved.
 */

import java.util.*;

import geditcom.JNFramework.JNUtilities;

public class Cracks
{
	private CmdViewer doc;
	
	private StringBuffer settings = null;
	
	private StringBuffer crackList = null;
	private StringBuffer currentCrack = null;
	private boolean crackFixed;
	private String crackFriction;
	private String crackTractProp;
	private String crackThickness;
	private String movePlane;
	private double cx,cy;
	private String Friction;
	private String altPropagate;
	private String propagate;
	private String propLength;
	
	// for 3D cracks
	private ArrayList<Double> coords = new ArrayList<Double>(30);
	private Hashtable<String,Integer> ptIDs = new Hashtable<String,Integer>();
	private ArrayList<Integer> facetKeys = new ArrayList<Integer>(30);
	private ArrayList<Double> facetLengths = new ArrayList<Double>(30);
	private ArrayList<Integer> facetTraction = new ArrayList<Integer>(30);
	boolean currentPlane;
	
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
			
			settings.append(" terms='"+terms+"'");
		}
		
		// use grid for energies
		int useGrid=0;
		if(args.size()>3)
		{	// options
			HashMap<String,Integer> options = new HashMap<String,Integer>(10);
			options.put("no", new Integer(0));
			options.put("yes", new Integer(1));
			options.put("false", new Integer(0));
			options.put("true", new Integer(1));
			useGrid = doc.readIntOption(args.get(3),options,"J grid energy option");
			
			settings.append(" gridenergy='"+useGrid+"'");
		}
		
		// finish the command
		settings.append("/>\n");
	}
	
	// NewCrack x,y,<#3>,<#4>,<#5>
	// #3 = material (tip=), exterior (tip=-2), or fixed (type='fixed'), or free (no tip)
	// #4 = friction setting (friction=) or "traction"
	// #5 - traction material (mat=)
	// For 3D is is #1,#2,#3,#4 for x,y,z of first point and label (all required)
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
		crackTractProp = null;
		currentPlane = false;
		
		// 3D cracks simpler
		if(doc.isMPM3D())
		{	if(args.size()<5)
				throw new Exception("'NewCrack' has too few parameters:\n"+args);
			double ptx = doc.readDoubleArg(args.get(1));
			double pty = doc.readDoubleArg(args.get(2));
			double ptz = doc.readDoubleArg(args.get(3));
			String key1 = doc.readStringArg(args.get(4));
			
			// optional contact law
			if(args.size()>5)
			{	String law = doc.readStringArg(args.get(5));
				int lawMat = doc.mats.getMatID(law);
				if(lawMat>0)
					crackFriction=" law='"+lawMat+"'";
				else
				{	throw new Exception("'NewCrack' has invalid contact law material:\n"+args);
				}
			}
			
			// start coordinates list
			coords.clear();
			coords.add(new Double(ptx));
			coords.add(new Double(pty));
			coords.add(new Double(ptz));
			
			// list of key points labels
			ptIDs.clear();
			ptIDs.put(key1,new Integer(1));
			
			// clear facets
			facetKeys.clear();
			facetLengths.clear();
			facetTraction.clear();
			
			// all done
			return;
		}
		
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
		
		// contact law ID, Deprecated friction setting, or "traction"
		int mat=-1;
		if(args.size()>4)
		{	String frict = doc.readStringArg(args.get(4));
		
			// look for traction law
			if(frict.toLowerCase().equals("traction"))
			{	// get traction material (but can't check if traction law)
				if(args.size()>5)
				{	mat = doc.mats.getMatID(doc.readStringArg(args.get(5)));
					if(mat==-1)
						throw new Exception("'NewCrack' has unknown traction law material:\n"+args);
				}
			}
			
			// did not set traction, so look for contact law
			if(mat==-1)
			{	// look for contact law material
				int lawMat = doc.mats.getMatID(frict);
				if(lawMat>0)
					crackFriction=" law='"+lawMat+"'";
				else
				{	// look for custom friction setting
					try
					{	ArrayList<String> fargs = new ArrayList<String>();
						fargs.add("Friction");
						fargs.add(frict);
						crackFriction=doc.doFriction(fargs,3);
					}
					catch(Exception e)
					{	throw new Exception("'NewCrack' has invalid contact law material:\n"+args);
					}
				}
			}
		}
		
		// start current crack XML
		currentCrack.append("    <pt x='"+doc.formatDble(ptx)+"' y='"+doc.formatDble(pty)+"'");
		if(tip!=-1) currentCrack.append(" tip='"+tip+"'");
		if(mat!=-1) currentCrack.append(" mat='"+mat+"'");
		currentCrack.append("/>\n");
		
		cx = ptx;
		cy = pty;
	}
	
	// Add keypoint to a 3D crack
	public void AddCrackKeypoint(ArrayList<String> args) throws Exception
	{
		// check analysis type and verify unneste
		if(currentCrack==null)
		{	throw new Exception("'CrackKeypoint' commands must come after a 'NewCrack' command:\n"+args);
		}
		if(!doc.isMPM3D())
		{	throw new Exception("'CrackKeypoint' command only alowed for 3D cracks:\n"+args);
		}

		// read x,y,z, label in #1 to #2 (required)
		if(args.size()<5)
			throw new Exception("'CrackKeypoint' has too few parameters:\n"+args);
		double ptx = doc.readDoubleArg(args.get(1));
		double pty = doc.readDoubleArg(args.get(2));
		double ptz = doc.readDoubleArg(args.get(3));
		String keyn = doc.readStringArg(args.get(4));
		
		// add three coordinates
		coords.add(new Double(ptx));
		coords.add(new Double(pty));
		coords.add(new Double(ptz));
		
		// add a number
		int ptnum = coords.size()/3;
		ptIDs.put(keyn,new Integer(ptnum));
	}

	
	// For 3D crack facet is #1,#2,#3,<#4>,<#5> for ptIDs and maxLength and traction law
	public void AddCrackFacet(ArrayList<String> args,boolean isPlane) throws Exception
	{
		String cmd = isPlane ? "CrackPlane" : "CrackFacet" ;
		
		// check analysis type and verify unneste
		if(currentCrack==null)
		{	throw new Exception("'"+cmd+"' commands must come after a 'NewCrack' command:\n"+args);
		}
		if(!doc.isMPM3D())
		{	throw new Exception("'"+cmd+"' command only alowed for 3D cracks:\n"+args);
		}
		
		// a plane must be the only element
		if(isPlane && facetKeys.size()>0)
		{	throw new Exception("'CrackPlane' command must be only facet for the crack:\n"+args);
		}
		currentPlane = isPlane;
		
		// read three labels (required)
		// read x,y,z, label in #1 to #2 (required)
		if(args.size()<4)
			throw new Exception("'"+cmd+"' has too few parameters:\n"+args);
		String nd1 = doc.readStringArg(args.get(1));
		Integer kp1 = ptIDs.get(nd1);
		String nd2 = doc.readStringArg(args.get(2));
		Integer kp2 = ptIDs.get(nd2);
		String nd3 = doc.readStringArg(args.get(3));
		Integer kp3 = ptIDs.get(nd3);
		if(kp1==null || kp2==null ||kp3==null)
			throw new Exception("'"+cmd+"' command has an invalid keypoint IDs:\n"+args);
			
		// optional length
		double maxLength = -1.;
		if(args.size()>4)
			maxLength = doc.readDoubleArg(args.get(4));
		
		// look for traction law material
		int mat = -1;
		if(args.size()>5)
		{	mat = doc.mats.getMatID(doc.readStringArg(args.get(5)));
			if(mat==-1)
				throw new Exception("'"+cmd+"' has unknown traction law material:\n"+args);
		}
		
		// keypoints
		facetKeys.add(kp1);
		facetKeys.add(kp2);
		facetKeys.add(kp3);
		
		// length
		facetLengths.add(new Double(maxLength));
		
		// traction law
		facetTraction.add(new Integer(mat));
	}

	// 0: GrowCrack x,y,<tip>,<mat>
	// 1: GrowCrackLine x,y,segs,<tip>,<mat>
	// 2: GrowCrackArc x1,y1,x2,y2,segs,ang1,ang2,<tip>,<mat>
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
		
		// Arc gets two more coordinates
		double ptx2=0.,pty2=0.;
		if(option==2)
		{	if(args.size()<5)
				throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
			ptx2 = doc.readDoubleArg(args.get(3));
			pty2 = doc.readDoubleArg(args.get(4));
			param = 5;
		}
		
		// segments (GrowCrackLine and GrowCrackArc)
		int segs=1;
		if(option>0)
		{	if(args.size()>param)
			{	segs = doc.readIntArg(args.get(param));
				param++;
			}
		}
		
		// angles (GrowCrackArc)
		double ang1=0.,ang2=0.;
		if(option==2)
		{	if(args.size()<param+2)
				throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
			ang1 = doc.readDoubleArg(args.get(param));
			ang2 = doc.readDoubleArg(args.get(param+1));
			param += 2;
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
		{	currentCrack.append("    <pt x='"+doc.formatDble(ptx)+"' y='"+doc.formatDble(pty)+"'");
			if(tip!=-1) currentCrack.append(" tip='"+tip+"'");
		}
		else if(option==1)
		{	currentCrack.append("    <Line xmin='"+doc.formatDble(cx)+"' ymin='"+doc.formatDble(cy)+"'");
			currentCrack.append(" xmax='"+doc.formatDble(ptx)+"' ymax='"+doc.formatDble(pty)+"' resolution='"+segs+"'");
			if(tip!=-1) currentCrack.append(" end_tip='"+tip+"'");
		}
		else
		{	currentCrack.append("    <Circle xmin='"+doc.formatDble(ptx)+"' ymin='"+doc.formatDble(pty)+"'");
			currentCrack.append(" xmax='"+doc.formatDble(ptx2)+"' ymax='"+doc.formatDble(pty2)+"' resolution='"+segs+"'");
			currentCrack.append(" start_angle='"+doc.formatDble(ang1)+"' end_angle='"+doc.formatDble(ang2)+"'");
			if(tip!=-1) currentCrack.append(" end_tip='"+tip+"'");
		}
		if(mat!=-1) currentCrack.append(" mat='"+mat+"'");
		currentCrack.append("/>\n");
		
		// save location
		if(option==2)
		{	// find the final point on the elipse
			//double xmid = 0.5*(ptx+ptx2);
			//double ymid = 0.5*(pty+pty2);
			//double ea = 0.5*(ptx2-ptx);
			//double eb = 0.5*(pty2-pty);
			cx = ptx2;
			cy = ptx2;
		}
		else
		{	cx = ptx;
			cy = pty;
		}
	}
	
	// set global friction command
	public void setFriction(String cmd) { Friction = cmd; }

	// set global friction command
	public void setCrackFriction(ArrayList<String> args,String cmd) throws Exception
	{	// must be in a crack
		if(currentCrack==null)
			throw new Exception("'CrackInterface' command must be in an active crack:\n"+args);
		crackFriction = cmd;
	}

	// CrackThickness #1
	public void doProagateLength(ArrayList<String> args) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);
		
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' command with too few arguments:\n"+args);
		
		double pl = doc.readDoubleArg(args.get(1));
		propLength = "      <PropagateLength>"+doc.formatDble(pl)+"</PropagateLength>\n";
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
		
		double ct = doc.readDoubleArg(args.get(1));
		crackThickness = "    <Thickness>"+doc.formatDble(ct)+"</Thickness>\n";
	}
	
	// CrackThickness #1
	public void doCrackTractionProp(ArrayList<String> args) throws Exception
	{	// MPM Only
		doc.requiresMPM(args);
		
		// must be in a crack
		if(currentCrack==null)
			throw new Exception("'"+args.get(0)+"' command must be in an active crack:\n"+args);
		
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' command with too few arguments:\n"+args);
		
		// get law ID
		int lawnum = doc.mats.getMatID(doc.readStringArg(args.get(1)));
		if(lawnum <= 0)
			throw new Exception("'" + args.get(0) + "' traction law has unknown material ID:\n" + args);

		crackTractProp = " Tprop='"+lawnum+"'";
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
		if(crackFriction!=null) crackList.append(crackFriction);
		
		if(!doc.isMPM3D())
		{	// 2D cracks
			if(crackFixed) crackList.append(" type='fixed'");
			if(crackTractProp!=null) crackList.append(crackTractProp);
		
			// finish up
			crackList.append(">\n"+currentCrack.toString());
			if(crackThickness!=null) crackList.append(crackThickness);
			crackList.append("  </CrackList>\n\n");
		}
		else if(currentPlane)
		{	// a single plane element for a 3D crack
			crackList.append(">\n     <Plane");
			
			System.out.println(ptIDs);
			System.out.println(coords);
			System.out.println(facetKeys);
			
			for(int p=0;p<3;p++)
			{	int pt1 = facetKeys.get(p).intValue();
				int index = 3*(pt1-1);
				System.out.println(p+","+pt1+","+facetKeys.get(p)+","+index);
				double ptX = coords.get(index).doubleValue();
				double ptY = coords.get(index+1).doubleValue();
				double ptZ = coords.get(index+2).doubleValue();
				crackList.append(" V"+p+"x='"+JNUtilities.formatDouble(ptX)+"'");
				crackList.append(" V"+p+"y='"+JNUtilities.formatDouble(ptY)+"'");
				crackList.append(" V"+p+"z='"+JNUtilities.formatDouble(ptZ)+"'");
			}
			
			// length
			double length = facetLengths.get(0).doubleValue();
			if(length>0.)
				crackList.append(" length='"+JNUtilities.formatDouble(length)+"'");
			
			// traction law
			int tnum = facetTraction.get(0).intValue();
			if(tnum>0)
				crackList.append(" mat='"+tnum+"'");
			
			// finish up
			crackList.append("/>\n");
			crackList.append("  </CrackList>\n\n");
		}
		else
		{	// 3D crack with keypoints and facets
			crackList.append("  <CrackList>\n    <Mesh>\n");
			
			// Node list
			crackList.append("      <NodeList>\n");
			int j = 0;
			double ptX,ptY,ptZ;
			while(j<coords.size())
			{	ptX=coords.get(j).doubleValue();
				ptY=coords.get(j+1).doubleValue();
				ptZ=coords.get(j+2).doubleValue();
				j+=3;
				crackList.append("        <pt");
				crackList.append(" x='"+JNUtilities.formatDouble(ptX)+"'");
				crackList.append(" y='"+JNUtilities.formatDouble(ptY)+"'");
				crackList.append(" z='"+JNUtilities.formatDouble(ptZ)+"'");
				crackList.append("/>\n");
			}
			crackList.append("      </NodeList>\n");
			
			// elements in blocks of 4 in facets, traction in tractions
			crackList.append("      <ElementList>\n");
			int jt = 0;
			j = 0;
			int pt1,pt2,pt3,tnum;
			double length;
			while(j<facetKeys.size())
			{	pt1 = facetKeys.get(j).intValue();
				pt2 = facetKeys.get(j+1).intValue();
				pt3 = facetKeys.get(j+2).intValue();
				j += 3;
				
				length = facetLengths.get(jt).doubleValue();
				tnum = facetTraction.get(jt).intValue();
				jt++;
				
				crackList.append("        <elem type='1'");
				if(length>0)
					crackList.append(" length='"+JNUtilities.formatDouble(length)+"'");
				if(tnum>0)
					crackList.append(" czm='"+tnum+"'");
				crackList.append(">"+pt1+","+pt2+","+pt3+"</elem>\n");
			}
			crackList.append("      </ElementList>\n");
			
			// finish up
			crackList.append("    </Mesh>\n");
			crackList.append("  </CrackList>\n\n");
		}

	}

	// when done or start new crack, append current one
	public void appendXMLCrack(String xmlData)
	{
		// output current crack first
		appendCurrentCrack();
		
		crackList.append("  <CrackList>\n");
		crackList.append(xmlData);
		crackList.append("  </CrackList>\n\n");
		
		currentCrack = null;
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
