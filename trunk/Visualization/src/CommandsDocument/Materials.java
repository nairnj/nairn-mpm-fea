/*
 * Materials.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 17 Dec 2012.
 * Copyright (c) 2012 RSAC Software. All rights reserved.
 */

import java.io.*;
import java.util.*;

import javax.swing.*;

import geditcom.JNFramework.*;

public class Materials
{
	private HashMap<String,Integer> matIDs;
	private int numMats;
	private StringBuffer xmldata;
	private boolean inMaterial;
	private CmdViewer doc;
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public Materials(CmdViewer cmdDoc)
	{	// save parent CmdViewer
		doc = cmdDoc;
	}
	
	public void initRunSettings()
	{
		matIDs = new HashMap<String,Integer>(10);
		numMats = 0;
		xmldata = new StringBuffer("");
		inMaterial = false;
	}
	
	//----------------------------------------------------------------------------
	// Methods
	//----------------------------------------------------------------------------
	
	// scripted material
	// needs Material ID, name, type
	public void StartMaterial(ArrayList<String> args) throws Exception
	{	// needs all args
		if(args.size()<4)
			throw new Exception("'Material' command with two few arguments: "+args);
		
		// get ID
		String newID = doc.readStringArg(args.get(1));
		if(matIDs.get(newID) != null)
			throw new Exception("Duplicate material ID ("+newID+") was used: "+args);
		
		// add to list
		numMats++;
		matIDs.put(newID, new Integer(numMats));

		// get name
		String matName = doc.readStringArg(args.get(2));
		
		// and type
		HashMap<String,Integer> options = new HashMap<String,Integer>(10);
		options.put("isotropic", new Integer(1));
		options.put("transverse 1", new Integer(2));
		options.put("transverse 2", new Integer(3));
		options.put("orthotropic", new Integer(4));
		options.put("interface", new Integer(5));
		int matInt = doc.readIntOption(args.get(3),options,null);
		if(matInt<0)
			throw new Exception("'Material' type not supported in scripting commands.\nUse XML method instead: "+args);
		
		// start the command
		xmldata.append("  <Material Type='"+matInt+"' Name='"+matName+"'>\n");
		inMaterial = true;
	}
	
	// material defined using XML commands
	public void StartXMLMaterial(String newID,String xml) throws Exception
	{
		if(matIDs.get(newID) != null)
			throw new Exception("Duplicate material ID ("+newID+") was used.");
		
		// add to list
		numMats++;
		matIDs.put(newID, new Integer(numMats));
		
		// if xml material, just add to xml data now
		if(xml!=null)
		{	xmldata.append(xml+"\n");
		}
	}
	
	// in material mode
	public void doMaterialProperty(String theCmd,ArrayList<String> args) throws Exception
	{
		// is it done
		if(theCmd.equals("done"))
		{	xmldata.append("  </Material>\n\n");
			inMaterial = false;
			return;
		}
		
		// perhaps trap some special commands
		
		// the property
		if(args.size()<2)
			throw new Exception("A material property has no value: "+args);
		
		String prop = args.get(0);
				
		// these commands require an integer
		
		// the rest are "Property double" but some need special cases
		// to account for difference in NairnFEAMPM
		if(prop.equals("a"))
			prop = "alpha";
		else if(prop.equals("aA"))
			prop = "alphaA";
		else if(prop.equals("aT"))
			prop = "alphaT";
		else if(prop.equals("ax"))
			prop = "alphax";
		else if(prop.equals("ay"))
			prop = "alphay";
		else if(prop.equals("az"))
			prop = "alphaz";
		else if(prop.equals("kA"))
			prop = "kCondA";
		else if(prop.equals("kT"))
			prop = "kCondT";
		else if(prop.equals("kx"))
			prop = "kCondx";
		else if(prop.equals("ky"))
			prop = "kCondy";
		else if(prop.equals("kz"))
			prop = "kCondz";
		
		// now add it
		xmldata.append("    <"+prop+">"+doc.readDoubleArg(args.get(1))+"</"+prop+">\n");
	}
	
	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	public String toXMLString() { return xmldata.toString(); }
	
	// numeric id for material (or <0 if not found)
	public int getMatID(String theID)
	{	Integer matnum = matIDs.get(theID);
		if(matnum == null) return -1;
		return matnum.intValue();
	}
	
	// is material active
	public boolean isInMaterial() { return inMaterial; }
}
