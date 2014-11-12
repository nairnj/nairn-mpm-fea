/*
 * Materials.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 17 Dec 2012.
 * Copyright (c) 2012 RSAC Software. All rights reserved.
 */

import java.util.*;

public class Materials
{
	private HashMap<String,Integer> matIDs;
	private int numMats;
	private StringBuffer xmldata;
	private boolean inMaterial;
	private int matType;
	private CmdViewer doc;
	private int criterion,direction,traction;
	private int altCriterion,altDirection,altTraction;
	
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
		criterion = -1;
		direction = -1;
		traction = -1;
		altCriterion = -1;
		altDirection = -1;
		altTraction = -1;
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
		HashMap<String,Integer> options = new HashMap<String,Integer>(25);
		options.put("isotropic", new Integer(1));
		options.put("transverse 1", new Integer(2));
		options.put("transverse 2", new Integer(3));
		options.put("orthotropic", new Integer(4));
		options.put("interface", new Integer(5));
		options.put("viscoelastic", new Integer(6));
		options.put("mooney", new Integer(8));
		options.put("vonmises", new Integer(9));
		options.put("isoplasticity", new Integer(9));
		options.put("bistable", new Integer(10));
		options.put("rigid", new Integer(11));
		options.put("triangulartraction", new Integer(12));
		options.put("lineartraction", new Integer(13));
		options.put("cubictraction", new Integer(14));
		options.put("hillplastic", new Integer(15));
		options.put("johnsoncook", new Integer(16));
		options.put("mgscglmaterial", new Integer(17));
		options.put("mgeosmaterial", new Integer(17));
		options.put("trilineartraction", new Integer(20));
		options.put("heanisotropic", new Integer(21));
		options.put("idealgas", new Integer(22));
		options.put("coupledsawtooth", new Integer(23));
		options.put("heisotropic", new Integer(24));
		options.put("hemgeosmaterial", new Integer(25));
		options.put("pressuretraction", new Integer(26));
		options.put("taitliquid", new Integer(27));
		options.put("neohookean", new Integer(28));
		matType = doc.readIntOption(args.get(3),options,null);
		if(matType<0)
			throw new Exception("'Material' type not yet supported in scripting commands.\nUse XML method instead: "+args);
		
		// start the command
		xmldata.append("  <Material Type='"+matType+"' Name='"+matName+"'>\n");
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
	public void doMaterialProperty(String theCmd,ArrayList<String> args,CmdViewer parent) throws Exception
	{
		// is it done
		if(theCmd.equals("done"))
		{	if(criterion>=0)
			{	xmldata.append("    <Propagate criterion='"+criterion+"'");
				if(direction>=0) xmldata.append(" direction='"+direction+"'");
				if(traction>=0) xmldata.append(" traction='"+traction+"'");
				xmldata.append("/>\n");
			}
			if(altCriterion>=0)
			{	xmldata.append("    <AltPropagate criterion='"+altCriterion+"'");
				if(altDirection>=0) xmldata.append(" direction='"+altDirection+"'");
				if(altTraction>=0) xmldata.append(" traction='"+altTraction+"'");
				xmldata.append("/>\n");
			}
			xmldata.append("  </Material>\n\n");
			inMaterial = false;
			return;
		}
		
		// perhaps trap some special commands
		
		// the property
		if(args.size()<2)
			throw new Exception("A material property has no value: "+args);
		
		String prop = args.get(0);
				
		// These are "Property double" but need special case for property name
		// to account for difference in NairnFEAMPM
		if(prop.equals("a"))
			prop = "alpha";
		else if(prop.equals("a0"))
			prop = "alpha0";
		else if(prop.equals("ad"))
			prop = "alphad";
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
		
		// these commands require an integer
		else if(prop.toLowerCase().equals("ujoption"))
		{	prop = "UJOption";
			// 0, 1, or 2
			int value = doc.readIntArg(args.get(1));
			if(value<0 || value>2)
				throw new Exception("'UJOption' must be integer 0, 1, or 2:\n"+args);
			xmldata.append("    <"+prop+">"+value+"</"+prop+">\n");
			return;
		}
		else if(prop.toLowerCase().equals("idealrubber"))
		{	prop = "IdealRubber";
			// 0 or 1
			int value = doc.readIntArg(args.get(1));
			if(value<0 || value>1)
				throw new Exception("'IdealRubber' property must be integer 0 or 1:\n"+args);
			xmldata.append("    <"+prop+">"+value+"</"+prop+">\n");
			return;
		}
		else if(prop.toLowerCase().equals("transition"))
		{	prop = "transition";
			// 1, 2, or 3
			int value = doc.readIntArg(args.get(1));
			if(value<1 || value>3)
				throw new Exception("'transition' property must be integer 1, 2, or 3:\n"+args);
			xmldata.append("    <"+prop+">"+value+"</"+prop+">\n");
			return;
		}
		else if(prop.toLowerCase().equals("reversible"))
		{	// 0 or 1 (or no or yes)
			HashMap<String,Integer> options = new HashMap<String,Integer>(2);
			options.put("no", new Integer(0));
			options.put("yes", new Integer(1));
			int value = doc.readIntOption(args.get(1),options,"reversible property");
			if(value==1)
				xmldata.append("    <reversible/>\n");
			else
				xmldata.append("    <irreversible/>\n");
			return;
		}
		else if(prop.toLowerCase().equals("direction"))
		{	if(matType==11)
			{	// in rigid material
				prop = "SetDirection";
				// 0 to 8
				int value = doc.readIntArg(args.get(1));
				if(value<0 || value>8)
					throw new Exception("Rigid 'direction' property must be integer 0 to 8:\n"+args);
				xmldata.append("    <"+prop+">"+value+"</"+prop+">\n");
			}
			else
			{	// crack growth direction
				direction = doc.cracks.decodeDirection(args.get(1));
			}
			return;
		}
		else if(prop.toLowerCase().equals("altdirection"))
		{	// alt crack growth direction
			altDirection = doc.cracks.decodeDirection(args.get(1));
			return;
		}
		else if(prop.toLowerCase().equals("criterion"))
		{	// propagation criterion
			criterion = doc.cracks.decodeCriterion(args.get(1));
			return;
		}
		else if(prop.toLowerCase().equals("altcriterion"))
		{	// propagation criterion
			altCriterion = doc.cracks.decodeCriterion(args.get(1));
			return;
		}
		else if(prop.toLowerCase().equals("traction"))
		{	String tract = doc.readStringArg(args.get(1));
			traction = doc.mats.getMatID(tract);
			if(traction==-1) traction = doc.readIntArg(args.get(1));
			if(traction<=0)
				throw new Exception("'"+args.get(0)+"' material property has unknown traction law material:\n"+args);
			return;
		}
		else if(prop.toLowerCase().equals("alttraction"))
		{	String tract = doc.readStringArg(args.get(1));
			altTraction = doc.mats.getMatID(tract);
			if(altTraction==-1)
				throw new Exception("'"+args.get(0)+"' material property has unknown traction law material:\n"+args);
			return;
		}
		else if(prop.toLowerCase().equals("temperature"))
		{	prop = "SetTemperature";
			// 0 or 1
			int value = doc.readIntArg(args.get(1));
			if(value<0 || value>1)
				throw new Exception("Rigid 'SetTemperature' property must be integer 0 or 1:\n"+args);
			xmldata.append("    <"+prop+">"+value+"</"+prop+">\n");
			return;
		}
		else if(prop.toLowerCase().equals("concentration"))
		{	prop = "SetConcentration";
			// 0 or 1
			int value = doc.readIntArg(args.get(1));
			if(value<0 || value>1)
				throw new Exception("Rigid 'SetConcentration' property must be integer 0 or 1:\n"+args);
			xmldata.append("    <"+prop+">"+value+"</"+prop+">\n");
			return;
		}
		else if(prop.toLowerCase().equals("mirrored"))
		{	prop = "mirrored";
			//-1, 0, or 1
			int value = doc.readIntArg(args.get(1));
			if(value<-1) value = -1;
			if(value>1) value = 1;
			xmldata.append("    <"+prop+">"+value+"</"+prop+">\n");
			return;
		}
		else if(prop.toLowerCase().equals("allowscracks"))
		{	prop = "allowsCracks";
			// 0 or 1
			int value = doc.readIntArg(args.get(1));
			if(value<0 || value>1)
				throw new Exception("The 'allowsCracks' property must be integer 0 or 1:\n"+args);
			xmldata.append("    <"+prop+">"+value+"</"+prop+">\n");
			return;
		}
		else if(prop.toLowerCase().equals("sharematfield"))
		{	String smf = doc.readStringArg(args.get(1));
			int smat = doc.mats.getMatID(smf);
			if(smat==-1)
				throw new Exception("'"+args.get(0)+"' material property has unknown material ID:\n"+args);
			xmldata.append("    <shareMatField>"+smat+"</shareMatField>\n");
			return;
		}
		
		// These require special handling
		else if(prop.toLowerCase().equals("friction"))
		{	xmldata.append(parent.doFriction(args,2));
			return;
		}
		else if(prop.toLowerCase().equals("interface"))
		{	xmldata.append(parent.doImperfectInterface(args,2));
			return;
		}
		else if(prop.toLowerCase().equals("color"))
		{	double red = doc.readDoubleArg(args.get(1));
			double blue=red,green=red,alpha=1.;
			if(args.size()>2)
			{	if(args.size()<4)
					throw new Exception("Color material needs, 1, 3, or 4 values.");
				green = doc.readDoubleArg(args.get(3));
				blue = doc.readDoubleArg(args.get(2));
				if(args.size()>4) alpha = doc.readDoubleArg(args.get(4));
			}
			xmldata.append("    <color red='"+red+"' green='"+green+
								"' blue='"+blue+"' alpha='"+alpha+"'/>\n");
			return;
		}
		else if(prop.toLowerCase().equals("artificialvisc"))
		{	String value = doc.readStringArg(args.get(1)).toLowerCase();
			if(!value.equals("on") && !value.equals("off"))
				throw new Exception("ArtificialVisc must be on or off.");
			if(value.equals("on"))
				xmldata.append("    <ArtificialVisc/>\n");
			return;
		}
		else if(prop.toLowerCase().equals("settingfunction") || 
				prop.toLowerCase().equals("settingfunction1") ||
				prop.toLowerCase().equals("settingfunctionx"))
		{	xmldata.append("    <SettingFunction>"+doc.readStringArg(args.get(1))+"</SettingFunction>\n");
			return;
		}
		else if(prop.toLowerCase().equals("settingfunction2") ||
				prop.toLowerCase().equals("settingfunctiony"))
		{	xmldata.append("    <SettingFunction2>"+doc.readStringArg(args.get(1))+"</SettingFunction2>\n");
			return;
		}
		else if(prop.toLowerCase().equals("settingfunction3") ||
				prop.toLowerCase().equals("settingfunctionz"))
		{	xmldata.append("    <SettingFunction3>"+doc.readStringArg(args.get(1))+"</SettingFunction3>\n");
			return;
		}
		else if(prop.toLowerCase().equals("valuefunction"))
		{	xmldata.append("    <ValueFunction>"+doc.readStringArg(args.get(1))+"</ValueFunction>\n");
			return;
		}
		else if(prop.toLowerCase().equals("function"))
		{	xmldata.append("    <function>"+doc.readStringArg(args.get(1))+"</function>\n");
			return;
		}
		else if(prop.toLowerCase().equals("hardening"))
		{	xmldata.append("    <Hardening>"+doc.readStringArg(args.get(1))+"</Hardening>\n");
			return;
		}
		
		// now add it (if not done already)
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
