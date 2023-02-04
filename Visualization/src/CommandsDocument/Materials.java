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
	private StringBuffer taukGk;
	private StringBuffer taukKk;
	private boolean inMaterial;
	private int matType;
	private int ntaus,nGs;
	private int ntausK,nKs;
	private CmdViewer doc;
	private int criterion,direction,traction;
	private int altCriterion,altDirection,altTraction;
	private double matDamping,matPIC;
	private double GT0,GA0,KT0,en0,ell0;
	private StringBuffer taukGT;
	private StringBuffer taukGA;
	private StringBuffer taukKT;
	private StringBuffer taukn;
	private StringBuffer taukell;
	private int ntauGT,nGTs;
	private int ntauGA,nGAs;
	private int ntauKT,nKTs;
	private int ntaun,nns;
	private int ntauell,nells;
	private int whichOne;

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
		taukGk = null;
		inMaterial = false;
		criterion = -1;
		direction = -1;
		traction = -1;
		altCriterion = -1;
		altDirection = -1;
		altTraction = -1;
		matDamping = -1.1e12;
		matPIC = -1.;
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
		HashMap<String,Integer> options = new HashMap<String,Integer>(60);
		options.put("isotropic", new Integer(1));
		options.put("transverse", new Integer(2));
		options.put("orthotropic", new Integer(4));
		if(doc.isFEA())
		{	options.put("interface", new Integer(5));
		}
		else
		{	options.put("tiviscoelastic", new Integer(5));
			options.put("viscoelastic", new Integer(7));
			options.put("mooney", new Integer(8));
			options.put("vonmises", new Integer(9));
			options.put("isoplasticity", new Integer(9));
			options.put("bistable", new Integer(10));
			options.put("rigid", new Integer(11));
			options.put("rigidbc", new Integer(11));
			options.put("triangulartraction", new Integer(12));
			options.put("lineartraction", new Integer(13));
			options.put("cubictraction", new Integer(14));
			options.put("hillplastic", new Integer(15));
			options.put("johnsoncook", new Integer(16));			// historical only
			options.put("mgscglmaterial", new Integer(17));			// historical only
			options.put("mgeosmaterial", new Integer(17));			// historical only
			options.put("trilineartraction", new Integer(20));
			options.put("heanisotropic", new Integer(21));
			options.put("idealgas", new Integer(22));
			options.put("coupledsawtooth", new Integer(23));
			options.put("coupledtraction", new Integer(23));
			options.put("heisotropic", new Integer(24));
			options.put("hemgeosmaterial", new Integer(25));
			options.put("pressuretraction", new Integer(26));
			options.put("taitliquid", new Integer(27));
			options.put("neohookean", new Integer(28));
			options.put("clampedneohookean", new Integer(29));
			options.put("phasetransition", new Integer(30));
			options.put("reactionphase", new Integer(31));
			options.put("jwlplusplus", new Integer(32));
			options.put("mixedmodetraction", new Integer(33));
			options.put("exponentialtraction", new Integer(34));
			options.put("rigidcontact", new Integer(35));
			options.put("rigidblock", new Integer(36));
			options.put("mooneymembrane", new Integer(40));
			options.put("isosoftening", new Integer(50));
			options.put("transisosoftening", new Integer(51));
			options.put("isoplasticsoftening", new Integer(53));
			options.put("orthosoftening", new Integer(54));
			options.put("isoplasticinterface", new Integer(55));
			options.put("orthoplasticsoftening", new Integer(56));
			options.put("isophasefieldsoftening", new Integer(57));
			options.put("isodamagemechanics", new Integer(58));
			options.put("phasetransition", new Integer(30));
			options.put("ignorecontact", new Integer(60));
			options.put("coulombfriction", new Integer(61));
			options.put("adhesivefriction", new Integer(63));
			options.put("linearinterface", new Integer(62));
			options.put("liquidcontact", new Integer(64));
			options.put("nonlinearinterface", new Integer(65));
			options.put("debondinginterface", new Integer(66));
		}
		matType = doc.readIntOption(args.get(3),options,null);
		
		// deprecated materials
		if(matType<0)
		{	options = new HashMap<String,Integer>(6);
			options.put("transverse 1", new Integer(2));
			options.put("transverse 2", new Integer(3));
			if(doc.isMPM())
			{	options.put("tiviscoelastic 1", new Integer(5));
				options.put("tiviscoelastic 2", new Integer(6));
				options.put("transisosoftening 1", new Integer(51));
				options.put("transisosoftening 2", new Integer(52));
			}
			matType = doc.readIntOption(args.get(3),options,null);
		}
		
		// if not found then error
		if(matType<0)
		{	if(doc.isFEA())
				throw new Exception("'Material' type not supported in FEA commands.\nUse XML method if possible: "+args);
			throw new Exception("'Material' type not supported in MPM commands.\nUse XML method if possible: "+args);
		}
		
		// start the command
		xmldata.append("  <Material Type='"+matType+"' Name='"+matName+"'>\n");
		inMaterial = true;
		
		// start viscoelastic
		taukGT = null;
		taukGA = null;
		taukKT = null;
		taukn = null;
		taukell = null;
		taukGk = null;
		taukKk = null;
		if(matType==7)
		{	taukGk = new StringBuffer("");
			ntaus = 0;
			nGs = 0;
			whichOne = -1;
			taukKk = new StringBuffer("");
			ntausK = 0;
			nKs = 0;
		}
		else if(matType==5 || matType==6)
		{	whichOne = 0;
		}
		else whichOne = -2;
			
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
			
			// damping
			if(matDamping>-1.e12 || matPIC>=0.)
			{	xmldata.append("    <PDamping");
				if(matPIC>=0.)
					xmldata.append(" PIC='"+doc.formatDble(matPIC)+"'");
				if(matDamping>-1.e12)
					xmldata.append(">"+doc.formatDble(matDamping)+"</PDamping>\n");
				else
					xmldata.append("/>\n");
			}
			
			// viscoelastic
			if(taukGk!= null)
			{	if(ntaus!=nGs)
					throw new Exception("A viscoelastic material property does not have same number of tauks and Gks");
				xmldata.append("    <ntaus>"+ntaus+"</ntaus>\n");
				xmldata.append(taukGk);
			}
			if(taukGT!=null)
			{	if(ntauGT!=nGTs)
					throw new Exception("A viscoelastic GT property does not have same number of tauks and Pks");
				xmldata.append("    <GT0>"+doc.formatDble(GT0)+"</GT0>\n");
				xmldata.append("    <ntaus>"+ntauGT+"</ntaus>\n");
				xmldata.append(taukGT);
			}
			if(taukGA!=null)
			{	if(ntauGA!=nGAs)
					throw new Exception("A viscoelastic GA property does not have same number of tauks and Pks");
				xmldata.append("    <GA0>"+doc.formatDble(GA0)+"</GA0>\n");
				xmldata.append("    <ntaus>"+ntauGA+"</ntaus>\n");
				xmldata.append(taukGA);
			}
			if(taukKT!=null)
			{	if(ntauKT!=nKTs)
					throw new Exception("A viscoelastic KT property does not have same number of tauks and Pks");
				xmldata.append("    <KT0>"+doc.formatDble(KT0)+"</KT0>\n");
				xmldata.append("    <ntaus>"+ntauKT+"</ntaus>\n");
				xmldata.append(taukKT);
			}
			if(taukn!=null)
			{	if(ntaun!=nns)
					throw new Exception("A viscoelastic en property does not have same number of tauks and Pks");
				xmldata.append("    <en0>"+doc.formatDble(en0)+"</en0>\n");
				xmldata.append("    <ntaus>"+ntaun+"</ntaus>\n");
				xmldata.append(taukn);
			}
			if(taukell!=null)
			{	if(ntauell!=nells)
					throw new Exception("A viscoelastic ell property does not have same number of tauks and Pks");
				xmldata.append("    <ell0>"+doc.formatDble(ell0)+"</ell0>\n");
				xmldata.append("    <ntaus>"+ntauell+"</ntaus>\n");
				xmldata.append(taukell);
			}
			if(taukKk!= null && ntausK>0)
			{	if(ntausK!=nKs)
					throw new Exception("A viscoelastic material property does not have same number of tauKks and Kks");
				xmldata.append("    <ntausK>"+ntausK+"</ntausK>\n");
				xmldata.append(taukKk);
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
		else if(prop.toLowerCase().equals("pressurelaw"))
		{	prop = "pressureLaw";
			int value = doc.readIntArg(args.get(1));
			xmldata.append("    <"+prop+">"+value+"</"+prop+">\n");
			return;
		}
		else if(prop.toLowerCase().equals("transition"))
		{	prop = "transition";
			// 1, 2, or 3 (or dilation, distortion, or vonmises)
			HashMap<String,Integer> options = new HashMap<String,Integer>(2);
			options.put("dilation", new Integer(1));
			options.put("distortion", new Integer(2));
			options.put("vonmises", new Integer(3));
			int value = doc.readIntOption(args.get(1),options,"transition property");
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
		else if(prop.toLowerCase().equals("solidphase"))
		{	String phase = doc.readStringArg(args.get(1));
			int phaseNum = doc.mats.getMatID(phase);
			if(phaseNum==-1) phaseNum = doc.readIntArg(args.get(1));
			if(phaseNum<=0)
				throw new Exception("'"+args.get(0)+"' material property has unknown solid phase material:\n"+args);
			xmldata.append("    <SolidPhase>"+phaseNum+"</SolidPhase>\n");
			return;
		}
		else if(prop.toLowerCase().equals("liquidphase"))
		{	String phase = doc.readStringArg(args.get(1));
			int phaseNum = doc.mats.getMatID(phase);
			if(phaseNum==-1) phaseNum = doc.readIntArg(args.get(1));
			if(phaseNum<=0)
				throw new Exception("'"+args.get(0)+"' material property has unknown liquid phase material:\n"+args);
			xmldata.append("    <LiquidPhase>"+phaseNum+"</LiquidPhase>\n");
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
				throw new Exception("Rigid 'concentration' property must be integer 0 or 1:\n"+args);
			xmldata.append("    <"+prop+">"+value+"</"+prop+">\n");
			return;
		}
		else if(prop.toLowerCase().equals("porepressure"))
		{	prop = "SetConcentration";
			// 0 or 1
			int value = doc.readIntArg(args.get(1));
			if(value<0 || value>1)
				throw new Exception("Rigid 'porepressure' property must be integer 0 or 1:\n"+args);
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
		else if(prop.toLowerCase().equals("contact"))
		{	xmldata.append(parent.doContactLaw(args,2));
			return;
		}
		else if(prop.toLowerCase().equals("friction"))
		{	// Deprecated - use contact instead
			xmldata.append(parent.doFriction(args,2));
			return;
		}
		else if(prop.toLowerCase().equals("interface"))
		{	// Deprecated - use contact instead
			xmldata.append(parent.doImperfectInterface(args,2));
			return;
		}
		else if(prop.toLowerCase().equals("color"))
		{	double red = doc.readDoubleArg(args.get(1));
			double blue=red,green=red,alpha=1.;
			if(args.size()>2)
			{	if(args.size()<4)
					throw new Exception("Color material needs, 1, 3, or 4 values.");
				green = doc.readDoubleArg(args.get(2));
				blue = doc.readDoubleArg(args.get(3));
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
		else if(prop.toLowerCase().equals("initialpressure"))
		{	xmldata.append("    <InitialPressure>"+doc.readStringArg(args.get(1))+"</InitialPressure>\n");
			return;
		}
		else if(prop.toLowerCase().equals("hardening"))
		{	xmldata.append("    <Hardening>"+doc.readStringArg(args.get(1))+"</Hardening>\n");
			return;
		}
		else if(prop.toLowerCase().equals("initiation"))
		{	xmldata.append("    <Initiation>"+doc.readStringArg(args.get(1))+"</Initiation>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningi"))
		{	xmldata.append("    <SofteningI>"+doc.readStringArg(args.get(1))+"</SofteningI>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningii"))
		{	xmldata.append("    <SofteningII>"+doc.readStringArg(args.get(1))+"</SofteningII>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningai"))
		{	xmldata.append("    <SofteningAI>"+doc.readStringArg(args.get(1))+"</SofteningAI>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningtii"))
		{	xmldata.append("    <SofteningTII>"+doc.readStringArg(args.get(1))+"</SofteningTII>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningaii"))
		{	xmldata.append("    <SofteningAII>"+doc.readStringArg(args.get(1))+"</SofteningAII>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningxx"))
		{	xmldata.append("    <SofteningXX>"+doc.readStringArg(args.get(1))+"</SofteningXX>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningyy"))
		{	xmldata.append("    <SofteningYY>"+doc.readStringArg(args.get(1))+"</SofteningYY>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningzz"))
		{	xmldata.append("    <SofteningZZ>"+doc.readStringArg(args.get(1))+"</SofteningZZ>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningxyx"))
		{	xmldata.append("    <SofteningXYX>"+doc.readStringArg(args.get(1))+"</SofteningXYX>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningxyy"))
		{	xmldata.append("    <SofteningXYY>"+doc.readStringArg(args.get(1))+"</SofteningXYY>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningxzx"))
		{	xmldata.append("    <SofteningXZX>"+doc.readStringArg(args.get(1))+"</SofteningXZX>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningxzz"))
		{	xmldata.append("    <SofteningXZZ>"+doc.readStringArg(args.get(1))+"</SofteningXZZ>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningyzy"))
		{	xmldata.append("    <SofteningYZY>"+doc.readStringArg(args.get(1))+"</SofteningYZY>\n");
			return;
		}
		else if(prop.toLowerCase().equals("softeningyzz"))
		{	xmldata.append("    <SofteningYZZ>"+doc.readStringArg(args.get(1))+"</SofteningYZZ>\n");
			return;
		}
		else if(prop.toLowerCase().equals("tauk"))
		{	double tk = doc.readDoubleArg(args.get(1));
			String taustr = "    <tauk>"+doc.formatDble(tk)+"</tauk>\n";
			switch(whichOne)
			{	case -1:
					taukGk.append(taustr);
					ntaus++;
					break;
				case 1:
					taukGT.append(taustr);
					ntauGT++;
					break;
				case 2:
					taukGA.append(taustr);
					ntauGA++;
					break;
				case 3:
					taukKT.append(taustr);
					ntauKT++;
					break;
				case 4:
					taukn.append(taustr);
					ntaun++;
					break;
				case 5:
					taukell.append(taustr);
					ntauell++;
					break;
				default:
					throw new Exception("tauk found at invalid location");
			}
			return;
		}
		else if(prop.toLowerCase().equals("taukk"))
		{	double tk = doc.readDoubleArg(args.get(1));
			taukKk.append("    <tauKk>"+doc.formatDble(tk)+"</tauKk>\n");
			ntausK++;
			return;
		}
		else if(prop.toLowerCase().equals("gk"))
		{	if(taukGk==null)
				throw new Exception("Gk found at invalid location");
			double gk = doc.readDoubleArg(args.get(1));
			taukGk.append("    <Gk>"+doc.formatDble(gk)+"</Gk>\n");
			nGs++;
			return;
		}
		else if(prop.toLowerCase().equals("kk"))
		{	if(taukKk==null)
				throw new Exception("Kk found at invalid location");
			double kk = doc.readDoubleArg(args.get(1));
			taukKk.append("    <Kk>"+doc.formatDble(kk)+"</Kk>\n");
			nKs++;
			return;
		}
		else if(prop.toLowerCase().equals("pk"))
		{	double pk = doc.readDoubleArg(args.get(1));
			String pkstr = "    <Pk>"+doc.formatDble(pk)+"</Pk>\n";
			switch(whichOne)
			{	case 1:
					taukGT.append(pkstr);
					nGTs++;
					break;
				case 2:
					taukGA.append(pkstr);
					nGAs++;
					break;
				case 3:
					taukKT.append(pkstr);
					nKTs++;
					break;
				case 4:
					taukn.append(pkstr);
					nns++;
					break;
				case 5:
					taukell.append(pkstr);
					nells++;
					break;
				default:
					throw new Exception("Pk found at invalid location");
			}
			return;
		}
		else if(prop.toLowerCase().equals("gt0"))
		{	if(taukGT!=null)
				throw new Exception("Found two entries for GT0");
			if(whichOne<0 || whichOne>5)
				throw new Exception("Found invalid GT0 property");
			GT0 = doc.readDoubleArg(args.get(1));
			taukGT = new StringBuffer("");
			ntauGT = 0;
			nGTs = 0;
			whichOne = 1;
			return;
		}
		else if(prop.toLowerCase().equals("ga0"))
		{	if(taukGA!=null)
				throw new Exception("Found two entries for GA0");
			if(whichOne<0 || whichOne>5)
				throw new Exception("Found invalid GA0 property");
			GA0 = doc.readDoubleArg(args.get(1));
			taukGA = new StringBuffer("");
			ntauGA = 0;
			nGAs = 0;
			whichOne = 2;
			return;
		}
		else if(prop.toLowerCase().equals("kt0"))
		{	if(taukKT!=null)
				throw new Exception("Found two entries for KT0");
			if(whichOne<0 || whichOne>5)
				throw new Exception("Found invalid KT0 property");
			KT0 = doc.readDoubleArg(args.get(1));
			taukKT = new StringBuffer("");
			ntauKT = 0;
			nKTs = 0;
			whichOne = 3;
			return;
		}
		else if(prop.toLowerCase().equals("en0"))
		{	if(taukn!=null)
				throw new Exception("Found two entries for en0");
			if(whichOne<0 || whichOne>5)
				throw new Exception("Found invalid en0 property");
			en0 = doc.readDoubleArg(args.get(1));
			taukn = new StringBuffer("");
			ntaun = 0;
			nns = 0;
			whichOne = 4;
			return;
		}
		else if(prop.toLowerCase().equals("ell0"))
		{	if(taukell!=null)
				throw new Exception("Found two entries for ell0");
			if(whichOne<0 || whichOne>5)
				throw new Exception("Found invalid ell0 property");
			ell0 = doc.readDoubleArg(args.get(1));
			taukell = new StringBuffer("");
			ntauell = 0;
			nells = 0;
			whichOne = 5;
			return;
		}
		else if(prop.toLowerCase().equals("matdamping"))
		{	matDamping = doc.readDoubleArg(args.get(1));
			return;
		}
		else if(prop.toLowerCase().equals("matpic"))
		{	matPIC = doc.readDoubleArg(args.get(1));
			return;
		}
		
		// now add it (if not done already)
		double mprop = doc.readDoubleArg(args.get(1));
		xmldata.append("    <"+prop+">"+doc.formatDble(mprop)+"</"+prop+">\n");
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
