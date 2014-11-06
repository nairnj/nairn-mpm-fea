/*
 * CmdViewer
 * NairnFEAMPMViz Application
 * 
 * Created 
 */

import java.awt.Toolkit;
import java.awt.event.*;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.Set;

import javax.swing.*;

import geditcom.JNFramework.*;

public class CmdViewer extends JNCmdTextDocument
{
	private static final long serialVersionUID = 1L;
	
	private ConsolePane soutConsole;
	private DocViewer linkedResults=null;
	
	// Actions
	private BgAnalysisAction bgAnalysisCammand = new BgAnalysisAction();
	private CheckAnalysisAction checkAnalysisCammand = new CheckAnalysisAction();
	private InterpretCommandsAction interpretCammand = new InterpretCommandsAction();
	private ShowPartnerAction showPartnerCommand = new ShowPartnerAction();
	private ExportXMLAction exportXMLCommand = new ExportXMLAction();
	private StopCurrentModeAction stopModeCommand = new StopCurrentModeAction();
	
	// analysis runner
	public NFMAnalysis nfmAnalysis = null;
	private int openMesh;
	private boolean useBackground;

	// global settings
	private String title;
	private String username;
	private StringBuffer header;
	private int np;
	private int processors=1;
	private int lnameEl;
	private HashMap<String,String> xmldata = null;
	private HashMap<String,String> entities = null;
	public Materials mats = null;
	public Areas areas = null;
	public Regions regions = null;
	public MPMGrid gridinfo = null;
	private FEABCs feaBCs = null;
	private MPMGridBCs mpmGridBCs = null;
	public MPMParticleBCs mpmParticleBCs = null;
	public Cracks cracks = null;
	private StringBuffer outFlags;
	private int mpmMethod;
	private String shapeMethod;
	private String archiveRoot;
	private String archiveTime;
	private String timeStep;
	private String maxTime;
	private String globalArchive;
	private String damping;
	private String pdamping;
	private String fbDamping;
	private String pfbDamping;
	private String leaveLimit;
	private String ptsPerElement;
	private String diffusion;
	private String conduction;
	private String gravity;
	private StringBuffer mpmOrder;
	private StringBuffer crackOrder;
	private boolean mpmMeshToFile;
	private String feaTemp;
	private double stressFreeTemp;
	private double rampTime;
	private double rampStart;
	private double rampDiff;
	private boolean stopCommand;
	private double MMVmin;
	private int MMDcheck;
	private int MMNormals;		// <0 means no multimaterial mode
	private double MMRigidBias;
	private double MMAzimuth;
	private double MMPolar;
	private String ContactPosition;
	private String FrictionMM;
	private StringBuffer customTasks;
	private String currentCustomTask = null;
	
	// scripting
	private boolean runningScript;
	private HashMap<String,Object> objs = null;
	private File scriptOutput;
	private CmdViewer theScript;

	// constants
	public static final int PLANE_STRAIN=0;
	public static final int PLANE_STRESS=1;
	public static final int AXI_SYM=2;
	public static final int THREE_DIM=3;
	public static final int BEGIN_MPM_TYPES=9;
	public static final int PLANE_STRAIN_MPM=10;
	public static final int PLANE_STRESS_MPM=11;
	public static final int THREED_MPM=12;
	public static final int AXI_SYM_MPM=13;
	
	public static final int NO_ELEMENT=0;
	
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public CmdViewer(String aType)
	{	super(aType,null,new ConsolePane());
		soutConsole=(ConsolePane)console;
	
		setFramePrefs("Commands Window Width",600,"Commands Window Height",800);
		
		makeMenuBar();
		
		// tool bar icons
		addDefaultToolBar(this);
		addToolBarIcon(null,"saveDocument","Save this document to a file.",this);
		
		Class<?> baseClass=JNApplication.main.getClass();
		addToolBarBar();
		ImageIcon goNext=new ImageIcon(baseClass.getResource("Resources/go-next.png"));
		addToolBarIcon(goNext,null,"Run FEA or MPM Analysis.",getRunAnalysisAction());
		ImageIcon goLast=new ImageIcon(baseClass.getResource("Resources/go-last.png"));
		addToolBarIcon(goLast,null,"Check mesh for FEA or MPM Analysis.",checkAnalysisCammand);
		ImageIcon scriptIcon=new ImageIcon(baseClass.getResource("Resources/scripticon.png"));
		addToolBarIcon(scriptIcon,null,"Interpret calculation or control script commands.",interpretCammand);
		ImageIcon doStop=new ImageIcon(baseClass.getResource("Resources/process-stop.png"));
		addToolBarIcon(doStop,null,"Stop currently running FEA or MPM Analysis.",stopModeCommand);

		addToolBarBar();
		ImageIcon showRes=new ImageIcon(baseClass.getResource("Resources/image-x-generic.png"));
		addToolBarIcon(showRes,null,"Show associated simulation results (if available).",showPartnerCommand);

		finishFrameworkWindow(true);
		
		// create persistent objects
		mats = new Materials(this);
		areas = new Areas(this);
		regions = new Regions(this);
		feaBCs = new FEABCs(this);
		gridinfo = new MPMGrid(this);
		mpmGridBCs = new MPMGridBCs(this);
		mpmParticleBCs = new MPMParticleBCs(this);
		cracks = new Cracks(this);
	}
	
	// make menu bar on launch
	protected void makeMenuBar()
	{
		// Create the menu bar
		JMenuBar menuBar = new JMenuBar();
		if(!JNApplication.isMacLNF())
			menuBar.add(defaultApplicationMenu());		// Application menu
		JMenu fileMenu = defaultFileMenu(this);
		fileMenu.add(exportXMLCommand);
		menuBar.add(fileMenu);				// File menu
		
		// Edit menu
		JMenu menu = defaultEditMenu(true);
		menuBar.add(menu);
		menu.addSeparator();
		menu.add(getGoToLineAction());
		
		// Analyze menu
		menu = new JMenu("Analyze");
		menuBar.add(menu);
		menu.add(getRunAnalysisAction("Run FEA/MPM Analysis"));
		menu.add(bgAnalysisCammand);
		menu.add(checkAnalysisCammand);
		menu.add(interpretCammand);
		menu.addSeparator();
		menu.add(stopModeCommand);
		
		// Window
		menu = new JMenu("Window");
		menuBar.add(menu);
		if(JNApplication.isMacLNF())
		{	menu.add(JNApplication.main.getOpenHelpAction());
		}
		menu.add(showPartnerCommand);
		menu.addSeparator();
		setWindowMenu(menu);

		// add the menu bar
		setJMenuBar(menuBar);
	}
	
	// default file menu referring to document target
	public static JMenu defaultFileMenu(JNDocument target)
	{	JMenu fileMenu=target.defaultFileMenu();
		JMenuItem newMPM=fileMenu.getItem(1);
		newMPM.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N,
				JNApplication.menuKeyMask()+ActionEvent.SHIFT_MASK));
		return fileMenu;
	}
	
	// default tool bar icons
	public static void addDefaultToolBar(JNDocument target)
	{
		// tool bar icons
		target.addToolBarBar();
		target.addToolBarIcon(null,null,"Open the preferences window.",JNApplication.main.getOpenPreferencesAction());
		target.addToolBarIcon(null,null,"Open the help information window.",JNApplication.main.getOpenHelpAction());
		
		target.addToolBarBar();
		target.addToolBarIcon(null,"openDocument","Open a saved document file.",JNApplication.main);
		Class<?> baseClass=JNApplication.main.getClass();
		ImageIcon newMPM=new ImageIcon(baseClass.getResource("Resources/document-new.png"));
		target.addToolBarIcon(newMPM,"newDocumentMPMCmd","Create a new MPM commands document.",JNApplication.main);
		ImageIcon newFEA=new ImageIcon(baseClass.getResource("Resources/document-newfea.png"));
		target.addToolBarIcon(newFEA,"newDocumentFEACmd","Create a new FEA command document.",JNApplication.main);
	}
	
	//----------------------------------------------------------------------------
	// Interpret commands to XML input commands
	//----------------------------------------------------------------------------
	
	// default command run method passed to custom one
	public void runAnalysis() { runNFMAnalysis(false,NFMAnalysis.FULL_ANALYSIS,null); }
	
	public void runNFMAnalysis(boolean doBackground,int runType,CmdViewer scriptDoc)
	{
		// only allowed if the commands have been saved
		if(getFile()==null)
		{	JNApplication.appBeep();
			JOptionPane.showMessageDialog(this,"The input commands have to be saved to a file before running an analysis.");
			return;
		}
		
		// self call sets scriptDoc to self
		if(scriptDoc==null)
			theScript = this;
		else
			theScript = scriptDoc;
		
		// create once
		if(nfmAnalysis == null)
			nfmAnalysis = new NFMAnalysis(this);
		
		// what if process is current running?
		if(nfmAnalysis.isRunning())
		{	JNApplication.appBeep();
			String message="An FEA or MPM process is currently running.\nDo you want stop it and start a new process?";
			int result=JOptionPane.showConfirmDialog(this, message, "Continue?", JOptionPane.OK_CANCEL_OPTION);
			if(result==JOptionPane.CANCEL_OPTION) return;
			nfmAnalysis.stopRunning();
		}
		
		// save this document (commands saved before each run - preference would be better)
		//if(!saveDocument()) return;
		
		// check if XML file
		soutConsole.clear();
		int offset=cmdField.getCommands().indexOf("<?xml ");
		if(offset<0 || offset>10)
		{	// interpret commands
			useBackground = doBackground;
			openMesh = runType;
			// call in super class initiates command interpertation
			super.runAnalysis();
			
			// when interpretaiont done, will launch the analysis in analysisFinished()
			return;
		}
		else
		{	// look for processors command in XML commands
			processors = 1;
			offset = cmdField.getCommands().indexOf("<!--processors ");
			if(offset>0)
			{	int endoffset = cmdField.getCommands().indexOf("-->",offset);
				if(endoffset>0)
				{	String procs = cmdField.getCommands().substring(offset+15,endoffset);
					Scanner getProcs=new Scanner(procs);
					if(getProcs.hasNextInt()) processors =  getProcs.nextInt();
					if(processors<1) processors = 1;
					getProcs.close();
				}
			}
		}
		
		// launch analysis with DTD commands in the field
		nfmAnalysis.runNFMAnalysis(doBackground,runType,cmdField.getCommands(),
					soutConsole,processors,theScript.getOutputFile());
	}
	
	// when analysis is done, proceed with calculations (if OKO)
	public void analysisFinished(boolean status)
	{	// give up on error
		if(status==false || stopCommand==true) return;
		
		if(runningScript)
		{	toFront();
			return;
		}
		
		// launch analysis with DTD commands in the field
		nfmAnalysis.runNFMAnalysis(useBackground,openMesh,buildXMLCommands(),
					soutConsole,processors,theScript.getOutputFile());
	}
	
	// initialize when interpreting commands
	public void initRunSettings() throws Exception
	{
		title = "NairnFEAMPMViz Calculations";
		username = null;
		header = new StringBuffer("");
		np = -1;
		lnameEl = NO_ELEMENT;
		xmldata = new HashMap<String,String>(10);
		entities = new HashMap<String,String>(10);
		mats.initRunSettings();
		areas.initRunSettings();
		regions.initRunSettings();
		feaBCs.initRunSettings();
		gridinfo.initRunSettings();
		mpmGridBCs.initRunSettings();
		mpmParticleBCs.initRunSettings();
		cracks.initRunSettings();
		mpmMeshToFile = true;
		outFlags = null;
		mpmOrder = null;
		crackOrder = null;
		mpmMethod = 0;
		shapeMethod = "uGIMP";
		archiveRoot = "    <ArchiveRoot>Results/data.</ArchiveRoot>\n";
		archiveTime = "";
		timeStep = "    <TimeStep units='ms'>1e15</TimeStep>\n";
		maxTime = "";
		globalArchive="";
		feaTemp = null;
		stressFreeTemp = 0.;
		rampStart = -1.;
		rampTime = -2.;
		stopCommand = false;
		damping = null;
		pdamping = null;
		fbDamping = null;
		pfbDamping = null;
		leaveLimit = null;
		ptsPerElement = null;
		diffusion = null;
		conduction = null;
		gravity = null;
		MMVmin = 0.0;
		MMDcheck = 0;
		MMNormals = -1;
		MMRigidBias = 1.0;
		MMAzimuth = 0.0;
		MMPolar = 0.0;
		ContactPosition = null;
		FrictionMM = null;
		customTasks = new StringBuffer("");
		
		runningScript = false;
		objs = new HashMap<String,Object>(10);
		scriptOutput = null;
		
		// is it called from a script?
		if(theScript!=this)
		{	variables.put("_ScriptMode_",new Double(1.));
			variables.putAll(theScript.getVariables());
			variablesStrs.putAll(theScript.getVariablesStrs());
		}
	}
	
	// handle commands
	public void doCommand(String theCmd,ArrayList<String> args) throws Exception
	{	
		// if script, switch to script commands
		if(runningScript)
		{	doScriptCommand(theCmd,args);
		}
		
		else if(theCmd.equals("script"))
		{	if(openMesh!=NFMAnalysis.SCRIPT_ONLY && openMesh!=NFMAnalysis.INTERPRET_ONLY)
				throw new Exception("Scripts can only run by using the 'Interpret Commands...'");
			runningScript = true;
			openMesh = NFMAnalysis.SCRIPT_ONLY;
		}
		
		else if(mats.isInMaterial())
		{	// commands go to material class when material (keep this option first)
			
			// but first see if language control command
			// (which means cannot match any material property)
			try
			{	super.doCommand(theCmd, args);
			}
			catch(Exception e)
			{	mats.doMaterialProperty(theCmd,args,this);
			}
		}
			
		else if(theCmd.equals("title"))
		{	// set analysis title
			if(args.size()<2)
				throw new Exception("'Title' command does not have a title");
			title = readStringArg(args.get(1));
		}
		
		else if(theCmd.equals("name"))
		{	// set analysis title
			if(args.size()<2)
				throw new Exception("'Name' command does not have a name");
			username = readStringArg(args.get(1));
		}
		
		else if(theCmd.equals("header"))
			header.append(readVerbatim("endheader"));
		
		else if(theCmd.equals("comment"))
		{	int i;
			header.append("Comment: ");
			for(i=1;i<args.size();i++)
			{	if(i>1) header.append(", ");
				header.append(readStringArg(args.get(i)));
			}
			header.append("\n");
		}
		
		else if(theCmd.equals("processors"))
		{	processors = readIntArg(args.get(1));
			if(processors<1) processors = 1;
		}
		
		else if(theCmd.equals("analysis"))
			doAnalysis(args);
		
		else if(theCmd.equals("mpmmethod"))
			doMPMMethod(args);
		
		else if(theCmd.equals("archive"))
			doArchive(args,false);
		
		else if(theCmd.equals("globalarchive"))
			doGlobalArchive(args);
		
		else if(theCmd.equals("toarchive"))
			doToArchive(args);
		
		else if(theCmd.equals("archiveunique"))
			doArchive(args,true);
		
		else if(theCmd.equals("archivetime"))
			doArchiveTime(args);
		
		else if(theCmd.equals("globalarchivetime"))
			doGlobalArchiveTime(args);

		else if(theCmd.equals("timestep"))
			doTimeStep(args);
		
		else if(theCmd.equals("maximumtime"))
			doMaxTime(args);
		
		else if(theCmd.equals("element"))
			doElement(args);
		
		else if(theCmd.equals("xmldata"))
			doXmldata(args);
		
		else if(theCmd.equals("entity"))
			doEntity(args);
		
		else if(theCmd.equals("area"))
			areas.StartArea(args);
		
		else if(theCmd.equals("endarea"))
			areas.EndArea(args);
		
		else if(theCmd.equals("path"))
			areas.StartPath(args);
		
		else if(theCmd.equals("endpath"))
			areas.EndPath(args);
		
		else if(theCmd.equals("paths"))
			areas.AddPaths(args);
		
		else if(theCmd.equals("keypoint"))
			areas.AddKeypoint(args);
		
		else if(theCmd.equals("keypoints"))
			areas.AddKeypoints(args);
		
		else if(theCmd.equals("fliptriangles"))
			areas.setFlipTriangles(args);
		
		else if(theCmd.equals("fixline"))
			feaBCs.StartFixLine(args);
		
		else if(theCmd.equals("selectline"))
		{	feaBCs.StartFixLine(args);
			feaBCs.AddSelect(args);
			feaBCs.EndFixLine(args);
		}
		
		else if(theCmd.equals("endfixline"))
			feaBCs.EndFixLine(args);
		
		else if(theCmd.equals("fixpoint"))
			feaBCs.StartFixPoint(args);
		
		else if(theCmd.equals("selectpoint"))
		{	feaBCs.StartFixPoint(args);
			feaBCs.AddSelect(args);
			feaBCs.EndFixPoint(args);
		}
		
		else if(theCmd.equals("endfixpoint"))
			feaBCs.EndFixPoint(args);
		
		else if(theCmd.equals("displacement"))
			feaBCs.AddDisplacement(args);
		
		else if(theCmd.equals("load"))
		{	if(isFEA())
				feaBCs.AddLoad(args);
			else
				mpmParticleBCs.AddCondition(args, MPMParticleBCs.ADD_LOAD);
		}
		
		else if(theCmd.equals("rotate"))
		{	if(isFEA())
				feaBCs.AddRotate(args);
			else
				regions.AddRotate(args);
		}
		
		else if(theCmd.equals("stress"))
			feaBCs.AddStress(args);
		
		else if(theCmd.equals("periodic"))
			feaBCs.AddPeriodic(args);
		
		else if(theCmd.equals("cracktip"))
			feaBCs.CrackTip(args);
		
		else if(theCmd.equals("resequence"))
			feaBCs.Resequence(args);
		
		else if(theCmd.equals("select"))
			feaBCs.AddSelect(args);
		
		else if(theCmd.equals("moveline"))
			mpmGridBCs.StartMoveLine(args,MPMGridBCs.MOVELINE_BC);
		
		else if(theCmd.equals("movearc"))
			mpmGridBCs.StartMoveLine(args,MPMGridBCs.MOVEARC_BC);
		
		else if(theCmd.equals("movebox"))
			mpmGridBCs.StartMoveBox(args,MPMGridBCs.MOVEBOX_BC);
		
		else if(theCmd.equals("endmoveline"))
			mpmGridBCs.EndMoveBlock(args,MPMGridBCs.MOVELINE_BC);
		
		else if(theCmd.equals("endmovearc"))
			mpmGridBCs.EndMoveBlock(args,MPMGridBCs.MOVEARC_BC);
		
		else if(theCmd.equals("endmovebox"))
			mpmGridBCs.EndMoveBlock(args,MPMGridBCs.MOVEBOX_BC);
		
		else if(theCmd.equals("boundaryid"))
			mpmGridBCs.SetBoundaryID(args);
		
		else if(theCmd.equals("velocity"))
			mpmGridBCs.AddVelocity(args);
		
		else if(theCmd.equals("concentration"))
			mpmGridBCs.AddTempConc(args,MPMGridBCs.ADD_CONCENTRATION);

		else if(theCmd.equals("temperature"))
		{	if(isFEA())
				doTemperature(args);
			else
				mpmGridBCs.AddTempConc(args,MPMGridBCs.ADD_TEMPERATURE);
		}

		else if(theCmd.equals("loadline"))
			mpmGridBCs.StartMoveLine(args,MPMGridBCs.LOADLINE_BC);
		
		else if(theCmd.equals("loadarc"))
			mpmGridBCs.StartMoveLine(args,MPMGridBCs.LOADARC_BC);
		
		else if(theCmd.equals("loadbox"))
			mpmGridBCs.StartMoveBox(args,MPMGridBCs.LOADBOX_BC);
		
		else if(theCmd.equals("loadrect"))
			mpmParticleBCs.StartLoadRect(args);
		
		else if(theCmd.equals("loadtype"))
			mpmParticleBCs.doLoadType(args);
		
		else if(theCmd.equals("endloadline"))
			mpmParticleBCs.EndLoadBlock(args,MPMParticleBCs.LOADLINE_BC);
		
		else if(theCmd.equals("endloadarc"))
			mpmParticleBCs.EndLoadBlock(args,MPMParticleBCs.LOADARC_BC);
		
		else if(theCmd.equals("endloadrect"))
			mpmParticleBCs.EndLoadBlock(args,MPMParticleBCs.LOADRECT_BC);
		
		else if(theCmd.equals("endloadbox"))
			mpmParticleBCs.EndLoadBlock(args,MPMParticleBCs.LOADBOX_BC);
		
		else if(theCmd.equals("traction"))
			mpmParticleBCs.AddCondition(args,MPMParticleBCs.ADD_TRACTION);
		
		else if(theCmd.equals("heatflux"))
			mpmParticleBCs.AddCondition(args,MPMParticleBCs.ADD_HEATFLUX);
		
		else if(theCmd.equals("concentrationflux"))
			mpmParticleBCs.AddCondition(args,MPMParticleBCs.ADD_CONCENTRATIONFLUX);
		
		else if(theCmd.equals("origin"))
		{	if(regions.isInBMPRegion())
				regions.setOrigin(args);
			else
				areas.setOrigin(args);
		}
		
		else if(theCmd.equals("intensity"))
			regions.AddIntensity(args);
		
		else if(theCmd.equals("material"))
			mats.StartMaterial(args);
		
		else if(theCmd.equals("output"))
			doOutput(args);
		
		else if(theCmd.equals("region"))
			regions.StartRegion(args);
		
		else if(theCmd.equals("bmpregion"))
			regions.StartBMPRegion(args);
		
		else if(theCmd.equals("endregion"))
			regions.EndRegion(args);
		
		else if(theCmd.equals("hole"))
			regions.StartHole(args);
		
		else if(theCmd.equals("endhole"))
			regions.EndHole(args);
		
		else if(theCmd.equals("rect"))
			regions.AddRectOrOval(args,"Rect");
		
		else if(theCmd.equals("oval"))
			regions.AddRectOrOval(args,"Oval");
		
		else if(theCmd.equals("polypt"))
			regions.AddPolypoint(args);
		
		else if(theCmd.equals("box"))
			regions.AddBox(args,"Box");
		
		else if(theCmd.equals("cylinder"))
			regions.AddBox(args,"Cylinder");
		
		else if(theCmd.equals("sphere"))
			regions.AddBox(args,"Sphere");
		
		else if(theCmd.equals("gridhoriz"))
			gridinfo.doGridAxis(args,0);
		
		else if(theCmd.equals("gridvert"))
			gridinfo.doGridAxis(args,1);
		
		else if(theCmd.equals("griddepth"))
			gridinfo.doGridAxis(args,2);
		
		else if(theCmd.equals("gridrect"))
			gridinfo.doGridRect(args);
		
		else if(theCmd.equals("gridthickness"))
			gridinfo.doGridThickness(args);
		
		else if(theCmd.equals("stressfreetemp"))
			doStressFreeTemp(args);
		
		else if(theCmd.equals("thermalramp"))
			doThermalRamp(args);
		
		else if(theCmd.equals("rampstart"))
			doRampStart(args);
		
		else if(theCmd.equals("damping"))
			doDamping(args,"Damping");
		
		else if(theCmd.equals("pdamping"))
			doDamping(args,"PDamping");
		
		else if(theCmd.equals("feedbackdamping"))
			doFBDamping(args,"FeedbackDamping");
		
		else if(theCmd.equals("pfeedbackdamping"))
			doFBDamping(args,"PFeedbackDamping");
		
		else if(theCmd.equals("leavelimit"))
			doLeaveLimit(args);
		
		else if(theCmd.equals("diffusion"))
			doDiffusion(args);
		
		else if(theCmd.equals("conduction"))
			doConduction(args);
		
		else if(theCmd.equals("gravity"))
			doGravity(args);
		
		else if(theCmd.equals("ptsperelement"))
			doPtsPerElement(args);
		
		else if(theCmd.equals("multimaterialmode"))
			doMultimaterialMode(args);
		
		else if(theCmd.equals("contactposition"))
			doContactPosition(args);
		
		else if(theCmd.equals("friction"))
			doFriction(args,0);
		
		else if(theCmd.equals("imperfectinterface"))
			doImperfectInterface(args,0);
		
		else if(theCmd.equals("frictionmm"))
			doFriction(args,1);
		
		else if(theCmd.equals("imperfectinterfacemm"))
			doImperfectInterface(args,1);

		else if(theCmd.equals("crackinterface"))
			doImperfectInterface(args,3);

		else if(theCmd.equals("jcontour"))
			cracks.doJContour(args);
		
		else if(theCmd.equals("newcrack"))
			cracks.StartCrack(args);

		else if(theCmd.equals("growcrack"))
			cracks.GrowCrack(args,0);

		else if(theCmd.equals("growcrackline"))
			cracks.GrowCrack(args,1);
		
		else if(theCmd.equals("growcrackarc"))
			cracks.GrowCrack(args,2);
		
		else if(theCmd.equals("crackthickness"))
			cracks.doCrackThickness(args);
		
		else if(theCmd.equals("propagate"))
			cracks.doPropagate(args,false);
		
		else if(theCmd.equals("altpropagate"))
			cracks.doPropagate(args,true);
		
		else if(theCmd.equals("moveplane"))
			cracks.doMovePlane(args);
		
		else if(theCmd.equals("propagatelength"))
			cracks.doProagateLength(args);
		
		else if(theCmd.equals("customtask"))
			doCustomTask(args);

		else if(theCmd.equals("parameter"))
			doParameter(args);

		else if(theCmd.equals("stop"))
		{	super.doCommand(theCmd,args);
			stopCommand = true;
		}
		
		else
		{	//System.out.println(args);
			super.doCommand(theCmd, args);
		}
	}
	
	// handle commands
	public void doScriptCommand(String theCmd,ArrayList<String> args) throws Exception
	{
		if(theCmd.equals("open"))
		{	// Open objName,path (omit path for dialog, can be relative path)
			if(args.size()<2)
				throw new Exception("The first argument in an 'Open' command must provide an object variable.\n"+args);
			String objectVar = args.get(1);
			if(!validObjectName(objectVar))
				throw new Exception("The first argument in an 'Open' command must be valid object name.\n"+args);
			
			// file by path or null
			File oneDoc = null;
			if(args.size()>2)
			{	oneDoc = scriptPath(readStringArg(args.get(2)),args,false);
				oneDoc = new File(oneDoc.getCanonicalPath());
				
				// see if already open
				JNDocument currentDoc = NairnFEAMPMViz.main.findDocument(oneDoc);
				if(currentDoc!=null)
				{	currentDoc.setVisible(true);
					currentDoc.toFront();
					objs.put(objectVar,currentDoc);
					return;
				}
			}
			
			// open now (exit if cancelled or error)
			JNDocument currentDoc = NairnFEAMPMViz.main.frontDocument();
			NairnFEAMPMViz.main.openDocument(oneDoc);
			if(currentDoc == NairnFEAMPMViz.main.frontDocument())
			{	// open failed or was canceled
				running = false;
				return;
			}
			objs.put(objectVar,NairnFEAMPMViz.main.frontDocument());
		}
		
		else if(theCmd.equals("openfolder"))
		{	// openFolder - string var name,title
			if(args.size()<2)
				throw new Exception("The first argument in an 'OpenFolder' command must be string variable name.\n"+args);
			
			String varName = args.get(1);
			if(!JNEvaluatorStrs.validStrVariableName(varName))
				throw new Exception("The first argument in an 'OpenFolder' command must be  a valid string variable name.\n"+args);
			
			// optional dialog title
			String fldrTitle = "Select a folder";
			if(args.size()>2)
			{	String userTitle = readStringArg(args.get(2));
				if(userTitle.length()>0) fldrTitle = userTitle;
			}
			
			// if path use it, otherwise dialog box
			String fldrPath = "";
			if(args.size()>3)
			{	File oneFldr = scriptPath(readStringArg(args.get(3)),args,true);
			
				// only need to create if does not exist
				if(!oneFldr.exists())
				{	if(!oneFldr.mkdirs())
						throw new Exception("File error creating the folder(s).\n"+args);
				}
			
				// get file path
				fldrPath = oneFldr.getCanonicalPath();
			}
			else
			{	JFileChooser fldrChooser=new JFileChooser();
				fldrChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				NFMVPrefs.setWorkspace(fldrChooser);
				fldrChooser.setDialogTitle(fldrTitle);
				int result = fldrChooser.showSaveDialog(this);
				if(result == JFileChooser.APPROVE_OPTION)
				{	fldrPath = fldrChooser.getSelectedFile().getPath();
				}
			}
			
			// save in variable with terminal path delimiter
			if(fldrPath.length()>0)
			{	if(JNApplication.isWindowsOS())
					fldrPath += "\\";
				else
					fldrPath += "/";
			}
			variablesStrs.put(varName,fldrPath);
		}
		
		else
		{	// look for object command
			String objCmd = args.get(0);
			int dot = objCmd.indexOf(".");
			if(dot>0)
			{	String objName = objCmd.substring(0,dot);
				Object obj = objs.get(objName);
				if(obj!=null)
				{	doObjectCommand(obj,objCmd.substring(dot+1).toLowerCase(),args);
					return;
				}
			}
			
			//System.out.println(args);
			super.doCommand(theCmd, args);
		}
	}
	
	// handle commands to an object
	public void doObjectCommand(Object obj,String theCmd,ArrayList<String> args) throws Exception
	{
		if(theCmd.equals("interpret"))
		{	// interpret the commands (no arguments)
			
			if(!obj.getClass().equals(CmdViewer.class))
				throw new Exception("The 'interpret' command can only by used on commands documents.\n"+args);
			
			scriptOutput = null;
			((CmdViewer)obj).runNFMAnalysis(false,NFMAnalysis.INTERPRET_ONLY,this);
			
			// wait for interpret to be done
			while(true)
			{	Thread.sleep(100);
				if(!((CmdViewer)obj).isRunning()) break;
			}
		}
		
		else if(theCmd.equals("run"))
		{	// run obj,outpath
			if(!obj.getClass().equals(CmdViewer.class))
				throw new Exception("The 'run' command can only by used on command documents.\n"+args);
						
			// need to provide path to save the file
			if(args.size()<3)
				throw new Exception("'Run' command needs object name and output file path.\n"+args);
			
			String objectVar = args.get(1);
			if(!validObjectName(objectVar))
				throw new Exception("The first argument in an 'Open' command must be valid object name.\n"+args);
			
			// get path
			File outDoc = scriptPath(readStringArg(args.get(2)),args,false);
			scriptOutput = new File(outDoc.getCanonicalPath());
			if(!scriptOutput.getParentFile().exists())
				throw new Exception("The folder selected for output does not exist.\n"+args);

			// start analysis
			((CmdViewer)obj).runNFMAnalysis(false,NFMAnalysis.FULL_ANALYSIS,this);
			
			// wait for interpret to be done
			while(true)
			{	Thread.sleep(1000);
				if(!((CmdViewer)obj).isRunning()) break;
			}
			
			// set obj to output document
			objs.put(objectVar,((DocViewer)NairnFEAMPMViz.main.frontDocument()).resDoc);
		}
		
		else if(theCmd.equals("export"))
		{	// run obj,outpath
			if(!obj.getClass().equals(CmdViewer.class))
				throw new Exception("The 'export' command can only by used on commands documents.\n"+args);
						
			// file by path or null
			File oneDoc = null;
			if(args.size()>1)
			{	String fPath = readStringArg(args.get(1));
				if(fPath.length()<2)
					throw new Exception("'export' command has empty path name.\n"+args);
				if(fPath.charAt(0)!='/' && fPath.charAt(1)!=':')
					oneDoc = new File(getFile().getParent(),fPath);
				else
					oneDoc = new File(fPath);
			}
			
			if(!((CmdViewer)obj).exportOutput(oneDoc,null))
				throw new Exception("The 'export' command failed.\n"+args);					
		}

		else
			throw new Exception("An unrecognized object command.\n"+args);
	}
	
	// object names begin in letter (not '#')
	// Rest letters, numbers, and underscore
	public static boolean validObjectName(String v)
	{	// need at least letter
		if(v.length()<1) return false;
		// other letters letter or number
		for(int i=0;i<v.length()-1;i++)
		{	char c = v.charAt(i);
			if ((c > 'z' || c < 'a')  && (c > 'Z' || c < 'A'))
			{	// first must be letter
				if(i==0) return false;
				// others can be numbers of underscore
				if((c > '9' || c < '0') && c!='_') return false;
			}
		}
		return true;
	}
	
	// decode argument to path for a script
	// allows relative or full path and allows Mac/Linux or Windows
	// if file exists, it must be folder or file is wantFolder is true or false
	public File scriptPath(String fPath,ArrayList<String> args,boolean wantFolder) throws Exception
	{	// empty is not allowed
		if(fPath.length()==0)
			throw new Exception("'"+args.get(0)+"' command has empty path name.\n"+args);
		
		// Mac/Linux full path begins in / at at lease 1 more character
		File oneDoc = null;
		if(fPath.charAt(0)=='/')
		{	// needs at least on more letter
			if(fPath.length()<2)
				throw new Exception("'"+args.get(0)+"' command has incomplete full path.\n"+args);
			oneDoc = new File(fPath);
		}
		else if(fPath.length()>3)
		{	// full Windows full path needs "c:\a" or at least 4 letters with : an \ in 2nd and 3rd
			if(fPath.charAt(1)==':' && fPath.charAt(2)=='\\')
				oneDoc = new File(fPath);
		}
		
		// it is a relative path
		if(oneDoc==null) oneDoc = new File(getFile().getParent(),fPath);
		
		// if already exists, it better be a folder
		if(oneDoc.exists())
		{	if(wantFolder)
			{	if(!oneDoc.isDirectory())
					throw new Exception("A specified folder name already exists but is not a folder.\n"+args);
			}
			else
			{	if(oneDoc.isDirectory())
					throw new Exception("A specified file name already exists but it is a folder.\n"+args);
			}
		}
		
		// return it
		return oneDoc;

	}

	// Analysis (type),(element)
	public void doAnalysis(ArrayList<String> args) throws Exception
	{
		if(np>=0)
			throw new Exception("Only one 'Analysis' command is allowed.");
		
		if(args.size()<2)
			throw new Exception("'Analysis' command has no argument.");
		
		// options
		HashMap<String,Integer> options = new HashMap<String,Integer>(10);
		options.put("plane strain", new Integer(PLANE_STRAIN));
		options.put("plane strain fea", new Integer(PLANE_STRAIN));
		options.put("plane stress", new Integer(PLANE_STRESS));
		options.put("plane stress fea", new Integer(PLANE_STRESS));
		options.put("axisymmetric", new Integer(AXI_SYM));
		options.put("axisymmetric fea", new Integer(AXI_SYM));
		options.put("plane strain mpm", new Integer(PLANE_STRAIN_MPM));
		options.put("plane stress mpm", new Integer(PLANE_STRESS_MPM));
		options.put("axisymmetric mpm", new Integer(AXI_SYM_MPM));
		options.put("3d mpm", new Integer(THREED_MPM));
		
		// read it
		np = readIntOption(args.get(1),options,"Analysis type");
		
		// optional element second
		if(args.size()>2)
		{	args.remove(1);
			doElement(args);
		}
		else if(lnameEl==NO_ELEMENT)
        {   if(np>BEGIN_MPM_TYPES)
            	lnameEl=ElementBase.FOUR_NODE_ISO;
        }
	}
	
	// MPMMethd #1,#2
	public void doMPMMethod(ArrayList<String> args) throws Exception
	{
	    // MPM Only
		requiresMPM(args);

	    // read analysis type
		if(args.size()<2)
			throw new Exception("'MPMMethod' has too few parameters:\n"+args);
		
		// options
		HashMap<String,Integer> options = new HashMap<String,Integer>(10);
		options.put("usf", new Integer(0));
		options.put("usavg", new Integer(2));
		options.put("szs", new Integer(3));
		mpmMethod = readIntOption(args.get(1),options,"MPM update method");
		
		// shape functions
		if(args.size()>2)
		{	String shape = readStringArg(args.get(2)).toLowerCase();
			if(shape.equals("gimp") || shape.equals("ugimp"))
				shapeMethod = "uGIMP";
			else if(shape.equals("cpdi") || shape.equals("lcpdi"))
				shapeMethod = "lCPDI";
			else if(shape.equals("qcpdi"))
				shapeMethod = "qCPDI";
			else if(shape.equals("classic") || shape.equals("dirac"))
				shapeMethod = "Dirac";
			else
				throw new Exception("The selected MPM shape function method was not recognized:\n"+args);
		
		}
	}
	
	// Archive #1 (if #2 and #3 give, passed to ArchiveTime command)
	public void doArchive(ArrayList<String> args,boolean makeUnique) throws Exception
	{
		// MPM Only
		requiresMPM(args);
		
	    // read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		String relPath = readStringArg(args.get(1));
		if(relPath.length()==0)
			throw new Exception("'"+args.get(0)+"' path has zero length:\n"+args);
		
		if(makeUnique)
			archiveRoot = "    <ArchiveRoot unique='1'>"+relPath+"</ArchiveRoot>\n";
		else
			archiveRoot = "    <ArchiveRoot>"+relPath+"</ArchiveRoot>\n";
		
		// optional element second
		if(args.size()>2)
		{	args.remove(1);
			doArchiveTime(args);
		}
	}

	// GlobalArchive #1,#2 for type and optional material ID
	public void doGlobalArchive(ArrayList<String> args) throws Exception
	{
		// MPM Only
		requiresMPM(args);
		
	    // read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		String type = readStringArg(args.get(1));
		if(type.length()==0)
			throw new Exception("'"+args.get(0)+"' quantity to archive has zero length:\n"+args);
		
		// optional material ID
		int matnum=0;
		if(args.size()>2)
		{	matnum = mats.getMatID(readStringArg(args.get(2)));
			if(matnum<=0)
			{	if(type.equals("reactionx") || type.equals("reactionx") ||
							type.equals("reactionz") || type.equals("reactionr"))
				{	matnum = readIntArg(args.get(2));
				}
				else
					throw new Exception("'"+args.get(0)+"' command has unknown material ID:\n"+args);
			}

			globalArchive = globalArchive + "    <GlobalArchive type='"+type+
							"' material='"+matnum+"'/>\n";
		}
		else
			globalArchive = globalArchive + "    <GlobalArchive type='"+type+"'/>\n";
	}

	// ArchiveTime #1,#2 (archive time and optional first archive time)
	public void doArchiveTime(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'ArchiveTime' has too few parameters:\n"+args);
		
		// archive time
		Object aTime = readNumberOrEntityArg(args.get(1),false);
		
		// optional max props
		int props = 0;
		if(args.size()>3) props = readIntArg(args.get(3));
		
		// get archiveTime
		if(props>0)
			archiveTime = "    <ArchiveTime units='ms' maxProps='"+props+"'>"+aTime+"</ArchiveTime>\n";
		else
			archiveTime = "    <ArchiveTime units='ms'>"+aTime+"</ArchiveTime>\n";
				
		// optional first archive time
		if(args.size()>2)
		{	Object firstArchiveTime = readNumberOrEntityArg(args.get(2),false);
			archiveTime = archiveTime + "    <FirstArchiveTime units='ms'>"+firstArchiveTime+"</FirstArchiveTime>\n";
		}
	}
	
	// ArchiveTime #1,#2 (archive time and optional first archive time)
	public void doGlobalArchiveTime(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'GlobalArchiveTime' has too few parameters:\n"+args);
		
		// archive time
		Object aTime = readNumberOrEntityArg(args.get(1),false);
		globalArchive = globalArchive+"    <GlobalArchiveTime units='ms'>"+aTime+"</GlobalArchiveTime>\n";
	}
	
	// TimeStep #1,#2,#3 (time step and optional max time and Courant factor)
	public void doTimeStep(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'TimeStep' has too few parameters:\n"+args);
		
		// archive time
		double aTime = readDoubleArg(args.get(1));
		timeStep = "    <TimeStep units='ms'>"+aTime+"</TimeStep>\n";
		
		// max time
		if(args.size()>2)
		{	aTime = readDoubleArg(args.get(2));
			maxTime = "    <MaxTime units='ms'>"+aTime+"</MaxTime>\n";
		}
		
		// Courant time
		if(args.size()>3)
		{	aTime = readDoubleArg(args.get(3));
			timeStep = timeStep + "    <TimeFactor>"+aTime+"</TimeFactor>\n";
		}
	}
		
	// TimeStep #1,#2,#3 (time step and optional max time and courant factor)
	public void doMaxTime(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'MaximumTime' has too few parameters:\n"+args);
		
		// archive time
		Object aTime = readNumberOrEntityArg(args.get(1),false);
		maxTime = "    <MaxTime units='ms'>"+aTime+"</MaxTime>\n";
	}
	
	// ToArchive #1,...
	public void doToArchive(ArrayList<String> args) throws Exception
	{
	    // MPM Only
	    requiresMPM(args);
	    
	    // needs at least one
	    if(args.size()<2)
	    	throw new Exception("'ToArchive' has too few parameters:\n"+args);
	    
	    // first time
    	int i;
	    if(mpmOrder == null)
	    {	// default settings
	    	mpmOrder = new StringBuffer("iY");
	    	for(i=2;i<ReadArchive.ARCH_MAXMPMITEMS;i++)
	    		mpmOrder.append('N');
	    	
	    	crackOrder = new StringBuffer("iY");
	    	for(i=2;i<ReadArchive.ARCH_MAXCRACKITEMS;i++)
	    		crackOrder.append('N');
	    }
		
		// initial history (integer starting at 0x30 to 0x3F)
		char historyChar = mpmOrder.charAt(ReadArchive.ARCH_History);
		int history;
		if(historyChar=='N')
			history=0x30;
		else if(historyChar=='Y')
			history=0x31;
		else
			history = (int)historyChar;
		int origHistory = history;
	    
		// set all options in this command
	    for(i=1;i<args.size();i++)
	    {	String archive = readStringArg(args.get(i)).toLowerCase();
	    	int loc = -1;
	    	int cloc = -1;
	        if(archive.equals("velocity"))
	        	loc = ReadArchive.ARCH_Velocity;
	        else if(archive.equals("stress"))
	        	loc = ReadArchive.ARCH_Stress;
	        else if(archive.equals("strain"))
	        	loc = ReadArchive.ARCH_Strain;
	        else if(archive.equals("plasticstrain"))
	        	loc = ReadArchive.ARCH_PlasticStrain;
	        else if(archive.equals("workenergy"))
	        	loc = ReadArchive.ARCH_WorkEnergy;
	        else if(archive.equals("temperature"))
	        	loc = ReadArchive.ARCH_DeltaTemp;
	        else if(archive.equals("plasticenergy"))
	        	loc = ReadArchive.ARCH_PlasticEnergy;
	        else if(archive.equals("shearcomponents"))
	        	loc = ReadArchive.ARCH_ShearComponents;
	        else if(archive.equals("strainenergy"))
	        	loc = ReadArchive.ARCH_StrainEnergy;
	        else if(archive.equals("jintegral"))
	        	cloc = ReadArchive.ARCH_JIntegral;
	        else if(archive.equals("stressintensity"))
	        	cloc = ReadArchive.ARCH_StressIntensity;
	        else if(archive.equals("history1"))
			{	history = history | 1;
				loc = 0;
			}
	        else if(archive.equals("history2"))
	        {	history = history | 2;
				loc = 0;
			}
	        else if(archive.equals("history3"))
	        {	history = history | 4;
				loc = 0;
			}
	        else if(archive.equals("history4"))
	        {	history = history | 8;
				loc = 0;
	        }
	        else if(archive.equals("heatenergy"))
	        	loc = ReadArchive.ARCH_HeatEnergy;
	        else if(archive.equals("concentration"))
	        	loc = ReadArchive.ARCH_Concentration;
	        else if(archive.equals("energybalance"))
	        	cloc = ReadArchive.ARCH_BalanceResults;
	        else if(archive.equals("elementcrossings"))
	        	loc = ReadArchive.ARCH_ElementCrossings;
	        else if(archive.equals("rotstrain"))
	        	loc = ReadArchive.ARCH_RotStrain;
	        
	        if(loc<0 && cloc<0)
	        	throw new Exception("'"+archive+"' is not a valid archiving option:\n"+args);
	        
	        if(loc>0)
	        	mpmOrder.replace(loc,loc+1,"Y");
	        if(cloc>0)
	        	crackOrder.replace(cloc,cloc+1,"Y");
	    }
		
		// replace the history character
		if(history != origHistory)
		{	char hchr = history==0x31 ? 'Y' : (char)history ;
			mpmOrder.setCharAt(ReadArchive.ARCH_History,hchr);
		}
	}

	// Element (element type)
	public void doElement(ArrayList<String> args) throws Exception
	{
		if(np<0)
			throw new Exception("The 'Element' command must come after the 'Analysis' command");
		
		if(args.size()<2)
			throw new Exception("'Element' command has no argument.");
		
		// options
		HashMap<String,Integer> options = new HashMap<String,Integer>(10);
		options.put("3 node triangle", new Integer(ElementBase.CST));
		options.put("4 node quadrilateral", new Integer(ElementBase.FOUR_NODE_ISO));
		options.put("8 node quadrilateral", new Integer(ElementBase.EIGHT_NODE_ISO));
		options.put("6 node triangle", new Integer(ElementBase.ISO_TRIANGLE));
		options.put("4 node interface", new Integer(ElementBase.LINEAR_INTERFACE));
		options.put("6 node interface", new Integer(ElementBase.QUAD_INTERFACE));
		options.put("8 node brick", new Integer(ElementBase.EIGHT_NODE_ISO_BRICK));
		options.put("9 node lagrange", new Integer(ElementBase.LAGRANGE_2D));
		
		int oldnameEl = lnameEl;
		lnameEl = readIntOption(args.get(1),options,"Element type");
		
		if(!ElementBase.CompatibleElements(lnameEl,oldnameEl,np))
		{	throw new Exception("Element type ("+args.get(1)+") not allowed or\nincompatible with other elements.");
		}
		
		// pass to FEA areas
		areas.setElementType(lnameEl);
	}
	
	// XMLData (section),(material ID)
	// Warning: does not check that section is valid
	// Valid are: Header, Mesh, MPMHeader, MaterialPoints, CrackList, Material (must have ID)
	//		GridBCs, ParticleBCs, Thermal, Gravity, CustomTasks, end (just append to end)
	public void doXmldata(ArrayList<String> args) throws Exception
	{	
		String section = "end";
		if(args.size()>1)
			section = readStringArg(args.get(1));
		
		// grab text
		String newXML = readVerbatim("endxmldata");
		
		// check for material section
		if(section.equals("Material"))
		{	if(args.size()<3)
				throw new Exception("XMLData command for a material needs to specify a material ID");
			String matID = readStringArg(args.get(2));
			mats.StartXMLMaterial(matID,newXML);
			return;
			
		}
		
		// check GridBCs block and intersperse
		else if(section.equals("GridBCs"))
		{	if(isFEA())
				feaBCs.AddXML(newXML);
			else
				mpmGridBCs.AddXML(newXML);
			return;
		}
		
		// check MaterialPoints block and intersperse
		else if(section.equals("MaterialPoints"))
		{	if(isFEA())
				regions.AddXML(newXML);
			else
				regions.AddXML(newXML);
			return;
		}
		
		// check previous option
		String currentXML = xmldata.get(section);
		if(currentXML != null) newXML = currentXML+newXML;
		
		// set value
		xmldata.put(section,newXML);
	}
	
	// Add an entity
	public void doEntity(ArrayList<String> args) throws Exception
	{	// read entity and value
		if(args.size()<3)
			throw new Exception("'Entity' command has too few arguments:\n"+args);
		String ent = readStringArg(args.get(1));
		String val = readStringArg(args.get(2));
		entities.put(ent, val);
	}
	
	// Output command for FEA analysis
	public void doOutput(ArrayList<String> args) throws Exception
	{
	    // FEA Only
	    requiresFEA(args);
		
		// first time set the defaults
		if(outFlags == null)
			outFlags = new StringBuffer("YYYYNY");
	    
	    // read quantity and option material
		if(args.size()<3)
			throw new Exception("'Output' command has too few arguments:\n"+args);
		
		String quant = readStringArg(args.get(1)).toLowerCase();
		String option = readStringArg(args.get(2)).toLowerCase();
			
		// add to flags
		// enum { DISPLACEMENT_OUT=0,FORCE_OUT,ELEMSTRESS_OUT,AVGSTRESS_OUT,
	    //                               REACT_OUT,ENERGY_OUT,NUMBER_OUT };
		if(option.equals("selected"))
			option="C";
		else if(option.equals("no"))
			option="N";
		else if(option.equals("yes"))
			option="Y";
		else
			throw new Exception("'Output' option must be 'yes', 'no', or 'selected':\n"+args);
		
		int offset;
		if(quant.equals("displacements"))
			offset=0;
		else if(quant.equals("forces"))
			offset=1;
		else if(quant.equals("elementstresses"))
			offset=2;
		else if(quant.equals("nodalstresses"))
			offset=3;
		else if(quant.equals("reactivities"))
			offset=4;
		else if(quant.equals("energy"))
			offset=5;
		else
	    	throw new Exception("Unrecognized 'Output' option:\n"+args);
		
		// make the change
		outFlags.deleteCharAt(offset);
		outFlags.insert(offset,option);
		
		// is there an archive time?
		if(args.size()>3)
		{	args.remove(1);
			args.remove(1);
			doOutput(args);
		}
	}

	// Temperature (FEA) only
	public void doTemperature(ArrayList<String> args) throws Exception
	{	// Temperature #1 which is a function
		if(args.size()<2)
			throw new Exception("'Temperature' command with too few arguments:\n"+args);
			
		feaTemp = readStringArg(args.get(1));
	}

	// Stress Free Temperature
	public void doStressFreeTemp(ArrayList<String> args) throws Exception
	{	if(args.size()<2)
			throw new Exception("'StressFreeTemp' command with too few arguments:\n"+args);
			
		stressFreeTemp = readDoubleArg(args.get(1));
	}

	// ThermalRamp (diff),<(time)>
	public void doThermalRamp(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
		
		if(args.size()<2)
			throw new Exception("'ThermalRamp' command with too few arguments:\n"+args);
			
		rampDiff = readDoubleArg(args.get(1));
		
		rampTime = -1.;
		if(args.size()>2)
		{	rampTime = readDoubleArg(args.get(2));
			if(rampTime<=0.) rampTime = -1.;
		}
	}
	
	// RampStart (start time)
	public void doRampStart(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
		
		if(args.size()<2)
			throw new Exception("'RampStart' command with too few arguments:\n"+args);
			
		rampStart = readDoubleArg(args.get(1));
	}
	
	// Damping #1 (number or function),#2 (0 to 1 for PIC)
	// also does PDamping command
	public void doDamping(ArrayList<String> args,String dcmd) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+dcmd+"' has too few parameters:\n"+args);
		
		// damping factor (required)
		double damp=0.;
		String dampcmd;
		Object dampArg = readStringOrDoubleArg(args.get(1));
		if(dampArg.getClass().equals(Double.class))
		{	damp = ((Double)dampArg).doubleValue();
			dampcmd = "    <"+dcmd;
		}
		else
		{	dampcmd = "    <"+dcmd+ " function='"+dampArg+"'";
		}
		
		// PIC fraction (optional)
		if(args.size()>2)
		{	double pic = readDoubleArg(args.get(2));
			if(pic<0 || pic>1)
				throw new Exception("PIC damping in '"+dcmd+"' must be from 0 to 1:\n"+args);
			dampcmd = dampcmd+" PIC='"+pic+"'>"+damp+"</"+dcmd+">\n";
		}
		else
			dampcmd = dampcmd+">"+damp+"</"+dcmd+">\n";
		
		if(dcmd.equals("Damping"))
			damping = dampcmd;
		else
			pdamping = dampcmd;
	}
	
	// MultimaterialMode Vmin,Dcheck,Normals,RigidBias
	public void doMultimaterialMode(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
		
		// turn it on
		MMNormals = 2;		// avggrad default
		
		// Vmin
		if(args.size()>1)
		{	MMVmin = readDoubleArg(args.get(1));
			if(MMVmin<0.) MMVmin = 0.;
		}
		
		// Dcheck
		if(args.size()>2)
		{	HashMap<String,Integer> options = new HashMap<String,Integer>(4);
			options.put("enabled", new Integer(1));
			options.put("yes", new Integer(1));
			options.put("disabled", new Integer(0));
			options.put("no", new Integer(0));
			MMDcheck = readIntOption(args.get(2),options,"Displacement check option");
		}
		
		// Normals
		if(args.size()>3)
		{	HashMap<String,Integer> options = new HashMap<String,Integer>(4);
			options.put("maxgrad", new Integer(0));
			options.put("maxvol", new Integer(1));
			options.put("avggrad", new Integer(2));
			options.put("owngrad", new Integer(3));
			options.put("specify", new Integer(4));
			MMNormals = readIntOption(args.get(3),options,"Normals option");
		}
		
		if(MMNormals==4)
		{	// polar angles
			if(args.size()>4)
				MMAzimuth = readDoubleArg(args.get(4));
			if(args.size()>5)
				MMPolar = readDoubleArg(args.get(5));
		}
		else if(args.size()>4)
		{	//Rigid Bias
			MMRigidBias = readDoubleArg(args.get(4));
			if(MMRigidBias<0.) MMRigidBias = 0.;
		}
	}
	
	// ContactPosition Value
	public void doContactPosition(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
		
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		double cp = readDoubleArg(args.get(1));
		ContactPosition = "      <ContactPosition>"+cp+"</ContactPosition>\n";
	}
	
	// Friction (number or stick, single (ignore), none),<material ID (only as material prop)>
	// MMMode = 0 (cracks), 1 (multimaterial), 2 (material property), 3 an attribute
	public String doFriction(ArrayList<String> args,int MMMode) throws Exception
	{	// MPM Only
		requiresMPM(args);
		
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		// see if nonnegative number
		double frict = 0.;
		try
		{	frict = readDoubleArg(args.get(1));
			if(frict<0)
				throw new Exception("The friction coefficient must be positive:\n"+args);
		}
		catch(Exception e)
		{	HashMap<String,Integer> options = new HashMap<String,Integer>(4);
			options.put("stick", new Integer(0));
			options.put("single", new Integer(1));
			options.put("ignore", new Integer(1));
			options.put("none", new Integer(2));
			int foption = readIntOption(args.get(1),options,"Friction setting");
			if(foption==0)
				frict = -5.;			// number between -1 and -9
			else if(foption==1)
				frict = -11.;			// number <-10
			else
				frict = 0.0;			// frictionless
		}
		
		// material property needs material ID
		if(MMMode==2)
		{	if(args.size()<3)
				throw new Exception("'"+args.get(0)+"' as material property has too few parameters:\n"+args);
		
			int matnum = mats.getMatID(readStringArg(args.get(2)));
			if(matnum<=0)
				throw new Exception("'"+args.get(0)+"' as material property has unknown material ID:\n"+args);
			
			String cmd = "    <Friction mat='"+matnum+"'>"+frict+"</Friction>\n";
			return cmd;
		}
		
		// Friction for cracks or multimaterial mode
		String cmd = "      <Friction>"+frict+"</Friction>\n";
		if(MMMode==1)
			FrictionMM = cmd;
		else if(MMMode==0)
			cracks.setFriction(cmd);
		else
			return " frict='"+frict+"'";
		return null;
	}
	
	// ImperfectInterface Dt,Dn,<Dnc>
	// MMMode = 0 (cracks), 1 (multimaterial), 2 (material property), 3 (CrackInterface commmand)
	public String doImperfectInterface(ArrayList<String> args,int MMMode) throws Exception
	{	// MPM Only
		requiresMPM(args);
		
		// read analysis type
		if(args.size()<3 || (MMMode==2 && args.size()<5))
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		// read doubles
		double Dt = readDoubleArg(args.get(1));
		double Dnt = readDoubleArg(args.get(2));
		double Dnc=0.;
		String cmd;
		if(args.size()>3)
			Dnc = readDoubleArg(args.get(3));
		
		if(MMMode==2)
		{	// get material ID
			int matnum = mats.getMatID(readStringArg(args.get(4)));
			if(matnum<=0)
				throw new Exception("'"+args.get(0)+"' as material property has unknown material ID:\n"+args);
			
			cmd = "    <Friction Dt='"+Dt+"' Dnt='"+Dnt+"' Dnc='"+Dnc+
					"' mat='"+matnum+"'>11</Friction>\n";
		}
		else if(MMMode==3)
		{	if(args.size()>3)
				cracks.setCrackFriction(args," Dt='"+Dt+"' Dnt='"+Dnt+"' Dnc='"+Dnc+"'");
			else
				cracks.setCrackFriction(args," Dt='"+Dt+"' Dnt='"+Dnt+"'");
			return null;
		}
		else
		{	if(args.size()>3)
				cmd = "      <Friction Dt='"+Dt+"' Dnt='"+Dnt+"' Dnc='"+Dnc+"'>11</Friction>\n";
			else
				cmd = "      <Friction Dt='"+Dt+"' Dn='"+Dnt+"'>11</Friction>\n";
		
			if(MMMode==1)
				FrictionMM = cmd;
			else
				cracks.setFriction(cmd);
		}
		
		return cmd;
	}
		
	// FeedbackDamping #1,#2,#3 (number,function,number)
	// PFeedbackDamping #1,#2,#3 (same for particle damping)
	public void doFBDamping(ArrayList<String> args,String dfbcmd) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+dfbcmd+"' has too few parameters:\n"+args);
		
		// read gain (required)
		double damp = readDoubleArg(args.get(1));
		
		// optional target and max alpha
		String target = null;
		double maxdamp = -1.;
		if(args.size()>2)
			target = readStringArg(args.get(2));
		if(args.size()>3)
			maxdamp = readDoubleArg(args.get(3));
		
		String fb;
		if(target==null)
			fb = "    <"+dfbcmd+">"+damp+"</"+dfbcmd+">\n";
		else if(maxdamp<0.)
			fb = "    <"+dfbcmd+" target='"+target+"'>"+damp+"</"+dfbcmd+">\n";
		else
		{	fb = "    <"+dfbcmd+" target='"+target+"' max='"+maxdamp+
								"'>"+damp+"</"+dfbcmd+">\n";
		}
		
		if(dfbcmd.equals("FeedbackDamping"))
			fbDamping = fb;
		else
			pfbDamping = fb;
	}
	
	// LeaveLimit #1 (integer)
	public void doLeaveLimit(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		// leav limit (required)
		int leave = readIntArg(args.get(1));
		leaveLimit = "    <LeaveLimit>"+leave+"</LeaveLimit>\n";
	}
	
	// Diffusion #1,<#2>
	public void doDiffusion(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		// diffusion on or off (required)
		String option = readStringArg(args.get(1));
		if(option.toLowerCase().equals("no"))
		{	diffusion = null;
			return;
		}
		else if(!option.toLowerCase().equals("yes"))
			throw new Exception("'"+args.get(0)+"' first parameter must be yes or no:\n"+args);
		
		double ref=0.;
		if(args.size()>2)
		{	ref = readDoubleArg(args.get(2));
			if(ref<0. || ref>1.)
				throw new Exception("'"+args.get(0)+"' second parameter must be 0 to 1:\n"+args);
		}

		// the command
		diffusion = "    <Diffusion reference='"+ref+"/>\n";
	}
	
	// Conduction (yes or no),<adibatic (or mechanical energy) or isothermal or "Crack Tips">
	public void doConduction(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		// yes or no
		boolean hasConduction;
		String option = readStringArg(args.get(1));
		if(option.toLowerCase().equals("no"))
			hasConduction = false;
		else if(option.toLowerCase().equals("yes"))
			hasConduction = true;
		else
			throw new Exception("'"+args.get(0)+"' first parameter must be yes or no:\n"+args);
		
		// options
		boolean hasCoupling = false;
		boolean hasTips = false;
		boolean hasFriction = false;
		boolean hasCrackFriction = false;
		HashMap<String,Integer> options = new HashMap<String,Integer>(4);
		options.put("adiabatic", new Integer(1));
		options.put("mechanical energy", new Integer(1));
		options.put("isothermal", new Integer(2));
		options.put("crack tips", new Integer(3));
		options.put("friction", new Integer(4));
		options.put("crack friction", new Integer(5));
		
		// each one
		int arg = 2;
		while(args.size()>arg)
		{	int opt = readIntOption(args.get(arg),options,"Conduction option");
			if(opt==1)
				hasCoupling = true;
			else if(opt==2)
				hasCoupling = false;
			else if(opt==3)
				hasTips = true;
			else if(opt==4)
				hasFriction = true;
			else if(opt==5)
				hasCrackFriction = true;
			arg++;
		}
		
		// <Conduction/>, <CrackTipHeating/>, <EnergyCoupling/>
		if(hasConduction)
			conduction = "    <Conduction/>\n";
		else
			conduction = "";
		if(hasCoupling) conduction = conduction + "    <EnergyCoupling/>\n";
		if(hasTips) conduction = conduction + "    <CrackTipHeating/>\n";
		if(hasFriction) conduction = conduction + "    <ContactHeating/>\n";
		if(hasCrackFriction) conduction = conduction + "    <CrackContactHeating/>\n";
	}
	
	// CustomTask name
	public void doCustomTask(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
		
		// read task name
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' is missing a custom task name:\n"+args);
		currentCustomTask = readStringArg(args.get(1));
		
		// finish last task
		if(customTasks.length()>0)
		{	customTasks.append("    </Schedule>\n");
		}
		
		// start new custom task
		customTasks.append("    <Schedule name='"+currentCustomTask+"'>\n");
	}
	
	// Parameter #1,<#2>
	public void doParameter(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
		
		// must be in custom task
		if(customTasks.length()==0)
			throw new Exception("'"+args.get(0)+"' must can after CustomTask command:\n"+args);
		
		// read parameter name
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' is missing the parameter name:\n"+args);
		String paramName = readStringArg(args.get(1));
		
		// handle special commands in some tasks
		if(currentCustomTask.equals("ReverseLoad"))
		{	// quantity combines into a single command
			if(paramName.equals("quantity") && args.size()>2)
			{	String value = readStringArg(args.get(2));
				customTasks.append("      <Parameter name='global "+value+"'/>\n");
				return;
			}
			else if(paramName.equals("material") && args.size()>2)
			{	// material looks for material ID
				int matnum = mats.getMatID(readStringArg(args.get(2)));
				if(matnum<=0)
				{	// negative is allowed for reaction forces
					matnum = readIntArg(args.get(2));
					if(matnum>=0)
						throw new Exception("'"+args.get(0)+"' command has unknown material ID or invalid BC ID:\n"+args);
				}
				customTasks.append("      <Parameter name='mat'>"+matnum+"</Parameter>\n");
				return;
			}
		}
		
		// volume gradient requires a material ID
		else if(currentCustomTask.equals("VTKArchive"))
		{	if(paramName.equals("volumegradient"))
			{	if(args.size()<3)
					throw new Exception("'volumegradient' quantity requires a material ID:\n"+args);
				int matnum = mats.getMatID(readStringArg(args.get(2)));
				if(matnum<=0)
					throw new Exception("'volumegradient' quantity requires a valid material ID:\n"+args);
				customTasks.append("      <Parameter name='volumegradient'>"+matnum+"</Parameter>\n");
				return;
			}
		}
		
		// single parameter
		if(args.size()==2)
		{	customTasks.append("      <Parameter name='"+paramName+"'/>\n");
		}
		else
		{	String value = readStringArg(args.get(2));
			customTasks.append("      <Parameter name='"+paramName+"'>"+value+"</Parameter>\n");
		}
	}
	
	// Gravity <#1>,<#2>,<#3>
	public void doGravity(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// defaults
		double gx=0.,gy=-9806.65,gz=0.;
		
		if(args.size()>1) gx = 1000.*readDoubleArg(args.get(1));
		if(args.size()>2) gy = 1000.*readDoubleArg(args.get(2));
		if(args.size()>3) gz = 1000.*readDoubleArg(args.get(3));

		// the command
		gravity = "    <BodyXForce>"+gx+"</BodyXForce>\n"+
					"    <BodyYForce>"+gy+"</BodyYForce>\n"+
					"    <BodyZForce>"+gz+"</BodyZForce>\n";
	}
	// PtsPerElement #1 (integer)
	public void doPtsPerElement(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters:\n"+args);
		
		//  points per element
		int pts = readIntArg(args.get(1));
		if(pts<0 || pts>5 || (pts>3 && isMPM3D()))
			throw new Exception("'"+args.get(0)+"' has unsupported number of points per element:\n"+args);
		
		int numCell = pts*pts;
		if(isMPM3D()) numCell *= pts;
		ptsPerElement = "    <MatlPtsPerElement>"+numCell+"</MatlPtsPerElement>\n";
	}
	
	// convert @ expression to Double
	public Double getAtDouble(String s)
	{	// get as string and see if a number
		String expr = getAtString(s);
		if(expr==null) return null;
		try
		{	return new Double(Double.parseDouble(expr));
		}
		catch(Exception e)
		{
		}
		return null;
	}
	
	// convert @ expression to String
	public String getAtString(String s)
	{	// split at periods
		String[] atoms = s.substring(1).split("[.]");
		
		// process them
		try
		{	int i=0;
			while(i<atoms.length)
			{	String nextAtom = atoms[i];
		
				// read key points
				if(nextAtom.equals("key"))
				{	// make sure has data
					i+=2;
					return areas.getKeypointProperty(readStringArg(atoms[i-1]),readStringArg(atoms[i]));
				}
				else if(nextAtom.equals("path"))
				{	// make sure has data
					i+=2;
					return areas.getPathProperty(readStringArg(atoms[i-1]),readStringArg(atoms[i]));
				}
				
				else if(runningScript)
				{	// look for obj.property
					Object obj = objs.get(nextAtom);
					if(obj==null || i+1>=atoms.length) return null;
					i++;
					nextAtom = atoms[i];
					
					// object properties (string properties end in '$')
					if(nextAtom.equals("energy"))
					{	if(!obj.getClass().equals(ResultsDocument.class)) return null;
						return ((ResultsDocument)obj).getEnergy();
					}
					
					else if(nextAtom.equals("get"))
					{	i++;
						if(i>=atoms.length) return null;
						if(!obj.getClass().equals(CmdViewer.class)) return null;
						if(!obj.getClass().equals(CmdViewer.class)) return null;
						return ((CmdViewer)obj).getVariable(atoms[i]);
					}
					
					else if(nextAtom.equals("section"))
					{	i++;
						if(i>=atoms.length) return null;
						int alen = atoms[i].length();
						if(alen<2) return null;
						if(atoms[i].charAt(0)=='"' && atoms[i].charAt(alen-1)=='"')
							atoms[i] = atoms[i].substring(1,alen-1);
						else
						{	// a string variable is allowed
							String strAtom = variablesStrs.get(atoms[i]);
							if(strAtom!=null) atoms[i] = strAtom;
						}
						if(!obj.getClass().equals(ResultsDocument.class)) return null;
						return ((ResultsDocument)obj).section(atoms[i]);
					}
					
				}
			}
		}
		catch(Exception e)
		{
		}
		
		// if here, than bad expression
		return null;
	}

	// when analysis is done create XML commands
	public String buildXMLCommands()
	{	// start buffer for XML commands
		String more;
		StringBuffer xml = new StringBuffer("<?xml version='1.0'?>\n");
		xml.append("<!DOCTYPE JANFEAInput SYSTEM 'pathto.dtd'");
		if(entities.size()>0)
		{	xml.append("\n[\n");
			Set<String> keys = entities.keySet();
			String [] allkeys = new String [entities.size()];
			allkeys = keys.toArray(allkeys);
			int i;
			for(i=0;i<entities.size();i++)
			{	String value = entities.get(allkeys[i]);
				xml.append("  <!ENTITY "+allkeys[i]+" \""+value+"\">\n");
			}
			xml.append("]>\n");
		}
		else
			xml.append(">\n");
		xml.append("<JANFEAInput version='3'>\n\n");
		
		// Header
		//-----------------------------------------------------------
		xml.append("  <Header>\n    <Description>\n");
		xml.append("Title: "+title+"\n");
		if(username != null)
			xml.append("User Name: "+username+"\n");
		if(header != null)
			xml.append(header);
		xml.append("    </Description>\n");
		xml.append("    <Analysis>"+np+"</Analysis>\n");
		if(outFlags!=null) xml.append("    <Output>"+outFlags+"</Output>\n");
		more = xmldata.get("Header");
		if(more != null) xml.append(more);
		xml.append("  </Header>\n\n");
		
		// FEA section: Mesh
		//-----------------------------------------------------------
		if(isFEA())
		{	xml.append("  <Mesh>\n"+areas.toXMLString());
		
			// check added xml
			more = xmldata.get("Mesh");
			if(more != null) xml.append(more);
			
			// done
			xml.append("  </Mesh>\n\n");
			
			// BMPRegion, Body, and Hole blocks
			xml.append(regions.toXMLString());
		}
		
		// MPM sections: MPMHeader, Mesh, MaterialPoints, CrackList
		//-----------------------------------------------------------
		if(isMPM())
		{	// MPM Header
			//-----------------------------------------------------------
			xml.append("  <MPMHeader>\n");
			
			// MPM method and GIMP
			xml.append("    <MPMMethod>"+mpmMethod+"</MPMMethod>\n");
			xml.append("    <GIMP type='"+shapeMethod+"'/>\n");
			if(ptsPerElement!=null) xml.append(ptsPerElement);
			xml.append(timeStep);
			xml.append(maxTime);
			xml.append(archiveRoot);
			xml.append(archiveTime);
			if(mpmOrder == null) mpmOrder = new StringBuffer("iYYYYYNYYYNNNNNNNN");
			xml.append("    <MPMArchiveOrder>"+mpmOrder+"</MPMArchiveOrder>\n");
			if(crackOrder == null) crackOrder = new StringBuffer("iYNNN");
			xml.append("    <CrackArchiveOrder>"+crackOrder+"</CrackArchiveOrder>\n");
			
			// global archive
			if(globalArchive.length()>0)
				xml.append(globalArchive);
			
			// damping, leave limit, diffusion
			if(damping!=null) xml.append(damping);
			if(pdamping!=null) xml.append(pdamping);
			if(fbDamping!=null) xml.append(fbDamping);
			if(pfbDamping!=null) xml.append(pfbDamping);
			if(leaveLimit!=null) xml.append(leaveLimit);
			if(diffusion!=null) xml.append(diffusion);
			
			// cracks
			more = cracks.getSettings(MMNormals,ContactPosition);
			if(more != null) xml.append(more);
			
			// Multimaterial mode <MultiMaterialMode Vmin='0.0' Dcheck='0' Normals='0' RigidBias='100'>
			// Subordinate friction and contact position
			if(MMNormals>=0)
			{	if(MMNormals==4)
				{	xml.append("    <MultiMaterialMode Vmin='"+MMVmin+"' Dcheck='"+MMDcheck+
						"' Normals='"+MMNormals+"' Azimuth='"+MMAzimuth+"' Polar='"+MMPolar+"'>\n");
				}
				else
				{	xml.append("    <MultiMaterialMode Vmin='"+MMVmin+"' Dcheck='"+MMDcheck+
								"' Normals='"+MMNormals+"' RigidBias='"+MMRigidBias+"'>\n");
				}
				if(FrictionMM!=null) xml.append(FrictionMM);
				if(ContactPosition!=null) xml.append(ContactPosition);
				xml.append("    </MultiMaterialMode>\n");
			}
			
			// stress free temperature
			if(stressFreeTemp!=0.)
				xml.append("    <StressFreeTemp>"+stressFreeTemp+"</StressFreeTemp>\n");
			
			// check added xml
			more = xmldata.get("MPMHeader");
			if(more != null) xml.append(more);
			
			// done
			xml.append("  </MPMHeader>\n\n");
			
			// MPM Mesh
			//-----------------------------------------------------------
			if(mpmMeshToFile)
				xml.append("  <Mesh output='file'>\n");
			else
				xml.append("  <Mesh>\n");
			
			xml.append(gridinfo.toXMLString());
			
			// check added xml
			more = xmldata.get("Mesh");
			if(more != null) xml.append(more);
			
			// done
			xml.append("  </Mesh>\n\n");
			
			// MPM Material Points (XMLData was already added, if any)
			//-----------------------------------------------------------
			xml.append("  <MaterialPoints>\n"+regions.toXMLString());
			xml.append("  </MaterialPoints>\n\n");
			
			// MPM Cracks
			more = cracks.getCrackList();
			if(more!=null) xml.append(more);
		}
		
		// Materials (XMLData was already added, if any)
		//-----------------------------------------------------------
		xml.append(mats.toXMLString());
		
		// GridBCs
		//-----------------------------------------------------------
		String gridXml = null;
		if(isFEA())
			gridXml = feaBCs.toXMLString();
		else
			gridXml = mpmGridBCs.toXMLString();
		
		if(gridXml.length()>0)
			xml.append("  <GridBCs>\n"+gridXml+"  </GridBCs>\n\n");
		
		// ParticleBCs
		//-----------------------------------------------------------
		if(isMPM())
		{	String partXml = mpmParticleBCs.toXMLString();
			more = xmldata.get("ParticleBCs");
			if(partXml.length()>0 || more!=null)
			{	xml.append("  <ParticleBCs>\n"+partXml);
				if(more != null) xml.append(more);
				xml.append("  </ParticleBCs>\n\n");
			}
		}
	
		// FEA: Thermal, MPM: Thermal, Gravity, CustomTasks
		//-----------------------------------------------------------
		if(isFEA())
		{	// FEA: Thermal
			//-----------------------------------------------------------
			more = xmldata.get("Thermal");
			if(more!=null || feaTemp!=null || stressFreeTemp!=0.)
			{	xml.append("  <Thermal>\n");
				
				if(feaTemp!=null)
					xml.append("    <Temperature>"+feaTemp+"</Temperature>\n");
				
				if(stressFreeTemp!=0.)
					xml.append("    <StressFreeTemp>"+stressFreeTemp+"</StressFreeTemp>\n");
				
				// check added xml
				if(more != null) xml.append(more);

				// done
				xml.append("  </Thermal>\n\n");
			}
		}
		
		else
		{	// MPM: Thermal
			//-----------------------------------------------------------
			more = xmldata.get("Thermal");
			if(more!=null || conduction!=null || rampTime>-1.5)
			{	xml.append("  <Thermal>\n");

				// conduction
				if(conduction!=null) xml.append(conduction);
				
				// <Isothermal time="(time)" start="(start time)">(diff)</Isothermal>
				if(rampTime>-1.5)
				{	xml.append("    <Isothermal");
					if(rampTime>0.) xml.append(" time='"+rampTime+"'");
					if(rampStart>0.) xml.append(" start='"+rampStart+"'");
					xml.append(">"+rampDiff+"</Isothermal>\n");
				}
					
			
				// check added xml
				if(more != null) xml.append(more);

				// done
				xml.append("  </Thermal>\n\n");
			}
			
			// MPM: Gravity
			//-----------------------------------------------------------
			more = xmldata.get("Gravity");
			if(more!=null || gravity!=null)
			{	xml.append("  <Gravity>\n");

				// check added xml
				if(gravity != null) xml.append(gravity);
				if(more != null) xml.append(more);

				// done
				xml.append("  </Gravity>\n\n");
			}
			
			// MPM: CustomTasks
			//-----------------------------------------------------------
			more = xmldata.get("CustomTasks");
			if(customTasks.length()>0 || more!=null)
			{	xml.append("  <CustomTasks>\n");
			
				// add tasks
				if(customTasks.length()>0)
				{	xml.append(customTasks);
					xml.append("    </Schedule>\n");
				}

				// check added xml
				if(more != null) xml.append(more);

				// done
				xml.append("  </CustomTasks>\n\n");
			}
		}
		
		// check added xml
		more = xmldata.get("end");
		if(more != null) xml.append(more);
		
		// convert to string and return
		xml.append("</JANFEAInput>\n");
		return xml.toString();
	}
		
	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	// empty if no text or if has welcome message
	public boolean isEmptyDocument()
	{	String cmds=cmdField.getCommands();
		if(cmds.length()==0) return true;
		if(cmds.equals("Welcome to "+NairnFEAMPMViz.appNameReadable+", "+NairnFEAMPMViz.versionReadable)) return true;
		return false;
	}
	
	// call by results when it closes
	public void setLinkedResults(DocViewer someResults) { linkedResults=someResults; }
	
	// tell linked results your are closing
	public void windowClosed(WindowEvent e)
	{	if(linkedResults!=null) linkedResults.setCommandsWindow(null);
		super.windowClosed(e);
	}

	// called when analysis is done and should link to new results in console
	public DocViewer linkToResults()
	{	if(linkedResults!=null)
		{	linkedResults.windowClosing(null);
		}
		NairnFEAMPMViz.main.openDocument(soutConsole.getFile());
		linkedResults=(DocViewer)NairnFEAMPMViz.main.findDocument(soutConsole.getFile());
		linkedResults.setCommandsWindow(this);
		return linkedResults;
	}
	
	// type of analysis
	public boolean isFEA() { return np>=0 && np<BEGIN_MPM_TYPES ; }
	public boolean isMPM() { return np>BEGIN_MPM_TYPES ; }
	public boolean isMPM3D() { return np==THREED_MPM; }
	
	// verify FEA or MPM
	public void requiresFEA(ArrayList<String> args) throws Exception
	{	if(isFEA()) return;
		if(args != null)
		{	if(args.size()>1)
				throw new Exception("The command '"+args.get(0)+"' is only allowed in FEA calculations:\n"+args);
		}
		throw new Exception("Some unknown command is only allowed in FEA calculations.");
	}
	public void requiresMPM(ArrayList<String> args) throws Exception
	{	if(isMPM()) return;
		if(args != null)
		{	if(args.size()>1)
				throw new Exception("The command '"+args.get(0)+"' is only allowed in MPM calculations:\n"+args);
		}
		throw new Exception("Some unknown command is only allowed in MPM calculations.");
	}
	
	// return Double object or look for entity
	public Object readNumberOrEntityArg(String text,boolean isInt) throws Exception
	{	Object arg = readStringOrDoubleArg(text);
		if(arg.getClass().equals(Double.class))
		{	if(isInt)
			{	Integer intarg = new Integer(((Double)arg).intValue());
				return intarg;
			}
			return arg;
		}
		
		// Strip & and ; if there
		String ent = (String)arg;
		if(ent.startsWith("&") && ent.endsWith(";"))
			ent = ent.substring(1, ent.length()-1);
		System.out.println(ent);
		
		// look for valid entity
		if(entities.get(ent)==null)
			throw new Exception("The argument '"+text+"'\nis neither a number nor a valid entity");
		return "&"+ent+";";
	}
	
	// override to check commands or anlaysis running
	public boolean isRunning()
	{	if(super.isRunning()) return true;
		if(nfmAnalysis==null) return false;
		if(nfmAnalysis.isRunning()) return true;
		return false;
	}
	
	// scripts return an output file to use for output
	// otherwise return null to select default output
	public File getOutputFile()
	{	if(!runningScript) return null;
		return scriptOutput;
	}
	
	// return variable value (or null if none) as string
	public String getVariable(String varName)
	{	Double nvar = variables.get(varName);
		if(nvar!=null)
			return JNUtilities.formatDouble(nvar.doubleValue());
		String svar = variablesStrs.get(varName);
		return svar;
	}
	
	// export the output file
	// return true is done or false if error or if cancelled
	public boolean exportOutput(File exportFile,String etitle)
	{
		if(etitle==null) etitle = "Export contents of output panel";
	
		if(exportFile==null)
		{	JFileChooser expChooser=new JFileChooser();
			expChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			NFMVPrefs.setWorkspace(expChooser);
			expChooser.setDialogTitle(etitle);
			int result = expChooser.showSaveDialog(this);
			if(result != JFileChooser.APPROVE_OPTION) return false;
			exportFile = expChooser.getSelectedFile();
		}
		
		// save output text to exportFile
		try
		{	FileWriter theFile=new FileWriter(exportFile);
			theFile.write(soutConsole.console.getText());
			theFile.flush();
			theFile.close();
		}
		catch (Exception fe)
		{	Toolkit.getDefaultToolkit().beep();
			JOptionPane.showMessageDialog(this,"Error exporting output results: " + fe);
			return false;
		}
				
		return true;
	}
	
	//----------------------------------------------------------------------------
	// Actions as inner classes
	//----------------------------------------------------------------------------
	
	// action for stop analysis menu command
	protected class BgAnalysisAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public BgAnalysisAction()
		{	super("Background FEA/MPM Analysis...",KeyEvent.VK_B);
		}
 
		public void actionPerformed(ActionEvent e) { runNFMAnalysis(true,NFMAnalysis.FULL_ANALYSIS,null); }
	}

	// action for stop analysis menu command
	protected class CheckAnalysisAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public CheckAnalysisAction()
		{	super("Test FEA/MPM Mesh...",KeyEvent.VK_T);
		}
 
		public void actionPerformed(ActionEvent e) { runNFMAnalysis(false,NFMAnalysis.RUN_CHECK_MESH,null); }
	}
	
	// action for stop analysis menu command
	protected class InterpretCommandsAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public InterpretCommandsAction()
		{	super("Interpret Commands...",KeyEvent.VK_I);
		}
 
		public void actionPerformed(ActionEvent e) { runNFMAnalysis(false,NFMAnalysis.INTERPRET_ONLY,null); }
	}
	
	// action to shaw partner menu command
	protected class ShowPartnerAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public ShowPartnerAction()
		{	super("Show Results");
		}
 
		public void actionPerformed(ActionEvent e)
		{	if(linkedResults!=null)
				linkedResults.toFront();
			else
				JNApplication.appBeep();
		}
	}
	
	// action for stop analysis menu command
	protected class ExportXMLAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public ExportXMLAction()
		{	super("Export Output...",KeyEvent.VK_S,true);
		}
 
		public void actionPerformed(ActionEvent e) { exportOutput(null,null); }
	}
	
	// action for stop analysis menu command
	protected class StopCurrentModeAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public StopCurrentModeAction()
		{	super("Stop Analysis",KeyEvent.VK_PERIOD);
		}
 
		public void actionPerformed(ActionEvent e)
		{	System.out.println(running+","+nfmAnalysis.isRunning());
			if(running)
				running = false;
			else if(nfmAnalysis!=null)
			{	if(nfmAnalysis.isRunning())
					nfmAnalysis.stopRunning();
			}
		}
	}


}
