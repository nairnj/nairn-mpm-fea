
/*
 * CmdViewer
 * NairnFEAMPMViz Application
 * 
 * Created 
 */

import java.awt.Font;
import java.awt.Toolkit;
import java.awt.event.*;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;
import java.util.Scanner;
import java.util.Set;

import javax.swing.*;

import geditcom.JNFramework.*;

public class CmdViewer extends JNCmdTextDocument
{
	private static final long serialVersionUID = 1L;

	private ConsolePane soutConsole;
	private DocViewer linkedResults = null;

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
	private String consistentUnits;
	private boolean plusSpin;
	private int processors = 1;
	private int lnameEl;
	private HashMap<String, String> xmldata = null;
	private HashMap<String, String> entities = null;
	public Materials mats = null;
	public Areas areas = null;
	public Regions regions = null;
	public MPMGrid gridinfo = null;
	private FEABCs feaBCs = null;
	public MPMGridBCs mpmGridBCs = null;
	public MPMParticleBCs mpmParticleBCs = null;
	public Cracks cracks = null;
	private StringBuffer outFlags;
	private int mpmMethod;
	private boolean skipExtrap;
	private String shapeMethod;
	private String archiveRoot;
	private String archiveTime;
	private String timeStep;
	private String maxTime;
	private String cflFactor;
	private String CPDIrcrit;
	private String exactTractions;
	private String transCflFactor;
	private String globalArchive;
	private String damping;
	private String pdamping;
	private String xpic;
	private String fbDamping;
	private String pfbDamping;
	private String leaveLimit;
	private String deleteLimit;
	private String extrapolateRigid;
	private String ptsPerElement;
	private String diffusion;
	private StringBuffer otherDiffusion;
	private String conduction;
	private String gravity;
	private StringBuffer mpmOrder;
	private StringBuffer crackOrder;
	private boolean mpmMeshToFile;
	private String feaTemp;
	private double stressFreeTemp;
	private boolean stopCommand;
	private int MMLump;
	private int MMNormals; // <0 means no multimaterial mode
	private double MMRigidBias;
	private double MMAzimuth;
	private double MMPolar;
	private String ContactPosition;
	private String ContactPositionCracks;
	private String FrictionMM;
	private StringBuffer customTasks;
	private String currentCustomTask = null;

	// scripting
	private boolean runningScript;
	private HashMap<String, Object> objs = null;
	private ArrayList<String> scriptParams;
	private CmdViewer theScript;

	// constants
	public static final int PLANE_STRAIN = 0;
	public static final int PLANE_STRESS = 1;
	public static final int AXI_SYM = 2;
	public static final int THREE_DIM = 3;
	public static final int BEGIN_MPM_TYPES = 9;
	public static final int PLANE_STRAIN_MPM = 10;
	public static final int PLANE_STRESS_MPM = 11;
	public static final int THREED_MPM = 12;
	public static final int AXI_SYM_MPM = 13;

	public static final int NO_ELEMENT = 0;

	// ----------------------------------------------------------------------------
	// Initialize
	// ----------------------------------------------------------------------------

	public CmdViewer(String aType)
	{ // font is for output pane on lower half of the window
		super(aType, null,
				new ConsolePane(new Font(NFMVPrefs.prefs.get(NFMVPrefs.OutputFontKey, NFMVPrefs.OutputFontDef),
						Font.PLAIN, NFMVPrefs.prefs.getInt(NFMVPrefs.OutputFontSizeKey, NFMVPrefs.OutputFontSizeDef))));
		soutConsole = (ConsolePane) console;

		setFramePrefs("Commands Window Width", 600, "Commands Window Height", 800);

		makeMenuBar();

		// tool bar icons
		addDefaultToolBar(this);
		addToolBarIcon(null, "saveDocument", "Save this document to a file.", this);

		Class<?> baseClass = JNApplication.main.getClass();
		addToolBarBar();
		ImageIcon goNext = new ImageIcon(baseClass.getResource("Resources/go-next.png"));
		addToolBarIcon(goNext, null, "Run FEA or MPM Analysis.", getRunAnalysisAction());
		ImageIcon goLast = new ImageIcon(baseClass.getResource("Resources/go-last.png"));
		addToolBarIcon(goLast, null, "Check mesh for FEA or MPM Analysis.", checkAnalysisCammand);
		ImageIcon scriptIcon = new ImageIcon(baseClass.getResource("Resources/scripticon.png"));
		addToolBarIcon(scriptIcon, null, "Interpret calculation or control script commands.", interpretCammand);
		ImageIcon doStop = new ImageIcon(baseClass.getResource("Resources/process-stop.png"));
		addToolBarIcon(doStop, null, "Stop currently running FEA or MPM Analysis.", stopModeCommand);

		addToolBarBar();
		ImageIcon showRes = new ImageIcon(baseClass.getResource("Resources/image-x-generic.png"));
		addToolBarIcon(showRes, null, "Show associated simulation results (if available).", showPartnerCommand);

		// change font on editing field on top of window
		Font theFont = new Font(NFMVPrefs.prefs.get(NFMVPrefs.ScriptFontKey, NFMVPrefs.ScriptFontDef), Font.PLAIN,
				NFMVPrefs.prefs.getInt(NFMVPrefs.ScriptFontSizeKey, NFMVPrefs.ScriptFontSizeDef));
		cmdField.textPane.setFont(theFont);

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
			menuBar.add(defaultApplicationMenu()); // Application menu
		JMenu fileMenu = defaultFileMenu(this);
		fileMenu.add(exportXMLCommand);
		menuBar.add(fileMenu); // File menu
		NairnFEAMPMViz.addExamplesMenu(menuBar,"File");

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
		{
			menu.add(JNApplication.main.getOpenHelpAction());
		}
		menu.add(showPartnerCommand);
		menu.addSeparator();
		setWindowMenu(menu);

		// add the menu bar
		setJMenuBar(menuBar);
	}

	// default file menu referring to document target
	public static JMenu defaultFileMenu(JNDocument target)
	{
		JMenu fileMenu = target.defaultFileMenu();
		JMenuItem newMPM = fileMenu.getItem(1);
		newMPM.setAccelerator(
				KeyStroke.getKeyStroke(KeyEvent.VK_N, JNApplication.menuKeyMask() + ActionEvent.SHIFT_MASK));
		return fileMenu;
	}

	// default tool bar icons
	public static void addDefaultToolBar(JNDocument target)
	{
		// tool bar icons
		target.addToolBarBar();
		target.addToolBarIcon(null, null, "Open the preferences window.",
				JNApplication.main.getOpenPreferencesAction());
		target.addToolBarIcon(null, null, "Open the help information window.", JNApplication.main.getOpenHelpAction());

		target.addToolBarBar();
		target.addToolBarIcon(null, "openDocument", "Open a saved document file.", JNApplication.main);
		Class<?> baseClass = JNApplication.main.getClass();
		ImageIcon newMPM = new ImageIcon(baseClass.getResource("Resources/document-new.png"));
		target.addToolBarIcon(newMPM, "newDocumentMPMCmd", "Create a new MPM commands document.", JNApplication.main);
		ImageIcon newFEA = new ImageIcon(baseClass.getResource("Resources/document-newfea.png"));
		target.addToolBarIcon(newFEA, "newDocumentFEACmd", "Create a new FEA command document.", JNApplication.main);
	}

	// ----------------------------------------------------------------------------
	// Interpret commands to XML input commands
	// ----------------------------------------------------------------------------

	// default command run method passed to custom one
	public void runAnalysis()
	{
		runNFMAnalysis(false, NFMAnalysis.FULL_ANALYSIS, null);
	}

	public void runNFMAnalysis(boolean doBackground,int runType,CmdViewer scriptDoc)
	{
		// System.out.println("Check real file");
		// only allowed if the commands have been saved
		if(getFile() == null)
		{	JNApplication.appBeep();
			JNUtilities.showMessage(this,
					"The input commands have to be saved to a file before running an analysis.");
			return;
		}

		// self call sets scriptDoc to self
		if(scriptDoc == null)
			theScript = this;
		else
			theScript = scriptDoc;

		// System.out.println("Check NFMAnalysis task");
		// create once
		if(nfmAnalysis == null)
			nfmAnalysis = new NFMAnalysis(this);

		// what if process is currently running?
		// System.out.println("Check nothing running");
		if(nfmAnalysis.isRunning())
		{	JNApplication.appBeep();
			String message = "An FEA or MPM process is currently running.\nDo you want stop it and start a new process?";
			int result = JOptionPane.showConfirmDialog(this, message, "Continue?", JOptionPane.OK_CANCEL_OPTION);
			if(result == JOptionPane.CANCEL_OPTION)
				return;
			nfmAnalysis.stopRunning();
		}

		// save this document (commands saved before each run - preference would
		// be better)
		if(!saveDocument())
			return;

		// check if XML file
		soutConsole.clear();
		int offset = cmdField.getCommands().indexOf("<?xml ");
		if(offset < 0 || offset > 10)
		{	// interpret commands
			useBackground = doBackground;
			openMesh = runType;
			// call in super class initiates command interpretation
			// System.out.println("Interpreting Commands");
			super.runAnalysis();

			// when interpretation done, will launch the analysis in
			// analysisFinished()
			return;
		}
		else
		{	// look for processors command in XML commands
			processors = 1;
			offset = cmdField.getCommands().indexOf("<!--processors ");
			if(offset > 0)
			{	int endoffset = cmdField.getCommands().indexOf("-->", offset);
				if(endoffset > 0)
				{	String procs = cmdField.getCommands().substring(offset + 15, endoffset);
					Scanner getProcs = new Scanner(procs);
					if(getProcs.hasNextInt())
						processors = getProcs.nextInt();
					if(processors < 1)
						processors = 1;
					getProcs.close();
				}
			}
		}

		// launch analysis with DTD commands in the field
		nfmAnalysis.launchNFMAnalysis(doBackground, runType, cmdField.getCommands(), soutConsole, processors,
				theScript.getScriptParams());
	}

	// when analysis is done, proceed with calculations (if OK)
	public void analysisFinished(boolean status)
	{ // give up on error
		if(status == false || stopCommand == true)
			return;

		if(runningScript)
		{	toFront();
			return;
		}

		// launch analysis with DTD commands in the field
		nfmAnalysis.launchNFMAnalysis(useBackground, openMesh, buildXMLCommands(), soutConsole, processors,
				theScript.getScriptParams());
	}

	// initialize when interpreting commands
	public void initRunSettings() throws Exception
	{
		title = JNApplication.appNameReadable + " Calculations";
		username = null;
		header = new StringBuffer("");
		np = -1;
		consistentUnits = null;
		processors = 1;
		plusSpin = false;
		lnameEl = NO_ELEMENT;
		xmldata = new HashMap<String, String>(10);
		entities = new HashMap<String, String>(10);
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
		skipExtrap = true;
		shapeMethod = "uGIMP";
		archiveRoot = "    <ArchiveRoot>Results/data.</ArchiveRoot>\n";
		archiveTime = "";
		timeStep = "    <TimeStep>1e15</TimeStep>\n";
		maxTime = "";
		cflFactor = "";
		transCflFactor = "";
		CPDIrcrit = "";
		exactTractions = "";
		globalArchive = "";
		feaTemp = null;
		stressFreeTemp = 0.;
		stopCommand = false;
		damping = null;
		pdamping = null;
		xpic = null;
		fbDamping = null;
		pfbDamping = null;
		leaveLimit = null;
		deleteLimit = null;
		extrapolateRigid = null;
		ptsPerElement = null;
		diffusion = null;
		otherDiffusion = new StringBuffer();
		conduction = null;
		gravity = null;
		MMLump = -1;
		MMNormals = -1;
		MMRigidBias = 1.0;
		MMAzimuth = 0.0;
		MMPolar = 0.0;
		ContactPosition = null;
		ContactPositionCracks = null;
		FrictionMM = null;
		customTasks = new StringBuffer("");

		runningScript = false;
		// dictionary of script object variables
		objs = new HashMap<String, Object>(10);
		objs.put("nfmapp",JNApplication.main);
		scriptParams = null;

		// is it called from a script?
		if(theScript != this)
		{	variablesStrs.put("_ScriptMode_","1");
			variablesStrs.putAll(theScript.getVariablesStrs());
		}
	}

	// handle commands
	public void doCommand(String theCmd,ArrayList<String> args) throws Exception
	{	//System.out.println(theCmd+":"+args);

		// if script, switch to script commands
		if(runningScript)
		{	doScriptCommand(theCmd, args);
		}

		else if(theCmd.equals("script"))
		{
			if(openMesh != NFMAnalysis.SCRIPT_ONLY && openMesh != NFMAnalysis.INTERPRET_ONLY)
				throw new Exception("Scripts can only run by using the 'Interpret Commands...'");
			runningScript = true;
			openMesh = NFMAnalysis.SCRIPT_ONLY;
		}

		else if(mats.isInMaterial())
		{	// commands go to material class when
			// material (keep this option first)

			// but first see if language control command
			// (which means cannot match any material property)
			if(theCmd.toLowerCase().equals("color"))
			{	mats.doMaterialProperty(theCmd, args, this);
			}
			else
			{	try
				{	super.doCommand(theCmd, args);
				}
				catch (Exception e)
				{	mats.doMaterialProperty(theCmd, args, this);
				}
			}
		}

		else if(theCmd.equals("title"))
		{	// set analysis title
			if(args.size() < 2)
				throw new Exception("'Title' command does not have a title");
			title = readStringArg(args.get(1));
		}

		else if(theCmd.equals("name"))
		{	// set analysis title
			if(args.size() < 2)
				throw new Exception("'Name' command does not have a name");
			username = readStringArg(args.get(1));
		}

		else if(theCmd.equals("header"))
			header.append(readVerbatim("endheader"));

		else if(theCmd.equals("comment"))
		{	header.append("Comment: ");
			for(int i=1; i<args.size(); i++)
			{
				if(i > 1)
					header.append(", ");
				header.append(readStringArg(args.get(i)));
			}
			header.append("\n");
		}

		else if(theCmd.equals("processors"))
		{	processors = readIntArg(args.get(1));
			if(processors < 1)
				processors = 1;
		}

		else if(theCmd.equals("analysis"))
			doAnalysis(args);

		else if(theCmd.equals("consistentunits"))
			doConsistentUnits(args);

		else if(theCmd.equals("mpmmethod"))
			doMPMMethod(args);

		else if(theCmd.equals("archive"))
			doArchive(args, false);

		else if(theCmd.equals("globalarchive"))
			doGlobalArchive(args);

		else if(theCmd.equals("toarchive"))
			doToArchive(args);

		else if(theCmd.equals("archiveunique"))
			doArchive(args, true);

		else if(theCmd.equals("archivetime"))
			doArchiveTime(args);

		else if(theCmd.equals("globalarchivetime"))
			doGlobalArchiveTime(args);

		else if(theCmd.equals("timestep"))
			doTimeStep(args);

		else if(theCmd.equals("maximumtime"))
			doMaxTime(args);

		else if(theCmd.equals("cflfactor"))
			doCFLFactor(args);

		else if(theCmd.equals("cpdircrit"))
			doCPDIrcrit(args);

		else if(theCmd.equals("exacttractions"))
			doExactTractions(args);

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

		else if(theCmd.equals("transform"))
			regions.AddTransform(args);

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

		else if(theCmd.equals("gridbc"))
			mpmGridBCs.StartMoveLine(args, MPMGridBCs.GRID_BC);

		else if(theCmd.equals("endgridbc"))
			mpmGridBCs.EndMoveBlock(args, MPMGridBCs.GRID_BC);

		else if(theCmd.equals("moveline"))
			mpmGridBCs.StartMoveLine(args, MPMGridBCs.MOVELINE_BC);

		else if(theCmd.equals("movearc"))
			mpmGridBCs.StartMoveLine(args, MPMGridBCs.MOVEARC_BC);

		else if(theCmd.equals("movebox"))
			mpmGridBCs.StartMoveBox(args, MPMGridBCs.MOVEBOX_BC);

		else if(theCmd.equals("endmoveline"))
			mpmGridBCs.EndMoveBlock(args, MPMGridBCs.MOVELINE_BC);

		else if(theCmd.equals("endmovearc"))
			mpmGridBCs.EndMoveBlock(args, MPMGridBCs.MOVEARC_BC);

		else if(theCmd.equals("endmovebox"))
			mpmGridBCs.EndMoveBlock(args, MPMGridBCs.MOVEBOX_BC);

		else if(theCmd.equals("boundaryid"))
			mpmGridBCs.SetBoundaryID(args);

		else if(theCmd.equals("velocity"))
			mpmGridBCs.AddVelocity(args);

		else if(theCmd.equals("concentration") || theCmd.equals("porepressure"))
			mpmGridBCs.AddTempConc(args, MPMGridBCs.ADD_CONCENTRATION);

		else if(theCmd.equals("temperature"))
		{	if(isFEA())
				doTemperature(args);
			else
				mpmGridBCs.AddTempConc(args, MPMGridBCs.ADD_TEMPERATURE);
		}

		else if(theCmd.equals("particlebc"))
			mpmGridBCs.StartMoveLine(args, MPMGridBCs.PARTICLE_BC);

		else if(theCmd.equals("endparticlebc"))
			mpmParticleBCs.EndLoadBlock(args, MPMGridBCs.PARTICLE_BC);

		else if(theCmd.equals("loadline"))
			mpmGridBCs.StartMoveLine(args, MPMGridBCs.LOADLINE_BC);

		else if(theCmd.equals("loadarc"))
			mpmGridBCs.StartMoveLine(args, MPMGridBCs.LOADARC_BC);

		else if(theCmd.equals("loadbox"))
			mpmGridBCs.StartMoveBox(args, MPMGridBCs.LOADBOX_BC);

		else if(theCmd.equals("loadrect"))
			mpmParticleBCs.StartLoadRect(args);

		else if(theCmd.equals("loadtype"))
			mpmParticleBCs.doLoadType(args);

		else if(theCmd.equals("endloadline"))
			mpmParticleBCs.EndLoadBlock(args, MPMGridBCs.LOADLINE_BC);

		else if(theCmd.equals("endloadarc"))
			mpmParticleBCs.EndLoadBlock(args, MPMGridBCs.LOADARC_BC);

		else if(theCmd.equals("endloadrect"))
			mpmParticleBCs.EndLoadBlock(args, MPMGridBCs.LOADRECT_BC);

		else if(theCmd.equals("endloadbox"))
			mpmParticleBCs.EndLoadBlock(args, MPMGridBCs.LOADBOX_BC);

		else if(theCmd.equals("traction"))
			mpmParticleBCs.AddCondition(args, MPMParticleBCs.ADD_TRACTION);

		else if(theCmd.equals("damage"))
			mpmParticleBCs.AddCondition(args, MPMParticleBCs.ADD_DAMAGE);

		else if(theCmd.equals("phasefield"))
			mpmParticleBCs.AddCondition(args, MPMParticleBCs.ADD_PHASEFIELD);

		else if(theCmd.equals("heatflux"))
			mpmParticleBCs.AddCondition(args, MPMParticleBCs.ADD_HEATFLUX);

		else if(theCmd.equals("concentrationflux") || theCmd.equals("porepressureflux"))
			mpmParticleBCs.AddCondition(args, MPMParticleBCs.ADD_CONCENTRATIONFLUX);

		else if(theCmd.equals("origin"))
		{
			if(regions.isInBMPRegion())
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

		else if(theCmd.equals("fill"))
			regions.FillReservoir(args);

		else if(theCmd.equals("bmpregion"))
			regions.StartBMPRegion(args);

		else if(theCmd.equals("angularmom0"))
			regions.AddAngularMom0(args);

		else if(theCmd.equals("endregion"))
			regions.EndRegion(args, theCmd);

		else if(theCmd.equals("hole"))
			regions.StartHole(args);

		else if(theCmd.equals("endhole"))
			regions.EndRegion(args, theCmd);

		else if(theCmd.equals("rect"))
			regions.AddRectOrOval(args, "Rect", 0);

		else if(theCmd.equals("oval"))
			regions.AddRectOrOval(args, "Oval", 0);

		else if(theCmd.equals("polypt"))
			regions.AddPolypoint(args, 0);

		else if(theCmd.equals("arc"))
			regions.AddRectOrOval(args, "Arc", 0);

		else if(theCmd.equals("line"))
		{
			if(isMPM3D())
				regions.AddBox(args, "Line", 0);
			else
				regions.AddRectOrOval(args, "Line", 0);
		}

		else if(theCmd.equals("box"))
			regions.AddBox(args, "Box", 0);

		else if(theCmd.equals("cylinder"))
			regions.AddBox(args, "Cylinder", 0);

		else if(theCmd.equals("torus"))
			regions.AddBox(args, "Torus", 0);

		else if(theCmd.equals("sphere"))
			regions.AddBox(args, "Sphere", 0);

		else if(theCmd.equals("cut"))
			regions.AddCutShape(args);

		else if(theCmd.equals("gridhoriz"))
			gridinfo.doGridAxis(args, 0);

		else if(theCmd.equals("gridvert"))
			gridinfo.doGridAxis(args, 1);

		else if(theCmd.equals("griddepth"))
			gridinfo.doGridAxis(args, 2);

		else if(theCmd.equals("gridrect"))
			gridinfo.doGridRect(args);

		else if(theCmd.equals("gridthickness"))
			gridinfo.doGridThickness(args);
		
		else if(theCmd.equals("tartanborder"))
			gridinfo.doTartanBorder(args);
		
		else if(theCmd.equals("tartanaoi"))
			gridinfo.doTartanAOI(args); 

		else if(theCmd.equals("stressfreetemp"))
			doStressFreeTemp(args);

		else if(theCmd.equals("damping"))
			doDamping(args, "Damping");

		else if(theCmd.equals("pdamping"))
			doDamping(args, "PDamping");

		else if(theCmd.equals("feedbackdamping"))
			doFBDamping(args, "FeedbackDamping");

		else if(theCmd.equals("pfeedbackdamping"))
			doFBDamping(args, "PFeedbackDamping");

		else if(theCmd.equals("leavelimit"))
			doLeaveLimit(args);

		else if(theCmd.equals("deletelimit"))
			doDeleteLimit(args);

		else if(theCmd.equals("extrapolaterigid"))
			doExtrapolateRigid(args);

		else if(theCmd.equals("diffusion"))
			doDiffusion(args, true);

		else if(theCmd.equals("poroelasticity"))
			doDiffusion(args, false);

		else if(theCmd.equals("conduction"))
			doConduction(args);

		else if(theCmd.equals("gravity"))
			doGravity(args);

		else if(theCmd.equals("ptsperelement"))
			doPtsPerElement(args);

		else if(theCmd.equals("multimaterialmode"))
			doMultimaterialMode(args);

		else if(theCmd.equals("contactposition"))
			doContactPosition(args,false);

		else if(theCmd.equals("contactpositioncracks"))
			doContactPosition(args,true);

		else if(theCmd.equals("contactcracks"))
			doContactLaw(args, 0);

		else if(theCmd.equals("friction"))
		{ // Deprecated - Use ContactCracks
			doFriction(args, 0);
		}

		else if(theCmd.equals("imperfectinterface"))
		{ // Deprecated - Use ContactCracks
			doImperfectInterface(args, 0);
		}

		else if(theCmd.equals("contactmm"))
			doContactLaw(args, 1);

		else if(theCmd.equals("frictionmm"))
		{ // deprecated - use ContactMM
			doFriction(args, 1);
		}

		else if(theCmd.equals("imperfectinterfacemm"))
		{ // deprecated - use
				// ContactMM
			doImperfectInterface(args, 1);
		}

		else if(theCmd.equals("crackinterface"))
		{ // Deprecated - used (frict) parameter for newCrack
			doImperfectInterface(args, 3);
		}

		else if(theCmd.equals("jcontour"))
			cracks.doJContour(args);

		else if(theCmd.equals("newcrack"))
			cracks.StartCrack(args);

		else if(theCmd.equals("crackkeypoint"))
			cracks.AddCrackKeypoint(args);

		else if(theCmd.equals("crackfacet"))
			cracks.AddCrackFacet(args,false);

		else if(theCmd.equals("crackplane"))
			cracks.AddCrackFacet(args,true);

		else if(theCmd.equals("growcrack"))
			cracks.GrowCrack(args, 0);

		else if(theCmd.equals("growcrackline"))
			cracks.GrowCrack(args, 1);

		else if(theCmd.equals("growcrackarc"))
			cracks.GrowCrack(args, 2);

		else if(theCmd.equals("crackthickness"))
			cracks.doCrackThickness(args);

		else if(theCmd.equals("cracktractionprop"))
			cracks.doCrackTractionProp(args);

		else if(theCmd.equals("propagate"))
			cracks.doPropagate(args, false);

		else if(theCmd.equals("altpropagate"))
			cracks.doPropagate(args, true);

		else if(theCmd.equals("moveplane"))
			cracks.doMovePlane(args);

		else if(theCmd.equals("propagatelength"))
			cracks.doProagateLength(args);

		else if(theCmd.equals("customtask"))
			doCustomTask(args);

		else if(theCmd.equals("parameter"))
			doParameter(args);

		else if(theCmd.equals("stop"))
		{
			super.doCommand(theCmd, args);
			stopCommand = true;
		}

		else
		{	//System.out.println(args);
			super.doCommand(theCmd, args);
		}
	}

	// Handle all script commands
	//  1. heck direct commands (in alphabetical order) as "commmand arg1,arg2,..."
	//  2. If not found look for object command "obj[.#]].command arg1,arg2,..."
	//		and if found pass to doObjectCommand()
	//  3. If not handle pass to super.doDommand()
	public void doScriptCommand(String theCmd,ArrayList<String> args) throws Exception
	{
		// ----------- CreateDictionary objVar
		if(theCmd.equals("createdictionary"))
		{	// first parameter must be an object name
			String objectVar = getObjVarName(1,args,"CreateDictionary",true);
			
			// create list and add to object variables
			ISDictType dict = new ISDictType();
			objs.put(objectVar,dict);
		}

		// ----------- CreateList objVar,item1,item2,...
		else if(theCmd.equals("createlist"))
		{	// first parameter must be an object name
			String objectVar = getObjVarName(1,args,"CreateList",true);
			
			// create list and add to object variables
			ISListType list = new ISListType(null);
			
			// any more added as strings or objects
			for(int i=2;i<args.size();i++)
			{	String expr = readStringArg(args.get(i));
				Object obj = scriptObjectForKey(expr);
				if(obj!=null)
				{	// add an object
					list.gcis_addObject(obj);
				}
				else
				{	// if not an object, just add the text
					list.gcis_addObject(expr);
				}
			}
			
			// add list to object varibles
			objs.put(objectVar,list);
		}

		// ---------- export path
		// export script console to file by path or null
		else if(theCmd.equals("export"))
		{	File oneDoc = null;
			if(args.size() > 1)
			{	oneDoc = scriptPath(readStringArg(args.get(1)), args, false);
			}
			if(!exportOutput(oneDoc, null))
				throw new Exception("The 'export' command failed.\n" + args);
		}
		
		// ----------- GetObject (objName),(obj).[#i]
		// set objName to and existing object
		else if(theCmd.equals("getobject"))
		{	// object name to return
			String objectVar = getObjVarName(1,args,"GetObject",true);
			
			// second is dereferenced object variable name
			if(args.size()<3)
				throw new Exception("The GetObject needs two parameters.\n" + args);
			String objVar = readStringArg(args.get(2));
			Object objValue = findExistingObject(objVar);
			
			// set to found object or to none type
			if(objValue!=null)
			{	// set new variable to the object
				objs.put(objectVar,objValue);
			}
			else
			{	// set new variable to None object
				objs.put(objectVar,new ISNoneType());
			}
		}
		
		// ---------- openfolder objName,path (omit path for dialog, can be relative path)
		// path can be "_results_" or "_commands_" for front document of that type
		else if(theCmd.equals("open"))
		{	// object name to return
			String objectVar = getObjVarName(1,args,"open",true);

			// file by path or null or _results_ or _commands_
			JNDocument currentDoc = NairnFEAMPMViz.main.frontDocument();
			JNDocument openedDoc = null;
			File docToOpen = null;
			
			if(args.size() > 2)
			{	String openArg = readStringArg(args.get(2));
				if(openArg.toLowerCase().equals("_results_"))
				{	// set openedDoc to first DocViewer, exit if none
					ArrayList<JNDocument> docs = NairnFEAMPMViz.main.getDocuments();
					for(int i=0;i<docs.size();i++)
					{	JNDocument adoc = docs.get(i);
						if(adoc.getClass().equals(DocViewer.class))
						{	openedDoc = adoc;
							break;
						}
					}
					if(openedDoc==null)
						throw new Exception("An 'Open' command found no such results document.\n" + args);
					
				}
				else if(openArg.toLowerCase().equals("_commands_"))
				{	// set openedDoc to first CmdViewer, exit if none
					ArrayList<JNDocument> docs = NairnFEAMPMViz.main.getDocuments();
					for(int i=0;i<docs.size();i++)
					{	JNDocument adoc = docs.get(i);
						if(adoc.getClass().equals(CmdViewer.class) && (adoc!=this))
						{	openedDoc = adoc;
							break;
						}
					}
					if(openedDoc==null)
						throw new Exception("An 'Open' command found no such commands document.\n" + args);
					
				}
				else
				{	// get path and set openedDoc if opened
					docToOpen = scriptPath(readStringArg(openArg), args, false);
					docToOpen = new File(docToOpen.getCanonicalPath());
					openedDoc = NairnFEAMPMViz.main.findDocument(docToOpen);
				}
				
				if(openedDoc!=null)
				{	openedDoc.setVisible(true);
					openedDoc.toFront();
				}
			}

			// if no existing document, try to open file (with path or null)
			if(openedDoc==null)
				NairnFEAMPMViz.main.openDocument(docToOpen);
			
			// Hack - wait for it to be front window, give up after some tries
			int maxLoops = 10;
			while(currentDoc == NairnFEAMPMViz.main.frontDocument() && maxLoops>0)
			{	Thread.sleep(20);
				maxLoops--;
			}
			if(maxLoops==0)
			{	running = false;
				return;
			}
			objs.put(objectVar, NairnFEAMPMViz.main.frontDocument());
		}

		// ---------- openfolder varName,title
		else if(theCmd.equals("openfolder"))
		{	// variable to return the path
			String varName = getObjVarName(1,args,"openFolder",false);

			// optional dialog title
			String fldrTitle = "Select a folder";
			if(args.size() > 2)
			{	String userTitle = readStringArg(args.get(2));
				if(userTitle.length() > 0)
					fldrTitle = userTitle;
			}

			// if path use it, otherwise dialog box
			String fldrPath = "";
			if(args.size() > 3)
			{	File oneFldr = scriptPath(readStringArg(args.get(3)), args, true);

				// only need to create if does not exist
				if(!oneFldr.exists())
				{	if(!oneFldr.mkdirs())
						throw new Exception("File error creating the folder(s).\n" + args);
				}

				// get file path
				fldrPath = oneFldr.getCanonicalPath();
			}
			else
			{	JFileChooser fldrChooser = new JFileChooser();
				fldrChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				NFMVPrefs.setWorkspace(fldrChooser);
				fldrChooser.setDialogTitle(fldrTitle);
				int result = fldrChooser.showSaveDialog(this);
				if(result == JFileChooser.APPROVE_OPTION)
				{	fldrPath = fldrChooser.getSelectedFile().getPath();
				}
			}

			// save in variable with terminal path delimiter
			if(fldrPath.length() > 0)
			{	if(JNApplication.isWindowsOS())
					fldrPath += "\\";
				else
					fldrPath += "/";
			}
			variablesStrs.put(varName, fldrPath);
		}
		
		// ---------- userInput "#var",(title),(msg),(initial)
		else if(theCmd.contentEquals("userinput"))
		{	// variable to return the path
			String varName = getObjVarName(1,args,"userinput",false);

			// dialog box options
			String dtitle = "";
			if(args.size() > 2)
			{	dtitle = readStringArg(args.get(2));
			}
			String prompt = "Enter some text";
			if(args.size() > 3)
			{	prompt = readStringArg(args.get(3));
			}
			String initText = "";
			if(args.size() > 4)
			{	initText = readStringArg(args.get(4));
			}
			String newExpr = (String)JOptionPane.showInputDialog(null,prompt,dtitle,
					JOptionPane.PLAIN_MESSAGE,null,null,initText);
			
			if(newExpr!=null)
			{	variablesStrs.put(varName+"[1]", "OK");
				variablesStrs.put(varName+"[2]",newExpr);
			}
			else
			{	variablesStrs.put(varName+"[1]", "Cancel");
				variablesStrs.put(varName+"[2]","");
			}
			variablesStrs.put(varName+"[0]", "2");
		}
		
		// ---------- object commands
		// if no direct command found look for obj.[#i.]command format object command
		else
		{	// look for object command
			String objCmd = args.get(0);
			int dot = objCmd.indexOf(".");
			if(dot > 0)
			{	String objName = objCmd.substring(0, dot);
				Object obj = scriptObjectForKey(objName);
				if(obj != null)
				{	doObjectCommand(obj, objCmd.substring(dot+1), args);
					return;
				}
			}

			// if not an object command, pass to super to handle it
			// System.out.println(args);
			super.doCommand(theCmd, args);
		}
	}

	// handle commands to object in obj
	// theCmd is period delimit list of atoms
	// args are the parameters to the command
	public void doObjectCommand(Object obj,String theCmd,ArrayList<String> args) throws Exception
	{
		// obj is object
		// if obj is a list, the command may be #i.theCmd or just theCmd
		
		// will have 1 or 2
		String [] atoms = theCmd.trim().split("[.]");
		
		// is there a list selector?
		int i=0;
		if(obj.getClass().equals(ISListType.class))
		{	String nextAtom = grabAtom(atoms,i);
			int index = decodeIndex(nextAtom);
			if(index>=0)
			{	if(index>=((ISListType)obj).count())
					throw new Exception("Command targeting a list item has index that is out of bounds.");
				// get list element, get rest as new command
				obj = ((ISListType)obj).objectAtIndex(index);
				i++;
				theCmd = grabAtom(atoms,i);
			}
			else if(!ISListType.supportsCommand(nextAtom))
			{	// [#i] was not a number AND obj.cmd is not a command for a list
				throw new Exception("Command targeting a list must be followed by valid locator or list command");
			}
		}
		
		// Object command for obj is theCmd
		// Possible commands follow in alphabetical order
		theCmd = theCmd.toLowerCase();
		
		// ----------- addList item1,item2,...
		if(theCmd.equals("addobject"))
		{	if(!obj.getClass().equals(ISListType.class))
				throw new Exception("The addObject command can only by used on lists.\n" + args);

			ISListType list = (ISListType)obj;
			
			// any more added as strings or objects
			for(i=1;i<args.size();i++)
			{	String objVar = readStringArg(args.get(i));
				Object objValue = findExistingObject(objVar);
				if(objValue!=null)
					list.gcis_addObject(objValue);
			}
		}

		// ----------- addList item1,item2,...
		else if(theCmd.equals("addstring"))
		{	if(!obj.getClass().equals(ISListType.class))
				throw new Exception("The addString command can only by used on lists.\n" + args);

			ISListType list = (ISListType)obj;
			
			// any more added as strings or objects
			for(i=1;i<args.size();i++)
			{	String expr = readStringArg(args.get(i));
				list.gcis_addObject(expr);
			}
		}
		
		// ---------- export path
		// export to a file
		else if(theCmd.equals("export"))
		{	if(!obj.getClass().equals(CmdViewer.class))
				throw new Exception("The export command can only by used on commands documents.\n" + args);

			// file by path or null
			File oneDoc = null;
			if(args.size() > 1)
			{	String fPath = readStringArg(args.get(1));
				if(fPath.length() < 2)
					throw new Exception("The export command has empty path name.\n" + args);
				if(fPath.charAt(0) != '/' && fPath.charAt(1) != ':')
					oneDoc = new File(getFile().getParent(), fPath);
				else
					oneDoc = new File(fPath);
			}

			if(!((CmdViewer) obj).exportOutput(oneDoc, null))
				throw new Exception("The export command failed.\n" + args);
		}
		
		// ---------- obj[#i].get ("#var" or objName)
		else if(theCmd.equals("get"))
		{   // get the object or variable name
			if(args.size()<2)
				throw new Exception("get command missing its variable or object name.\n"+args);
			String objVar = readStringArg(args.get(1));

			// if target object is a string or number extracted from list then done
			if(obj.getClass().equals(String.class))
			{	if(!JNExpression.validVariableName(objVar))
				{	throw new Exception(
						"First parameter in get command for a string must be a valid variable name.\n"+args);
				}
				variablesStrs.put(objVar,(String)obj);
				return;
			}
			else if(obj.getClass().equals(Double.class))
			{	if(!JNExpression.validVariableName(objVar))
				{	throw new Exception(
						"First parameter in get command for a number must be a valid variable name.\n"+args);
				}
				variablesStrs.put(objVar,JNUtilities.formatDouble(((Double)obj).doubleValue()));
				return;
			}

			// rest are objects
			if(!validObjectName(objVar))
			{	throw new Exception(
					"First parameter in get command for an object must be a valid object variable name.\n"+args);
			}

			// If not any more arguments, set that object variable to obj
			if(args.size()==2)
			{	// if no more arguments, set to "self"
				objs.put(objVar,obj);
				return;
			}
			
			// look for object attribute of current object
			// Only allow object attributes (to-many or to-one) and not properties
			String attr = readStringArg(args.get(2));
			Object objValue = null;
			
			// handle each object type here
			if(obj.getClass().equals(DocViewer.class))
				objValue = ((DocViewer)obj).gcis_getObjectAttribute(attr,this);
			
			else if(obj.getClass().equals(ISDictType.class))
				objValue = ((ISDictType)obj).gcis_getObjectAttribute(attr,this);
			
			else if(obj.getClass().equals(NairnFEAMPMViz.class))
				objValue = ((NairnFEAMPMViz)obj).gcis_getObjectAttribute(attr,this);

			// set variable or error
			if(objValue != null)
				objs.put(objVar,objValue);
			else
			{	throw new Exception(
					"Class "+obj.getClass().toString()+" does not have that object attributes.\n"+args);
			}
		}

		// ----------- insertObject obj,index,...
		else if(theCmd.equals("insertobject"))
		{	if(!obj.getClass().equals(ISListType.class))
				throw new Exception("The insertObject command can only by used on lists.\n" + args);
			if(args.size()<3)
				throw new Exception("The insertObject command needs object and insert location.\n" + args);

			String objVar = readStringArg(args.get(1));
			Object objValue = findExistingObject(objVar);
			if(objValue==null)
				throw new Exception("The insertObject command needs an existing.\n" + args);
			int index = readIntArg(args.get(2));
			if(!((ISListType)obj).gcis_insertObject(objValue,index))
			{	throw new Exception("Insert index for a list is out of bounds.\n" + args);
			}
		}

		// ----------- insertString expr,index,...
		else if(theCmd.equals("insertstring"))
		{	if(!obj.getClass().equals(ISListType.class))
				throw new Exception("The insertString command can only by used on lists.\n" + args);
			if(args.size()<3)
				throw new Exception("The insertString command needs string and insert location.\n" + args);

			String expr = readStringArg(args.get(1));
			int index = readIntArg(args.get(2));
			if(!((ISListType)obj).gcis_insertObject(expr,index))
			{	throw new Exception("Insert index for a list is out of bounds.\n" + args);
			}
		}
		
		// ---------- interpret
		// interpret commands, no arguments
		else if(theCmd.equals("interpret"))
		{	if(!obj.getClass().equals(CmdViewer.class))
				throw new Exception("The interpret command can only by used on commands documents.\n" + args);

			scriptParams = null;
			((CmdViewer) obj).runNFMAnalysis(false, NFMAnalysis.INTERPRET_ONLY, this);

			// wait for interpret to be done
			while (runningScript)
			{	Thread.sleep(100);
				if(!((CmdViewer) obj).isRunning())
					break;
			}
			if(((CmdViewer) obj).isRunning())
			{	((CmdViewer) obj).stopRunning();
				throw new Exception("The script was stopped.\n" + args);
			}
		}

		// ---------- plottable (table)
		// plot table of data in a new plot window
		else if(theCmd.equals("plottable"))
		{	// run obj,plottable data
			if(!obj.getClass().equals(DocViewer.class))
				throw new Exception("The plottable command can only by used on results documents.\n" + args);
			
			if(args.size() > 1)
			{	String tableData = readStringArg(args.get(1));
				XYPlotWindow plotWindow = new XYPlotWindow((DocViewer)obj);
				plotWindow.plotPPTable(tableData);
			}
			else
				throw new Exception("The plottable command had no data.\n" + args);
		}

		// -----------  remove index
		else if(theCmd.equals("remove"))
		{	if(!obj.getClass().equals(ISListType.class))
				throw new Exception("The remove command can only by used on lists.\n" + args);
			if(args.size()<2)
				throw new Exception("The remove command needs an index to remove.\n" + args);

			int index = readIntArg(args.get(1));
			if(!((ISListType)obj).gcis_removeObjectAtIndex(index))
			{	throw new Exception("Remove index for a list is out of bounds.\n" + args);
			}
		}

		// ---------- removeValueForKey key
		else if(theCmd.equals("removevalueforkey"))
		{	if(!obj.getClass().equals(ISDictType.class))
				throw new Exception("The removeValueForKey command can only by for dictionaries.\n" + args);
			if(args.size()<2)
				throw new Exception("The removeValueForKey command needs key.\n" + args);
			
			// first parameter must a string expression for the key
			String key = readStringArg(args.get(1));
			((ISDictType)obj).gcis_removeValueForKey(key);
		}
			
		// ---------- run (resDoc),(path)
		// run analysis and set resDoc to results document
		else if(theCmd.equals("run"))
		{	// run obj,outpath
			if(!obj.getClass().equals(CmdViewer.class))
				throw new Exception("The run command can only by used on commands documents.\n" + args);

			// need to provide path to save the file
			if(args.size() < 3)
				throw new Exception("Run command needs object name and output file path.\n" + args);

			// return object name
			String objectVar = getObjVarName(1,args,"run",true);

			// get path
			File outDoc = scriptPath(readStringArg(args.get(2)), args, false);
			scriptParams = new ArrayList<String>(1);
			String scriptPath = outDoc.getCanonicalPath();
			File scriptOutput = new File(scriptPath);
			if(!scriptOutput.getParentFile().exists())
				throw new Exception("The folder selected for output does not exist.\n" + args);
			scriptParams.add(scriptPath);

			// start analysis
			NFMVPrefs.setRemoteMode(false);
			((CmdViewer) obj).setVisible(true);
			((CmdViewer) obj).toFront();
			((CmdViewer) obj).runNFMAnalysis(false, NFMAnalysis.FULL_ANALYSIS, this);

			// wait for run to be done
			while(runningScript)
			{	Thread.sleep(1000);
				if(!((CmdViewer) obj).isRunning())
					break;
			}
			if(((CmdViewer) obj).isRunning())
			{	((CmdViewer) obj).stopRunning();
				throw new Exception("A script process was stopped.\n" + args);
			}
			NFMVPrefs.restoreRemoteMode();

			// set obj to output document
			objs.put(objectVar, ((CmdViewer) obj).getLinkedResults());
			
			// exit if was aborted with partial results
			// but this only works to abort the first calculations
			if(((CmdViewer) obj).wasAborted())
			{	throw new Exception("Script calculations were stopped.\n" + args);
			}
		}

		// ---------- runremote obj,outpath,outoption,localpath,localoption
		// run remotely
		else if(theCmd.equals("runremote"))
		{	if(!obj.getClass().equals(CmdViewer.class))
				throw new Exception("The runRemote command can only by used on commands documents.\n" + args);

			// need to provide remote path on the server
			// folder/.../name - requires slashes and at least one parent folder
			if(args.size() < 3)
				throw new Exception("RunRemote command needs object name and remote file path.\n" + args);

			// get object variable
			String objectVar = getObjVarName(1,args,"runRemote",true);

			// get remote path (require internal /, but check for .mpm or .fea
			// when run)
			String remotePath = readStringArg(args.get(2));
			int offset = remotePath.lastIndexOf("/");
			if(offset < 1 || offset >= remotePath.length() - 1)
			{	throw new Exception(
					"The remote folder must have at least a folder and a name (e.g., 'folder/name').\n" + args);
			}

			// get output option (default is overwrite or unique or clear)
			String remoteOption = "overwrite";
			if(args.size() > 3)
			{	remoteOption = readStringArg(args.get(3)).toLowerCase();
				if(!remoteOption.equals("overwrite") && !remoteOption.equals("unique") && !remoteOption.equals("clear"))
					throw new Exception("Invalid remote output option (" + remoteOption + "?)\n" + args);
			}

			// get local save option first
			String localOption = "download";
			if(args.size() > 5)
			{	localOption = readStringArg(args.get(5)).toLowerCase();
				if(!localOption.equals("download") && !localOption.equals("nodownload") && !localOption.equals("home"))
					throw new Exception("Invalid local folder option (" + localOption + "?)\n" + args);
			}

			// get local save path
			String outFolder = "";
			if(args.size() > 4)
			{	outFolder = readStringArg(args.get(4));
				if(!localOption.equals("nodownload"))
				{	File outDoc = scriptPath(outFolder, args, true);
					if(!outDoc.exists())
						throw new Exception("The local download folder must already exist\n" + args);
					outFolder = outDoc.getCanonicalPath();
				}
				else
				{	outFolder = "";
				}
			}
			else
				localOption = "nodownload";

			// pack script parameters
			scriptParams = new ArrayList<String>(4);
			scriptParams.add(remotePath);
			scriptParams.add(remoteOption);
			scriptParams.add(outFolder);
			scriptParams.add(localOption);

			// start analysis
			NFMVPrefs.setRemoteMode(true);
			((CmdViewer) obj).runNFMAnalysis(false, NFMAnalysis.FULL_ANALYSIS, this);

			// wait for run to be done
			while (runningScript)
			{	Thread.sleep(1000);
				if(!((CmdViewer) obj).isRunning())
					break;
			}
			if(((CmdViewer) obj).isRunning())
			{	((CmdViewer) obj).stopRunning();
				throw new Exception("The script was stopped.\n" + args);
			}
			NFMVPrefs.restoreRemoteMode();

			// set obj to output document
			if(!localOption.equals("nodownload"))
				objs.put(objectVar, ((DocViewer) NairnFEAMPMViz.main.frontDocument()).resDoc);
		}
		
		// ---------- setNumberForKey number,key
		// set object in a dictionary
		else if(theCmd.equals("setnumberforkey"))
		{	if(!obj.getClass().equals(ISDictType.class))
				throw new Exception("The setNumberForKey command can only by for dictionaries.\n" + args);
			if(args.size()<3)
				throw new Exception("The setNumberForKey needs two arguments.\n" + args);

			// first parameter must a numeric expression
			double adble = readDoubleArg(args.get(1));
		
			// second is key
			String key = readStringArg(args.get(2));
			((ISDictType)obj).gcis_setObjectforKey(new Double(adble),key);
		}
		
		// ---------- setObjectForKey object,key
		// set object in a dictionary
		else if(theCmd.equals("setobjectforkey"))
		{	if(!obj.getClass().equals(ISDictType.class))
				throw new Exception("The setObjectForKey command can only by for dictionaries.\n" + args);
			if(args.size()<3)
				throw new Exception("The setObjectForKey needs two arguments.\n" + args);

			// first parameter must be an existing object name
			String objVar = readStringArg(args.get(1));
			Object objValue = findExistingObject(objVar);
			if(objValue==null)
				throw new Exception("First parameter to setObjectForKey command must be a valid object.\n" + args);
				
			// second is key
			String key = readStringArg(args.get(2));
			((ISDictType)obj).gcis_setObjectforKey(objValue,key);
		}

		// ---------- setStringForKey string,key
		// set object in a dictionary
		else if(theCmd.equals("setstringforkey"))
		{	if(!obj.getClass().equals(ISDictType.class))
				throw new Exception("The setStringForKey command can only by for dictionaries.\n" + args);
			if(args.size()<3)
				throw new Exception("The setStringForKey needs two arguments.\n" + args);

			// first parameter must a string expression
			String objValue = readStringArg(args.get(1));
		
			// second is key
			String key = readStringArg(args.get(2));
			((ISDictType)obj).gcis_setObjectforKey(objValue,key);
		}
		
		// ---------- timeplot objName,settings
		// plot time data and return list of two lists for x and y data
		else if(theCmd.equals("timeplot"))
		{	if(!obj.getClass().equals(DocViewer.class))
				throw new Exception("The timeplot command can only be used for results documents.\n" + args);
		
			// first parameter must object name
			String objname = getObjVarName(1,args,"timeplot",true);
	
			// second parameter must be an existing dictionary with settings
			if(args.size()<3)
				throw new Exception("The timeplot command missing plot settings.\n" + args);
			String objVar = readStringArg(args.get(2));
			Object objValue = findExistingObject(objVar);
			if(objValue==null)
				throw new Exception("The timeplot command missing plot settings.\n" + args);
			
			// generate the plot
			ISListType plotResults = ((DocViewer)obj).scriptTimeplot((ISDictType)objValue);
			objs.put(objname, plotResults);
		}
		
		// ---------- valueForKey (#var or objName),key
		else if(theCmd.equals("valueforkey"))
		{	if(!obj.getClass().equals(ISDictType.class))
				throw new Exception("The valueForKey command can only by for dictionaries.\n" + args);
			if(args.size()<3)
				throw new Exception("The valueForKey command needs variable and key.\n" + args);
			
			// get object or variable name
			String objVar = readStringArg(args.get(1));
			
			// second parameter must a string expression for the key
			String key = readStringArg(args.get(2));
			
			Object getObj = ((ISDictType)obj).gcis_objectForKey(key);
			
			// set object from the dictionary
			if(getObj!=null)
			{	// if fetched is a string, then set a string variable
				if(getObj.getClass().equals(String.class))
				{	if(!JNExpression.validVariableName(objVar))
					{	throw new Exception(
							"First parameter in valueForKey command for a string must be a valid variable name.\n"+args);
					}
					variablesStrs.put(objVar,(String)getObj);
					return;
				}
				else if(getObj.getClass().equals(Double.class))
				{	if(!JNExpression.validVariableName(objVar))
					{	throw new Exception(
							"First parameter in valueForKey command for a number must be a valid variable name.\n"+args);
					}
					variablesStrs.put(objVar,JNUtilities.formatDouble(((Double)getObj).doubleValue()));
					return;
				}

				// rest are objects
				if(!validObjectName(objVar))
				{	throw new Exception(
						"First parameter in valueForKey command for an object must be a valid object name.\n"+args);
				}
				objs.put(objVar,getObj);
			}
			else
			{	// is it a variable beginning in '#'
				if(JNExpression.validVariableName(objVar))
				{	// undefine the object
					if(variablesStrs.get(objVar)!=null)
						variablesStrs.remove(objVar);
					return;
				}
				else if(validObjectName(objVar))
				{	// set to NoneType
					objs.put(objVar,new ISNoneType());
					return;
				}
				
				throw new Exception(
					"First parameter in valueForKey command must be valid variable or object name.\n"+args);
			}
		}

		// ---------- xypplot objName,settings
		// plot mesh 2D data and return list of two lists for x and y data
		else if(theCmd.equals("xyplot"))
		{	if(!obj.getClass().equals(DocViewer.class))
				throw new Exception("The xyplot command can only be used for results documents.\n" + args);
		
			// first parameter must object name
			String objname = getObjVarName(1,args,"xyplot",true);
	
			// second parameter must be an existing dictionary with settings
			if(args.size()<3)
				throw new Exception("The xyplot command missing plot settings.\n" + args);
			String objVar = readStringArg(args.get(2));
			Object objValue = findExistingObject(objVar);
			if(objValue==null)
				throw new Exception("The xyplot command missing plot settings.\n" + args);
			
			// generate the plot
			ISListType plotResults = ((DocViewer)obj).scriptXYplot((ISDictType)objValue);
			objs.put(objname, plotResults);
		}
		// ---------- invalid object command is an error
		else
			throw new Exception("An unrecognized object command.\n" + args);
	}

	// Look for existing object with optional deference to a list
	// 1. If just object name, return its obj
	// 2. If not look for obj.(locator)
	//		a. if obj not a valid list, return nil
	//		b. if valid object return it (but if gets string return nil)
	//		c. If not found (by number, id or name) return nil
	public Object findExistingObject(String objName)
	{
		// check name first
		Object obj = scriptObjectForKey(objName);
		if(obj!=null) return obj;
		
		// break up by periods
		String [] atoms = objName.trim().split("[.]");

		// can only be 2 for list.index
		if(atoms.length!=2) return null;
		
		// first must be object variable (and cannot be an expression, but can be a variable)
		String objectVar = grabAtom(atoms,0);
		char firstChar = objectVar.charAt(0);
		if(firstChar=='#')
		{	String deref = variablesStrs.get(objectVar);
			if(deref!=null) objectVar = deref;
		}
		else if(firstChar=='^' && objectVar.length()>1)
		{	objectVar = objectVar.substring(1);
			String deref = variablesStrs.get(objectVar);
			if(deref!=null) objectVar = deref;
		}

		// see if ended up with object variable for a list
		obj = objs.get(objectVar);
		if(obj==null) return null;
		if(!obj.getClass().equals(ISListType.class)) return null;

		// look for index into the array
		String nextAtom = grabAtom(atoms,1);
		int index = decodeIndex(nextAtom);
		if(index>=0)
		{	if(index>=((ISListType)obj).count()) return null;
				
			// switch obj to be an object in the list
			obj = ((ISListType)obj).objectAtIndex(index);
				
			// if list element is a string or number, then not a valid object
			if(obj.getClass().equals(String.class) || obj.getClass().equals(Double.class))
				return null;
			
			// return this object
			return obj;
		}
		
		// not found in the list
		return null;
	}

	// get argument that should be object variable name
	public String getObjVarName(int argNum,ArrayList<String> args,String acmd,boolean forObj) throws Exception
	{	if(args.size() < argNum+1)
			throw new Exception("The "+acmd+" command must provide an object or variable name.\n"+args);
		String objectVar = readStringArg(args.get(argNum));
		if(forObj)
		{	if(!validObjectName(objectVar))
				throw new Exception("The "+acmd+" command must provide a valid object name.\n"+args);
		}
		else
		{	if(!JNExpression.validVariableName(objectVar))
				throw new Exception("The "+acmd+" command must provide a valid variable name.\n"+args);
		}
		return objectVar;
	}

	// check next item in list of atoms
	// if begins in # see if matches a variable (can be #x[3], but not #x[#i])
	// note that variables without # are not interpreted here
	public String grabAtom(String [] atoms,int i)
	{
		if(i >= atoms.length) return "";
		
		// grab and exit if zero or 1 characters
		String atom = atoms[i];
		if(atom.length()<2) return atom;
		
		char firstChar = atom.charAt(0);
		if(firstChar=='#')
		{	String deref = variablesStrs.get(atom);
			if(deref!=null) atom = deref;
		}
		
		return atom;
	}
	
	// Decode atom in @ expression or in object command. The next
	// atom might be a number to pick an element from the list in obj
	// Requirements: a number (all digitis)
	public int decodeIndex(String atom)
	{
		// empty or null not a number
		if(atom==null) return -1;
		if(atom.length()==0) return -1;
		
		// all can accept a number
		boolean allDigits = true;
		for(int i=0;i<atom.length();i++)
		{	char c = atom.charAt(i);
			if(c<'0' || c>'9')
			{	allDigits = false;
				break;
			};
		}
		if(allDigits) return Integer.parseInt(atom);
		return -1;
	}
	
	// check for object variable when in a script
	public boolean isAltVariable(String testVar)
	{	if(!runningScript) return false;
		Object obj = scriptObjectForKey(testVar);
		if(obj==null) return false;
		if(obj.getClass().equals(ISNoneType.class)) return false;
		return true;
	}

	
	// get script object with name in objectVar
	// If objectVar begins in '^', check rest for valid variable
	// return nil if all fails
	public Object scriptObjectForKey(String objectVar)
	{
		Object obj = objs.get(objectVar);
		
		// if not found, check for &var where var is variable containing an object name
		if(obj==null && objectVar.length()>1)
		{	if(objectVar.charAt(0)=='^')
			{	String var = objectVar.substring(1);
				objectVar = variablesStrs.get(var);
				if(objectVar!=null)
					obj = objs.get(objectVar);
			}
		}
		
		return obj;
	}

	// object names begin in letter (not '#')
	// Rest letters, numbers, and underscore
	public static boolean validObjectName(String v)
	{	// need at least letter
		if(v.length() < 1)
			return false;
		// other letters letter or number
		for(int i = 0; i < v.length(); i++)
		{	char c = v.charAt(i);
			if((c > 'z' || c < 'a') && (c > 'Z' || c < 'A'))
			{	// first must be letter
				if(i == 0)
					return false;
				// others can be numbers of underscore
				if((c > '9' || c < '0') && c != '_')
					return false;
			}
		}
		return true;
	}

	// decode argument to path for a script
	// allows relative or full path and allows Mac/Linux or Windows
	// if file exists, it must be folder or file if wantFolder is true or false
	public File scriptPath(String fPath,ArrayList<String> args,boolean wantFolder) throws Exception
	{	// empty is not allowed
		if(fPath.length() == 0)
			throw new Exception("'" + args.get(0) + "' command has empty path name.\n" + args);

		// Mac/Linux full path begins in / at at lease 1 more character
		File oneDoc = null;
		if(fPath.charAt(0) == '/')
		{ // needs at least on more letter
			if(fPath.length() < 2)
				throw new Exception("'" + args.get(0) + "' command has incomplete full path.\n" + args);
			oneDoc = new File(fPath);
		}
		else if(fPath.length() > 3)
		{	// full Windows full path needs "c:\a" or at least 4 letters with : an \
			// in 2nd and 3rd
			if(fPath.charAt(1) == ':' && fPath.charAt(2) == '\\')
				oneDoc = new File(fPath);
		}

		// it is a relative path
		if(oneDoc == null)
			oneDoc = new File(getFile().getParent(), fPath);

		// if already exists, it better be a folder
		if(oneDoc.exists())
		{	if(wantFolder)
			{	if(!oneDoc.isDirectory())
					throw new Exception("A specified folder name already exists but is not a folder.\n" + args);
			}
			else
			{	if(oneDoc.isDirectory())
					throw new Exception("A specified file name already exists but it is a folder.\n" + args);
			}
		}

		// return it
		return oneDoc;
	}
	
	// convert @ expression to String
	public String getScriptAtString(String s)
	{	
		String badResult = "ERROR: expression error";
		
	    // break up by periods
		String [] atoms = s.trim().split("[.]");
	    if(atoms.length<2)
	    {   return "ERROR: bar @ expression: "+s;
	    }
	    
		// first must be object variable (and cannot be an expression, but can be a variable)
	    String objectVar = atoms[0];
		// remove the "@"
	    if(objectVar.length()>1) objectVar = objectVar.substring(1);
		// check for variable
		char firstChar = objectVar.charAt(0);
		if(firstChar=='#')
		{	String deref = variablesStrs.get(objectVar);
			if(deref!=null) objectVar = deref;
		}
		else if(firstChar=='^' && objectVar.length()>1)
		{	objectVar = objectVar.substring(1);
			String deref = variablesStrs.get(objectVar);
			if(deref!=null) objectVar = deref;
		}
		
		// see if ended up with object variable
	    Object obj = objs.get(objectVar);
	    if(obj==null)
	    {   return "ERROR: @ expression with invalid object: "+s;
	    }
	    
		// if obj is an array, get index into the array
		String nextAtom;
	    int i=1;
	    if(obj.getClass().equals(ISListType.class))
		{	// get number (or other)
			nextAtom = grabAtom(atoms,i);
			
			// attributes of lists itself
			if(nextAtom.equals("count") || nextAtom.equals("length"))
				return Integer.toString(((ISListType)obj).count());
			else if(nextAtom.equals("class"))
			{	return "list";
			}
			
			// is it an integer?
			int index = decodeIndex(nextAtom);
			if(index>=0)
			{	if(index>=((ISListType)obj).count())
					return "ERROR: list index in @ expression out of bounds: "+s;
				
				// switch obj to be an object in the list
				obj = ((ISListType)obj).objectAtIndex(index);
				i++;
				
				// if list element is a string or number, just return that string
				// subsequent atoms ignored, except for class
				if(obj.getClass().equals(String.class))
				{	nextAtom = grabAtom(atoms,i);
					if(nextAtom.equals("class"))
						return "string-number";
					return (String)obj;
				}
				else if(obj.getClass().equals(Double.class))
				{	nextAtom = grabAtom(atoms,i);
					if(nextAtom.equals("class"))
						return "string-number";
					return ((Double)obj).toString();
				}
			}
			else
			{	return "ERROR: list in @ must be followed by a valid index: "+s;
			}
		}

		// In string expression, we are at an object and need atom to get an attribute
		if(i>=atoms.length)
		{	return "ERROR: @ expression is lacking a property: "+s;
		}

		// case sensitive
		nextAtom = grabAtom(atoms,i);
		String expr = null;
		
		// handle each object type here
		if(obj.getClass().equals(DocViewer.class))
			expr = ((DocViewer)obj).gcis_getAttribute(atoms,i,this);
		
		else if(obj.getClass().equals(CmdViewer.class))
			expr = ((CmdViewer)obj).gcis_getAttribute(atoms,i,this);
		
		else if(obj.getClass().equals(ISDictType.class))
			expr = ((ISDictType)obj).gcis_getAttribute(atoms,i,this);
		
		else if(obj.getClass().equals(ISListType.class))
			expr = ((ISListType)obj).gcis_getAttribute(atoms,i,this);
			
		else if(obj.getClass().equals(ISNoneType.class))
			expr = ((ISNoneType)obj).gcis_getAttribute(atoms,i,this);
			
		else if(obj.getClass().equals(NairnFEAMPMViz.class))
			expr = ((NairnFEAMPMViz)obj).gcis_getAttribute(atoms,i,this);
			
		else if(obj.getClass().equals(MPMArchive.class))
			expr = ((MPMArchive)obj).gcis_getAttribute(atoms,i,this);
			
		if(expr==null)
		{	return "ERROR: class "+obj.getClass().toString()+" does not have any properties: "+s;
		}
		
		// return the result
		return expr;
	}

	// Scripting attributes for internal scripts for commands document
	public String gcis_getAttribute(String [] atoms,int i,CmdViewer server)
	{
	    String attr = server.grabAtom(atoms,i);
	    
		if(attr.equals("get"))
		{	i++;
			if(i >= atoms.length)
				return "ERROR: The variable name is missing";
			attr = server.grabAtom(atoms,i);
			return variablesStrs.get(attr);
		}

		else if(attr.equals("name"))
		{	File theFile = getFile();
			if(theFile!=null)
				return theFile.getName();
			return "Untitled";
		}

		else if(attr.equals("path"))
		{	File theFile = getFile();
			if(theFile!=null)
				return theFile.getPath();
			return "Untitled";
		}

		else if(attr.equals("folder"))
		{	File theFile = getFile();
			if(theFile!=null)
				return theFile.getParent();
			return "";
		}
	    
		else if(attr.equals("class"))
			return "CommandDocument";

		return null;
	}

	// Analysis (type),(element)
	public void doAnalysis(ArrayList<String> args) throws Exception
	{
		if(np >= 0)
			throw new Exception("Only one 'Analysis' command is allowed.");

		if(args.size() < 2)
			throw new Exception("'Analysis' command has no argument.");

		// options
		HashMap<String, Integer> options = new HashMap<String, Integer>(10);
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
		options.put("plane strain mpm+ps", new Integer(PLANE_STRAIN_MPM + 100));
		options.put("plane stress mpm+ps", new Integer(PLANE_STRESS_MPM + 100));
		options.put("axisymmetric mpm+ps", new Integer(AXI_SYM_MPM + 100));
		options.put("3d mpm+ps", new Integer(THREED_MPM + 100));

		// read it
		np = readIntOption(args.get(1), options, "Analysis type");

		// decode particle spin
		if(np > 100)
		{
			np -= 100;
			plusSpin = true;
		}

		// optional element second
		if(args.size() > 2)
		{
			args.remove(1);
			doElement(args);
		}
		else if(lnameEl == NO_ELEMENT)
		{
			if(np > BEGIN_MPM_TYPES)
			{
				if(np == THREED_MPM)
					lnameEl = ElementBase.EIGHT_NODE_ISO_BRICK;
				else
					lnameEl = ElementBase.FOUR_NODE_ISO;
			}
		}
	}

	// ConsistentUnits Length,mass,time
	public void doConsistentUnits(ArrayList<String> args) throws Exception
	{
		if(np >= 0)
			throw new Exception("'ConsistentUnits' command should come before analysis command.");

		if(consistentUnits != null)
			throw new Exception("Only one 'ConsistentUnits' command is allowed.");

		// length, mass, and time (they have to be nil if got here)
		String lengthU = "", massU = "", timeU = "";
		if(args.size() < 2)
		{
			lengthU = "m";
			massU = "kg";
			timeU = "s";
		}
		else
		{
			if(args.size() < 4)
				throw new Exception("ConsistentUnits command must have 0 or 3 arguments");

			lengthU = readStringArg(args.get(1)).toLowerCase();
			massU = readStringArg(args.get(2)).toLowerCase();
			timeU = readStringArg(args.get(3)).toLowerCase();

			// km, m, dm cm, mm um (or microns), or nm
			if(!lengthU.equals("km") && !lengthU.equals("m") && !lengthU.equals("dm") && !lengthU.equals("cm")
					&& !lengthU.equals("mm") && !lengthU.equals("um") && !lengthU.equals("micron")
					&& !lengthU.equals("nm"))
			{
				throw new Exception(
						"Length (#1) in ConsistentUnits command must be km, m, dm cm, mm um (or microns), or nm.");
			}

			// kg, g, mg, or ug
			if(!massU.equals("kg") && !massU.equals("g") && !massU.equals("mg") && !massU.equals("ug"))
			{
				throw new Exception("Mass (#1) in ConsistentUnits command must be kg, g, mg, or ug.");
			}

			// s (or sec), ms (or msec), or us
			if(!timeU.equals("s") && !timeU.equals("sec") && !timeU.equals("ms") && !timeU.equals("msec")
					&& !timeU.equals("us"))
			{
				throw new Exception("Time (#3) in ConsistentUnits command must be s (or sec), ms (or msec), or us.");
			}
		}

		consistentUnits = "    <ConsistentUnits length='" + lengthU + "' mass='" + massU + "' time='" + timeU + "'/>\n";
	}

	// MPMMethd #1,#2
	public void doMPMMethod(ArrayList<String> args) throws Exception
	{
		// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'MPMMethod' has too few parameters:\n" + args);

		// options
		HashMap<String, Integer> options = new HashMap<String, Integer>(10);
		options.put("usf", new Integer(0));
		options.put("usavg", new Integer(2));
		options.put("usavg+", new Integer(2));
		options.put("usavg-", new Integer(2));
		options.put("usl", new Integer(3));
		options.put("usl+", new Integer(3));
		options.put("usl-", new Integer(3));
		mpmMethod = readIntOption(args.get(1), options, "MPM update method");

		// look for skipping extrapolations
		String uoption = readStringArg(args.get(1));
		int nu = uoption.length();
		if(nu>1 && mpmMethod!=0)
		{	if(uoption.charAt(nu - 1) == '-')
				skipExtrap = true;
			else if(uoption.charAt(nu - 1) == '+')
				skipExtrap = false;
			else
				throw new Exception("The USAVG and USL methods must now specify '+' or '-':\n" + args);
		}
		else
			skipExtrap = false;

		// shape functions
		if(args.size() > 2)
		{	String shape = readStringArg(args.get(2)).toLowerCase();
			if(shape.equals("gimp") || shape.equals("ugimp") || shape.equals("1"))
				shapeMethod = "uGIMP";
			else if(shape.equals("cpdi") || shape.equals("lcpdi") || shape.equals("2"))
				shapeMethod = "lCPDI";
			else if(shape.equals("qcpdi") || shape.equals("3"))
				shapeMethod = "qCPDI";
			else if(shape.equals("finite") || shape.equals("4"))
				shapeMethod = "Finite";
			else if(shape.equals("classic") || shape.equals("dirac") || shape.equals("0"))
				shapeMethod = "Dirac";
			else if(shape.equals("b2gimp") || shape.equals("5"))
				shapeMethod = "B2GIMP";
			else if(shape.equals("b2spline") || shape.equals("6"))
				shapeMethod = "B2SPLINE";
			else if(shape.equals("b2cpdi") || shape.equals("7"))
				shapeMethod = "B2CPDI";
			else
				throw new Exception("The selected MPM shape function method was not recognized:\n" + args);

		}
	}

	// Archive #1 (if #2 and #3 give, passed to ArchiveTime command)
	public void doArchive(ArrayList<String> args,boolean makeUnique) throws Exception
	{
		// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' has too few parameters:\n" + args);

		String relPath = readStringArg(args.get(1));
		if(relPath.length() == 0)
			throw new Exception("'" + args.get(0) + "' path has zero length:\n" + args);

		if(makeUnique)
			archiveRoot = "    <ArchiveRoot unique='1'>" + relPath + "</ArchiveRoot>\n";
		else
			archiveRoot = "    <ArchiveRoot>" + relPath + "</ArchiveRoot>\n";

		// optional element second
		if(args.size() > 2)
		{
			args.remove(1);
			doArchiveTime(args);
		}
	}

	// GlobalArchive #1,#2 for type and optional material ID
	public void doGlobalArchive(ArrayList<String> args) throws Exception
	{
		// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' has too few parameters:\n" + args);

		String type = readStringArg(args.get(1));
		if(type.length() == 0)
			throw new Exception("'" + args.get(0) + "' quantity to archive has zero length:\n" + args);

		// optional material ID or tracer particle position
		int matnum = 0;
		if(args.size() > 3)
		{	// look for x, y, z
			double xpos = readDoubleArg(args.get(2));
			double ypos = readDoubleArg(args.get(3));
			String pos = formatDble(xpos)+","+formatDble(ypos);
			if(args.size() > 4)
			{	double zpos = readDoubleArg(args.get(4));
				pos = pos + "," + formatDble(zpos);
			}
			globalArchive = globalArchive + "    <GlobalArchive type='" + type + "' pt='" + pos + "'/>\n";
		}
		else if(args.size() > 2)
		{	// look for material IS
			matnum = mats.getMatID(readStringArg(args.get(2)));
			if(matnum <= 0)
			{	// look for boundaryID
				if(type.equals("reactionx") || type.equals("reactiony") || type.equals("reactionz")
						|| type.equals("reactionR") || type.equals("reactionZ") || type.equals("heatWatts"))
				{
					matnum = readIntArg(args.get(2));
				}
				else
					throw new Exception("'" + args.get(0) + "' command has unknown material ID:\n" + args);
			}

			globalArchive = globalArchive + "    <GlobalArchive type='" + type + "' material='" + matnum + "'/>\n";
		}
		else
			globalArchive = globalArchive + "    <GlobalArchive type='" + type + "'/>\n";
	}

	// ArchiveTime #1,#2,#3 (archive time and optional first archive time and
	// optional max props)
	public void doArchiveTime(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'ArchiveTime' has too few parameters:\n" + args);

		// archive time
		Object aTime = readNumberOrEntityArg(args.get(1), false, 1.e-3);

		// optional max props
		int props = 0;
		if(args.size() > 3)
			props = readIntArg(args.get(3));

		// get archiveTime
		if(props > 0)
			archiveTime = archiveTime + "    <ArchiveTime maxProps='" + props + "'>" + aTime + "</ArchiveTime>\n";
		else
			archiveTime = archiveTime + "    <ArchiveTime>" + aTime + "</ArchiveTime>\n";

		// optional first archive time
		if(args.size() > 2)
		{
			Object firstArchiveTime = readNumberOrEntityArg(args.get(2), false, 1.e-3);
			archiveTime = archiveTime + "    <FirstArchiveTime>" + firstArchiveTime + "</FirstArchiveTime>\n";
		}
	}

	// GlobalArchiveTime #1
	public void doGlobalArchiveTime(ArrayList<String> args) throws Exception
	{ 	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'GlobalArchiveTime' has too few parameters:\n" + args);

		// archive time
		Object aTime = readNumberOrEntityArg(args.get(1), false, 1.e-3);
		globalArchive = globalArchive + "    <GlobalArchiveTime>" + aTime + "</GlobalArchiveTime>\n";
	}

	// TimeStep #1,#2,#3 (time step and optional max time and Courant factor)
	public void doTimeStep(ArrayList<String> args) throws Exception
	{ 	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'TimeStep' has too few parameters:\n" + args);

		// archive time (in sec)
		double aTime = readDoubleArg(args.get(1)) * legacyUnitScaling(1.e-3);
		timeStep = "    <TimeStep>" + formatDble(aTime) + "</TimeStep>\n";

		// max time (in sec)
		if(args.size() > 2)
		{
			aTime = readDoubleArg(args.get(2)) * legacyUnitScaling(1.e-3);
			maxTime = "    <MaxTime>" + formatDble(aTime) + "</MaxTime>\n";
		}

		// Courant time
		if(args.size() > 3)
		{
			aTime = readDoubleArg(args.get(3));
			cflFactor = "    <TimeFactor>" + formatDble(aTime) + "</TimeFactor>\n";
		}
	}

	// CFLFactor #1,<#2>
	public void doCFLFactor(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'CFLFactor' has too few parameters:\n" + args);

		double aCFL = readDoubleArg(args.get(1));
		cflFactor = "    <TimeFactor>" + formatDble(aCFL) + "</TimeFactor>\n";

		// Transport Courant time
		if(args.size() > 2)
		{
			aCFL = readDoubleArg(args.get(2));
			transCflFactor = "    <TransTimeFactor>" + formatDble(aCFL) + "</TransTimeFactor>\n";
		}
	}
	
	// CPDIrcrit #1
	public void doCPDIrcrit(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'CFLFactor' has too few parameters:\n" + args);

		double rcrit = readDoubleArg(args.get(1));
		CPDIrcrit = "    <CPDIrcrit>" + formatDble(rcrit) + "</CPDIrcrit>\n";
	}

	// TimeStep #1,#2,#3 (time step and optional max time and courant factor)
	public void doMaxTime(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'MaximumTime' has too few parameters:\n" + args);

		// archive time
		Object aTime = readNumberOrEntityArg(args.get(1), false, 1.e-3);
		maxTime = "    <MaxTime>" + aTime + "</MaxTime>\n";
	}

	// ToArchive #1,...
	public void doToArchive(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// needs at least one
		if(args.size() < 2)
			throw new Exception("'ToArchive' has too few parameters:\n" + args);

		// first time
		int i;
		if(mpmOrder == null)
		{ // default settings
			mpmOrder = new StringBuffer("iY");
			for(i = 2; i < ReadArchive.ARCH_MAXMPMITEMS; i++)
				mpmOrder.append('N');

			crackOrder = new StringBuffer("iY");
			for(i = 2; i < ReadArchive.ARCH_MAXCRACKITEMS; i++)
				crackOrder.append('N');
		}

		// initial history options
		char historyChar = mpmOrder.charAt(ReadArchive.ARCH_History);
		int history;
		if(historyChar == 'N')
			history = 0x30;
		else if(historyChar == 'Y')
			history = 0x31;
		else
			history = (int) historyChar;
		int origHistory = history;

		historyChar = mpmOrder.charAt(ReadArchive.ARCH_History59);
		int history59;
		if(historyChar == 'N')
			history59 = 0x20;
		else if(historyChar == 'Y')
			history59 = 0x21;
		else
			history59 = (int) historyChar;
		int origHistory59 = history59;

		historyChar = mpmOrder.charAt(ReadArchive.ARCH_History1014);
		int history1014;
		if(historyChar == 'N')
			history1014 = 0x20;
		else if(historyChar == 'Y')
			history1014 = 0x21;
		else
			history1014 = (int) historyChar;
		int origHistory1014 = history1014;

		historyChar = mpmOrder.charAt(ReadArchive.ARCH_History1519);
		int history1519;
		if(historyChar == 'N')
			history1519 = 0x20;
		else if(historyChar == 'Y')
			history1519 = 0x21;
		else
			history1519 = (int) historyChar;
		int origHistory1519 = history1519;

		historyChar = crackOrder.charAt(ReadArchive.ARCH_Traction15);
		int traction15;
		if(historyChar == 'N')
			traction15 = 0x20;
		else if(historyChar == 'Y')
			traction15 = 0x21;
		else
			traction15 = (int) historyChar;
		int origTraction15 = traction15;

		historyChar = crackOrder.charAt(ReadArchive.ARCH_Traction610);
		int traction610;
		if(historyChar == 'N')
			traction610 = 0x20;
		else if(historyChar == 'Y')
			traction610 = 0x21;
		else
			traction610 = (int) historyChar;
		int origTraction610 = traction610;
		
		// set all options in this command
		for(i = 1; i < args.size(); i++)
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
			else if(archive.equals("history5"))
			{	history59 = history59 | 1;
				loc = 0;
			}
			else if(archive.equals("history6"))
			{	history59 = history59 | 2;
				loc = 0;
			}
			else if(archive.equals("history7"))
			{	history59 = history59 | 4;
				loc = 0;
			}
			else if(archive.equals("history8"))
			{	history59 = history59 | 8;
				loc = 0;
			}
			else if(archive.equals("history9"))
			{	history59 = history59 | 16;
				loc = 0;
			}
			else if(archive.equals("history10"))
			{	history1014 = history1014 | 1;
				loc = 0;
			}
			else if(archive.equals("history11"))
			{	history1014 = history1014 | 2;
				loc = 0;
			}
			else if(archive.equals("history12"))
			{	history1014 = history1014 | 4;
				loc = 0;
			}
			else if(archive.equals("history13"))
			{	history1014 = history1014 | 8;
				loc = 0;
			}
			else if(archive.equals("history14"))
			{	history1014 = history1014 | 16;
				loc = 0;
			}
			else if(archive.equals("history15"))
			{	history1519 = history1519 | 1;
				loc = 0;
			}
			else if(archive.equals("history16"))
			{	history1519 = history1519 | 2;
				loc = 0;
			}
			else if(archive.equals("history17"))
			{	history1519 = history1519 | 4;
				loc = 0;
			}
			else if(archive.equals("history18"))
			{	history1519 = history1519 | 8;
				loc = 0;
			}
			else if(archive.equals("history19"))
			{	history1519 = history1519 | 16;
				loc = 0;
			}
			else if(archive.equals("traction1"))
			{	traction15 = traction15 | 1;
				loc = 0;
			}
			else if(archive.equals("traction2"))
			{	traction15 = traction15 | 2;
				loc = 0;
			}
			else if(archive.equals("traction3"))
			{	traction15 = traction15 | 4;
				loc = 0;
			}
			else if(archive.equals("traction4"))
			{	traction15 = traction15 | 8;
				loc = 0;
			}
			else if(archive.equals("traction5"))
			{	traction15 = traction15 | 16;
				loc = 0;
			}
			else if(archive.equals("traction6"))
			{	traction610 = traction610 | 1;
				loc = 0;
			}
			else if(archive.equals("traction7"))
			{	traction610 = traction610 | 2;
				loc = 0;
			}
			else if(archive.equals("traction8"))
			{	traction610 = traction610 | 4;
				loc = 0;
			}
			else if(archive.equals("traction9"))
			{	traction610 = traction610 | 8;
				loc = 0;
			}
			else if(archive.equals("traction10"))
			{	traction610 = traction610 | 16;
				loc = 0;
			}
			else if(archive.equals("heatenergy"))
				loc = ReadArchive.ARCH_HeatEnergy;
			else if(archive.equals("concentration") || archive.equals("porepressure"))
				loc = ReadArchive.ARCH_Concentration;
			else if(archive.equals("czmdisp"))
				cloc = ReadArchive.ARCH_CZMDeltaG;
			else if(archive.equals("elementcrossings"))
				loc = ReadArchive.ARCH_ElementCrossings;
			else if(archive.equals("rotstrain"))
				loc = ReadArchive.ARCH_RotStrain;
			else if(archive.equals("damagenormal"))
				loc = ReadArchive.ARCH_DamageNormal;
			else if(archive.equals("lp"))
				loc = ReadArchive.ARCH_SpinMomentum;
			else if(archive.equals("wp"))
				loc = ReadArchive.ARCH_SpinVelocity;
			else if(archive.equals("size"))
				loc = ReadArchive.ARCH_Size;

			if(loc < 0 && cloc < 0)
				throw new Exception("'" + archive + "' is not a valid archiving option:\n" + args);

			if(loc > 0)
				mpmOrder.replace(loc, loc + 1, "Y");
			if(cloc > 0)
				crackOrder.replace(cloc, cloc + 1, "Y");
		}

		// replace the history character
		if(history != origHistory)
		{	char hchr = history == 0x31 ? 'Y' : (char) history;
			mpmOrder.setCharAt(ReadArchive.ARCH_History, hchr);
		}
		if(history59 != origHistory59)
		{	char hchr = history59 == 0x21 ? 'Y' : (char) history59;
			mpmOrder.setCharAt(ReadArchive.ARCH_History59, hchr);
		}
		if(history1014 != origHistory1014)
		{	char hchr = history1014 == 0x21 ? 'Y' : (char) history1014;
			mpmOrder.setCharAt(ReadArchive.ARCH_History1014, hchr);
		}
		if(history1519 != origHistory1519)
		{	char hchr = history1519 == 0x21 ? 'Y' : (char) history1519;
			mpmOrder.setCharAt(ReadArchive.ARCH_History1519, hchr);
		}
		if(traction15 != origTraction15)
		{	char hchr = traction15 == 0x21 ? 'Y' : (char) traction15;
			crackOrder.setCharAt(ReadArchive.ARCH_Traction15, hchr);
		}
		if(traction610 != origTraction610)
		{	char hchr = traction610 == 0x21 ? 'Y' : (char) traction610;
			crackOrder.setCharAt(ReadArchive.ARCH_Traction610, hchr);
		}
	}

	// Element (element type)
	public void doElement(ArrayList<String> args) throws Exception
	{
		if(np < 0)
			throw new Exception("The 'Element' command must come after the 'Analysis' command");

		if(args.size() < 2)
			throw new Exception("'Element' command has no argument.");

		// options
		HashMap<String, Integer> options = new HashMap<String, Integer>(10);
		options.put("3 node triangle", new Integer(ElementBase.CST));
		options.put("4 node quadrilateral", new Integer(ElementBase.FOUR_NODE_ISO));
		options.put("8 node quadrilateral", new Integer(ElementBase.EIGHT_NODE_ISO));
		options.put("6 node triangle", new Integer(ElementBase.ISO_TRIANGLE));
		options.put("4 node interface", new Integer(ElementBase.LINEAR_INTERFACE));
		options.put("6 node interface", new Integer(ElementBase.QUAD_INTERFACE));
		options.put("8 node brick", new Integer(ElementBase.EIGHT_NODE_ISO_BRICK));
		options.put("9 node lagrange", new Integer(ElementBase.LAGRANGE_2D));

		int oldnameEl = lnameEl;
		lnameEl = readIntOption(args.get(1), options, "Element type");

		if(!ElementBase.CompatibleElements(lnameEl, oldnameEl, np))
		{	throw new Exception("Element type (" + args.get(1) + ") not allowed or\nincompatible with other elements.");
		}

		// pass to FEA areas
		areas.setElementType(lnameEl);
	}

	// XMLData (section),(material ID)
	// Warning: does not check that section is valid
	// Valid are: Header, Mesh, MPMHeader, MaterialPoints, CrackList, Material
	// (must have ID)
	// GridBCs, ParticleBCs, Thermal, Gravity, CustomTasks, end (just append to
	// end)
	public void doXmldata(ArrayList<String> args) throws Exception
	{
		String section = "end";
		if(args.size() > 1)
			section = readStringArg(args.get(1)).toLowerCase();

		// grab text
		String newXML = readVerbatim("endxmldata");

		// check for material section
		if(section.equals("material"))
		{	if(args.size() < 3)
				throw new Exception("XMLData command for a material needs to specify a material ID");
			String matID = readStringArg(args.get(2));
			mats.StartXMLMaterial(matID, newXML);
			return;

		}

		// check GridBCs block and intersperse
		else if(section.equals("gridbcs"))
		{	if(isFEA())
				feaBCs.AddXML(newXML);
			else
				mpmGridBCs.AddXML(newXML);
			return;
		}

		// check MaterialPoints block and intersperse
		else if(section.equals("materialpoints"))
		{	if(regions.isInRegion())
				throw new Exception("XMLData insert in 'materialpoints' must be between regions");
			if(isFEA())
				regions.AddXML(newXML);
			else
				regions.AddXML(newXML);
			return;
		}

		else if(section.equals("cracklist"))
		{
			cracks.appendXMLCrack(newXML);
			return;
		}

		// check previous option
		String currentXML = xmldata.get(section);
		if(currentXML != null)
			newXML = currentXML + newXML;

		// set value
		xmldata.put(section, newXML);
	}

	// Add an entity
	public void doEntity(ArrayList<String> args) throws Exception
	{ // read
			// entity
		// and value
		if(args.size() < 3)
			throw new Exception("'Entity' command has too few arguments:\n" + args);
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
		if(args.size() < 3)
			throw new Exception("'Output' command has too few arguments:\n" + args);

		String quant = readStringArg(args.get(1)).toLowerCase();
		String option = readStringArg(args.get(2)).toLowerCase();

		// add to flags
		// enum { DISPLACEMENT_OUT=0,FORCE_OUT,ELEMSTRESS_OUT,AVGSTRESS_OUT,
		// REACT_OUT,ENERGY_OUT,NUMBER_OUT };
		if(option.equals("selected"))
			option = "C";
		else if(option.equals("no"))
			option = "N";
		else if(option.equals("yes"))
			option = "Y";
		else
			throw new Exception("'Output' option must be 'yes', 'no', or 'selected':\n" + args);

		int offset;
		if(quant.equals("displacements"))
			offset = 0;
		else if(quant.equals("forces"))
			offset = 1;
		else if(quant.equals("elementstresses"))
			offset = 2;
		else if(quant.equals("nodalstresses"))
			offset = 3;
		else if(quant.equals("reactivities"))
			offset = 4;
		else if(quant.equals("energy"))
			offset = 5;
		else
			throw new Exception("Unrecognized 'Output' option:\n" + args);

		// make the change
		outFlags.deleteCharAt(offset);
		outFlags.insert(offset, option);

		// is there an archive time?
		if(args.size() > 3)
		{
			args.remove(1);
			args.remove(1);
			doOutput(args);
		}
	}

	// Temperature (FEA) only
	public void doTemperature(ArrayList<String> args) throws Exception
	{ // Temperature
			// #1
		// which
		// is
		// a
		// function
		if(args.size() < 2)
			throw new Exception("'Temperature' command with too few arguments:\n" + args);

		feaTemp = readStringArg(args.get(1));
	}

	// Stress Free Temperature
	public void doStressFreeTemp(ArrayList<String> args) throws Exception
	{
		if(args.size() < 2)
			throw new Exception("'StressFreeTemp' command with too few arguments:\n" + args);

		stressFreeTemp = readDoubleArg(args.get(1));
	}

	// Damping #1 (number or function),#2 (0 to 1 for PIC),#3 (>0 int for XPIC)
	// also does PDamping command
	public void doDamping(ArrayList<String> args,String dcmd) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'" + dcmd + "' has too few parameters:\n" + args);

		// damping factor (required)
		double damp = 0.;
		String dampcmd;
		Object dampArg = readStringOrDoubleArg(args.get(1));
		if(dampArg.getClass().equals(Double.class))
		{	damp = ((Double) dampArg).doubleValue();
			dampcmd = "    <" + dcmd;
		}
		else
		{	dampcmd = "    <" + dcmd + " function='" + dampArg + "'";
		}

		// PIC fraction (optional)
		double pic = -1.;
		if(args.size() > 2)
		{	pic = readDoubleArg(args.get(2));
			if(pic < 0 || pic > 1)
				throw new Exception("PIC damping in '" + dcmd + "' must be from 0 to 1:\n" + args);
			dampcmd = dampcmd + " PIC='" + pic + "'>" + formatDble(damp) + "</" + dcmd + ">\n";
		}
		else
			dampcmd = dampcmd + ">" + formatDble(damp) + "</" + dcmd + ">\n";

		if(dcmd.equals("Damping"))
			damping = dampcmd;
		else
			pdamping = dampcmd;

		// XPIC (optional)
		if(args.size() > 3)
		{	int xpicOrder = readIntArg(args.get(3));
			if(xpicOrder < 1)
				throw new Exception("XPIC order in '" + dcmd + "' must be integer > 0:\n" + args);
			xpic = "    <XPIC order='" + xpicOrder + "'/>\n";
		}

		// but delete if no PIC
		if(pic <= 0.)
			xpic = null;
	}

	// MultimaterialMode Lumping,Dcheck,Normals,RigidBias
	public void doMultimaterialMode(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// turn it on
		MMLump = -1; // not used by default
		MMNormals = 2; // avggrad default (>=0 means multimaterial mode)
		MMRigidBias = 1.0;
		MMAzimuth = 0.0;
		MMPolar = 0.0;

		// Lumping (used to be Vmin), but current not used even if set
		if(args.size() > 1)
		{	MMLump = readIntArg(args.get(1));
		}

		// Dcheck (but no longer used)
		/*
		 * if(args.size() > 2) { HashMap<String, Integer> options = new HashMap<String, Integer>(4);
		 * options.put("enabled", new Integer(1)); options.put("yes", new Integer(1)); options.put("disabled", new
		 * Integer(0)); options.put("no", new Integer(0)); MMDcheck = readIntOption(args.get(2), options,
		 * "Displacement check option"); }
		 */

		// Normals
		if(args.size() > 3)
		{
			HashMap<String, Integer> options = new HashMap<String, Integer>(4);
			options.put("maxgrad", new Integer(0));
			options.put("maxvol", new Integer(1));
			options.put("avggrad", new Integer(2));
			options.put("owngrad", new Integer(3));
			options.put("specify", new Integer(4));
			options.put("linreg", new Integer(5));
			options.put("logreg", new Integer(6));
			MMNormals = readIntOption(args.get(3), options, "Normals option");
		}

		if(MMNormals == 4)
		{ // polar angles
			if(args.size() > 4)
				MMAzimuth = readDoubleArg(args.get(4));
			if(args.size() > 5)
				MMPolar = readDoubleArg(args.get(5));
		}
		else if(args.size() > 4)
		{ // Rigid Bias
			MMRigidBias = readDoubleArg(args.get(4));
			if(MMRigidBias < 0.)
				MMRigidBias = 0.;
		}
	}

	// ContactPosition Value
	public void doContactPosition(ArrayList<String> args,boolean forCracks) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' has too few parameters:\n" + args);

		double cp = readDoubleArg(args.get(1));
		if(forCracks)
			ContactPositionCracks = "      <ContactPosition>" + formatDble(cp) + "</ContactPosition>\n";
		else
			ContactPosition = "      <ContactPosition>" + formatDble(cp) + "</ContactPosition>\n";
	}

	// ContactMM (LawID),<material ID (only as material prop)>
	// MMMode = 0 (cracks), 1 (multimaterial), 2 (material property), 3 an
	// attribute
	public String doContactLaw(ArrayList<String> args,int MMMode) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' has too few parameters:\n" + args);

		// get law ID
		int lawnum = mats.getMatID(readStringArg(args.get(1)));
		if(lawnum <= 0)
			throw new Exception("'" + args.get(0) + "' as contact law property has unknown material ID:\n" + args);

		// material property needs material ID
		if(MMMode == 2)
		{
			if(args.size() < 3)
				throw new Exception("'" + args.get(0) + "' as material property has too few parameters:\n" + args);

			int matnum = mats.getMatID(readStringArg(args.get(2)));
			if(matnum <= 0)
				throw new Exception("'" + args.get(0) + "' as material property has unknown material ID:\n" + args);

			String cmd = "    <Friction law='" + lawnum + "' mat='" + matnum + "'/>\n";
			return cmd;
		}

		// Friction for cracks or multimaterial mode
		String cmd = "      <Friction law='" + lawnum + "'/>\n";
		if(MMMode == 1)
			FrictionMM = cmd;
		else if(MMMode == 0)
			cracks.setFriction(cmd);
		return null;
	}

	// Deprecated - use doContactLaw instead
	// Friction (number or stick, single (ignore), none),<material ID (only as
	// material prop)>
	// MMMode = 0 (cracks), 1 (multimaterial), 2 (material property), 3 an
	// attribute for CrackList
	public String doFriction(ArrayList<String> args,int MMMode) throws Exception
	{ // MPM  Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' has too few parameters:\n" + args);

		// see if nonnegative number
		double frict = 0.;
		try
		{
			frict = readDoubleArg(args.get(1));
			if(frict < 0)
				throw new Exception("The friction coefficient must be positive:\n" + args);
		}
		catch (Exception e)
		{
			HashMap<String, Integer> options = new HashMap<String, Integer>(4);
			options.put("stick", new Integer(0));
			options.put("single", new Integer(1));
			options.put("ignore", new Integer(1));
			options.put("none", new Integer(2));
			int foption = readIntOption(args.get(1), options, "Friction setting");
			if(foption == 0)
				frict = -5.; // number between -1 and -9
			else if(foption == 1)
				frict = -11.; // number <-10
			else
				frict = 0.0; // frictionless
		}

		// material property needs material ID
		if(MMMode == 2)
		{
			if(args.size() < 3)
				throw new Exception("'" + args.get(0) + "' as material property has too few parameters:\n" + args);

			int matnum = mats.getMatID(readStringArg(args.get(2)));
			if(matnum <= 0)
				throw new Exception("'" + args.get(0) + "' as material property has unknown material ID:\n" + args);

			String cmd = "    <Friction mat='" + matnum + "'>" + formatDble(frict) + "</Friction>\n";
			return cmd;
		}

		// Friction for cracks or multimaterial mode
		String cmd = "      <Friction>" + formatDble(frict) + "</Friction>\n";
		if(MMMode == 1)
			FrictionMM = cmd;
		else if(MMMode == 0)
			cracks.setFriction(cmd);
		else
			return " frict='" + frict + "'";
		return null;
	}

	// Deprecated - use doContactLaw instead (3 uses NewCrack command (frict) ar
	// instead)
	// ImperfectInterface Dt,Dn,<Dnc>
	// MMMode = 0 (cracks), 1 (multimaterial), 2 (material property), 3
	// (CrackInterface commmand)
	public String doImperfectInterface(ArrayList<String> args,int MMMode) throws Exception
	{ // MPMOnly
		requiresMPM(args);

		// read analysis type
		if(args.size() < 3 || (MMMode == 2 && args.size() < 5))
			throw new Exception("'" + args.get(0) + "' has too few parameters:\n" + args);

		// read doubles
		double Dt = readDoubleArg(args.get(1));
		double Dnt = readDoubleArg(args.get(2));
		double Dnc = 0.;
		String cmd;
		if(args.size() > 3)
			Dnc = readDoubleArg(args.get(3));

		if(MMMode == 2)
		{ // get material ID
			int matnum = mats.getMatID(readStringArg(args.get(4)));
			if(matnum <= 0)
				throw new Exception("'" + args.get(0) + "' as material property has unknown material ID:\n" + args);

			cmd = "    <Friction Dt='" + formatDble(Dt) + "' Dnt='" + formatDble(Dnt) + "' Dnc='" + formatDble(Dnc)
					+ "' mat='" + matnum + "'>11</Friction>\n";
		}
		else if(MMMode == 3)
		{
			if(args.size() > 3)
				cracks.setCrackFriction(args,
						" Dt='" + formatDble(Dt) + "' Dnt='" + formatDble(Dnt) + "' Dnc='" + formatDble(Dnc) + "'");
			else
				cracks.setCrackFriction(args, " Dt='" + formatDble(Dt) + "' Dnt='" + formatDble(Dnt) + "'");
			return null;
		}
		else
		{
			if(args.size() > 3)
				cmd = "      <Friction Dt='" + formatDble(Dt) + "' Dnt='" + formatDble(Dnt) + "' Dnc='"
						+ formatDble(Dnc) + "'>11</Friction>\n";
			else
				cmd = "      <Friction Dt='" + formatDble(Dt) + "' Dn='" + formatDble(Dnt) + "'>11</Friction>\n";

			if(MMMode == 1)
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
		if(args.size() < 2)
			throw new Exception("'" + dfbcmd + "' has too few parameters:\n" + args);

		// read gain (required)
		double damp = readDoubleArg(args.get(1));

		// optional target and max alpha
		String target = null;
		double maxdamp = -1.;
		if(args.size() > 2)
			target = readStringArg(args.get(2));
		if(args.size() > 3)
			maxdamp = readDoubleArg(args.get(3));

		String fb;
		if(target == null)
			fb = "    <" + dfbcmd + ">" + formatDble(damp) + "</" + dfbcmd + ">\n";
		else if(maxdamp < 0.)
			fb = "    <" + dfbcmd + " target='" + target + "'>" + formatDble(damp) + "</" + dfbcmd + ">\n";
		else
		{
			fb = "    <" + dfbcmd + " target='" + target + "' max='" + formatDble(maxdamp) + "'>" + formatDble(damp)
					+ "</" + dfbcmd + ">\n";
		}

		if(dfbcmd.equals("FeedbackDamping"))
			fbDamping = fb;
		else
			pfbDamping = fb;
	}

	// ExtrapolateRigid
	public void doExtrapolateRigid(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// activate with no argument
		if(args.size() < 2)
			extrapolateRigid = "    <ExtrapolateRigid/>\n";
		else
		{	String option = readStringArg(args.get(1));
			if(option.toLowerCase().equals("yes"))
				extrapolateRigid = "    <ExtrapolateRigid/>\n";
			else
				extrapolateRigid = null;
		}
	}

	// ExactTraction (yes or no)
	public void doExactTractions(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// activate with no argument
		if(args.size() < 2)
			exactTractions = "    <ExactTractions/>\n";
		else
		{	String option = readStringArg(args.get(1));
			if(option.toLowerCase().equals("yes"))
				exactTractions = "    <ExactTractions/>\n";
			else
				exactTractions = "";
		}
	}

	// LeaveLimit #1 (integer)
	public void doLeaveLimit(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' has too few parameters:\n" + args);

		// leave limit (required)
		int leave = readIntArg(args.get(1));
		leaveLimit = "    <LeaveLimit>" + leave + "</LeaveLimit>\n";
	}

	// DeleteLimit #1 (integer)
	public void doDeleteLimit(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' has too few parameters:\n" + args);

		// delete limit (required)
		int del = readIntArg(args.get(1));
		deleteLimit = "    <DeleteLimit>" + del + "</DeleteLimit>\n";
	}

	// Diffusion #1,<#2> or Poroelasticity #1,<#2>,<#3> or Diffsion,extra
	//		where extra = fracture (or 3), battery (or 4), or conduction (or 5)
	public void doDiffusion(ArrayList<String> args,boolean isMoisture) throws Exception
	{ // MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' has too few parameters:\n" + args);

		// diffusion or poroelasticity on or off (required)
		String option = readStringArg(args.get(1));
		int diffMode=-1;
		if(!isMoisture)
		{	if(option.toLowerCase().equals("yes") || option.equals("2"))
				diffMode = 2;
		}
		else if(option.toLowerCase().equals("yes") || option.toLowerCase().equals("solvent")
				|| option.toLowerCase().equals("moisture") || option.equals("1"))
		{	diffMode = 1;
		}
		else if(option.toLowerCase().equals("fracture") || option.equals("3"))
		{	diffMode = 3;
		}
		else if(option.toLowerCase().equals("battery") || option.equals("4"))
		{	diffMode = 4;
		}
		else if(option.toLowerCase().equals("conduction") || option.equals("4"))
		{	diffMode = 4;
		}
		else if(option.toLowerCase().equals("no"))
		{	diffusion = null;
			return;
		}
		else
			throw new Exception("'" + args.get(0) + "' first parameter no valid:\n" + args);
		
		if(diffusion!=null && diffMode<=2)
			throw new Exception("Cannot have diffusion and poroelastity or two of either:\n" + args);

		// Moisture and Poroelasticy has a references (optionally)
		double ref = 0.;
		if(diffMode<=2 && args.size() > 2)
		{	// require ref>=0, but allow any positive number
			ref = readDoubleArg(args.get(2));
			if(ref < 0.)
				throw new Exception("'" + args.get(0) + "' second parameter must >= 0:\n" + args);
		}

		// viscosity for poroelasticity
		double visc = 1.;
		if(diffMode==2)
		{
			if(args.size() > 3)
			{
				visc = readDoubleArg(args.get(3));
				if(visc <= 0.)
					throw new Exception("'" + args.get(0) + "' third parameter must positive:\n" + args);
			}
		}

		// the command
		if(diffMode==1)
			diffusion = "    <Diffusion reference='" + formatDble(ref) + "'/>\n";
		else if(diffMode==2)
		{
			diffusion = "    <Poroelasticity reference='" + formatDble(ref) + "' viscosity='" + formatDble(visc)
					+ "'/>\n";
		}
		else
			otherDiffusion.append("    <Diffusion style='" + diffMode + "'/>\n");
	}

	// Conduction (yes or no),<adibatic (or mechanical energy) or isothermal or
	// "Crack Tips">
	public void doConduction(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' has too few parameters:\n" + args);

		// yes or no
		boolean hasConduction;
		String option = readStringArg(args.get(1));
		if(option.toLowerCase().equals("no"))
			hasConduction = false;
		else if(option.toLowerCase().equals("yes"))
			hasConduction = true;
		else
			throw new Exception("'" + args.get(0) + "' first parameter must be yes or no:\n" + args);

		// options
		boolean hasCoupling = false;
		boolean hasTips = false;
		boolean hasFriction = false;
		boolean hasCrackFriction = false;
		HashMap<String, Integer> options = new HashMap<String, Integer>(4);
		options.put("adiabatic", new Integer(1));
		options.put("mechanical energy", new Integer(1));
		options.put("isothermal", new Integer(2));
		options.put("crack tips", new Integer(3));
		options.put("friction", new Integer(4));
		options.put("crack friction", new Integer(5));

		// each one
		int arg = 2;
		while (args.size() > arg)
		{
			int opt = readIntOption(args.get(arg), options, "Conduction option");
			if(opt == 1)
				hasCoupling = true;
			else if(opt == 2)
				hasCoupling = false;
			else if(opt == 3)
				hasTips = true;
			else if(opt == 4)
				hasFriction = true;
			else if(opt == 5)
				hasCrackFriction = true;
			arg++;
		}

		// <Conduction/>, <CrackTipHeating/>, <EnergyCoupling/>
		if(hasConduction)
			conduction = "    <Conduction/>\n";
		else
			conduction = "";
		if(hasCoupling)
			conduction = conduction + "    <EnergyCoupling/>\n";
		if(hasTips)
			conduction = conduction + "    <CrackTipHeating/>\n";
		if(hasFriction)
			conduction = conduction + "    <ContactHeating/>\n";
		if(hasCrackFriction)
			conduction = conduction + "    <CrackContactHeating/>\n";
	}

	// CustomTask name
	public void doCustomTask(ArrayList<String> args) throws Exception
	{
		// MPM Only
		requiresMPM(args);

		// read task name
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' is missing a custom task name:\n" + args);
		currentCustomTask = readStringArg(args.get(1));

		// finish last task
		if(customTasks.length() > 0)
		{
			customTasks.append("    </Schedule>\n");
		}

		// start new custom task
		customTasks.append("    <Schedule name='" + currentCustomTask + "'>\n");
	}

	// Parameter #1,<#2>
	public void doParameter(ArrayList<String> args) throws Exception
	{
		// MPM Only
		requiresMPM(args);

		// must be in custom task
		if(customTasks.length() == 0)
			throw new Exception("'" + args.get(0) + "' must can after CustomTask command:\n" + args);

		// read parameter name
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' is missing the parameter name:\n" + args);
		String paramName = readStringArg(args.get(1));

		// handle special commands in some tasks
		if(currentCustomTask.equals("ReverseLoad"))
		{	// quantity combines into a single command
			if(paramName.equals("quantity") && args.size() > 2)
			{
				String value = readStringArg(args.get(2));
				customTasks.append("      <Parameter name='global " + value + "'/>\n");
				return;
			}
			else if(paramName.equals("material") && args.size() > 2)
			{	// material looks for material ID
				int matnum = mats.getMatID(readStringArg(args.get(2)));
				if(matnum <= 0)
				{	// negative is allowed for reaction forces
					matnum = readIntArg(args.get(2));
					if(matnum >= 0)
						throw new Exception(
								"'" + args.get(0) + "' command has unknown material ID or invalid BC ID:\n" + args);
				}
				customTasks.append("      <Parameter name='mat'>" + matnum + "</Parameter>\n");
				return;
			}
			else if(paramName.equals("style") && args.size()>2)
			{	// allow text entery
				HashMap<String, Integer> options = new HashMap<String, Integer>(10);
				options.put("reverse", new Integer(0));
				options.put("hold", new Integer(1));
				options.put("continue", new Integer(2));
				options.put("abort", new Integer(3));
				int styleNum = readIntOption(args.get(2), options, "ReverseLoad style");
				if(styleNum<0 || styleNum>3)
				{	throw new Exception("ReverseLoad 'style' must be 0 to 3");
				}
				customTasks.append("      <Parameter name='style'>" + styleNum + "</Parameter>\n");
				return;
			}
		}

		// volume gradient requires a material ID
		else if(currentCustomTask.equals("VTKArchive"))
		{
			if(paramName.equals("volumegradient"))
			{
				if(args.size() < 3)
					throw new Exception("'volumegradient' quantity requires a material ID:\n" + args);
				int matnum = mats.getMatID(readStringArg(args.get(2)));
				if(matnum <= 0)
					throw new Exception("'volumegradient' quantity requires a valid material ID:\n" + args);
				customTasks.append("      <Parameter name='volumegradient'>" + matnum + "</Parameter>\n");
				return;
			}
		}

		// single parameter
		if(args.size() == 2)
		{
			customTasks.append("      <Parameter name='" + paramName + "'/>\n");
		}
		else
		{
			String value = readStringArg(args.get(2));
			customTasks.append("      <Parameter name='" + paramName + "'>" + value + "</Parameter>\n");
		}
	}

	// Gravity <#1>,<#2>,<#3>
	public void doGravity(ArrayList<String> args) throws Exception
	{ // MPM Only
		requiresMPM(args);

		// defaults are 0,-9806.65,0 otherwise

		// x value
		if(args.size() > 1)
		{
			Object gxarg = readStringOrDoubleArg(args.get(1));
			if(gxarg.getClass().equals(Double.class))
			{
				double gx = ((Double) gxarg).doubleValue();
				gravity = "    <BodyXForce>" + formatDble(gx) + "</BodyXForce>\n";
			}
			else
				gravity = "    <GridBodyXForce>" + (String) gxarg + "</GridBodyXForce>\n";
		}
		else
			gravity = "    <BodyXForce>0.0</BodyXForce>\n";

		// y value
		if(args.size() > 2)
		{
			Object gyarg = readStringOrDoubleArg(args.get(2));
			if(gyarg.getClass().equals(Double.class))
			{
				double gy = ((Double) gyarg).doubleValue();
				gravity = gravity + "    <BodyYForce>" + formatDble(gy) + "</BodyYForce>\n";
			}
			else
				gravity = gravity + "    <GridBodyYForce>" + (String) gyarg + "</GridBodyYForce>\n";
		}
		else if(consistentUnits != null)
			gravity = gravity + "    <BodyYForce>-9.80565</BodyYForce>\n";
		else
			gravity = gravity + "    <BodyYForce>-9805.65</BodyYForce>\n";

		// z value
		if(args.size() > 3)
		{
			Object gzarg = readStringOrDoubleArg(args.get(3));
			if(gzarg.getClass().equals(Double.class))
			{
				double gz = ((Double) gzarg).doubleValue();
				gravity = gravity + "    <BodyZForce>" + formatDble(gz) + "</BodyZForce>\n";
			}
			else
				gravity = gravity + "    <GridBodyZForce>" + (String) gzarg + "</GridBodyZForce>\n";
		}
		else
			gravity = gravity + "    <BodyZForce>0.0</BodyZForce>\n";
	}

	// PtsPerElement #1 (integer)
	public void doPtsPerElement(ArrayList<String> args) throws Exception
	{ // MPM Only
		requiresMPM(args);

		// read analysis type
		if(args.size() < 2)
			throw new Exception("'" + args.get(0) + "' has too few parameters:\n" + args);

		// points per element
		int pts = readIntArg(args.get(1));
		if(pts < 0 || pts > 5 || (pts > 3 && isMPM3D()))
			throw new Exception("'" + args.get(0) + "' has unsupported number of points per element:\n" + args);

		int numCell = pts * pts;
		if(isMPM3D())
			numCell *= pts;
		ptsPerElement = "    <MatlPtsPerElement>" + numCell + "</MatlPtsPerElement>\n";
	}

	// convert @ expression to String
	public String getAtString(String s)
	{	
		// scripts have their own att string method
		if(runningScript) return getScriptAtString(s);
		
		// rest when not running scripts
		
		// split at periods
		String[] atoms = s.substring(1).split("[.]");

		// process them
		boolean passException = false;
		try
		{
			int i = 0;
			while (i < atoms.length)
			{
				String nextAtom = atoms[i];

				// read key points
				if(nextAtom.equals("key"))
				{ // make sure has data
					i += 2;
					return areas.getKeypointProperty(readStringArg(atoms[i - 1]), readStringArg(atoms[i]));
				}
				else if(nextAtom.equals("path"))
				{ // make sure has data
					i += 2;
					return areas.getPathProperty(readStringArg(atoms[i - 1]), readStringArg(atoms[i]));
				}
			}
		}
		catch (Exception e)
		{
			if(passException)
				return "ERROR: " + e.getMessage();
		}

		// if here, than bad expression
		return null;
	}

	// when analysis is done create XML commands
	public String buildXMLCommands()
	{ // start buffer for XML commands
		String more;
		StringBuffer xml = new StringBuffer("<?xml version='1.0'?>\n");
		xml.append("<!DOCTYPE JANFEAInput SYSTEM 'pathto.dtd'");
		if(entities.size() > 0)
		{
			xml.append("\n[\n");
			Set<String> keys = entities.keySet();
			String[] allkeys = new String[entities.size()];
			allkeys = keys.toArray(allkeys);
			int i;
			for(i = 0; i < entities.size(); i++)
			{
				String value = entities.get(allkeys[i]);
				xml.append("  <!ENTITY " + allkeys[i] + " \"" + value + "\">\n");
			}
			xml.append("]>\n");
		}
		else
			xml.append(">\n");
		xml.append("<JANFEAInput version='3'>\n\n");

		// Header
		// -----------------------------------------------------------
		xml.append("  <Header>\n    <Description>\n");
		xml.append("Title: " + title + "\n");
		if(username != null)
			xml.append("User Name: " + username + "\n");
		if(header != null)
			xml.append(header);
		xml.append("    </Description>\n");
		if(consistentUnits != null)
			xml.append(consistentUnits);
		xml.append("    <Analysis>" + np + "</Analysis>\n");
		if(outFlags != null)
			xml.append("    <Output>" + outFlags + "</Output>\n");
		more = xmldata.get("header");
		if(more != null)
			xml.append(more);
		xml.append("  </Header>\n\n");

		// FEA section: Mesh
		// -----------------------------------------------------------
		if(isFEA())
		{	xml.append("  <Mesh>\n" + areas.toXMLString());

			// check added xml
			more = xmldata.get("mesh");
			if(more != null)
				xml.append(more);

			// done
			xml.append("  </Mesh>\n\n");

			// BMPRegion, Body, and Hole blocks
			xml.append(regions.toXMLString());
		}

		// MPM sections: MPMHeader, Mesh, MaterialPoints, CrackList
		// -----------------------------------------------------------
		if(isMPM())
		{ 	// MPM Header
			// -----------------------------------------------------------
			xml.append("  <MPMHeader>\n");

			// MPM method and GIMP
			xml.append("    <MPMMethod>" + mpmMethod + "</MPMMethod>\n");
			if(skipExtrap && mpmMethod != 0)
				xml.append("    <SkipPostExtrapolation/>\n");
			xml.append("    <GIMP type='" + shapeMethod + "'/>\n");
			if(shapeMethod.equals("lCPDI") || shapeMethod.equals("qCPDI") || shapeMethod.equals("B2CPDI"))
				xml.append(CPDIrcrit);
			if(plusSpin)
				xml.append("    <TrackParticleSpin/>\n");
			if(ptsPerElement != null)
				xml.append(ptsPerElement);
			xml.append(timeStep);
			xml.append(cflFactor);
			xml.append(transCflFactor);
			xml.append(maxTime);
			xml.append(archiveRoot);
			xml.append(archiveTime);
			if(mpmOrder == null)
				mpmOrder = new StringBuffer("iYYYYYNYYYNNNNNNNN");
			xml.append("    <MPMArchiveOrder>" + mpmOrder + "</MPMArchiveOrder>\n");
			if(crackOrder == null)
				crackOrder = new StringBuffer("iYNNN");
			xml.append("    <CrackArchiveOrder>" + crackOrder + "</CrackArchiveOrder>\n");

			// global archive
			if(globalArchive.length() > 0)
				xml.append(globalArchive);

			// damping, leave limit, diffusion
			if(damping != null)
				xml.append(damping);
			if(pdamping != null)
				xml.append(pdamping);
			if(xpic != null)
				xml.append(xpic);
			if(fbDamping != null)
				xml.append(fbDamping);
			if(pfbDamping != null)
				xml.append(pfbDamping);
			if(extrapolateRigid != null)
				xml.append(extrapolateRigid);
			if(leaveLimit != null)
				xml.append(leaveLimit);
			if(deleteLimit != null)
				xml.append(deleteLimit);
			if(diffusion != null)
				xml.append(diffusion);
			if(otherDiffusion.length()>0)
				xml.append(otherDiffusion);

			// cracks
			if(MMNormals<0)
			{	// for backward compatibility, single material mode without
				// ContactPositionCracks setting, will use its ContactPosition one
				if(ContactPositionCracks==null) ContactPositionCracks = ContactPosition;
			}
			more = cracks.getSettings(MMNormals, ContactPositionCracks);
			if(more != null)
				xml.append(more);

			// Multimaterial mode <MultiMaterialMode Lumping='1'
			// 							Normals='0' RigidBias='100'>
			// Add Aximuth and Polor for specify
			// Subordinate friction and contact position
			if(MMNormals >= 0)
			{
				xml.append("    <MultiMaterialMode");
				if(MMLump >= 0)
					xml.append(" Lumping='" + formatDble(MMLump) + "'");
				if(MMNormals == 4)
				{
					xml.append(" Normals='" + MMNormals + "' Azimuth='" + formatDble(MMAzimuth) + "' Polar='"
							+ formatDble(MMPolar) + "'>\n");
				}
				else
				{
					xml.append(" Normals='" + MMNormals + "' RigidBias='" + formatDble(MMRigidBias) + "'>\n");
				}
				if(FrictionMM != null)
					xml.append(FrictionMM);
				if(ContactPosition != null)
					xml.append(ContactPosition);
				xml.append("    </MultiMaterialMode>\n");
			}

			// stress free temperature
			if(stressFreeTemp != 0.)
				xml.append("    <StressFreeTemp>" + formatDble(stressFreeTemp) + "</StressFreeTemp>\n");
			
			// exact tractions
			xml.append(exactTractions);

			// check added xml
			more = xmldata.get("mpmheader");
			if(more != null)
				xml.append(more);

			// done
			xml.append("  </MPMHeader>\n\n");

			// MPM Mesh
			// -----------------------------------------------------------
			if(mpmMeshToFile)
				xml.append("  <Mesh output='file'>\n");
			else
				xml.append("  <Mesh>\n");

			xml.append(gridinfo.toXMLString());

			// check added xml
			more = xmldata.get("mesh");
			if(more != null)
				xml.append(more);

			// done
			xml.append("  </Mesh>\n\n");

			// MPM Material Points (XMLData was already added, if any)
			// -----------------------------------------------------------
			xml.append("  <MaterialPoints>\n" + regions.toXMLString());
			xml.append("  </MaterialPoints>\n\n");

			// MPM Cracks
			more = cracks.getCrackList();
			if(more != null)
				xml.append(more);
		}

		// Materials (XMLData was already added, if any)
		// -----------------------------------------------------------
		xml.append(mats.toXMLString());

		// GridBCs
		// -----------------------------------------------------------
		String gridXml = null;
		if(isFEA())
			gridXml = feaBCs.toXMLString();
		else
			gridXml = mpmGridBCs.toXMLString();

		if(gridXml.length() > 0)
			xml.append("  <GridBCs>\n" + gridXml + "  </GridBCs>\n\n");

		// ParticleBCs
		// -----------------------------------------------------------
		if(isMPM())
		{	String partXml = mpmParticleBCs.toXMLString();
			more = xmldata.get("particlebcs");
			if(partXml.length() > 0 || more != null)
			{
				xml.append("  <ParticleBCs>\n" + partXml);
				if(more != null)
					xml.append(more);
				xml.append("  </ParticleBCs>\n\n");
			}
		}

		// FEA: Thermal, MPM: Thermal, Gravity, CustomTasks
		// -----------------------------------------------------------
		if(isFEA())
		{ 	// FEA: Thermal
			// -----------------------------------------------------------
			more = xmldata.get("thermal");
			if(more != null || feaTemp != null || stressFreeTemp != 0.)
			{	xml.append("  <Thermal>\n");

				if(feaTemp != null)
					xml.append("    <Temperature>" + feaTemp + "</Temperature>\n");

				if(stressFreeTemp != 0.)
					xml.append("    <StressFreeTemp>" + formatDble(stressFreeTemp) + "</StressFreeTemp>\n");

				// check added xml
				if(more != null)
					xml.append(more);

				// done
				xml.append("  </Thermal>\n\n");
			}
		}

		else
		{ 	// MPM: Thermal
			// -----------------------------------------------------------
			more = xmldata.get("thermal");
			if(more != null || conduction != null)
			{	xml.append("  <Thermal>\n");

				// conduction
				if(conduction != null)
					xml.append(conduction);

				// check added xml
				if(more != null)
					xml.append(more);

				// done
				xml.append("  </Thermal>\n\n");
			}

			// MPM: Gravity
			// -----------------------------------------------------------
			more = xmldata.get("gravity");
			if(more != null || gravity != null)
			{	xml.append("  <Gravity>\n");

				// check added xml
				if(gravity != null)
					xml.append(gravity);
				if(more != null)
					xml.append(more);

				// done
				xml.append("  </Gravity>\n\n");
			}

			// MPM: CustomTasks
			// -----------------------------------------------------------
			more = xmldata.get("customtasks");
			if(customTasks.length() > 0 || more != null)
			{	xml.append("  <CustomTasks>\n");

				// add tasks
				if(customTasks.length() > 0)
				{	xml.append(customTasks);
					xml.append("    </Schedule>\n");
				}

				// check added xml
				if(more != null)
					xml.append(more);

				// done
				xml.append("  </CustomTasks>\n\n");
			}
		}

		// check added xml
		more = xmldata.get("end");
		if(more != null)
			xml.append(more);

		// convert to string and return
		xml.append("</JANFEAInput>\n");
		return xml.toString();
	}

	// ----------------------------------------------------------------------------
	// Accessors
	// ----------------------------------------------------------------------------

	// empty if no text or if has welcome message
	public boolean isEmptyDocument()
	{
		String cmds = cmdField.getCommands();
		if(cmds.length() == 0)
			return true;
		if(cmds.equals("Welcome to " + NairnFEAMPMViz.appNameReadable + ", " + NairnFEAMPMViz.versionReadable))
			return true;
		return false;
	}

	// call by results when it closes
	public void setLinkedResults(DocViewer someResults)
	{
		linkedResults = someResults;
	}
	public DocViewer getLinkedResults() { return linkedResults; }

	// tell linked results your are closing
	public void windowClosed(WindowEvent e)
	{
		if(linkedResults != null)
			linkedResults.setCommandsWindow(null);
		super.windowClosed(e);
	}

	// called when analysis is done and should link to new results in console
	public DocViewer linkToResults()
	{	if(linkedResults!=null)
		{	// reuse window if same file,but new results
			if(linkedResults.getFile().getPath().equals(soutConsole.getFile().getPath()))
			{	// new results for the same file
				linkedResults.loadNewTextFromFile();
				linkedResults.setVisible(true);
				linkedResults.toFront();
				return linkedResults;
			}
			// close previous linked results
			linkedResults.windowClosing(null);
			linkedResults = null;
		}
	
		// get new reuslts window
		NairnFEAMPMViz.main.openDocument(soutConsole.getFile());
		linkedResults = (DocViewer) NairnFEAMPMViz.main.findDocument(soutConsole.getFile());
		linkedResults.setCommandsWindow(this);
		return linkedResults;
	}

	// type of analysis
	public boolean isFEA()
	{
		return np >= 0 && np < BEGIN_MPM_TYPES;
	}

	public boolean isMPM()
	{
		return np > BEGIN_MPM_TYPES;
	}

	public boolean isMPM3D()
	{
		return np == THREED_MPM;
	}

	// verify FEA or MPM
	public void requiresFEA(ArrayList<String> args) throws Exception
	{
		if(isFEA())
			return;
		if(args != null)
		{
			if(args.size() > 1)
				throw new Exception("The command '" + args.get(0) + "' is only allowed in FEA calculations:\n" + args);
		}
		throw new Exception("Some unknown command is only allowed in FEA calculations.");
	}

	public void requiresMPM(ArrayList<String> args) throws Exception
	{
		if(isMPM())
			return;
		if(args != null)
		{
			if(args.size() > 1)
				throw new Exception("The command '" + args.get(0) + "' is only allowed in MPM calculations:\n" + args);
		}
		throw new Exception("Some unknown command is only allowed in MPM calculations.");
	}

	// return String of Double, String of entity, or Integer object
	// If number, scale by scaleNum if Legacy units
	public Object readNumberOrEntityArg(String text,boolean isInt,double scaleNum) throws Exception
	{
		Object arg = readStringOrDoubleArg(text);
		if(arg.getClass().equals(Double.class))
		{
			if(isInt)
			{
				Integer intarg = new Integer(((Double) arg).intValue() * ((int) scaleNum));
				return intarg;
			}

			// Double newarg = new Double(((Double)arg).doubleValue()*scaleNum);
			double newValue = ((Double) arg).doubleValue();
			return formatDble(newValue * legacyUnitScaling(scaleNum));
		}

		// Strip & and ; if there
		String ent = (String) arg;
		if(ent.startsWith("&") && ent.endsWith(";"))
			ent = ent.substring(1, ent.length() - 1);

		// look for valid entity
		if(entities.get(ent) == null)
			throw new Exception("The argument '" + text + "'\nis neither a number nor a valid entity");
		return "&" + ent + ";";
	}

	// return number if Legacy units or 1 if consistent units
	public double legacyUnitScaling(double scaling)
	{
		if(consistentUnits != null)
			return 1.;
		return scaling;
	}

	// format double and remove trailing zeros from number string (unless has e)
	public String formatDble(double dval)
	{
		String dstr = String.format(Locale.US, "%.8g", dval);
		int dloc = dstr.indexOf('.');
		int eloc = dstr.indexOf('e');
		int lastChar;
		if(eloc<0) eloc = dstr.indexOf('E');
		
		// exponential forms
		if(eloc>=0)
		{	// java might insert zeros after the decimal point
			if(dloc<0)	return dstr;
			
			// find last nonzero
			lastChar = eloc-1;
			while(lastChar>0 && dstr.charAt(lastChar)=='0')
				lastChar--;
			
			return dstr.substring(0,lastChar+1)+dstr.substring(eloc,dstr.length());
		}
		
		// done if no deciment
		if(dloc<0) return dstr;
		
		// remove trailing zeros (if has decimal point)
		lastChar = dstr.length() - 1;
		while (lastChar > 0 && dstr.charAt(lastChar) == '0')
			lastChar--;

		// remove decimal point
		if(dstr.charAt(lastChar) == '.')
			lastChar--;

		return dstr.substring(0, lastChar + 1);
	}

	// override to check commands or analysis running
	public boolean isRunning()
	{
		if(super.isRunning())
			return true;
		if(nfmAnalysis == null)
			return false;
		if(nfmAnalysis.isRunning())
			return true;
		return false;
	}
	
	// check if recent run was aborted
	public boolean wasAborted()
	{	if(nfmAnalysis == null) return false;
		return nfmAnalysis.wasAborted();
	}

	// stop analysis and command interpretation
	public void stopRunning()
	{
		if(nfmAnalysis != null)
		{
			if(nfmAnalysis.isRunning())
				nfmAnalysis.stopRunning();
		}
		if(running)
			running = false;
		if(runningScript)
			runningScript = false;
	}

	// Parameters are
	// 0 = path to output file
	public ArrayList<String> getScriptParams()
	{
		if(!runningScript || scriptParams == null)
			return null;
		return scriptParams;
	}

	// return variable value (or null if none) as string
	public String getVariable(String varName)
	{	String svar = variablesStrs.get(varName);
		return svar;
	}

	// xml text to insert
	public String getXMLData(String blockName)
	{
		return xmldata.get(blockName.toLowerCase());
	}

	// export the output file
	// return true is done or false if error or if cancelled
	public boolean exportOutput(File exportFile,String etitle)
	{
		if(etitle == null)
			etitle = "Export contents of output panel";

		if(exportFile == null)
		{	JFileChooser expChooser = new JFileChooser();
			expChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			NFMVPrefs.setWorkspace(expChooser);
			expChooser.setDialogTitle(etitle);
			int result = expChooser.showSaveDialog(this);
			if(result != JFileChooser.APPROVE_OPTION)
				return false;
			exportFile = expChooser.getSelectedFile();
			
			// check if overwriting
			exportFile = JNUtilities.CheckFileStatus(this, exportFile);
			if(exportFile == null) return false;
		}

		// save output text to exportFile
		try
		{	FileWriter theFile = new FileWriter(exportFile);
			theFile.write(soutConsole.console.getText());
			theFile.flush();
			theFile.close();
		}
		catch (Exception fe)
		{	Toolkit.getDefaultToolkit().beep();
			JNUtilities.showMessage(this, "Error exporting output results: " + fe);
			return false;
		}

		return true;
	}

	// ----------------------------------------------------------------------------
	// Actions as inner classes
	// ----------------------------------------------------------------------------

	// action for stop analysis menu command
	protected class BgAnalysisAction extends JNAction
	{
		private static final long serialVersionUID = 1L;

		public BgAnalysisAction()
		{
			super("Background FEA/MPM Analysis...", KeyEvent.VK_B);
		}

		public void actionPerformed(ActionEvent e)
		{
			runNFMAnalysis(true, NFMAnalysis.FULL_ANALYSIS, null);
		}
	}

	// action for stop analysis menu command
	protected class CheckAnalysisAction extends JNAction
	{
		private static final long serialVersionUID = 1L;

		public CheckAnalysisAction()
		{
			super("Test FEA/MPM Mesh...", KeyEvent.VK_T);
		}

		public void actionPerformed(ActionEvent e)
		{
			runNFMAnalysis(false, NFMAnalysis.RUN_CHECK_MESH, null);
		}
	}

	// action for stop analysis menu command
	protected class InterpretCommandsAction extends JNAction
	{
		private static final long serialVersionUID = 1L;

		public InterpretCommandsAction()
		{
			super("Interpret Commands...", KeyEvent.VK_I);
		}

		public void actionPerformed(ActionEvent e)
		{
			runNFMAnalysis(false, NFMAnalysis.INTERPRET_ONLY, null);
		}
	}

	// action to shaw partner menu command
	protected class ShowPartnerAction extends JNAction
	{
		private static final long serialVersionUID = 1L;

		public ShowPartnerAction()
		{
			super("Show Results");
		}

		public void actionPerformed(ActionEvent e)
		{
			if(linkedResults != null)
				linkedResults.toFront();
			else
				JNApplication.appBeep();
		}
	}

	// action for stop analysis menu command
	protected class ExportXMLAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public ExportXMLAction()
		{	super("Export Output...", KeyEvent.VK_S, true);
		}

		public void actionPerformed(ActionEvent e)
		{	exportOutput(null, null);
		}
	}

	// action for stop analysis menu command
	protected class StopCurrentModeAction extends JNAction
	{
		private static final long serialVersionUID = 1L;

		public StopCurrentModeAction()
		{
			super("Stop Analysis", KeyEvent.VK_PERIOD);
		}

		public void actionPerformed(ActionEvent e)
		{
			stopRunning();
		}
	}

}
