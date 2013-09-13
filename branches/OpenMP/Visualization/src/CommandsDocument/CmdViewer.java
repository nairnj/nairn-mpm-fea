/*
 * CmdViewer
 * NairnFEAMPMViz Application
 * 
 * Created 
 */

import java.awt.event.*;
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
	
	// analysis runner
	private NFMAnalysis nfmAnalysis = null;
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
	private StringBuffer outFlags;
	private int mpmMethod;
	private String shapeMethod;
	private String archiveRoot;
	private String archiveTime;
	private String timeStep;
	private String maxTime;
	private String globalArchive;
	private String damping;
	private String fbDamping;
	private String leaveLimit;
	private String ptsPerElement;
	private StringBuffer mpmOrder;
	private StringBuffer crackOrder;
	private boolean mpmMeshToFile;
	private String feaTemp;
	private double stressFreeTemp;
	private boolean stopCommand;
	private double MMVmin;
	private int MMDcheck;
	private int MMNormals;		// <0 means no multimaterial mode
	private double MMRigidBias;
	private String ContactPosition;
	private String FrictionMM;
	private String Friction;
	
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
		ImageIcon doStop=new ImageIcon(baseClass.getResource("Resources/process-stop.png"));
		addToolBarIcon(doStop,null,"Stop currently running FEA or MPM Analysis.",getStopAnalysisAction());

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
	}
	
	// make menu bar on launch
	protected void makeMenuBar()
	{
		// Create the menu bar
		JMenuBar menuBar = new JMenuBar();
		if(!JNApplication.isMacLNF())
			menuBar.add(defaultApplicationMenu());		// Application menu
		menuBar.add(defaultFileMenu(this));				// File menu
		
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
		menu.add(getStopAnalysisAction());
		
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
	
	// default file menu refering to document target
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
		target.addToolBarIcon(newMPM,"newDocumentMPMCmd","Create a new document.",JNApplication.main);
		ImageIcon newFEA=new ImageIcon(baseClass.getResource("Resources/document-newfea.png"));
		target.addToolBarIcon(newFEA,"newDocumentFEACmd","Create a new document.",JNApplication.main);
	}
	
	//----------------------------------------------------------------------------
	// Interpret commands to XML input commands
	//----------------------------------------------------------------------------
	
	// default command run method passed to custom one
	public void runAnalysis() { runNFMAnalysis(false,NFMAnalysis.FULL_ANALYSIS); }
	
	public void runNFMAnalysis(boolean doBackground,int runType)
	{
		// only allowed if the commands have been saved
		if(getFile()==null)
		{	JNApplication.appBeep();
			JOptionPane.showMessageDialog(this,"The input commands have to be saved to a file before running an analysis.");
			return;
		}
		
		// create once
		if(nfmAnalysis == null)
			nfmAnalysis = new NFMAnalysis(this);
		
		// what is process is current running?
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
			super.runAnalysis();
			
			// when done, will launch the analysis
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
				}
			}
		}
		
		// launch analysis with DTD commands in the field
		nfmAnalysis.runNFMAnalysis(doBackground,runType,cmdField.getCommands(),soutConsole,processors);
	}
	
	// when analysis is done, proceed with calculations (if OKO)
	public void analysisFinished(boolean status)
	{	// give up on error
		if(status==false || stopCommand==true) return;
		
		// launch analysis with DTD commands in the field
		nfmAnalysis.runNFMAnalysis(useBackground,openMesh,buildXMLCommands(),soutConsole,processors);
	}
	
	// initialize variables when intepreting commands
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
		stopCommand = false;
		damping = null;
		fbDamping = null;
		leaveLimit = null;
		ptsPerElement = null;
		MMVmin = 0.0;
		MMDcheck = 0;
		MMNormals = -1;
		MMRigidBias = 1.0;
		ContactPosition = null;
		FrictionMM = null;
		Friction = null;
	}
	
	// handle commands
	public void doCommand(String theCmd,ArrayList<String> args) throws Exception
	{	
		if(mats.isInMaterial())
		{	// commands go to material class when material (keep this option first)
			mats.doMaterialProperty(theCmd,args,this);
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
			feaBCs.AddLoad(args);
		
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
			mpmGridBCs.StartMoveLine(args);
		
		else if(theCmd.equals("movearc"))
			mpmGridBCs.StartMoveLine(args);
		
		else if(theCmd.equals("movebox"))
			mpmGridBCs.StartMoveBox(args);
		
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
			mpmGridBCs.AddConcentration(args);

		else if(theCmd.equals("temperature"))
			mpmGridBCs.AddTemperature(args);

		else if(theCmd.equals("origin"))
			areas.setOrigin(args);
		
		else if(theCmd.equals("material"))
			mats.StartMaterial(args);
		
		else if(theCmd.equals("output"))
			doOutput(args);
		
		else if(theCmd.equals("region"))
			regions.StartRegion(args);
		
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
			regions.AddBox(args);
		
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
		
		else if(theCmd.equals("temperature"))
			doTemperature(args);
		
		else if(theCmd.equals("stressfreetemp"))
			doStressFreeTemp(args);
		
		else if(theCmd.equals("damping"))
			doDamping(args);
		
		else if(theCmd.equals("feedbackdamping"))
			doFBDamping(args);
		
		else if(theCmd.equals("leavelimit"))
			doLeaveLimit(args);
		
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

		else if(theCmd.equals("stop"))
		{	super.doCommand(theCmd,args);
			stopCommand = true;
		}
		
		else
			super.doCommand(theCmd, args);
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
			throw new Exception("'MPMMethod' has too few parameters: "+args);
		
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
				throw new Exception("The selected MPM shape function method was not recognized: "+args);
		
		}
	}
	
	// Archive #1 (if #2 and #3 give, passed to ArchiveTime command)
	public void doArchive(ArrayList<String> args,boolean makeUnique) throws Exception
	{
		// MPM Only
		requiresMPM(args);
		
	    // read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters: "+args);
		
		String relPath = readStringArg(args.get(1));
		if(relPath.length()==0)
			throw new Exception("'"+args.get(0)+"' path has zero length: "+args);
		
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
			throw new Exception("'"+args.get(0)+"' has too few parameters: "+args);
		
		String type = readStringArg(args.get(1));
		if(type.length()==0)
			throw new Exception("'"+args.get(0)+"' quantity to archive has zero length: "+args);
		
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
					throw new Exception("'"+args.get(0)+"' command has unknown material ID: "+args);
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
			throw new Exception("'ArchiveTime' has too few parameters: "+args);
		
		// archive time
		double aTime = readDoubleArg(args.get(1));
		archiveTime = "    <ArchiveTime units='ms'>"+aTime+"</ArchiveTime>\n";
		
		// optional first archive time
		if(args.size()>2)
		{	double firstArchiveTime = readDoubleArg(args.get(2));
			archiveTime = archiveTime + "    <FirstArchiveTime units='ms'>"+firstArchiveTime+"</FirstArchiveTime>\n";
		}
	}
	
	// ArchiveTime #1,#2 (archive time and optional first archive time)
	public void doGlobalArchiveTime(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'GlobalArchiveTime' has too few parameters: "+args);
		
		// archive time
		double aTime = readDoubleArg(args.get(1));
		globalArchive = globalArchive+"    <GlobalArchiveTime units='ms'>"+aTime+"</GlobalArchiveTime>\n";
	}
	
	// TimeStep #1,#2,#3 (time step and optional max time and Courant factor)
	public void doTimeStep(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'TimeStep' has too few parameters: "+args);
		
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
			throw new Exception("'MaximumTime' has too few parameters: "+args);
		
		// archive time
		double aTime = readDoubleArg(args.get(1));
		maxTime = "    <MaxTime units='ms'>"+aTime+"</MaxTime>\n";
	}
	
	// ToArchive #1,...
	public void doToArchive(ArrayList<String> args) throws Exception
	{
	    // MPM Only
	    requiresMPM(args);
	    
	    // needs at least one
	    if(args.size()<2)
	    	throw new Exception("'ToArchive' has too few parameters: "+args);
	    
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
		
		// initial history
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
	        else if(archive.equals("externalwork"))
	        	loc = ReadArchive.ARCH_ExtWork;
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
	        else if(archive.equals("thermalenergy"))
	        	loc = ReadArchive.ARCH_ThermalEnergy;
	        else if(archive.equals("concentration"))
	        	loc = ReadArchive.ARCH_Concentration;
	        else if(archive.equals("energybalance"))
	        	cloc = ReadArchive.ARCH_BalanceResults;
	        else if(archive.equals("elementcrossings"))
	        	loc = ReadArchive.ARCH_ElementCrossings;
	        else if(archive.equals("rotstrain"))
	        	loc = ReadArchive.ARCH_RotStrain;
	        
	        if(loc<0 && cloc<0)
	        	throw new Exception("'"+archive+"' is not a valid archiving option: "+args);
	        
	        if(loc>0)
	        	mpmOrder.replace(loc,loc+1,"Y");
	        if(cloc>0)
	        	crackOrder.replace(cloc,cloc+1,"Y");
	    }
		
		// replace the history character
		if(history != origHistory)
		{	char[] hchr = new char[2];
			hchr[0] = (char)history;
			hchr[1] = 0;
			String hstr = history==0x31 ? "Y" : new String(hchr) ;
			history = ReadArchive.ARCH_History;
			mpmOrder.replace(history,history+1,hstr);
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
		{	throw new Exception("Element type ("+args.get(1)+") not allowed or incompatible with other elements.");
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
			throw new Exception("'Entity' command has too few arguments: "+args);
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
			throw new Exception("'Output' command has too few arguments: "+args);
		
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
			throw new Exception("'Output' option must be 'yes', 'no', or 'selected': "+args);
		
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
	    	throw new Exception("Unrecognized 'Output' option: "+args);
		
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

	// Temperature
	public void doTemperature(ArrayList<String> args) throws Exception
	{
		if(isFEA())
		{	// Temperature #1 which is a function
			if(args.size()<2)
				throw new Exception("'Temperature' command with too few arguments: "+args);
			
			feaTemp = readStringArg(args.get(1));
		}
		
		else
		{	throw new Exception("Temperature command not available for MPM yet: "+args);
		}
	}

	// Stress Free Temperature
	public void doStressFreeTemp(ArrayList<String> args) throws Exception
	{	if(args.size()<2)
			throw new Exception("'StressFreeTemp' command with too few arguments: "+args);
			
		stressFreeTemp = readDoubleArg(args.get(1));
	}

	// Damping #1 (number of function)
	public void doDamping(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'Damping' has too few parameters: "+args);
		
		// damping factor
		String damp = readStringArg(args.get(1));
		damping = "    <Damping>"+damp+"</Damping>\n";
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
			MMNormals = readIntOption(args.get(3),options,"Normals option");
		}
		
		// Rigid Bias
		if(args.size()>4)
		{	MMRigidBias = readDoubleArg(args.get(4));
			if(MMRigidBias<0.) MMRigidBias = 0.;
		}
	}
	
	// ContactPosition Value
	public void doContactPosition(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
		
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters: "+args);
		
		double cp = readDoubleArg(args.get(1));
		ContactPosition = "      <ContactPosition>"+cp+"</ContactPosition>\n";
	}
	
	// Friction (number or stick, single (ignore), none),<material ID (only as material prop)>
	// MMMode = 0 (cracks), 1 (multimaterial), 2 (material property)
	public String doFriction(ArrayList<String> args,int MMMode) throws Exception
	{	// MPM Only
		requiresMPM(args);
		
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters: "+args);
		
		// see if nonnegative number
		double frict = 0.;
		try
		{	frict = readDoubleArg(args.get(1));
			if(frict<0)
				throw new Exception("The friction coefficient must be positive: "+args);
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
				throw new Exception("'"+args.get(0)+"' as material property has too few parameters: "+args);
		
			int matnum = mats.getMatID(readStringArg(args.get(2)));
			if(matnum<=0)
				throw new Exception("'"+args.get(0)+"' as material property has unknown material ID: "+args);
			
			String cmd = "    <Friction mat='"+matnum+"'>"+frict+"</Friction>\n";
			return cmd;
		}
		
		// Friction for cracks or multimaterial mode
		String cmd = "      <Friction>"+frict+"</Friction>\n";
		if(MMMode==1)
			FrictionMM = cmd;
		else
			Friction = cmd;
		return null;
	}
	
	// ImperfectInterface Dt,Dn,<Dnc>
	// MMMode = 0 (cracks), 1 (multimaterial), 2 (material property)
	public String doImperfectInterface(ArrayList<String> args,int MMMode) throws Exception
	{	// MPM Only
		requiresMPM(args);
		
		// read analysis type
		if(args.size()<3 || (MMMode==2 && args.size()<5))
			throw new Exception("'"+args.get(0)+"' has too few parameters: "+args);
		
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
				throw new Exception("'"+args.get(0)+"' as material property has unknown material ID: "+args);
			
			cmd = "    <Friction Dt='"+Dt+"' Dnt='"+Dnt+"' Dnc='"+Dnc+
					"' mat='"+matnum+"'>11</Friction>\n";
		}
		else
		{	if(args.size()>3)
				cmd = "      <Friction Dt='"+Dt+"' Dnt='"+Dnt+"' Dnc='"+Dnc+"'>11</Friction>\n";
			else
				cmd = "      <Friction Dt='"+Dt+"' Dn='"+Dnt+"'>11</Friction>\n";
		
			if(MMMode==1)
				FrictionMM = cmd;
			else
				Friction = cmd;
		}
		
		return cmd;
	}
		
	// FeedbackDamping #1,#2,#3 (number,function,number)
	public void doFBDamping(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'FeedbackDamping' has too few parameters: "+args);
		
		// archive time
		double damp = readDoubleArg(args.get(1));
		String target = null;
		double maxdamp = -1.;
		if(args.size()>2)
			target = readStringArg(args.get(2));
		if(args.size()>3)
			maxdamp = readDoubleArg(args.get(3));
		
		if(target==null)
			fbDamping = "    <FeedbackDamping>"+damp+"</FeedbackDamping>\n";
		else if(maxdamp<0.)
			fbDamping = "    <FeedbackDamping target='"+target+"'>"+damp+"</FeedbackDamping>\n";
		else
		{	fbDamping = "    <FeedbackDamping target='"+target+"' max='"+maxdamp+
								"'>"+damp+"</FeedbackDamping>\n";
		}
	}
	
	// LeaveLimit #1 (integer)
	public void doLeaveLimit(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters: "+args);
		
		// damping factor
		int leave = readIntArg(args.get(1));
		leaveLimit = "    <LeaveLimit>"+leave+"</LeaveLimit>\n";
	}
	
	// PtsPerElement #1 (integer)
	public void doPtsPerElement(ArrayList<String> args) throws Exception
	{	// MPM Only
		requiresMPM(args);
	
		// read analysis type
		if(args.size()<2)
			throw new Exception("'"+args.get(0)+"' has too few parameters: "+args);
		
		// damping factor
		int pts = readIntArg(args.get(1));
		ptsPerElement = "    <MatlPtsPerElement>"+pts+"</MatlPtsPerElement>\n";
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
			
			// damping, leave limit
			if(damping!=null) xml.append(damping);
			if(fbDamping!=null) xml.append(fbDamping);
			if(leaveLimit!=null) xml.append(leaveLimit);
			
			// Multimaterial mode <MultiMaterialMode Vmin='0.0' Dcheck='0' Normals='0' RigidBias='100'>
			// Subordinate friction and contact position
			if(MMNormals>=0)
			{	xml.append("    <MultiMaterialMode Vmin='"+MMVmin+"' Dcheck='"+MMDcheck+
							"' Normals='"+MMNormals+"' RigidBias='"+MMRigidBias+"'>\n");
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
			
			// MPM Material Points
			//-----------------------------------------------------------
			xml.append("  <MaterialPoints>\n"+regions.toXMLString());
			xml.append("  </MaterialPoints>\n\n");
		}
		
		// Materials
		//-----------------------------------------------------------
		xml.append(mats.toXMLString());
		
		// GridBCs
		//-----------------------------------------------------------
		if(isFEA())
			xml.append("  <GridBCs>\n"+feaBCs.toXMLString());
		else
			xml.append("  <GridBCs>\n"+mpmGridBCs.toXMLString());
		
		// done
		xml.append("  </GridBCs>\n\n");
		
		// ParticleBCs
		//-----------------------------------------------------------
		if(isMPM())
		{	more = xmldata.get("ParticleBCs");
		
			if(more!=null)
			{	xml.append("  <ParticleBCs>\n");
	
				// check added xml
				if(more != null) xml.append(more);

				// done
				xml.append("  </ParticleBCs>\n\n");
			}
		}
	
		// FEA: Thermal, MPM: Thermal, Gravity, CustomTasks
		//-----------------------------------------------------------
		if(isFEA())
		{	// FEA: Thermal
			//-----------------------------------------------------------
			more = xmldata.get("Thermal");
			if(more!=null || feaTemp!=null)
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
			if(more!=null)
			{	xml.append("  <Thermal>\n");

				// check added xml
				if(more != null) xml.append(more);

				// done
				xml.append("  </Thermal>\n\n");
			}
			
			// MPM: Gravity
			//-----------------------------------------------------------
			more = xmldata.get("Gravity");
			if(more!=null)
			{	xml.append("  <Gravity>\n");

				// check added xml
				if(more != null) xml.append(more);

				// done
				xml.append("  </Gravity>\n\n");
			}
			
			// MPM: CustomTasks
			//-----------------------------------------------------------
			more = xmldata.get("CustomTasks");
			if(more!=null)
			{	xml.append("  <CustomTasks>\n");

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
		linkedResults=(DocViewer)NairnFEAMPMViz.main.frontDocument();
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
				throw new Exception("The command '"+args.get(0)+"' is only allowed in FEA calculations: "+args);
		}
		throw new Exception("Some unknown command is only allowed in FEA calculations.");
	}
	public void requiresMPM(ArrayList<String> args) throws Exception
	{	if(isMPM()) return;
		if(args != null)
		{	if(args.size()>1)
				throw new Exception("The command '"+args.get(0)+"' is only allowed in MPM calculations: "+args);
		}
		throw new Exception("Some unknown command is only allowed in MPM calculations.");
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
 
		public void actionPerformed(ActionEvent e) { runNFMAnalysis(true,NFMAnalysis.FULL_ANALYSIS); }
	}

	// action for stop analysis menu command
	protected class CheckAnalysisAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public CheckAnalysisAction()
		{	super("Test FEA/MPM Mesh...",KeyEvent.VK_T);
		}
 
		public void actionPerformed(ActionEvent e) { runNFMAnalysis(false,NFMAnalysis.RUN_CHECK_MESH); }
	}
	
	// action for stop analysis menu command
	protected class InterpretCommandsAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public InterpretCommandsAction()
		{	super("Interpret Commands...",KeyEvent.VK_E);
		}
 
		public void actionPerformed(ActionEvent e) { runNFMAnalysis(false,NFMAnalysis.INTERPRET_ONLY); }
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
}
