/*
 * CmdViewer
 * NairnFEAMPMViz Application
 * 
 * Created 
 */

import java.awt.event.*;
import java.util.ArrayList;
import java.util.HashMap;

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
	private int lnameEl;
	private HashMap<String,String> xmldata = null;
	public Materials mats = null;
	public Areas areas = null;
	public Regions regions = null;
	public MPMGrid gridinfo = null;
	private FEABCs feaBCs = null;
	private StringBuffer outFlags;
	private int mpmMethod;
	private String shapeMethod;
	private String archiveRoot;
	private String archiveTime;
	private String timeStep;
	private String maxTime;
	private StringBuffer mpmOrder;
	private StringBuffer crackOrder;
	private boolean mpmMeshToFile;
	private String feaTemp;
	private double stressFreeTemp;
	
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
		
		// launch analysis with DTD commands in the field
		nfmAnalysis.runNFMAnalysis(doBackground,runType,cmdField.getCommands(),soutConsole);
	}
	
	// when analysis is done, proceed with calculations (if OKO)
	public void analysisFinished(boolean status)
	{	// give up on error
		if(status==false) return;
		
		// launch analysis with DTD commands in the field
		nfmAnalysis.runNFMAnalysis(useBackground,openMesh,buildXMLCommands(),soutConsole);
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
		mats.initRunSettings();
		areas.initRunSettings();
		regions.initRunSettings();
		feaBCs.initRunSettings();
		gridinfo.initRunSettings();
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
		feaTemp = null;
		stressFreeTemp = 0.;
	}
	
	// handle commands
	public void doCommand(String theCmd,ArrayList<String> args) throws Exception
	{	
		if(mats.isInMaterial())
		{	// commands go to material class when material (keep this option first)
			mats.doMaterialProperty(theCmd,args);
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
		
		else if(theCmd.equals("analysis"))
			doAnalysis(args);
		
		else if(theCmd.equals("mpmmethod"))
			doMPMMethod(args);
		
		else if(theCmd.equals("archive"))
			doArchive(args,false);
		
		else if(theCmd.equals("toarchive"))
			doToArchive(args);
		
		else if(theCmd.equals("archiveunique"))
			doArchive(args,true);
		
		else if(theCmd.equals("archivetime"))
			doArchiveTime(args);
		
		else if(theCmd.equals("timestep"))
			doTimeStep(args);
		
		else if(theCmd.equals("maximumtime"))
			doMaxTime(args);
		
		else if(theCmd.equals("element"))
			doElement(args);
		
		else if(theCmd.equals("xmldata"))
			doXmldata(args);
		
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
		
		else if(theCmd.equals("resequence"))
			feaBCs.Resequence(args);
		
		else if(theCmd.equals("select"))
			feaBCs.AddSelect(args);
		
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
	
	// Archive #1 (if #2 and #3 give, passed to ArchiveTime command
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
	
	// TimeStep #1,#2,#3 (time step and optional max time and courant factor)
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
		
		// check previous option
		String currentXML = xmldata.get(section);
		if(currentXML != null) newXML = currentXML+newXML;
		
		// set value
		xmldata.put(section,newXML);
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

	// when analysis is done create XML commands
	public String buildXMLCommands()
	{	// start buffer for XML commands
		StringBuffer xml = new StringBuffer("<?xml version='1.0'?>\n");
		xml.append("<!DOCTYPE JANFEAInput SYSTEM 'pathto.dtd'>\n");
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
			String more = xmldata.get("Mesh");
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
			xml.append(timeStep);
			xml.append(maxTime);
			xml.append(archiveRoot);
			xml.append(archiveTime);
			if(mpmOrder == null) mpmOrder = new StringBuffer("iYYYYYNYYYNNNNNNNN");
			xml.append("    <MPMArchiveOrder>"+mpmOrder+"</MPMArchiveOrder>\n");
			if(crackOrder == null) crackOrder = new StringBuffer("iYNNN");
			xml.append("    <CrackArchiveOrder>"+crackOrder+"</CrackArchiveOrder>\n");
			
			// check added xml
			String more = xmldata.get("MPMHeader");
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
			
			// check added xml
			more = xmldata.get("MaterialPoints");
			if(more != null) xml.append(more);
			
			// done
			xml.append("  </MaterialPoints>\n\n");
		}
		
		// Materials
		//-----------------------------------------------------------
		xml.append(mats.toXMLString());
		
		// GridBCs
		//-----------------------------------------------------------
		if(isFEA())
		{	xml.append("  <GridBCs>\n"+feaBCs.toXMLString());
		
			// check added xml
			String more = xmldata.get("GridBCs");
			if(more != null) xml.append(more);
		
			// done
			xml.append("  </GridBCs>\n\n");
		}
		
		else
		{	String more = xmldata.get("GridBCs");
		
			if(more!=null)
			{	xml.append("  <GridBCs>\n");
		
				// check added xml
				if(more != null) xml.append(more);
	
				// done
				xml.append("  </GridBCs>\n\n");
			}
		}
		
		// ParticleBCs
		//-----------------------------------------------------------
		if(isMPM())
		{	String more = xmldata.get("ParticleBCs");
		
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
			String more = xmldata.get("Thermal");
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
			String more = xmldata.get("Thermal");
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
