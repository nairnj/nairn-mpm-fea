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
	private String header;
	private int np;
	private int lnameEl;
	private HashMap<String,String> xmldata = null;
	public Materials mats = null;
	public Areas areas = null;
	private FEABCs feaBCs = null;
	
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
		feaBCs = new FEABCs(this);
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
		
		// check if DTD file
		int offset=cmdField.getCommands().indexOf("<?xml ");
		if(offset<0 || offset>10)
		{	// interpret commands
			useBackground = doBackground;
			openMesh = runType;
			soutConsole.clear();
			super.runAnalysis();
			
			// when done, when launch the analysis
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
		header = null;
		np = -1;
		lnameEl = NO_ELEMENT;
		xmldata = new HashMap<String,String>(10);
		mats.initRunSettings();
		areas.initRunSettings();
		feaBCs.initRunSettings();
	}
	
	// handle commands
	public void doCommand(String theCmd,ArrayList<String> args) throws Exception
	{	
		if(theCmd.equals("title"))
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
			header = readVerbatim("endheader");
		
		else if(theCmd.equals("analysis"))
			doAnalysis(args);
		
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
		
		else if(theCmd.equals("endfixline"))
			feaBCs.EndFixLine(args);
		
		else if(theCmd.equals("fixpoint"))
			feaBCs.StartFixPoint(args);
		
		else if(theCmd.equals("endfixpoint"))
			feaBCs.EndFixPoint(args);
		
		else if(theCmd.equals("displacement"))
			feaBCs.AddDisplacement(args);
		
		else if(theCmd.equals("load"))
			feaBCs.AddLoad(args);
		
		else if(theCmd.equals("resequence"))
			feaBCs.Resequence(args);
		
		else if(theCmd.equals("origin"))
			areas.setOrigin(args);
		
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
			mats.StartNewMaterial(matID,newXML);
			return;
			
		}
		
		// check previous option
		String currentXML = xmldata.get(section);
		if(currentXML != null) newXML = currentXML+newXML;
		
		// set value
		xmldata.put(section,newXML);
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
		}
		
		// MPM sections: MPMHeader, Mesh, MaterialPoints, CrackList
		//-----------------------------------------------------------
		if(isMPM())
		{
			
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
		
		// ParticleBCs
		//-----------------------------------------------------------
		
		// FEA: Thermal
		//-----------------------------------------------------------
		
		// MPM: Thermal, Gravity, CustomTasks
		//-----------------------------------------------------------

		
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
