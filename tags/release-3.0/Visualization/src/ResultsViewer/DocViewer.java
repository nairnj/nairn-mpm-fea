/*
 * DocViewer.java
 * NairnFEAMPMViz Application
 * 
 * Created by John Nairn on Fri Mar 05 2004.
 * Copyright (c) 2004 RSAC Software. All rights reserved.
 */

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import geditcom.JNFramework.*;

public class DocViewer extends JNDocument
{
	private static final long serialVersionUID = 1L;

	// instance variables
	
	public ResultsDocument resDoc=new ResultsDocument();		// one results document in this viewer
	private TextDisplay display;		// split view with sections and text of sections
	public ControlPanel controls;		// controls to set plot options
	private JSplitPane full;			// window split view - TextDisplay on top, ControlPanel on bottom
	private boolean startWithMesh=false;
	private CmdViewer commandsWindow=null;

	public DocViewer(boolean openMesh)
	{   super("DocViewer");
	
		setFramePrefs("Results Window Width",600,"Results Window Height",800);
	
		makeMenuBar();
		
		Container content=getContentPane( );
		resDoc.setDocController(this);
		startWithMesh=openMesh;
		
		// tool bar icons
		CmdViewer.addDefaultToolBar(this);
		
		addToolBarBar();
		Class<?> baseClass=JNApplication.main.getClass();
		ImageIcon newMPM=new ImageIcon(baseClass.getResource("Resources/image-x-generic.png"));
		addToolBarIcon(newMPM,"Start Plot","Open plot using the currently selected options.",this);
		ImageIcon showCmds=new ImageIcon(baseClass.getResource("Resources/commands-editor.png"));
		addToolBarIcon(showCmds,"ShowPartner","Show associated commands window (if available).",this);

		// top displays sections and text of selected section
		display=new TextDisplay(resDoc);
		
		// control area
		controls=new ControlPanel(this);
		JScrollPane bottom=new JScrollPane(controls);
		
		// full window split pane
		full=new JSplitPane(JSplitPane.VERTICAL_SPLIT,display,bottom);
		
		// add to content
		content.add(full,BorderLayout.CENTER);
		
		// size and location
		finishFrameworkWindow(true);
	}
	
	// make menu bar on launch
	protected void makeMenuBar()
	{
		// Create the menu bar
		JMenuBar menuBar = new JMenuBar();
		if(!JNApplication.isMacLNF())
			menuBar.add(defaultApplicationMenu());		// Application menu
		menuBar.add(CmdViewer.defaultFileMenu(this));				// File menu
		
		// Edit menu
		JMenu menu = new JMenu("Edit");
		menuBar.add(menu);
		makeMenuItem(menu,"Find...","openFindPanel",JNApplication.main,KeyEvent.VK_F);
		makeMenuItem(menu,"Find Again","findAgain",JNApplication.main,KeyEvent.VK_G);
		
		// Analyze menu
		menu = new JMenu("Analyze");
		menuBar.add(menu);
		makeMenuItem(menu,"Plot Results","Start Plot",this,KeyEvent.VK_R);
		makeMenuItem(menu,"Scale Results...","Scale Results",this);

		// Window
		menu = new JMenu("Window");
		menuBar.add(menu);
		if(JNApplication.isMacLNF())
		{	makeMenuItem(menu,"Help","openHelp",JNApplication.main,0);
		}
		makeMenuItem(menu,"Show Commands","ShowPartner",this);
		menu.addSeparator();
		setWindowMenu(menu);

		// add the menu bar
		setJMenuBar(menuBar);
	}
	
	// read data, throw Exception error
	public void loadTextFromFile(String fileText) throws Exception
	{	display.readFile(getFile(),fileText);
		controls.fileHasLoaded();
		setTitle(resDoc.name);
	}
	
	//----------------------------------------------------------------------------
	// handle application commands
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{	String theCmd=e.getActionCommand();
	
		if(theCmd.equals("Start Plot"))
		{	startNewPlot(controls.getPlotType());
		}
		
		else if(theCmd.equals("ShowPartner"))
		{	if(commandsWindow!=null)
				commandsWindow.toFront();
			else
				JNApplication.appBeep();
		}
		
		else if(theCmd.equals("Scale Results"))
		{	scaleResults();	
		}
		
		else
			super.actionPerformed(e);
	}
	
	public void doFindReplaceAction(int frAction)
	{	display.textPane.doFindReplaceAction(frAction);
	}
	
	// start new plot using current settings
	public void startNewPlot(int newType)
	{	Class<?> currentClass=null;
		switch(newType)
		{	case LoadArchive.PARTICLE_PLOT:
				currentClass=MPMParticlePlotWindow.class;
			case LoadArchive.MESH_PLOT:
				if(newType==LoadArchive.MESH_PLOT)
				{	currentClass = resDoc.isMPMAnalysis() ?
						MeshPlotWindow.class : FEAMeshPlotWindow.class ;
				}
				MoviePlotWindow movieFrame=getMovieFrame();
				if(!currentMovieIsOK(currentClass,movieFrame))
				{	if(newType==LoadArchive.PARTICLE_PLOT)
						movieFrame=(MoviePlotWindow)(new MPMParticlePlotWindow(resDoc,this));
					else if(resDoc.isMPMAnalysis())
						movieFrame=(MoviePlotWindow)(new MeshPlotWindow(resDoc,this));
					else
						movieFrame=(MoviePlotWindow)(new FEAMeshPlotWindow(resDoc,this));
					
					// manually set to current setting in main window controls
					controls.plotOpened();
					int newComponent=controls.getPlotComponent();
					movieFrame.setTitle(resDoc.name+" ("+PlotQuantity.plotName(newComponent)+")");
					movieFrame.setVisible(true);
					movieFrame.toFront();
				}
				else
				{	// existing window
					int newComponent;
					if(newType==controls.getPlotType())
						newComponent=controls.getPlotComponent();
					else
						newComponent=movieFrame.movieControls.getPlotComponent();
					movieFrame.setTitle(resDoc.name+" ("+PlotQuantity.plotName(newComponent)+")");
					movieFrame.setVisible(true);
					movieFrame.toFront();
					movieFrame.beginNewIndexNewComponent(controls.getArchiveIndex(),newComponent);
					if(newType==controls.getPlotType())
						movieFrame.movieControls.syncPlotQuantityMenus();
				}
				break;
			
			case LoadArchive.TIME_PLOT:
				TimePlotWindow timeFrame=getTimeFrame();
				if(timeFrame==null)
				{	timeFrame=new TimePlotWindow(this);
					controls.plotOpened();
				}
				
				try
				{	timeFrame.addPlot(controls);
				}
				catch(Exception tpe)
				{	JNApplication.appBeep();
					JOptionPane.showMessageDialog(timeFrame,tpe.toString());
				}
				
				break;
			
			case LoadArchive.MESH2D_PLOT:
				XYPlotWindow xyPlotFrame=getXYPlotFrame();
				if(xyPlotFrame==null)
				{	xyPlotFrame=new XYPlotWindow(this);
					controls.plotOpened();
				}
				
				try
				{	xyPlotFrame.addPlot(controls);
				}
				catch(Exception tpe)
				{	JNApplication.appBeep();
					JOptionPane.showMessageDialog(xyPlotFrame,tpe.toString());
				}
				
				break;
			
			default:
				break;
		}
	}
	
	// see if current movie is same class and requested movie
	private boolean currentMovieIsOK(Class<?> theClass,MoviePlotWindow movieFrame)
	{	if(movieFrame==null) return false;
		if(theClass.equals(movieFrame.getClass())) return true;
		movieFrame.dispose();
		return false;
	}
	
	// draw check on check mesh run
	public void checkMeshNow()
	{	controls.setCheckMeshItem();
		if(resDoc.isMPMAnalysis())
			startNewPlot(LoadArchive.PARTICLE_PLOT);
		else
			startNewPlot(LoadArchive.MESH_PLOT);
	}
	
	// scale mesh dialog
	public void scaleResults()
	{	// run the dialog box
		ScaleResultsDialog sr=new ScaleResultsDialog(this);
		if(sr.getClickedButton()==JNDialog.CANCEL_BUTTON) return;
		
		boolean changeUnits=false;
		if(!sr.getLengthUnits().equals(resDoc.distU))
		{	changeUnits=true;
			resDoc.distU=sr.getLengthUnits();
			resDoc.lengthScale=sr.getLengthScale();
		}
		if(!sr.getTimeUnits().equals(resDoc.timeU))
		{	changeUnits=true;
			resDoc.timeU=sr.getTimeUnits();
			resDoc.timeScale=sr.getTimeScale();
		}
		
		// change of if new units were selected
		if(changeUnits)
		{	// close plot windows
			MoviePlotWindow mp=getMovieFrame();
			if(mp!=null) mp.windowClosing(null);
			TimePlotWindow tp=getTimeFrame();
			if(tp!=null) tp.windowClosing(null);
			XYPlotWindow xyp=getXYPlotFrame();
			if(xyp!=null) xyp.windowClosing(null);
			
			// reread the data
			try
			{	resDoc.clear(false);
				resDoc.DecodeFileSections(getFile());
				controls.updateTimeDisplay();
			}
			catch (Exception e)
			{	String emsg = e.getMessage();
				if(emsg == null)
					emsg = "Scanner error probably due to corrupted or misformatted data.";
				JOptionPane.showMessageDialog(this,"The analysis failed to read for rescaling:\n   "
												+ emsg);
				windowClosing(null);
			}
		}
	}

	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	public MoviePlotWindow getMovieFrame()
	{	return (MoviePlotWindow)getFirstChildOfType("movieFrame");
	}
	
	public TimePlotWindow getTimeFrame()
	{	return (TimePlotWindow)getFirstChildOfType("timeFrame");
	}
	
	public XYPlotWindow getXYPlotFrame()
	{	return (XYPlotWindow)getFirstChildOfType("xyPlotFrame");
	}
	
	// set command window when opened by simulation here
	public void setCommandsWindow(CmdViewer theCmds) { commandsWindow=theCmds; }

	//----------------------------------------------------------------------------
	// Window events
	//----------------------------------------------------------------------------
	
	public void windowOpened(WindowEvent e)
	{	double splitLoc=NFMVPrefs.prefs.getDouble(NFMVPrefs.ResultsTopSplitKey,
													NFMVPrefs.ResultsTopSplitDef);
		if(splitLoc<10.)
		{	Dimension winSize = JNPreferences.getFrameSize(this);
			splitLoc = (double)winSize.width/2.;
		}
	    display.setDividerLocation((int)splitLoc);
	    
		splitLoc=NFMVPrefs.prefs.getDouble(NFMVPrefs.ResultsSplitKey,
													NFMVPrefs.ResultsSplitDef);
		if(splitLoc<10.)
		{	Dimension winSize = JNPreferences.getFrameSize(this);
			splitLoc = (double)winSize.height/2.;
		}
		full.setDividerLocation((int)splitLoc);
		
		super.windowOpened(e);
		
		if(startWithMesh)
		{	controls.setCheckMeshItem();
			if(resDoc.isMPMAnalysis())
				startNewPlot(LoadArchive.PARTICLE_PLOT);
			else
				startNewPlot(LoadArchive.MESH_PLOT);
		}
	}
			
	public void windowClosing(WindowEvent e)
	{
		double loc=(double)display.getDividerLocation();
		NFMVPrefs.prefs.putDouble(NFMVPrefs.ResultsTopSplitKey,loc);
	    
	    loc=(double)full.getDividerLocation();
		NFMVPrefs.prefs.putDouble(NFMVPrefs.ResultsSplitKey,loc);

		super.windowClosing(e);
		
		// inform the commandsWindow
		if(commandsWindow!=null) commandsWindow.setLinkedResults(null);
	}
	
	public void childWindowDidClose(JNChildWindow child)
	{	super.childWindowDidClose(child);
		//might need to change button name
		controls.plotOpened();
	}

}
