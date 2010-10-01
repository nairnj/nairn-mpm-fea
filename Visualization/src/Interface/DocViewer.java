/*******************************************************************
	DocViewer.java
	NairnFEAMPMViz

	Created by John Nairn on Fri Mar 05 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import java.io.*;

public class DocViewer extends NFMVViewer
{
	static final long serialVersionUID=4L;
	
	//----------------------------------------------------------------------------
	// constants and variables
	//----------------------------------------------------------------------------
	
	public ResultsDocument resDoc=new ResultsDocument();		// one results document in this viewer
	private TextDisplay display;		// split view with sections and text of sections
	public ControlPanel controls;		// controls to set plot options
	private JSplitPane full;			// window split view - TextDisplay on top, ControlPanel on bottom
	private boolean startWithMesh=false;
	
	public MoviePlotWindow movieFrame=null;	// window to plot particles and meshes
	public TimePlotWindow timeFrame=null;	// window to plot data vs time
	public XYPlotWindow xyPlotFrame=null;	// window to plot along grid contour
	
	//----------------------------------------------------------------------------
	// constants and variables
	//----------------------------------------------------------------------------

	public DocViewer(boolean openMesh)
	{   super("DocViewer");
	
		Container content=getContentPane( );
		resDoc.setDocController(this);
		startWithMesh=openMesh;
		
		// more tools
		buttonBar.add(new JLabel("     "));
		addToolBarBtn(buttonBar,"image-x-generic.png","Start Plot",
				"Open plot using the currently selected options.",this);
		addToolBarBtn(buttonBar,"accessories-text-editor.png","ShowPartner",
				"Show associated commands window (if available).",this);

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
		setFramePrefs("Results Window Width",600,"Results Window Height",800);
		setFrameLocAndSize(this);
		
		// notifications
		//JNNotificationCenter.getInstance().addNameAndObjectForTarget("TimeSliderChanged",this,this);
	}
	

	// load the file
	public void loadTheFile(File file,String textContents) throws Exception
	{	display.readFile(file,textContents);
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
		{	if(partnerCtrl!=null)
				partnerCtrl.toFront();
			else
			{	Toolkit.getDefaultToolkit().beep();
				JOptionPane.showMessageDialog(this,"This window does not have an associated commands window.");
			}
		}
		
		else
			super.actionPerformed(e);
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
				if(!currentMovieIsOK(currentClass))
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
				if(timeFrame==null)
				{   timeFrame=new TimePlotWindow(resDoc);
					controls.plotOpened();
				}
				
				try
				{	timeFrame.addPlot(controls);
				}
				catch(Exception tpe)
				{	Toolkit.getDefaultToolkit().beep();
					JOptionPane.showMessageDialog(timeFrame,tpe.toString());
				}
				
				break;
			
			case LoadArchive.MESH2D_PLOT:
				if(xyPlotFrame==null)
				{   xyPlotFrame=new XYPlotWindow(resDoc);
					controls.plotOpened();
				}
				
				try
				{	xyPlotFrame.addPlot(controls);
				}
				catch(Exception tpe)
				{	Toolkit.getDefaultToolkit().beep();
					JOptionPane.showMessageDialog(timeFrame,tpe.toString());
				}
				
				break;
			
			default:
				break;
		}
	}
	
	// see if current movie is same class and requested movie
	private boolean currentMovieIsOK(Class<?> theClass)
	{	if(movieFrame==null) return false;
		if(theClass.equals(movieFrame.getClass())) return true;
		movieFrame.dispose();
		return false;
	}
	
	//----------------------------------------------------------------------------
	// Window events
	//----------------------------------------------------------------------------
	
	public void windowOpened(WindowEvent e)
	{	double splitLoc=NFMVPrefs.prefs.getDouble(NFMVPrefs.ResultsTopSplitKey,
													NFMVPrefs.ResultsTopSplitDef);
		if(splitLoc<0. || splitLoc>1.) splitLoc=NFMVPrefs.ResultsTopSplitDef;
	    display.setDividerLocation(splitLoc);
		splitLoc=NFMVPrefs.prefs.getDouble(NFMVPrefs.ResultsSplitKey,
													NFMVPrefs.ResultsSplitDef);
		if(splitLoc<0. || splitLoc>1.) splitLoc=NFMVPrefs.ResultsSplitDef;
		full.setDividerLocation(splitLoc);
		
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
		double loc=(double)display.getDividerLocation()/
					(double)(display.getMaximumDividerLocation()-display.getMinimumDividerLocation());
		NFMVPrefs.prefs.putDouble(NFMVPrefs.ResultsTopSplitKey,loc);
		loc=(double)full.getDividerLocation()/
					(double)(full.getMaximumDividerLocation()-full.getMinimumDividerLocation());
		NFMVPrefs.prefs.putDouble(NFMVPrefs.ResultsSplitKey,loc);

		super.windowClosing(e);
		if(movieFrame!=null) movieFrame.dispose();
		if(timeFrame!=null) timeFrame.dispose();
		if(xyPlotFrame!=null) xyPlotFrame.dispose();
	}
	
	public void closePlotFrame(NFMVFrame f)
	{	if(f==movieFrame)
			movieFrame=null;
		else if(f==timeFrame)
			timeFrame=null;
		else if(f==xyPlotFrame)
			xyPlotFrame=null;
		controls.plotOpened();
	}
	
	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	public String getFullPath() { return resDoc.fullPath; }
}

