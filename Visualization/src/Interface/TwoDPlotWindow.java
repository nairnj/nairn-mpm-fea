/*******************************************************************
	TwoDPlotWindow.java
	NairnFEAMPMViz

	Created by John Nairn on 1/22/08.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import javax.swing.*;
import java.awt.event.*;

import geditcom.plot2d.*;

public class TwoDPlotWindow extends NFMVFrame
{
	static final long serialVersionUID=24L;

	// variables
	protected ResultsDocument resDoc;
	protected JNPlotView plot2DView;
	
	// initialize
	public TwoDPlotWindow(ResultsDocument gResDoc)
	{
		super(gResDoc.name+" Time Plot");
		resDoc=gResDoc;
		setFramePrefs("Plot2D Window Width",800,"Plot2D Window Height",600);
		
		// size and location
		setFrameLocAndSize(this);
		
		// add plot view (entire window for now)
		Container content=getContentPane();
		content.add(plot2DView);
	}
	
	// make menu bar on launch
	protected void makeMenuBar()
	{	// need this object created first
		plot2DView=new JNPlotView(this);
		
		//Create the menu bar
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		
		// Plot menu
		JMenu menu = new JMenu("File");
		menuBar.add(menu);
		makeMenuItem(menu,"Open...","OpenFile",KeyEvent.VK_O,NairnFEAMPMViz.appCtrl);
		menu.addSeparator();
		makeMenuItem(menu,"Close","Close",KeyEvent.VK_W,this);
		menu.addSeparator();
		makeMenuItem(menu,"Export Graphics...","ExportImage",KeyEvent.VK_G,plot2DView);
		makeMenuItem(menu,"Export Data...","Export",KeyEvent.VK_E,plot2DView);
		makeMenuItem(menu,"Import...","Import",KeyEvent.VK_I,plot2DView);
		makeMenuItem(menu,"Print...","Print",KeyEvent.VK_P,plot2DView);
		
		// Edit menu
		menu = new JMenu("Edit");
		menuBar.add(menu);
		makeMenuItem(menu,"Cut","Cut",KeyEvent.VK_X,plot2DView);
		makeMenuItem(menu,"Copy","Copy",KeyEvent.VK_C,plot2DView);
		makeMenuItem(menu,"Paste","Paste",KeyEvent.VK_V,plot2DView);
		makeMenuItem(menu,"Delete","Delete",0,plot2DView);
		
		// Window menu
		menu = new JMenu("Window");
		menuBar.add(menu);
		makeMenuItem(menu,"Analysis Results","ShowResults",KeyEvent.VK_D,this);
	}
	
	//----------------------------------------------------------------------------
	// handle some commands
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{	String theCmd=e.getActionCommand();
	
		if(theCmd.equals("ShowResults"))
		{	resDoc.getDocController().toFront();
		}
		
		else
			super.actionPerformed(e);
	}
	
	//----------------------------------------------------------------------------
	// Window events
	//----------------------------------------------------------------------------
	
	public void	windowClosed(WindowEvent e)
	{	resDoc.getDocController().closePlotFrame(this);
	}

	
}
