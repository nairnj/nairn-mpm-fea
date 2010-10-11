/*
 * TwoDPlotWindow.java
 * NairnFEAMPMViz
 *
 * Created by John Nairn on 1/22/08.
 * Copyright 2007 RSAC Software. All rights reserved.
 */

import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import geditcom.JNFramework.*;
import geditcom.plot2d.*;

public class TwoDPlotWindow extends JNChildWindow
{
	static final long serialVersionUID=24L;

	// variables
	protected JNPlotView plot2DView;
	
	// initialize
	public TwoDPlotWindow(DocViewer parent)
	{	super(parent,parent.resDoc.name+" Time Plot");
		setFramePrefs("Plot2D Window Width",800,"Plot2D Window Height",600);
		
		makeMenuBar();
		
		// add plot view (entire window for now)
		Container content=getContentPane();
		content.add(plot2DView);
		
		finishFrameworkWindow(false);
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
		document.makeMenuItem(menu,"Open...","OpenFile",NairnFEAMPMViz.main,KeyEvent.VK_O);
		menu.addSeparator();
		document.makeMenuItem(menu,"Close","closeWindow",this,KeyEvent.VK_W);
		menu.addSeparator();
		document.makeMenuItem(menu,"Export Graphics...","ExportImage",plot2DView,KeyEvent.VK_G);
		document.makeMenuItem(menu,"Export Data...","Export",plot2DView,KeyEvent.VK_E);
		document.makeMenuItem(menu,"Import...","Import",plot2DView,KeyEvent.VK_I,true);
		document.makeMenuItem(menu,"Print...","Print",plot2DView,KeyEvent.VK_P);
		
		// Edit menu
		menuBar.add(document.defaultPlotEditMenu(plot2DView));
		
		// Window menu
		menu = new JMenu("Window");
		menuBar.add(menu);
		document.makeMenuItem(menu,"Analysis Results","ShowResults",this,KeyEvent.VK_D);
	}
	
	//----------------------------------------------------------------------------
	// handle some commands
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{	String theCmd=e.getActionCommand();
	
		if(theCmd.equals("ShowResults"))
		{	document.toFront();
		}
		
		else
			super.actionPerformed(e);
	}
	
}
