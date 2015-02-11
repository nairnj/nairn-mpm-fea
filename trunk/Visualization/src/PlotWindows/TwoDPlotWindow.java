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
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

import geditcom.JNFramework.*;
import geditcom.plot2d.*;

import de.erichseifert.vectorgraphics2d.PDFGraphics2D;

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
		
		// Plot File menu
		JMenu menu = new JMenu("File");
		JNDocument.makeMenuItem(menu,"Save PDF Graphics...","ExportPDF",this,KeyEvent.VK_S);
		JNDocument.makeMenuItem(menu,"Save JPG, BMP, or PNG Graphics...","ExportImage",plot2DView,KeyEvent.VK_S,true);
		JNDocument.makeMenuItem(menu,"Export PublishPlotJ or Table...","Export",plot2DView,KeyEvent.VK_E,true);
		JNDocument.makeMenuItem(menu,"Import...","Import",plot2DView,KeyEvent.VK_I,true);
		menu.addSeparator();
		JNDocument.makeMenuItem(menu,"Close","closeWindow",this,KeyEvent.VK_W);
		JNDocument.makeMenuItem(menu,"Print...","Print",plot2DView,KeyEvent.VK_P);
		menuBar.add(menu);
		
		// Edit menu
		menuBar.add(JNDocument.defaultPlotEditMenu(plot2DView));
		
		// Window menu
		menu = new JMenu("Window");
		menuBar.add(menu);
		if(JNApplication.isMacLNF())
		{	menu.add(JNApplication.main.getOpenHelpAction());
		}
		JNDocument.makeMenuItem(menu,"Analysis Results","ShowResults",this,KeyEvent.VK_D);
		menu.addSeparator();
		setWindowMenu(menu);
	}
	
	//----------------------------------------------------------------------------
	// handle some commands
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{	String theCmd=e.getActionCommand();
	
		if(theCmd.equals("ShowResults"))
		{	document.toFront();
		}
		
		if(theCmd.equals("ExportPDF"))
		{	PDFExport();
		}

		else
			super.actionPerformed(e);
	}
	
	public void PDFExport()
	{
		JFileChooser chooser=new JFileChooser();
		//JNFileFilter filter=getFilter(docType);
		int result = chooser.showSaveDialog(this);
		if(result != JFileChooser.APPROVE_OPTION) return;
		
		// check the extension
		ArrayList<String> exts = new ArrayList<String>(1);
		exts.add(".pdf");
		File tempFile=JNUtilities.CheckFileExtension(this,chooser.getSelectedFile(),
							exts,"an exported PDF file.");
		if(tempFile==null) return;
		
		// check if overwriting
		tempFile=JNUtilities.CheckFileStatus(this,tempFile);
		if(tempFile==null) return;
		
		Dimension psize = plot2DView.getSize();
		PDFGraphics2D g = new PDFGraphics2D(0.0, 0.0, (float)psize.width, (float)psize.height);
		JNPlotObject oldSelectedObject=plot2DView.getSelectedObject();
		plot2DView.setSelectedObject(null);
		plot2DView.paintComponent(g);
		plot2DView.setSelectedObject(oldSelectedObject);
		
		FileOutputStream file = null;
		try
		{	file = new FileOutputStream(tempFile);
			file.write(g.getBytes());
			file.close();
		}
		catch(IOException fe)
		{	JOptionPane.showMessageDialog(this,"Error writing PDF plot data: " + fe);
		}
		
	}

}
