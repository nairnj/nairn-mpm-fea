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
import java.awt.geom.Point2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.NumberFormat;
import java.util.ArrayList;

import geditcom.JNFramework.*;
import geditcom.plot2d.*;

import de.erichseifert.vectorgraphics2d.PDFGraphics2D;

public class TwoDPlotWindow extends JNChildWindow
{
	static final long serialVersionUID=24L;

	// variables
	protected JNPlotView plot2DView;
	
	private JPopupMenu shapePop;
	private JLabel xlabel=new JLabel("x:");
	private JLabel ylabel=new JLabel("y:");
	private JLabel wlabel=new JLabel("w:");
	private JLabel hlabel=new JLabel("h:");

	// initialize
	public TwoDPlotWindow(DocViewer parent)
	{	super(parent,parent.resDoc.name+" Time Plot");
		setFramePrefs("Plot2D Window Width",800,"Plot2D Window Height",600);
		
		// add plot view (entire window for now)
		plot2DView=new JNPlotView(this);
		Container content=getContentPane();
		content.add(plot2DView);
		
		makeMenuBar();
		
		// tool bar icons
		makeToolBar();
		toolBar.setLayout(null);
		JNUnderlineBorder lb = new JNUnderlineBorder(Color.gray);
		toolBar.setBorder(lb);
		toolBar.setPreferredSize(new Dimension(toolBar.getWidth(),44));
		
		int hpos = 6;
		int centerLoc = 4;
		Class<?> baseClass=JNApplication.main.getClass();

		ImageIcon addText=new ImageIcon(baseClass.getResource("Resources/text_tool.png"));
		JButton toolBtn = addToolBarIcon(addText,"AddLabel","Add text label to the plot",plot2DView,new Dimension(32,32));
		toolBtn.setSize(toolBtn.getPreferredSize());
		toolBtn.setLocation(hpos,centerLoc);
		hpos += toolBtn.getWidth()+6;
		
		ImageIcon shapesPNG=new ImageIcon(baseClass.getResource("Resources/Arrow_tool.png"));
		toolBtn = addToolBarIcon(shapesPNG,"","Click and hold to add a shape to the plot",null,new Dimension(32,32));
		toolBtn.setSize(toolBtn.getPreferredSize());
		toolBtn.setLocation(hpos,centerLoc);
		shapePop = new JPopupMenu();
		JMenuItem shapeItem = new JMenuItem("Add Arrow");
		shapeItem.setActionCommand("AddArrow");
		shapeItem.addActionListener(plot2DView);
		shapePop.add(shapeItem);
		shapeItem = new JMenuItem("Add Double Arrow");
		shapeItem.setActionCommand("AddDoubleArrow");
		shapeItem.addActionListener(plot2DView);
		shapePop.add(shapeItem);
		shapeItem = new JMenuItem("Add Rectangle");
		shapeItem.setActionCommand("AddRectangle");
		shapeItem.addActionListener(plot2DView);
		shapePop.add(shapeItem);
		shapeItem = new JMenuItem("Add Oval");
		shapeItem.setActionCommand("AddOval");
		shapeItem.addActionListener(plot2DView);
		shapePop.add(shapeItem);
		shapeItem = new JMenuItem("Add Line");
		shapeItem.setActionCommand("AddLine");
		shapeItem.addActionListener(plot2DView);
		shapePop.add(shapeItem);
		toolBtn.addMouseListener(new MouseAdapter()
		{	public void mousePressed(MouseEvent e)
			{	shapePop.show(e.getComponent(),0, 32);
			}
		});
		hpos += toolBtn.getWidth()+6;

		int infoWidth = 140;
		JPanel pinfo = new JPanel();
		pinfo.setSize(new Dimension(infoWidth,32));
		pinfo.setLocation(hpos,centerLoc);
		toolBar.add(pinfo);
		
		pinfo.setLayout(new GridLayout(2,2));
		Font labelFont=new Font("sanserif",Font.PLAIN,10);
		xlabel.setFont(labelFont);
		pinfo.add(xlabel);
		wlabel.setFont(labelFont);
		pinfo.add(wlabel);
		ylabel.setFont(labelFont);
		pinfo.add(ylabel);
		hlabel.setFont(labelFont);
		pinfo.add(hlabel);
		
		pinfo.setBorder(BorderFactory.createEmptyBorder());
		plot2DView.addMouseMotionListener(new MouseMotionAdapter()
		{	public void mouseMoved(MouseEvent e)
			{	// coordinates
				Point2D.Float pt = plot2DView.pixelsToPt(e.getPoint());
				xlabel.setText("x:"+JNUtilities.formatDouble(pt.x));
				ylabel.setText("y:"+JNUtilities.formatDouble(pt.y));
			}
		});
		
		plot2DView.addComponentListener(new ComponentAdapter()
		{	public void componentResized(ComponentEvent e)
			{	// panal size
				Dimension psize = plot2DView.getSize();
				NumberFormat nf=NumberFormat.getInstance();
				nf.setMaximumFractionDigits(2);
				nf.setMinimumFractionDigits(2);
				double pr = psize.getWidth()/psize.getHeight();
				if(pr<1.) pr = 1./pr;
				wlabel.setText("w:"+nf.format(psize.getWidth()/72.)+ " ("+nf.format(pr)+")");
				hlabel.setText("h:"+nf.format(psize.getHeight()/72.));
			}
		});
		
		finishFrameworkWindow(false);
	}
	
	// make menu bar on launch
	protected void makeMenuBar()
	{	//Create the menu bar
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		if(!JNApplication.isMacLNF())
			menuBar.add(JNDocument.defaultApplicationMenu());			// Application menu
		
		// File menu
		menuBar.add(document.defaultPlotFileMenu(plot2DView,this,this,this));	// File menu
		menuBar.add(JNDocument.defaultPlotEditMenu(plot2DView));	// Edit menu
		
		// Window
		JMenu menu = new JMenu("Window");
		menuBar.add(menu);
		if(JNApplication.isMacLNF())
		{	menu.add(JNApplication.main.getOpenHelpAction());
		}
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
		
		else if(theCmd.equals("ExportPDF"))
		{	PDFExport();
		}

		else if(theCmd.equals("SavePPLTJ"))
		{
			JFileChooser chooser=new JFileChooser();
			//JNFileFilter filter=getFilter(docType);
			int result = chooser.showSaveDialog(this);
			if(result != JFileChooser.APPROVE_OPTION) return;
			
			// check the extension
			ArrayList<String> exts = new ArrayList<String>(1);
			exts.add(".ppltj");
			File tempFile=JNUtilities.CheckFileExtension(this,chooser.getSelectedFile(),
								exts,"an exported PublishPlotJ file.");
			if(tempFile==null) return;
			
			// check if overwriting
			tempFile=JNUtilities.CheckFileStatus(this,tempFile);
			if(tempFile==null) return;
			
			// Get the file
			String ppltj = plot2DView.textRepresentation();
			
			// save if possible
			try
			{	FileOutputStream theFile = new FileOutputStream(tempFile);
				BufferedWriter out = new BufferedWriter(new OutputStreamWriter(
						theFile, "UTF-8"));
				out.write(ppltj);
				out.flush();
				out.close();
				theFile.flush();
				theFile.close();
			}
			catch (Exception fe)
			{	JNApplication.appBeep();
				JOptionPane.showMessageDialog(null, "Error writing plot file document: " + fe);
			}
			
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
