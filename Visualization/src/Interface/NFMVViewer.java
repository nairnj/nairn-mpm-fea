/*******************************************************************
	CmdViewer.java
	NairnFEAMPMViz

	Created by John Nairn on Feb 14 2008.
	Copyright (c) 2008 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.event.*;
import java.io.*;
import java.net.URL;
import java.awt.*;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;

public class NFMVViewer extends NFMVFrame
{
	static final long serialVersionUID=31L;
	
	//----------------------------------------------------------------------------
	// constants and variables
	//----------------------------------------------------------------------------
	
	NFMVViewer partnerCtrl=null;			// cmds->viz or viz->cmds
	protected JMenuItem showPartnerMenuItem;
	public JMenu windowMenu;
	protected JPanel buttonBar=new JPanel();
	
	//----------------------------------------------------------------------------
	// constants and variables
	//----------------------------------------------------------------------------

	public NFMVViewer(String viewerName)
	{   super(viewerName);
	
		Container content=getContentPane( );
		
		// button bar
		buttonBar.setLayout(new FlowLayout(FlowLayout.LEFT));
		addToolBarBtn(buttonBar,"document-open.png","OpenFile",
							"Open an existing commands or results file.",NairnFEAMPMViz.appCtrl);
		addToolBarBtn(buttonBar,"document-newfea.png","NewFEA","Create new FEA commands file.",
							NairnFEAMPMViz.appCtrl);
		addToolBarBtn(buttonBar,"document-new.png","NewMPM","Create new MPM commands file.",
							NairnFEAMPMViz.appCtrl);
		addToolBarBtn(buttonBar,"preferences-system.png","Preferences","Edit preferences.",
							NairnFEAMPMViz.appCtrl);
		addToolBarBtn(buttonBar,"help-browser.png","Help","Open the help window.",
							NairnFEAMPMViz.appCtrl);
	
		//buttonBar.setPreferredSize(new Dimension(600,22));
		content.add(buttonBar,BorderLayout.NORTH);
	}
	
	// add button to tool bar
	public void addToolBarBtn(JPanel buttonBar,String btnName,String btnCmd,String toolTip,
								ActionListener target)
	{	
		URL btnImage=Main.class.getResource(btnName);
		JButton toolBtn=new JButton(new ImageIcon(btnImage));
		toolBtn.setToolTipText(toolTip);
		toolBtn.setActionCommand(btnCmd);
		toolBtn.addActionListener(target);
		toolBtn.setContentAreaFilled(false);
		toolBtn.setFocusPainted(false);
		toolBtn.setBorderPainted(false);
		toolBtn.setPreferredSize(new Dimension(24,24));
		buttonBar.add(toolBtn);
	}
	
	// load the file
	public void loadTheFile(File file,String textContents) throws Exception
	{
	}

	// make menu bar on launch
	protected void makeMenuBar()
	{	//Create the menu bar
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		
		// NairnFEAMPMViz menu
		JMenu menu = new JMenu("NairnFEAMPMViz");
		menuBar.add(menu);
		makeMenuItem(menu,"About","About",0,NairnFEAMPMViz.appCtrl);
		makeMenuItem(menu,"Preferences...","Preferences",0,NairnFEAMPMViz.appCtrl);
		makeMenuItem(menu,"Help...","Help",0,NairnFEAMPMViz.appCtrl);
		menu.addSeparator();
		makeMenuItem(menu,NairnFEAMPMViz.appCtrl.quitText(),"Quit",KeyEvent.VK_Q,NairnFEAMPMViz.appCtrl);
		
		// File menu
		menu = new JMenu("File");
		menuBar.add(menu);
		makeMenuItem(menu,"Open...","OpenFile",KeyEvent.VK_O,NairnFEAMPMViz.appCtrl);
		makeMenuItem(menu,"New FEA Commands","NewFEA",KeyEvent.VK_N,NairnFEAMPMViz.appCtrl);
		makeMenuItem(menu,"New MPM Commands","NewMPM",KeyEvent.VK_M,NairnFEAMPMViz.appCtrl);
		menu.addSeparator();
		makeMenuItem(menu,"Close","Close",KeyEvent.VK_W,this);
		
		if(isCommandViewer())
		{	makeMenuItem(menu,"Save","Save",KeyEvent.VK_S,this);
			makeMenuItem(menu,"Save As...","SaveAs",0,this);
		}
		
		// Window menu
		windowMenu = new JMenu("Window");
		menuBar.add(windowMenu);
		if(isCommandViewer())
			showPartnerMenuItem=makeMenuItem(windowMenu,"Visualization","ShowPartner",KeyEvent.VK_D,this);
		else
			showPartnerMenuItem=makeMenuItem(windowMenu,"Input Commands","ShowPartner",KeyEvent.VK_D,this);
		showPartnerMenuItem.setEnabled(false);
		windowMenu.addSeparator();
		
		// add existing windows
		for(int i=0;i<NairnFEAMPMViz.appCtrl.windowMenu.getMenuComponentCount();i++)
		{	NFMVWindowMenuItem mItem=(NFMVWindowMenuItem)(NairnFEAMPMViz.appCtrl.windowMenu.getMenuComponent(i));
			windowMenu.add(new NFMVWindowMenuItem(mItem.getController(),mItem.getController().getTitle()));
		}

		// Analyze menu
		menu = new JMenu("Analyze");
		menuBar.add(menu);
		if(isCommandViewer())
		{	makeMenuItem(menu,"Run FEA/MPM Analysis","RunAnalysis",KeyEvent.VK_R,this);
			makeMenuItem(menu,"Background FEA/MPM Analysis...","BgAnalysis",KeyEvent.VK_B,this);
			makeMenuItem(menu,"Test FEA/MPM Mesh...", "CheckAnalysis",KeyEvent.VK_T,this);
			menu.addSeparator();
			makeMenuItem(menu,"Stop Analysis...","StopAnalysis",KeyEvent.VK_PERIOD,this);
		}
		else
			makeMenuItem(menu,"Plot Results","Start Plot",KeyEvent.VK_R,this);
		
	}
	
	//----------------------------------------------------------------------------
	// Window events
	//----------------------------------------------------------------------------
	
	public void windowClosing(WindowEvent e)
	{	dispose();
		NairnFEAMPMViz.appCtrl.closeViewer(this);
	}
	
	// subclass can block closing by return false
	public boolean confirmClosing() { return true; }
	
	// add window menu item
	public void addWindowMenuItem(NFMVViewer docCtrl)
	{	NFMVWindowMenuItem mItem=null;
		int i;
		for(i=2;i<windowMenu.getMenuComponentCount();i++)
		{	mItem=(NFMVWindowMenuItem)(windowMenu.getMenuComponent(i));
			if(mItem.getController()==docCtrl) break;
			mItem=null;
		}
		
		// create or change the menu item
		if(mItem==null)
		{	mItem=new NFMVWindowMenuItem(docCtrl,docCtrl.getTitle());
			windowMenu.add(mItem);
		}
		else
			mItem.setText(docCtrl.getTitle());
	}
	
	// when window closes, remove from controller list
	public void removeWindowMenuItem(NFMVViewer docCtrl)
	{	// remove from Window menu
		for(int i=2;i<windowMenu.getMenuComponentCount();i++)
		{	NFMVWindowMenuItem mItem=(NFMVWindowMenuItem)(windowMenu.getMenuComponent(i));
			if(mItem.getController()==docCtrl)
			{	windowMenu.remove(mItem);
				break;
			}
		}
	}
	
	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	// must override tnad return full path the the associated file
	public String getFullPath() { return "not_a_viewer"; }
	
	public void setPartner(NFMVViewer partner)
	{	partnerCtrl=partner;
		showPartnerMenuItem.setEnabled(partnerCtrl!=null);
	}
	public void removePartner(NFMVViewer partner)
	{	if(partner==partnerCtrl)
		{	partnerCtrl=null;
			showPartnerMenuItem.setEnabled(false);
		}
	}
	
	// type of viewer
	public boolean isCommandViewer() { return false; }
	
	public boolean isReuseable() { return false; }
}	
