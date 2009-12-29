/*******************************************************************
	NFMVFrame.java
	NairnFEAMPMViz

	Created by John Nairn on Feb 18 2008.
	Copyright (c) 2008 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class NFMVFrame extends JFrame implements ActionListener, WindowListener, ComponentListener
{
	static final long serialVersionUID=35L;
	
	//----------------------------------------------------------------------------
	// constants and variables
	//----------------------------------------------------------------------------
	
	public String frameWidthKey="Generic Window Width";
	public int frameWidthDef=400;
	public String frameHeightKey="Generic Window Height";
	public int frameHeightDef=400;
	
	//----------------------------------------------------------------------------
	// constants and variables
	//----------------------------------------------------------------------------

	public NFMVFrame(String frameName)
	{   super(frameName);
	
		addWindowListener(this);
		addComponentListener(this);
		
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		
		makeMenuBar();
	}
	
	// set size, offset, and make sure fits on screen. Return final dimensions
	public Dimension setFrameLocAndSize(NFMVFrame f)
	{	Dimension d=NFMVPrefs.getFrameSize(this);
		Rectangle bnds=NairnFEAMPMViz.appCtrl.getViewerBounds(d,getGraphicsConfiguration().getBounds());
		setLocation(bnds.x,bnds.y);
		d.width=bnds.width;
		d.height=bnds.height;
		setSize(d);
		return d;
	}
	
	// override to have a manu bar
	protected void makeMenuBar() {}

	// make a menu  item
	protected void makeMenuItem(JMenu menu,String menuTitle,String menuAction,int mKey,ActionListener target)
	{	JMenuItem menuItem = new JMenuItem(menuTitle);
		menuItem.setActionCommand(menuAction);
		menuItem.addActionListener(target);
		if(mKey!=0)
		{	menuItem.setAccelerator(KeyStroke.getKeyStroke(mKey,
					NairnFEAMPMViz.appCtrl.menuKeyMask()));
		}
		menu.add(menuItem);
	}
	
	//----------------------------------------------------------------------------
	// handle common commands
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{
		String theCmd=e.getActionCommand();
		
		if(theCmd.equals("Close"))
		{	windowClosing(null);
		}
	}
	
	//----------------------------------------------------------------------------
	// Window events
	//----------------------------------------------------------------------------
	
	public void windowClosing(WindowEvent e)
	{	dispose();
	}
	public void windowOpened(WindowEvent e) {}
	public void	windowActivated(WindowEvent e) {}
	public void	windowClosed(WindowEvent e) {}
	public void	windowDeactivated(WindowEvent e) {}
	public void	windowDeiconified(WindowEvent e) {}
	public void	windowIconified(WindowEvent e) {}
	
	//----------------------------------------------------------------------------
	// Component events
	//----------------------------------------------------------------------------

	public void componentResized(ComponentEvent e)
	{	NFMVPrefs.setFrameSize(this,getSize());
	}
	public void	componentHidden(ComponentEvent e) {}
	public void componentMoved(ComponentEvent e) {}
	public void componentShown(ComponentEvent e) {}
	
	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------

	// must call for all subclasses
	public void setFramePrefs(String wKey,int wDef,String hKey,int hDef)
	{	frameWidthKey=wKey;
		frameWidthDef=wDef;
		frameHeightKey=hKey;
		frameHeightDef=hDef;
	}

}
