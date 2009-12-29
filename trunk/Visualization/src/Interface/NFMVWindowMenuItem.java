/*******************************************************************
	NFMVWindowMenuItem.java
	NairnFEAMPMViz

	Created by John Nairn on Feb 18 2008.
	Copyright (c) 2008 RSAC Software. All rights reserved.
*******************************************************************/

import javax.swing.*;
import java.awt.event.*;

public class NFMVWindowMenuItem extends JMenuItem implements ActionListener
{
	static final long serialVersionUID=37L;
	
	//-----------------------------------------------------------------
	// Class variables and constants
	//-----------------------------------------------------------------
	
	private NFMVViewer docCtrl;	
  
	//-----------------------------------------------------------------
	// Initialize and respond
	//-----------------------------------------------------------------
	
	public NFMVWindowMenuItem(NFMVViewer dc,String title)
	{	super(title);
		docCtrl=dc;
		setActionCommand("DocToFront");
		addActionListener(this);
	}
	
	public void actionPerformed(ActionEvent e)
	{	docCtrl.toFront();
	}
		
	//-----------------------------------------------------------------
	// Accessors 
	//-----------------------------------------------------------------
	
	public NFMVViewer getController() { return docCtrl; }
}