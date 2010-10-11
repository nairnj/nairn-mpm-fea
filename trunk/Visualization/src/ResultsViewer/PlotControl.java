/*******************************************************************
	PlotControl.java
	NairnFEAMPMViz

	Created by John Nairn on Wed Mar 10 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
	
	Base class for control objects in the main window
*******************************************************************/

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;

public class PlotControl extends JPanel
{
	static final long serialVersionUID=15L;
	
	protected DocViewer docCtrl;
	
	// initialize
	PlotControl(int width,int height,DocViewer dc)
	{   super();
		docCtrl=dc;
		setSize(width,height);
		setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));			// of this JPanel
	}
	
	// override to enable and disable all sub controls
	public void setEnabled(boolean setting) { }

	// for connecting lines
	public Rectangle getControlRect()
	{	Point loc=getLocation();
		return new Rectangle(loc.x,loc.y,getWidth(),getHeight());
	}
	
}
