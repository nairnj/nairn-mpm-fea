/*******************************************************************
	PlotLaunch.java
	NairnFEAMPMViz

	Created by John Nairn on Tue Mar 09 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import javax.swing.*;
import java.awt.*;

public class PlotLaunch extends PlotControl
{
	static final long serialVersionUID=16L;
	
	// variables
	private JButton plotBtn=new JButton("Add Plot");		// longest text
	public JProgressBar progress=new JProgressBar(0,10);		// longest text
	
	// initialize
	PlotLaunch(DocViewer dc)
	{   super(ControlPanel.WIDTH,36,dc);
		setLayout(null);
	
		// start plot button
		plotBtn.setAction(docCtrl.getStartPlotAction());
		plotBtn.setActionCommand("Start Plot");
		plotBtn.setSize(plotBtn.getPreferredSize());
		int centerLoc=(36-plotBtn.getHeight())>>1;
		plotBtn.setLocation(10,centerLoc);
		//plotBtn.addActionListener(docCtrl);
		plotBtn.setText("Plot");
		add(plotBtn);
		
		Dimension size=progress.getPreferredSize();
		centerLoc=(36-size.height)>>1;
		progress.setSize(new Dimension(ControlPanel.WIDTH-30-plotBtn.getWidth(),progress.getHeight()));
		progress.setSize(new Dimension(ControlPanel.WIDTH-30-plotBtn.getWidth(),size.height));
		progress.setLocation(20+plotBtn.getWidth(),centerLoc);
		progress.setValue(0);
		progress.setEnabled(false);
		add(progress);
		
		setEnabled(true);
	}

	// call when plot window opened or close
	public void updateTitle(int plotType)
	{
		switch(plotType)
		{	case LoadArchive.PARTICLE_PLOT:
			case LoadArchive.MESH_PLOT:
				if(docCtrl.getMovieFrame()!=null)
					plotBtn.setText("Replot");
				else
					plotBtn.setText("Plot");
				break;
			
			case LoadArchive.TIME_PLOT:
				if(docCtrl.getTimeFrame()!=null)
					plotBtn.setText("Add Plot");
				else
					plotBtn.setText("Plot");
				break;
			
			default:
				if(docCtrl.getXYPlotFrame()!=null)
					plotBtn.setText("Add Plot");
				else
					plotBtn.setText("Plot");
				break;
		}
	}
}
