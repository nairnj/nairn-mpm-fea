/*******************************************************************
	TimeSelector.java
	NairnFEAMPMViz
	
	This same class is used for time selecter in both the results
	control panel and in the movie window. Its location is determined
	by whether movieCtrl==null (in results control panel) or
	movieCtrl!=null (in movie plot window). docCtrl is always the
	document controller.

	Created by John Nairn on Tue Mar 09 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import javax.swing.*;
import java.util.*;
import javax.swing.event.*;
import geditcom.JNFramework.*;

public class TimeSelector extends PlotControl
{
	static final long serialVersionUID=23L;
	
	// variables
	public JSlider select=new JSlider(0,50,0);
	private JLabel timeSelected=new JLabel();
	private ArrayList<Double> archiveTimes;
	private int value=-1;
	private MoviePlotWindow movieCtrl;
	private DocViewer docCtrl;
	private String prefix;
	
	// initialize
	TimeSelector(MoviePlotWindow mpc,DocViewer dc)
	{   super(ControlPanel.WIDTH,64,dc);
		setLayout(new GridLayout(2,1));

		// slider bar
		docCtrl=dc;
		movieCtrl=mpc;
		if(movieCtrl==null)
		{	select.setMajorTickSpacing(10);
			select.setPaintTicks(true);
			select.setToolTipText("Select the initial time results to be plotted");
		}
		else
			select.setToolTipText("Select different time results to be plotted");
		add(select);
		
		// the displayed value in a text field
		timeSelected.setHorizontalAlignment(JLabel.CENTER);
		add(timeSelected);
		
		// initialize
		if(movieCtrl!=null)
		{	setEnabled(LoadArchive.PARTICLE_PLOT);
			setBorder(BorderFactory.createEmptyBorder());		// of the PlotControl JPanel
			setBackground(Color.lightGray);
			select.setBackground(Color.lightGray);
			timeSelected.setFont(new Font("sanserif",Font.PLAIN,10));
			timeSelected.setText(archiveTimes.get(0) + " ms");
			prefix="";
		}
		else
		{	setEnabled(LoadArchive.NO_PLOT);
			timeSelected.setText("Time: 0.0 ms");
			prefix="Time: ";
		}
		select.setFocusable(false);
			
		// the slider listener
		select.addChangeListener(new ChangeListener()
		{   public void stateChanged(ChangeEvent e)
			{	timeSelected.setText(prefix + archiveTimes.get(select.getValue()) + " ms");
				if(!select.getValueIsAdjusting())
				{	if(value!=select.getValue())
					{	value=select.getValue();
						JNNotificationCenter.getInstance().postNotification("TimeSliderChanged",docCtrl,select);
					}
				}
			}
		});
	}
	
	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------
	
	// called when new file loaded
	public void setEnabled(int selected)
	{	boolean opened = (selected!=LoadArchive.TIME_PLOT) && (selected!=LoadArchive.NO_PLOT) && (docCtrl.resDoc.isMPMAnalysis());
	    if(opened)
		{   archiveTimes=docCtrl.resDoc.archiveTimes;
			select.setEnabled(true);
			select.setMaximum(archiveTimes.size()-1);
			int ticks=(archiveTimes.size()-1)/10;
			if(ticks<1) ticks=1;
			select.setMajorTickSpacing(ticks);
			//select.setValue(0);
			timeSelected.setEnabled(true);
		}
		else
		{   select.setEnabled(false);
			timeSelected.setEnabled(false);
		}
	}
	
	// return selected archive by index
	public int getArchiveIndex() { return select.getValue(); }
	
	// change value and return true or false if this change caused a reload
	// if it did not reload, will need to reload data if component has changed too
	public boolean setArchiveIndex(int newIndex)
	{	if(newIndex!=select.getValue())
		{	select.setValue(newIndex);
			return true;
		}
		else if(value<0)
		{	// only here on first call, so post notification to get it read if in movie controller
			value=select.getValue();
			JNNotificationCenter.getInstance().postNotification("TimeSliderChanged",docCtrl,select);
			return true;
		}
		return false;
	}
	
	// increment archive index if possible
	public boolean incrementArchiveIndex()
	{	int currentIndex=select.getValue();
		if(currentIndex==select.getMaximum()) return false;
		select.setValue(currentIndex+1);
		return true;
	}
}
