/*******************************************************************
	CrackSelector.java
	NairnFEAMPMViz

	Created by John Nairn on 8/2/07.
	Copyright (c) 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import javax.swing.*;

public class CrackSelector extends PlotControl
{
	static final long serialVersionUID=3L;

	public static final int CRACK_START=0;
	public static final int CRACK_END=1;

	private ButtonGroup whichEnd=new ButtonGroup();
	private JRadioButton startTip=new JRadioButton("Start");
	private JRadioButton endTip=new JRadioButton("End");
	private JLabel crackLabel=new JLabel("Crack Number:",JLabel.RIGHT);
	private JTextField crackNumberText=new JTextField("1");
	
	// initialize
	CrackSelector(DocViewer dc)
	{   super(ControlPanel.WIDTH,60,dc);
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		//setLayout(new GridLayout(2,2));
		setLayout(gridbag);

		c.fill = GridBagConstraints.BOTH;
		
		// Row 1: crack number text field
		c.insets=new Insets(6,3,1,0);
		c.gridwidth = 1;
		c.weightx = 1.0;
		gridbag.setConstraints(crackLabel, c);
		add(crackLabel);
		
		c.insets.set(6,3,1,3);
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 1.0;
		gridbag.setConstraints(crackNumberText, c);
		add(crackNumberText);
		
		// Row 2: which end radio group
		JLabel blank=new JLabel(" ");
		c.insets.set(0,3,6,0);
		c.gridwidth = 1;
		c.weightx = 1.0;
		add(blank);
		
		c.insets.set(0,3,6,0);
		c.weightx = 2.0;
		gridbag.setConstraints(startTip, c);
		startTip.setSelected(true);
		whichEnd.add(startTip);
		add(startTip);
		
		c.insets.set(0,0,6,3);
		gridbag.setConstraints(endTip, c);
		whichEnd.add(endTip);
		add(endTip);
		
		setEnabled(LoadArchive.NO_PLOT,0);
		
	}
	
	// set current state
	public void setEnabled(int plotType,int plotComponent)
	{
		if(plotType==LoadArchive.TIME_PLOT)
		{	switch(plotComponent)
			{	case PlotQuantity.MPMJ1:
				case PlotQuantity.MPMJ2:
				case PlotQuantity.MPMKI:
				case PlotQuantity.MPMKII:
				case PlotQuantity.MPMLENGTH:
				case PlotQuantity.MPMDEBONDLENGTH:
				case PlotQuantity.MPMCRACKRELEASE:
				case PlotQuantity.MPMCRACKABSORB:
				case PlotQuantity.MPMNORMALCTOD:
				case PlotQuantity.MPMSHEARCTOD:
				case PlotQuantity.MPMDEBONDNCTOD:
				case PlotQuantity.MPMDEBONDSCTOD:
					crackLabel.setEnabled(true);
					crackNumberText.setEnabled(true);
					startTip.setEnabled(true);
					endTip.setEnabled(true);
					break;
				
				default:
					setEnabled(LoadArchive.NO_PLOT,0);
					break;
			}
		}
		else if(plotType==LoadArchive.MESH2D_PLOT)
		{	switch(plotComponent)
			{	case PlotQuantity.MPMNORMALCTOD:
				case PlotQuantity.MPMSHEARCTOD:
				case PlotQuantity.MPMCRACKPROFILE:
				case PlotQuantity.MPMOPENINGFRACTION:
				case PlotQuantity.MPMSHEARFRACTION:
					crackLabel.setEnabled(true);
					crackNumberText.setEnabled(true);
					startTip.setEnabled(false);
					endTip.setEnabled(false);
					break;
				
				default:
					setEnabled(LoadArchive.NO_PLOT,0);
					break;
			}
		}
		else
		{	crackLabel.setEnabled(false);
			crackNumberText.setEnabled(false);
			startTip.setEnabled(false);
			endTip.setEnabled(false);
		}
	}
	
	// get crack number (must be >0)
	public int getCrackNumber() throws Exception
	{	int crack=ControlPanel.readInteger(crackNumberText,"crack number");
		if(crack<1)
			throw new Exception("The crack number must be a positive integer");
		return crack;
	}
	
	// get particle number option for time plots
	public int getCrackTip()
	{	if(startTip.isSelected())
			return CRACK_START;
		else
			return CRACK_END;
	}

}
