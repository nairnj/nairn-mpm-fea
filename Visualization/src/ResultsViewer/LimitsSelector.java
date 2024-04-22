/*******************************************************************
	LimitsSelector.java
	NairnFEAMPMViz

	Created by John Nairn on 8/2/07.
	Copyright (c) 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Point2D;

import javax.swing.*;

import geditcom.JNFramework.JNUtilities;

public class LimitsSelector extends PlotControl
{
	static final long serialVersionUID=27L;
	
	public static final int DYNAMIC_LIMITS=0;
	public static final int FIXED_LIMITS=1;
	public static final int GLOBAL_LIMITS=2;
	private int currentLimits=DYNAMIC_LIMITS;

	private ButtonGroup limitType=new ButtonGroup();
	private JRadioButton dynamicLimits=new JRadioButton("Dynamic Limits");
	private JRadioButton globalLimits=new JRadioButton("Global Dynamic Limits");
	private JRadioButton fixedLimits=new JRadioButton("Fixed Limits");
	private JLabel minLabel=new JLabel("Min:");
	private JTextField minText=new JTextField("0");
	private JLabel maxLabel=new JLabel("Max:");
	private JTextField maxText=new JTextField("1");
	
	private static int prevLimits = 0;
	private static String prevMinText = "0";
	private static String prevMaxText = "1";
	
	// initialize
	LimitsSelector(DocViewer dc)
	{   super(ControlPanel.WIDTH,104,dc);
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		//setLayout(new GridLayout(2,2));
		setLayout(gridbag);
		
		minText.setText(prevMinText);
		maxText.setText(prevMaxText);

		c.fill = GridBagConstraints.BOTH;
		
		// Row 1: First radio button
		c.insets=new Insets(0,0,0,0);
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 1.0;
		gridbag.setConstraints(dynamicLimits, c);
		if(prevLimits==0) dynamicLimits.setSelected(true);
		limitType.add(dynamicLimits);
		add(dynamicLimits);
		dynamicLimits.addActionListener(new ActionListener()
		{	public void actionPerformed(ActionEvent e)
			{	minText.setEnabled(false);
				maxText.setEnabled(false);
				minLabel.setEnabled(false);
				maxLabel.setEnabled(false);
				currentLimits=DYNAMIC_LIMITS;
			}
		});
		
		// Row 2: Second radio button
		gridbag.setConstraints(globalLimits, c);
		if(prevLimits==1) globalLimits.setSelected(true);
		limitType.add(globalLimits);
		add(globalLimits);
		globalLimits.addActionListener(new ActionListener()
		{	public void actionPerformed(ActionEvent e)
			{	minText.setEnabled(false);
				maxText.setEnabled(false);
				minLabel.setEnabled(false);
				maxLabel.setEnabled(false);
				if(currentLimits==FIXED_LIMITS)
				{	// when switching to global units from fixed units, must revalidate limits
					MoviePlotWindow movieFrame=docCtrl.getMovieFrame();
					movieFrame.plotView.dataLimitsSet=false;
				}
				currentLimits=GLOBAL_LIMITS;
			}
		});
		
		// Row 3: Third radio button
		gridbag.setConstraints(fixedLimits, c);
		if(prevLimits==2) fixedLimits.setSelected(true);
		limitType.add(fixedLimits);
		add(fixedLimits);
		fixedLimits.addActionListener(new ActionListener()
		{	public void actionPerformed(ActionEvent e)
			{	minText.setEnabled(true);
				maxText.setEnabled(true);
				minLabel.setEnabled(true);
				maxLabel.setEnabled(true);
				currentLimits=FIXED_LIMITS;
			}
		});
		
		// Row 4: settings
		minLabel.setHorizontalAlignment(JTextField.RIGHT);
		c.gridwidth = 1;
		c.weightx = 1.0;
		gridbag.setConstraints(minLabel, c);
		add(minLabel);
		
		c.weightx = 2.0;
		gridbag.setConstraints(minText, c);
		add(minText);
		
		maxLabel.setHorizontalAlignment(JTextField.RIGHT);
		c.weightx = 1.0;
		gridbag.setConstraints(maxLabel, c);
		add(maxLabel);
		
		c.insets=new Insets(0,0,0,6);
		c.weightx = 2.0;
		c.gridwidth=GridBagConstraints.REMAINDER;
		gridbag.setConstraints(maxText, c);
		add(maxText);
		
		setEnabled(LoadArchive.NO_PLOT,0);
	}
	
	// set current state
	public void setEnabled(int plotType,int plotComponent)
	{
		if(plotType==LoadArchive.MESH_PLOT || plotType==LoadArchive.PARTICLE_PLOT)
		{	this.setVisible(true);
			dynamicLimits.setEnabled(true);
			fixedLimits.setEnabled(true);
			if(dynamicLimits.isSelected())
			{	minText.setEnabled(false);
				maxText.setEnabled(false);
				minLabel.setEnabled(false);
				maxLabel.setEnabled(false);
			}
			else
			{	minText.setEnabled(true);
				maxText.setEnabled(true);
				minLabel.setEnabled(true);
				maxLabel.setEnabled(true);
			}
		}
		else
		{	this.setVisible(false);
			/*
			dynamicLimits.setEnabled(false);
			fixedLimits.setEnabled(false);
			minText.setEnabled(false);
			maxText.setEnabled(false);
			minLabel.setEnabled(false);
			maxLabel.setEnabled(false);
			*/
		}
	}
	
	// adjust limits if desired
	public Point2D.Double adjustLimits(double dmin,double dmax,double prevMin,double prevMax,boolean validLimits)
	{	try
		{	switch(getLimitType())
			{	case FIXED_LIMITS:
					dmin=ControlPanel.readDouble(minText,"minimum limit");
					dmax=ControlPanel.readDouble(maxText,"maximum limit");
					prevMinText = minText.getText();
					prevMaxText = maxText.getText();
					break;
				case GLOBAL_LIMITS:
					if(validLimits)
					{	dmin=Math.min(dmin,prevMin);
						dmax=Math.max(dmax,prevMax);
					}
					break;
				default:
					break;
					
			}
		}
		catch(Exception pe)
		{	Toolkit.getDefaultToolkit().beep();
			JNUtilities.showMessage(docCtrl,"Error: "+pe.getMessage()+" (plot will switch to dynamic limits)");
			dynamicLimits.setSelected(true);
			minText.setEnabled(false);
			maxText.setEnabled(false);
			minLabel.setEnabled(false);
			maxLabel.setEnabled(false);
		}
		return new Point2D.Double(dmin,dmax);
	}
	
	// get particle number option for time plots
	public int getLimitType()
	{	if(dynamicLimits.isSelected())
		{	prevLimits = 0;
			return DYNAMIC_LIMITS;
		}
		else if(globalLimits.isSelected())
		{	prevLimits = 1;
			return GLOBAL_LIMITS;
		}
		else
		{	prevLimits = 2;
			return FIXED_LIMITS;
		}
	}

}
