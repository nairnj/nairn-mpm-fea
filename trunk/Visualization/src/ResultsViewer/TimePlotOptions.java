/*
 * TimePlotOptions.java
 * NairnFEAMPMViz
 *
 * Created by John Nairn on 8/1/07.
 * Copyright (c) 2007 RSAC Software. All rights reserved.
 */

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;

public class TimePlotOptions extends PlotControl
{
	static final long serialVersionUID=21L;
	
	private ButtonGroup plotWhat=new ButtonGroup();
	private JTextField ptNumberText=new JTextField("1");
	private JTextField matNumberText=new JTextField("1");
	private JRadioButton plotPoint=new JRadioButton("Point Number");
	private JRadioButton plotAll=new JRadioButton("Total All Materials");
	private JRadioButton plot1Mat=new JRadioButton("Total Selected Material");
	
	private JComboBox xyContour=new JComboBox();
	private JTextField functionText=new JTextField("1");
	private JLabel plusMinus=new JLabel("+/-");
	private JTextField rangeText=new JTextField();
	
	private static final int left=3,mid=3,right=3;
	private static final int top=6,rows=1,bottom=6;
	
	// initialize
	TimePlotOptions(DocViewer dc)
	{   super(ControlPanel.WIDTH,top+bottom+5*26+4*rows,dc);
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		//setLayout(new GridLayout(5,2));
		setLayout(gridbag);
		
		c.fill = GridBagConstraints.BOTH;

		// pt number RB and text field
		c.insets=new Insets(top,left,rows,0);
		c.gridwidth = 2;
		c.weightx = 0.0;
		gridbag.setConstraints(plotPoint, c);
		plotPoint.setSelected(true);
		plotWhat.add(plotPoint);
		add(plotPoint);
		plotPoint.addActionListener(new ActionListener()
		{	public void actionPerformed(ActionEvent e)
			{	ptNumberText.setEnabled(true);
				matNumberText.setEnabled(false);
		}
		});
		
		c.insets.set(top,mid,rows,right);
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 1.0;
		gridbag.setConstraints(ptNumberText, c);
		add(ptNumberText);
		
		// all points and blank
		c.insets.set(0,left,rows,0);
		c.gridwidth = 2;
		c.weightx = 0.0;
		gridbag.setConstraints(plotAll, c);
		plotAll.addActionListener(new ActionListener()
		{	public void actionPerformed(ActionEvent e)
			{	ptNumberText.setEnabled(false);
				matNumberText.setEnabled(false);
			}
		});
		plotWhat.add(plotAll);
		add(plotAll);
		
		JLabel blank=new JLabel(" ");
		c.insets.set(0,mid,rows,right);
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 1.0;
		gridbag.setConstraints(blank, c);
		add(blank);
		
		// all points on material and field number
		c.insets.set(0,left,rows,0);
		c.gridwidth = 2;
		c.weightx = 0.0;
		gridbag.setConstraints(plot1Mat, c);
		plot1Mat.addActionListener(new ActionListener()
		{	public void actionPerformed(ActionEvent e)
			{	ptNumberText.setEnabled(false);
				matNumberText.setEnabled(true);
			}
		});
		plotWhat.add(plot1Mat);
		add(plot1Mat);
		
		c.insets.set(0,mid,rows,right);
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 1.0;
		gridbag.setConstraints(matNumberText, c);
		add(matNumberText);
		
		// contour plot options
		c.insets.set(0,left,rows,0);
		xyContour.addItem(new PlotMenuItem("x="));
		xyContour.addItem(new PlotMenuItem("y="));
		xyContour.addItem(new PlotMenuItem("D="));
		xyContour.addItem(new PlotMenuItem("T="));
		c.gridwidth = 1;
		c.weightx = 0.0;
		gridbag.setConstraints(xyContour, c);
		add(xyContour);
		
		c.insets.set(0,mid,rows,right);
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 1.0;
		gridbag.setConstraints(functionText, c);
		add(functionText);
		
		// +/- option
		c.insets.set(0,left,bottom,0);
		c.gridwidth = 2;
		c.weightx = 0.0;
		gridbag.setConstraints(plusMinus, c);
		plusMinus.setHorizontalAlignment(JTextField.RIGHT);
		add(plusMinus);
		
		c.insets.set(0,mid,bottom,right);
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 1.0;
		gridbag.setConstraints(rangeText, c);
		add(rangeText);
		
		// enable at first
		setEnabled(LoadArchive.NO_PLOT,0);
	}
	
	// set current state
	public void setEnabled(int plotType,int plotComponent)
	{
		if(plotType==LoadArchive.TIME_PLOT)
		{	switch(plotComponent)
			{	case PlotQuantity.MPMSTRENERGY:
				case PlotQuantity.MPMKINENERGY:
				case PlotQuantity.MPMENERGY:
				case PlotQuantity.MPMEXTWORK:
				case PlotQuantity.MPMPLASTICENERGY:
				case PlotQuantity.MPMTHERMALENERGY:
				case PlotQuantity.MPMELEMENTCROSSINGS:
					plotPoint.setEnabled(true);
					if(plotPoint.isSelected())
						ptNumberText.setEnabled(true);
					else
						ptNumberText.setEnabled(false);	
					plotAll.setEnabled(true);
					plot1Mat.setEnabled(true);
					if(plot1Mat.isSelected())
						matNumberText.setEnabled(true);
					else
						matNumberText.setEnabled(false);
					break;
				
				case PlotQuantity.MPMJ1:
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
					setEnabled(LoadArchive.NO_PLOT,0);
					break;
				
				default:
					plotPoint.setEnabled(true);
					if(plotPoint.isSelected())
						ptNumberText.setEnabled(true);
					else
						ptNumberText.setEnabled(false);	
					plotAll.setEnabled(true);
					plot1Mat.setEnabled(true);
					if(plot1Mat.isSelected())
						matNumberText.setEnabled(true);
					else
						matNumberText.setEnabled(false);
					break;
					/*
					plotPoint.setEnabled(true);
					ptNumberText.setEnabled(true);
					plotPoint.setSelected(true);
					plotAll.setEnabled(false);
					plot1Mat.setEnabled(false);
					matNumberText.setEnabled(false);
					break;
					*/
			}
			xyContour.setEnabled(false);
			functionText.setEnabled(false);
			plusMinus.setEnabled(false);
			rangeText.setEnabled(false);
		}
		else if(plotType==LoadArchive.MESH2D_PLOT)
		{	plotPoint.setEnabled(false);
			plotAll.setEnabled(false);
			plot1Mat.setEnabled(false);
			ptNumberText.setEnabled(false);			
			matNumberText.setEnabled(false);
			boolean popup=true;
			switch(plotComponent)
			{	case PlotQuantity.MPMNORMALCTOD:
				case PlotQuantity.MPMSHEARCTOD:
				case PlotQuantity.MPMCRACKPROFILE:
				case PlotQuantity.MPMOPENINGFRACTION:
				case PlotQuantity.MPMSHEARFRACTION:
					popup=false;
					break;
				default:
					break;
			}
			xyContour.setEnabled(popup);
			functionText.setEnabled(popup);
			plusMinus.setEnabled(popup);
			rangeText.setEnabled(popup);
		}
		else
		{	plotPoint.setEnabled(false);
			plotAll.setEnabled(false);
			plot1Mat.setEnabled(false);
			ptNumberText.setEnabled(false);			
			matNumberText.setEnabled(false);
			xyContour.setEnabled(false);
			functionText.setEnabled(false);
			plusMinus.setEnabled(false);
			rangeText.setEnabled(false);
		}
	}
	
	// get material point number option
	//		>0 to plot that material point
	//		==0 to plot all materials
	//		<0 to plot that material number (- the setting)
	public int getParticleNumber() throws Exception
	{	if(plotPoint.isSelected())
			return ControlPanel.readInteger(ptNumberText,"point number");
		else if(plotAll.isSelected())
			return 0;
		else
			return -ControlPanel.readInteger(matNumberText,"material number");
	}
	
	// get selected item for contour menu options
	public int getContour() { return xyContour.getSelectedIndex(); }
	
	// get contour expression
	public String getContourFunction() { return functionText.getText(); }
	
	// get contour +/- range (or zero if empty)
	public double getPlusMinus() throws Exception
	{
		String xyRange=rangeText.getText();
		if(xyRange.length()==0) return 0.d;
		Scanner validate=new Scanner(xyRange);		// will use user's locale
		if(validate.hasNextDouble()) return validate.nextDouble();
		throw new Exception("The optional '+/-' field must be a valid number or be empty");
	}

}
