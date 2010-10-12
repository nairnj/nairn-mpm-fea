/*
 * ScaleResultsDialog
 * NairnFEAMPMViz
 * 
 * Created 10/11/2010, John A. Nairn
 * Copyright 2010, RSAC Software, All Rights Reserved
 */

import geditcom.JNFramework.*;
import java.awt.*;
import javax.swing.*;

public class ScaleResultsDialog extends JNDialog
{
	private static final long serialVersionUID = 1L;
	
	private JPanel values=new JPanel();
	private JComboBox lengthUnits=null;
	private JComboBox timeUnits=null;
	
	// initialize
	ScaleResultsDialog(DocViewer docCtrl)
	{	super(docCtrl,"Scale Analysis Results","Select length and time units",
			" OK ","Cancel",null);
		
		// units menus
		values.setLayout(new GridLayout(2,2));
		
		values.add(new JLabel("Length Units:",JLabel.RIGHT));
		String[] lengthOptions={"nm","µm","mm","cm","m","mils","in","ft"};
		lengthUnits=new JComboBox(lengthOptions);
		for(int i=0;i<lengthOptions.length;i++)
		{	if(lengthOptions[i].equals(docCtrl.resDoc.distU))
			{	lengthUnits.setSelectedIndex(i);
				break;
			}
		}
		values.add(lengthUnits);
		
		values.add(new JLabel("Time Units:",JLabel.RIGHT));
		String[] timeOptions={"µs","ms","sec"};
		timeUnits=new JComboBox(timeOptions);
		for(int i=0;i<timeOptions.length;i++)
		{	if(timeOptions[i].equals(docCtrl.resDoc.timeU))
			{	timeUnits.setSelectedIndex(i);
				break;
			}
		}
		values.add(timeUnits);
		
		add(values,BorderLayout.CENTER);
		
		// finish up
		setSize(220,150," ","   ");
	}

	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------

	// get length units
	public String getLengthUnits()
	{	return (String)lengthUnits.getSelectedItem();
	}
	
	// convert scale factor from mm to new units
	public double getLengthScale()
	{	String lu=getLengthUnits();
		if(lu.equals("nm")) return 1.e6;
		if(lu.equals("µm")) return 1000.;
		if(lu.equals("mm")) return 1.;
		if(lu.equals("cm")) return .1;
		if(lu.equals("m")) return .001;
		if(lu.equals("mils")) return 1000./25.4;
		if(lu.equals("in")) return 1./25.4;
		return 1./(12.*25.4);
	}
	
	// get time units
	public String getTimeUnits()
	{	return (String)timeUnits.getSelectedItem();
	}

	// convert scale factor from ms to new units
	public double getTimeScale()
	{	String lu=getTimeUnits();
		if(lu.equals("µs")) return 1000.;
		if(lu.equals("ms")) return 1.;
		return 0.001;
	}

}
