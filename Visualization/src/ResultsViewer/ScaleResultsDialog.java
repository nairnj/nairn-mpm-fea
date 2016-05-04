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
	private JComboBox<String> lengthUnits=null;
	private JComboBox<String> timeUnits=null;
	
	// initialize
	ScaleResultsDialog(DocViewer docCtrl)
	{	super(docCtrl,"Scale Results","Select length and time units",
			" OK ","Cancel",null);
		
		// units menus
		values.setLayout(new GridLayout(2,2));
		
		values.add(new JLabel("Length Units:",JLabel.RIGHT));
		String[] lengthOptions={"nm","\u03BCm","mm","cm","m","mils","in","ft"};
		lengthUnits=new JComboBox<String>(lengthOptions);
		lengthUnits.setSelectedIndex(docCtrl.resDoc.units.lengthScaleIndex());
		values.add(lengthUnits);
		
		values.add(new JLabel("Time Units:",JLabel.RIGHT));
		String[] timeOptions={"\u03BCs","ms","sec"};
		timeUnits=new JComboBox<String>(timeOptions);
		timeUnits.setSelectedIndex(docCtrl.resDoc.units.timeScaleIndex());
		values.add(timeUnits);
		
		add(values,BorderLayout.CENTER);
		
		// finish up
		setSize(220,150," ","   ");
	}

	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------

	// get length units
	public int getLengthIndex()
	{	return lengthUnits.getSelectedIndex();
	}
	
	// get time units
	public int getTimeIndex()
	{	return timeUnits.getSelectedIndex();
	}

}
