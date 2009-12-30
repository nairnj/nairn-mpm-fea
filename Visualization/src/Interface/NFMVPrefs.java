/*******************************************************************
	NFMVPrefs.java
	NairnFEAMPMViz

	Created by John Nairn on 2/17/08.
	Copyright (c) 2008 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.io.File;
import java.util.prefs.*;
import javax.swing.*;

import java.awt.event.*;
import geditcom.plot2d.*;

public class NFMVPrefs extends JFrame implements ActionListener
{
	static final long serialVersionUID=26L;
	
	public static Preferences prefs=Preferences.userRoot().node("com/geditcom/NairnFEAMPMViz");
	
	// Results window split locations
	public static String ResultsTopSplitKey="Results Top Split";
	public static double ResultsTopSplitDef=0.35;
	public static String ResultsSplitKey="Results Split";
	public static double ResultsSplitDef=0.50;
	
	// Commands window split location
	public static String CommandsSplitKey="Commands Split";
	public static double CommandsSplitDef=0.60;
	
	// path to NairnMPM binary
	public static String NairnMPMKey="NairnMPM Path";
	public static String NairnMPMDef="/usr/local/bin/NairnMPM";
	
	// path to NairnMPM dtd
	public static String NairnMPMDTDKey="NairnMPM DTD File Path";
	public static String NairnMPMDTDDef="/usr/local/include/NairnMPM.dtd";
	public static String NairnMPMValidateKey="Validate MPM Calculations";
	public static boolean NairnMPMValidateDef=true;

	// path to NairnFEA binary
	public static String NairnFEAKey="NairnFEA Path";
	public static String NairnFEADef="/usr/local/bin/NairnFEA";
	
	// path to NairnFEA DTD
	public static String NairnFEADTDKey="NairnFEA DTD File Path";
	public static String NairnFEADTDDef="/usr/local/include/NairnFEA.dtd";
	public static String NairnFEAValidateKey="Validate MPM Calculations";
	public static boolean NairnFEAValidateDef=true;
	
	// path to shell command
	public static String ShellKey="Shell Path";
	public static String ShellDef="/bin/bash";
	
	// path to workspace folder
	public static String WorkSpaceKey="Workspace Folder";
	public static String WorkSpaceDef="";
	
	// plot spectrum
	public static String SpectrumKey="PlotSpectrum";
	public static int SpectrumDef=1;
	public static String NumContoursKey="NumberOfPlotContours";
	public static int NumContoursDef=0;
	
	// plot colors
	public static String backColorKey="Background Color";
	public static Color backColorDef=new Color(84,89,109);  		// 0.329f,0.349f,0.427f
	public static String meshLineColorKey="Mesh Line Color";
	public static Color meshLineColorDef=new Color(191,191,191);	// 0.75f,0.75f,0.75f
	public static String meshNodesColorKey="Nodal Point Color";
	public static Color meshNodesColorDef=new Color(212,212,212);	// 0.83f,0.83f,0.83f
	public static String textColorKey="Text Label Color";
	public static Color textColorDef=new Color(255,255,255);		// 1.f,1.f,1.f

	
	private JTextField mpmCodePath=new JTextField();
	private JTextField mpmDTDPath=new JTextField();
	private JTextField feaCodePath=new JTextField();
	private JTextField feaDTDPath=new JTextField();
	private JTextField shellPath=new JTextField();
	private JTextField workPath=new JTextField();
	private JTextField numContours=new JTextField();
	private JFileChooser chooser = new JFileChooser( );
	
	// prefs
	private JComboBox spectrumBox;
	
	private JCheckBox validateMPM=new JCheckBox("Validate",NairnMPMValidateDef);
	private JCheckBox validateFEA=new JCheckBox("Validate",NairnFEAValidateDef);
	
	//----------------------------------------------------------------------------
	// constructor
	//----------------------------------------------------------------------------

	public NFMVPrefs( )
	{   super("Preferences");
		
		JTabbedPane tabbedPane = new JTabbedPane();

		// code preferences pane -------------------------------------
		JPanel panel1 = new JPanel();
		
		panel1.add(filePathPanel("MPM Code",mpmCodePath,prefs.get(NairnMPMKey,NairnMPMDef),
									false,null));
		panel1.add(filePathPanel("MPM DTD File",mpmDTDPath,prefs.get(NairnMPMDTDKey,NairnMPMDTDDef),
									prefs.getBoolean(NairnMPMValidateKey,NairnMPMValidateDef),validateMPM));
		panel1.add(filePathPanel("FEA Code",feaCodePath,prefs.get(NairnFEAKey,NairnFEADef),
									false,null));
		panel1.add(filePathPanel("FEA DTD File",feaDTDPath,prefs.get(NairnFEADTDKey,NairnFEADTDDef),
									prefs.getBoolean(NairnFEAValidateKey,NairnFEAValidateDef),validateFEA));
		panel1.add(filePathPanel("Shell Command",shellPath,prefs.get(ShellKey,ShellDef),
									false,null));
		panel1.add(filePathPanel("Work Space Directory",workPath,prefs.get(WorkSpaceKey,WorkSpaceDef),
				false,null));
		
		tabbedPane.addTab("Code",null,panel1,"Settings for running calculations");

		// color preferences pane ----------------------------------------
		JPanel panel2 = new JPanel();
		tabbedPane.addTab("Colors",null,panel2,"Select colors used in plots");

		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints pc = new GridBagConstraints();		// for the pane
		panel2.setLayout(gridbag);
		
		// spectrum box
		GridBagLayout spectrumBag = new GridBagLayout();
		JNBoxPanel spectrumPan=new JNBoxPanel("Color Plot Spectrum",spectrumBag);
		GridBagConstraints c = new GridBagConstraints();		// for the box panel
		
		JLabel section=new JLabel("Type: ");
		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.WEST;
		c.gridx=0;
		c.gridwidth = 1;
		c.weightx = 0.0;
		c.insets=new Insets(6,6,1,0);			// tlbr
		spectrumBag.setConstraints(section,c);
		spectrumPan.add(section);
		
		String [] lines={"Rainbow 1 (Blue to Cyan)","Rainbow 2 (Blue to Red)","Rainbow 3 (Purple to Orange)","Grayscale (Black to White)","Cool Diverging (Blue to Red)"};
		spectrumBox=new JComboBox(lines);
		spectrumBox.setSelectedIndex(ColorPicker.getSpectrumType());
		spectrumBox.addActionListener(this);
		spectrumBox.setFocusable(false);
		spectrumBox.setToolTipText("Select one of the available spectra for coloring the plots.");
		c.insets=new Insets(6,3,1,6);			// tlbr
		c.gridx=1;
		c.weightx = 1.0;
		spectrumBag.setConstraints(spectrumBox,c);
		spectrumPan.add(spectrumBox);

		JLabel section1=new JLabel("Contours: ");
		c.insets=new Insets(6,6,1,0);			// tlbr
		c.gridx=0;
		c.weightx = 0.0;
		spectrumBag.setConstraints(section1,c);
		spectrumPan.add(section1);

		numContours.setText(""+ColorPicker.getNumberOfContours());
		numContours.setEditable(true);
		numContours.setColumns(5);
		numContours.setToolTipText("Enter 1 for continuous colors or >1 for number of discrete colors and hit return or enter");
		numContours.setActionCommand("change contours");
		numContours.addActionListener(this);
		c.insets=new Insets(6,3,1,6);			// tlbr
		c.gridx=1;
		c.weightx = 1.0;
		spectrumBag.setConstraints(numContours,c);
		spectrumPan.add(numContours);
		
		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx=0;
		c.gridwidth = 1;
		c.weightx = 1.0;

		pc.fill = GridBagConstraints.HORIZONTAL;
		pc.gridx=0;
		pc.gridwidth = 1;
		pc.weightx = 1.0;
		pc.insets=new Insets(1,6,1,6);			// tlbr
		pc.anchor = GridBagConstraints.CENTER;
		gridbag.setConstraints(spectrumPan,pc);
		panel2.add(spectrumPan);

		// plot element color wells box
		GridBagLayout wellsBag = new GridBagLayout();
		JNBoxPanel wellsPan=new JNBoxPanel("Plot Elements",wellsBag);
		GridBagConstraints ec = new GridBagConstraints();		// for the plot element box panel
		
		JLabel elab1=new JLabel("Background: ");
		ec.fill = GridBagConstraints.NONE;
		ec.anchor = GridBagConstraints.EAST;
		ec.insets=new Insets(10,0,0,0);			// tlbr
		ec.gridx=0;
		ec.weightx = 0.5;
		ec.gridwidth=1;
		wellsBag.setConstraints(elab1,ec);
		wellsPan.add(elab1);
		
		JNColorWell bgWell=new JNColorWell(getPrefColor(backColorKey,backColorDef));
		bgWell.addActionListener(this,"bgColor","Choose Background Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx=1;
		ec.weightx = 1.0;
		wellsBag.setConstraints(bgWell,ec);
		wellsPan.add(bgWell);
		
		JLabel elab2=new JLabel("Mesh Lines: ");
		ec.anchor = GridBagConstraints.EAST;
		ec.gridx=2;
		ec.weightx = 0.5;
		wellsBag.setConstraints(elab2,ec);
		wellsPan.add(elab2);
		
		JNColorWell mlWell=new JNColorWell(getPrefColor(meshLineColorKey,meshLineColorDef));
		mlWell.addActionListener(this,"mlColor","Choose Mesh Line Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx=3;
		ec.weightx = 1.0;
		wellsBag.setConstraints(mlWell,ec);
		wellsPan.add(mlWell);
		
		JLabel elab3=new JLabel("Labels: ");
		ec.insets=new Insets(10,0,6,0);			// tlbr
		ec.anchor = GridBagConstraints.EAST;
		ec.gridx=0;
		ec.weightx = 0.5;
		ec.gridwidth=1;
		wellsBag.setConstraints(elab3,ec);
		wellsPan.add(elab3);
		
		JNColorWell labWell=new JNColorWell(getPrefColor(textColorKey,textColorDef));
		labWell.addActionListener(this,"labColor","Choose Color for Text Labels");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx=1;
		ec.weightx = 1.0;
		wellsBag.setConstraints(labWell,ec);
		wellsPan.add(labWell);
		
		JLabel elab4=new JLabel("Nodal Points: ");
		ec.anchor = GridBagConstraints.EAST;
		ec.gridx=2;
		ec.weightx = 0.5;
		wellsBag.setConstraints(elab4,ec);
		wellsPan.add(elab4);
		
		JNColorWell npWell=new JNColorWell(getPrefColor(meshNodesColorKey,meshNodesColorDef));
		npWell.addActionListener(this,"npColor","Choose Nodal Point Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx=3;
		ec.weightx = 1.0;
		wellsBag.setConstraints(npWell,ec);
		wellsPan.add(npWell);

		gridbag.setConstraints(wellsPan,pc);
		panel2.add(wellsPan);
	
		// FEA and MPM preferences pane -----------------------------------------
		JPanel panel3 = new JPanel();
		tabbedPane.addTab("FEA and MPM",null,panel3,"FEA and MPM plot settings");
		
		add(tabbedPane);
		setSize(600,480);
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
	}
	
	public JPanel filePathPanel(String pathName,JTextField pathField,String pathText,boolean validPref,JCheckBox validate)
	{
		JPanel filePanel=new JPanel();
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		filePanel.setLayout(gridbag);

		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx=0;
		c.gridwidth = 1;
		c.weightx = 1.0;
		c.insets=new Insets(1,6,1,0);			// tlbr
		JLabel nameLabel=new JLabel(pathName+" Path:");
		gridbag.setConstraints(nameLabel,c);
		filePanel.add(nameLabel);
		
		int rowWidth=2;
		if(validate!=null)
		{	c.gridx++;
			validate.setSelected(validPref);
			validate.setActionCommand(pathName+" Use");
			validate.addActionListener(this);
			validate.setToolTipText("Check to validate input files using the provided DTD file");
			gridbag.setConstraints(validate, c);
			filePanel.add(validate);
			rowWidth++;
		}

		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.EAST;
		c.gridx++;
		c.insets=new Insets(1,0,1,6);			// tlbr
		JButton change=new JButton("Change...");
		change.setActionCommand(pathName);
		change.addActionListener(this);
		change.setToolTipText("Click to change the "+pathName+" path");
		gridbag.setConstraints(change,c);
		filePanel.add(change);
		
		c.anchor = GridBagConstraints.CENTER;
		c.gridx=0;
		c.gridwidth = rowWidth;
		c.insets=new Insets(0,6,1,9);			// tlbr
		pathField.setText(pathText);
		pathField.setEditable(false);
		pathField.setColumns(45);
		pathField.setToolTipText("The current "+pathName+" path; click 'Change...' to change it.");
		gridbag.setConstraints(pathField,c);
		filePanel.add(pathField);
		
		return filePanel;
	}

	//----------------------------------------------------------------------------
	// handle application commands
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{   String theCmd=e.getActionCommand();

		if(theCmd.equals("MPM Code"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			mpmCodePath.setText(newPath);
			prefs.put(NairnMPMKey,newPath);
		}
		
		else if(theCmd.equals("MPM DTD File Use"))
		{	prefs.putBoolean(NairnMPMValidateKey,validateMPM.isSelected());
		}
		
		else if(theCmd.equals("MPM DTD File"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			mpmDTDPath.setText(newPath);
			prefs.put(NairnMPMDTDKey,newPath);
		}
		
		else if(theCmd.equals("FEA Code"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			feaCodePath.setText(newPath);
			prefs.put(NairnFEAKey,newPath);
		}
		
		else if(theCmd.equals("FEA DTD File Use"))
		{	prefs.putBoolean(NairnFEAValidateKey,validateFEA.isSelected());
		}
		
		else if(theCmd.equals("bgColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(backColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("labColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(textColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("mlColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(meshLineColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("npColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(meshNodesColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("FEA DTD File"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			feaDTDPath.setText(newPath);
			prefs.put(NairnFEADTDKey,newPath);
		}
		
		else if(theCmd.equals("Shell Command"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			shellPath.setText(newPath);
			prefs.put(ShellKey,newPath);
		}
		
		else if(theCmd.equals("Work Space Directory"))
		{	chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			int result=chooser.showOpenDialog(this);
			String newPath=chooser.getSelectedFile().getPath();
			chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			if(result==JFileChooser.CANCEL_OPTION) return;
			workPath.setText(newPath);
			prefs.put(WorkSpaceKey,newPath);
			setWorkspace(NairnFEAMPMViz.appCtrl.chooser);
		}
		
		else if(theCmd.equals("change contours"))
		{	try
			{	int newContours=Integer.valueOf(numContours.getText());
				if(newContours<1) newContours=1;
				int oldNum=ColorPicker.getNumberOfContours();
				if(newContours!=oldNum)
				{	prefs.putInt(NumContoursKey,newContours);
					ColorPicker.setNumberOfContours();
				}
			}
			catch(Exception ie)
			{	Toolkit.getDefaultToolkit().beep();
			}
		}
		
		else
        {	JComboBox cb = (JComboBox)e.getSource();
			if(cb==spectrumBox)
			{	int oldType=ColorPicker.getSpectrumType();
				int newType=cb.getSelectedIndex();
				if(newType!=oldType)
				{	prefs.putInt(SpectrumKey,newType);
					ColorPicker.setSpectrumType();
				}
			}
		}
    }
		
	//----------------------------------------------------------------------------
	// static class methods
	//----------------------------------------------------------------------------
	
	// set/get NFMVFrame size
	public static void setFrameSize(NFMVFrame f,Dimension d)
	{	prefs.putInt(f.frameHeightKey,d.height);
		prefs.putInt(f.frameWidthKey,d.width);
	}
	public static Dimension getFrameSize(NFMVFrame f)
	{	return new Dimension(prefs.getInt(f.frameWidthKey,f.frameWidthDef),
				prefs.getInt(f.frameHeightKey,f.frameHeightDef));
	}
	
	// set/get color using three float preferences
	public void setPrefColor(String colorKey,Color setColor)
	{	prefs.putInt(colorKey+" Red", setColor.getRed());
		prefs.putInt(colorKey+" Green", setColor.getGreen());
		prefs.putInt(colorKey+" Blue", setColor.getBlue());
	}
	public static Color getPrefColor(String colorKey,Color colorDef)
	{	try
		{	return new Color(prefs.getInt(colorKey+" Red", colorDef.getRed()),
					prefs.getInt(colorKey+" Green", colorDef.getGreen()),
					prefs.getInt(colorKey+" Blue", colorDef.getBlue()));
		}
		catch(Exception ce)
		{	return new Color(128,128,128);
		}
	}

	// set a chooser to current work space
	public static void setWorkspace(JFileChooser chooser)
	{	String workPath=prefs.get(WorkSpaceKey,WorkSpaceDef);
		if(workPath.length()>0)
		{	try
			{	File workDirectory=new File(workPath);
				chooser.setCurrentDirectory(workDirectory);
			}
			catch(Exception e) {}
		}
	}
	
}
