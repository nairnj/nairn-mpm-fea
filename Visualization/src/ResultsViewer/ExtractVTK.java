/*
 * ExtractVTK
 * NairnFEAMPMViz
 * 
 * Created 12/16/20216, John A. Nairn
 * Copyright 2021, RSAC Software, All Rights Reserved
 */

import geditcom.JNFramework.*;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;
import javax.swing.*;

// The class provides a modal dialog box to run ExtractMPM
// and runs it too (all in modal window)
class ExtractVTK extends JNDialog implements Runnable
{
	private static final long serialVersionUID = 1L;
	
	private JNCommandEditor dataField;
	private JTextField nameText=new JTextField("Particles");
	private JComboBox<PlotMenuItem> quant=new JComboBox<PlotMenuItem>();
	private DocViewer docView;
	private ResultsDocument resDoc;

	private ProcessBuilder builder = null;
	private Thread extractThread;
	private Process process = null;
	private boolean running;
	private String lastCommands = "";
	private JButton resetBtn = new JButton("   Reset   ");
	private Dimension bs;
	
	private JButton defBtn;
	
	//----------------------------------------------------------------------------
	// initialize
	//----------------------------------------------------------------------------

	// Create an instance of ExtractVTK dialog to pick quantites to export to VTK files
	protected ExtractVTK(DocViewer docCtrl)
	{	super((JNFrame)docCtrl,"ExtractMPM","Export particle data to VTK Files","Save","Cancel",null);

		docView = docCtrl;
		resDoc = docCtrl.resDoc;
		running = false;
		
		JPanel vtkPanel=new JPanel();
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		vtkPanel.setLayout(gridbag);

		int top=6;
		int bottom = 3;
		int left = 2;
		int right = 2;
		double lweight = 0.,rweight = 8.;
		
		// -----------------------------
		// Row 1: Label, (edit fieldX2)
		
		// Root label
		c.insets=new Insets(top,left,bottom,right);			// tlbr
		c.gridx=0;
		c.gridwidth = 1;
		c.weightx = lweight;		
		c.weighty = 1;
		c.fill = GridBagConstraints.HORIZONTAL;
		JLabel nameLabel=new JLabel("Root Name:",JLabel.RIGHT);
		gridbag.setConstraints(nameLabel, c);
		vtkPanel.add(nameLabel);
		
		// Root editing field
		c.gridx=1;
		c.weightx = rweight;
		c.gridwidth = 2;
		gridbag.setConstraints(nameText, c);
		vtkPanel.add(nameText);
		
		// -------------------------------
		// Row 2: (menuX2), (data field 1 of 4)
		
		// -q menu to insert
		quant.removeAllItems();
		quant.addItem(new PlotMenuItem("Insert..."));
		byte [] arch=resDoc.archFormat.getBytes();
		quant.addItem(new PlotMenuItem("mat"));
		quant.addItem(new PlotMenuItem("mass"));
		if(arch[ReadArchive.ARCH_Size]=='Y')
			quant.addItem(new PlotMenuItem("meansize"));
		if(arch[ReadArchive.ARCH_Velocity]=='Y')
		{	quant.addItem(new PlotMenuItem("velocity"));
			quant.addItem(new PlotMenuItem("velx"));
			quant.addItem(new PlotMenuItem("vely"));
			quant.addItem(new PlotMenuItem("velz"));
		}
		quant.addItem(new PlotMenuItem("displacement"));
		quant.addItem(new PlotMenuItem("dispx"));
		quant.addItem(new PlotMenuItem("dispy"));
		quant.addItem(new PlotMenuItem("dispz"));
		if(arch[ReadArchive.ARCH_Stress]=='Y')
		{	quant.addItem(new PlotMenuItem("stress"));
			quant.addItem(new PlotMenuItem("sxx"));
			quant.addItem(new PlotMenuItem("syy"));
			quant.addItem(new PlotMenuItem("szz"));
			quant.addItem(new PlotMenuItem("sxy"));
			if(resDoc.is3D())
			{	quant.addItem(new PlotMenuItem("sxz"));
				quant.addItem(new PlotMenuItem("syz"));
			}
			quant.addItem(new PlotMenuItem("pressure"));
			quant.addItem(new PlotMenuItem("vonmises"));
		}
		if(arch[ReadArchive.ARCH_Strain]=='Y')
		{	quant.addItem(new PlotMenuItem("strain"));
			quant.addItem(new PlotMenuItem("exx"));
			quant.addItem(new PlotMenuItem("eyy"));
			quant.addItem(new PlotMenuItem("ezz"));
			quant.addItem(new PlotMenuItem("exy"));
			if(resDoc.is3D())
			{	quant.addItem(new PlotMenuItem("exz"));
				quant.addItem(new PlotMenuItem("eyz"));
			}
			quant.addItem(new PlotMenuItem("equivstrain"));
		}
		if(arch[ReadArchive.ARCH_PlasticStrain]=='Y')
		{	quant.addItem(new PlotMenuItem("plasticstrain"));
			quant.addItem(new PlotMenuItem("pxx"));
			quant.addItem(new PlotMenuItem("pyy"));
			quant.addItem(new PlotMenuItem("pzz"));
			quant.addItem(new PlotMenuItem("pxy"));
			if(resDoc.is3D())
			{	quant.addItem(new PlotMenuItem("pxz"));
				quant.addItem(new PlotMenuItem("pyz"));
			}
		}
		if(arch[ReadArchive.ARCH_WorkEnergy]=='Y')
			quant.addItem(new PlotMenuItem("work"));
		if(arch[ReadArchive.ARCH_HeatEnergy]=='Y')
			quant.addItem(new PlotMenuItem("heat"));
		if(arch[ReadArchive.ARCH_StrainEnergy]=='Y')
			quant.addItem(new PlotMenuItem("strerg"));
		if(arch[ReadArchive.ARCH_PlasticEnergy]=='Y')
			quant.addItem(new PlotMenuItem("plerg"));
		if(arch[ReadArchive.ARCH_DeltaTemp]=='Y')
			quant.addItem(new PlotMenuItem("temp"));
		if(arch[ReadArchive.ARCH_Concentration]=='Y')
			quant.addItem(new PlotMenuItem("conc"));
		if(resDoc.is3D())
		{	quant.addItem(new PlotMenuItem("angx"));
			quant.addItem(new PlotMenuItem("angy"));
		}
		quant.addItem(new PlotMenuItem("angz"));
		if(arch[ReadArchive.ARCH_History]=='Y')
			quant.addItem(new PlotMenuItem("hist"));
		else if(arch[ReadArchive.ARCH_History]!='N')
		{	int history=(int)arch[ReadArchive.ARCH_History];
			if((history & 0x01) !=0)
				quant.addItem(new PlotMenuItem("hist1"));
			if((history & 0x02) !=0)
				quant.addItem(new PlotMenuItem("hist2"));
			if((history & 0x04) !=0)
				quant.addItem(new PlotMenuItem("hist3"));
			if((history & 0x08) !=0)
				quant.addItem(new PlotMenuItem("hist4"));
		}
		if(arch[ReadArchive.ARCH_History59]=='Y')
			quant.addItem(new PlotMenuItem("hist5"));
		else if(arch[ReadArchive.ARCH_History59]!='N')
		{	int history=(int)arch[ReadArchive.ARCH_History59];
			if((history & 0x01) !=0)
				quant.addItem(new PlotMenuItem("hist5"));
			if((history & 0x02) !=0)
				quant.addItem(new PlotMenuItem("hist6"));
			if((history & 0x04) !=0)
				quant.addItem(new PlotMenuItem("hist7"));
			if((history & 0x08) !=0)
				quant.addItem(new PlotMenuItem("hist8"));
			if((history & 0x10) !=0)
				quant.addItem(new PlotMenuItem("hist9"));
		}
		if(arch[ReadArchive.ARCH_History1014]=='Y')
			quant.addItem(new PlotMenuItem("hist10"));
		else if(arch[ReadArchive.ARCH_History1014]!='N')
		{	int history=(int)arch[ReadArchive.ARCH_History1014];
			if((history & 0x01) !=0)
				quant.addItem(new PlotMenuItem("hist10"));
			if((history & 0x02) !=0)
				quant.addItem(new PlotMenuItem("hist11"));
			if((history & 0x04) !=0)
				quant.addItem(new PlotMenuItem("hist12"));
			if((history & 0x08) !=0)
				quant.addItem(new PlotMenuItem("hist13"));
			if((history & 0x10) !=0)
				quant.addItem(new PlotMenuItem("hist14"));
		}
		
		if(arch[ReadArchive.ARCH_History1519]=='Y')
			quant.addItem(new PlotMenuItem("History 15"));
		else if(arch[ReadArchive.ARCH_History1519]!='N')
		{	int history=(int)arch[ReadArchive.ARCH_History1519];
			if((history & 0x01) !=0)
				quant.addItem(new PlotMenuItem("hist15"));
			if((history & 0x02) !=0)
				quant.addItem(new PlotMenuItem("hist16"));
			if((history & 0x04) !=0)
				quant.addItem(new PlotMenuItem("hist17"));
			if((history & 0x08) !=0)
				quant.addItem(new PlotMenuItem("hist18"));
			if((history & 0x10) !=0)
				quant.addItem(new PlotMenuItem("hist19"));
		}
		
		quant.setToolTipText("Add quantity to be exported");
		quant.setFocusable(false);
		
		c.gridx=0;
		c.weightx = 1.0;
		c.weighty = 0.;
		c.fill = GridBagConstraints.NONE;
		gridbag.setConstraints(quant, c);
		vtkPanel.add(quant);
		
		// text field
		dataField=new JNCommandEditor();
		dataField.textPane.setLineWrap(false);
		
		c.gridx=2;
		c.weightx = rweight;
		c.weighty = 10.;
		c.gridheight = 4;
		c.gridwidth = 1;
		c.fill = GridBagConstraints.BOTH;
		gridbag.setConstraints(dataField, c);
		vtkPanel.add(dataField);
		
		// ------------------------------------
		// Row 3: (Default Button X2) (data field 2 below)
		
		// button
		c.gridx=0;
		c.weightx = 1.0;
		c.weighty = 0.;
		c.gridheight = 1;
		c.gridwidth = 2;
		c.insets=new Insets(0,left,3,right);			// tlbr
		c.fill = GridBagConstraints.NONE;
		defBtn = new JButton("Default Set");
		defBtn.setToolTipText("Set to default set of available quantities");
		defBtn.setFocusable(false);
		gridbag.setConstraints(defBtn, c);
		vtkPanel.add(defBtn);
		
		// ------------------------------------
		// Row 4: (Reset Button X2) (data field 2 below)
		
		// button
		resetBtn.setToolTipText("Reset text field to previous commands");
		resetBtn.setFocusable(false);
		resetBtn.setEnabled(false);
		gridbag.setConstraints(resetBtn, c);
		vtkPanel.add(resetBtn);
		
		// adjust size
		bs = defBtn.getPreferredSize();
		Dimension bs2 = resetBtn.getPreferredSize();
		if(bs2.width>bs.width) bs.width = bs2.width;
		defBtn.setPreferredSize(bs);
		resetBtn.setPreferredSize(bs);

		// ---------------------------------------
		// Row 5 (blank X 2) (data field 3 below)
		c.weighty = 8.0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.SOUTH;
		JLabel blank = new JLabel("See ExtractMPM info for more options");
		gridbag.setConstraints(blank, c);
		vtkPanel.add(blank);
		blank.setHorizontalTextPosition(JLabel.CENTER);

		// add to dialog
		add(vtkPanel,BorderLayout.CENTER);
		
		// when quantity changes, update component menu
		quant.addItemListener(new ItemListener()
		{	public void itemStateChanged(ItemEvent e)
			{	if(quant.getSelectedIndex()>0)
				{	PlotMenuItem pm=(PlotMenuItem)quant.getSelectedItem();
					String inserted = dataField.textPane.getText();
					String value = pm.toString();
					int offset = inserted.indexOf(value);
					if(offset>=0)
					{	// not a match if next character is not return
						offset += value.length();
						if(offset<inserted.length())
						{	if(inserted.charAt(offset)!='\n')
								offset = -1;
						}
					}
					if(offset<0)
					{	if(inserted.length()>0)
						{	if(inserted.charAt(inserted.length()-1)!='\n')
								dataField.textPane.append("\n");
						}
						dataField.appendLine(value);
					}
				}
			}
		});

		defBtn.addActionListener(new ActionListener()
		{	public void actionPerformed(ActionEvent e)
			{   String def = "mat\nmass\ndisplacement\n";
				if(arch[ReadArchive.ARCH_Size]=='Y') def=def+"meansize\n";
				if(arch[ReadArchive.ARCH_Velocity]=='Y') def=def+"velocity\n";
				if(arch[ReadArchive.ARCH_Stress]=='Y') def=def+"stress\n";
				if(arch[ReadArchive.ARCH_Strain]=='Y') def=def+"strain\n";
				if(arch[ReadArchive.ARCH_PlasticStrain]=='Y') def=def+"plasticstrain\n";
				if(arch[ReadArchive.ARCH_DeltaTemp]=='Y')  def=def+"temp\n";
				if(arch[ReadArchive.ARCH_Concentration]=='Y') def=def+"conc\n";
				if(arch[ReadArchive.ARCH_History]=='Y')
					def=def+"hist1\n";
				else if(arch[ReadArchive.ARCH_History]!='N')
				{	int history=(int)arch[ReadArchive.ARCH_History];
					if((history & 0x01) !=0)
						def=def+"hist1\n";
					if((history & 0x02) !=0)
						def=def+"hist2\n";
					if((history & 0x03) !=0)
						def=def+"hist3\n";
					if((history & 0x04) !=0)
						def=def+"hist4\n";
				}
				if(arch[ReadArchive.ARCH_Size]=='Y'&& arch[ReadArchive.ARCH_Strain]=='Y')
					def="-U\n"+def;
				dataField.textPane.setText(def);
			}
		});

		resetBtn.addActionListener(new ActionListener()
		{	public void actionPerformed(ActionEvent e)
			{	dataField.textPane.setText(lastCommands);
			}
		});
		// initial size then override to resizable
		add(new JLabel("  "),BorderLayout.WEST);
		add(new JLabel("  "),BorderLayout.EAST);
		int width=500;
		int height=400;
		super.setSize(width,height);
		try
		{	Point windowLoc=docView.getLocationOnScreen();
			Dimension psize=docView.getSize();
			setLocation(windowLoc.x+(psize.width-width)/2,windowLoc.y+40);
		}
		catch(Exception e)
		{	setLocation(40,40);
		}
		setResizable(true);
		setVisible(true);
	}
	
	// over ride while running
	public void actionPerformed(ActionEvent e)
	{	if(!running)
		{	super.actionPerformed(e);
			return;
		}
	
		String theCmd=e.getActionCommand();
		
		if(theCmd.equals("cancel"))
		{	running = false;
		}
		else if(theCmd.equals("ok"))
		{	JNApplication.appBeep();
		}
	}
	
	// When click "Save" change plot data set to the newly edited data.
	public boolean dialogDone(int result)
	{	// get the new text
		try
		{	runExtractMPM();
			return false;
		}
		catch(Exception e)
		{	Toolkit.getDefaultToolkit().beep();
			JNUtilities.showMessage(null,e.getMessage());
			return false;
		}
		
		//return true;
	}
	
	// run ExtractMPM
	public boolean runExtractMPM() throws Exception
	{
		// save the commands
		lastCommands = dataField.getCommands();
		
		// get shell command and set command style
		String pathDelim = "/";
		String exe = "";
		if(NairnFEAMPMViz.isWindowsOS())
		{	pathDelim = "\\";
			exe = ".exe";
		}
		
		// Get  myCmd = path to executable
		// get path to ExtractMPM required in same folder as NairnMPM
		String myCmd=NFMVPrefs.prefs.get(NFMVPrefs.NairnMPMKey,NFMVPrefs.NairnMPMDef);
		if(myCmd.indexOf("$(bundle)")>=0)
			myCmd=NairnFEAMPMViz.bundleFolder+"ExtractMPM"+exe;
		else
		{	int lastSep = myCmd.lastIndexOf(pathDelim);
			if(lastSep<0) throw new Exception("Invalid path name for ExtractMPM:\n"+myCmd);
			myCmd = myCmd.substring(0,lastSep+1) + "ExtractMPM"+exe;		
		}
					
		// get output file or remote settings for output files
		String fileRoot = nameText.getText();
		if(fileRoot.length()==0 || fileRoot.indexOf('\n')>=0 || fileRoot.indexOf('\r')>=0)
			throw new Exception("Invalid root file name for exporting");
		
		// Get working directory and input root
		String archiveFolder = resDoc.getArchiveRoot();
		
		// build shell commands for the working directory
		// Windows: "app.exe" (options) "(input file)"
		// Mac/Linux: 'app' (options) '(input file)'
		
		// pbcmds will be [(app),[options],(input file)]
		ArrayList<String> pbcmds=new ArrayList<String>(20);
				
		// as direct process builder command
		pbcmds.add(myCmd);
		
		// shell options
		pbcmds.add("-VPs");
		pbcmds.add("-o");
		pbcmds.add(fileRoot);
		
		// each quantity
		Scanner s=new Scanner(lastCommands);
		boolean hasCommands=false;
		s.useDelimiter("[\\n\\r]");
		while(s.hasNext())
		{	String quantity = s.next();
			if(quantity.length()>0)
			{	if(quantity.charAt(0)!='-')
				{	pbcmds.add("-q");
					pbcmds.add(quantity);
				}
				else
				{	// look for command
					String[] argWords = quantity.split("\\s+");
					if(argWords.length==1)
						pbcmds.add(argWords[0]);
					else if(argWords.length==2)
					{	pbcmds.add(argWords[0]);
						pbcmds.add(argWords[1]);
					}
					else
					{	s.close();
						throw new Exception("A dash (-) option has too many words");
					}
				
				}
				hasCommands = true;
			}
		}
				
		// form into a table
		s.close();
		
		// give up if now commands
		if(!hasCommands)
			throw new Exception("No valid ExtractMPM options were provided");
						
		// each input file
		for(int i=0;i<resDoc.mpmArchives.size();i++)
		{	File archive = resDoc.mpmArchives.get(i).getFile();
			pbcmds.add(archive.getName());
		}

		// create the process and set working directory
		builder = new ProcessBuilder(pbcmds);
		builder.directory(new File(archiveFolder));
		
		// change button before running
		JButton btn = getButton(OK_BUTTON);
		btn.setText("Running");
		btn.setEnabled(false);
		setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);

		// start the thread
		extractThread=new Thread(this);
		extractThread.start();
		
		return true;
	
	}
	
	//----------------------------------------------------------------------------
	// detachable to run the process
	//----------------------------------------------------------------------------
	
	public void run()
	{	
		// Start process first
		running = true;
		resetBtn.setEnabled(false);
		defBtn.setEnabled(false);
		resetBtn.setEnabled(false);
		quant.setEnabled(false);
		try
		{	process = builder.start();
		
			InputStream is = process.getInputStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
		
			InputStream es = process.getErrorStream();
			InputStreamReader esr = new InputStreamReader(es);
			BufferedReader ebr = new BufferedReader(esr);

			dataField.setCommands("");
			String line;
			String errMsg = "";
		
			// read standard output
			while((line = br.readLine()) != null)
			{	dataField.appendLine(line);
				if(!running) break;
			}
			
			// if still running, check for error message
			if(running)
			{	while((line = ebr.readLine()) != null)
				{	if(errMsg.length()>0) errMsg = errMsg + "\n";
					errMsg = errMsg+line;
					if(!running) break;
				}
			}
		
			// get results code or error message
			if(running)
			{	try
				{	process.waitFor();
					//result = process.exitValue();
				}
				catch(InterruptedException e)
				{	errMsg = e.getLocalizedMessage();
				}
			}
			
			if(!running)
				dataField.appendLine("Canceled");
		
			// close all
			is.close();
			es.close();
			
			if(errMsg.length()>0)
			{	//dataField.setCommands(lastCommands);
				throw new Exception(errMsg);
			}
		}
		catch(Exception tpe)
		{	JNApplication.appBeep();
			JNUtilities.showMessage(null,tpe.getLocalizedMessage());
		}
		
		JButton btn = getButton(OK_BUTTON);
		btn.setText("Save");
		btn.setEnabled(true);
		defBtn.setEnabled(true);
		setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		resetBtn.setEnabled(true);
		resetBtn.setSize(bs);
		quant.setEnabled(true);
		
		running = false;
	}
}
