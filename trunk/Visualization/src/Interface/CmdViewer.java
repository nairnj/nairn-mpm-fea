/*******************************************************************
	CmdViewer.java
	NairnFEAMPMViz

	Created by John Nairn on Feb 14 2008.
	Copyright (c) 2008 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import java.io.*;
import java.util.*;

public class CmdViewer extends NFMVViewer implements Runnable, FocusListener
{
	static final long serialVersionUID=30L;
	
	//----------------------------------------------------------------------------
	// constants and variables
	//----------------------------------------------------------------------------
	
	public CommandDocument cmdDoc=new CommandDocument();		// one results document in this viewer
	private JSplitPane full;			// window split view - TextDisplay on top, Console on bottom
	private CommandEdit cmdField;		// editing zone at top of window
	private ConsolePane console;		// console pane on the bottom
	private Component currentFocus=null;	// which text area has focus
	private ProcessBuilder builder;
	private boolean running=false;
	private boolean openMesh=false;
	private boolean wasSubmitted=false;
	private Thread runThread=null;
	
	//----------------------------------------------------------------------------
	// constants and variables
	//----------------------------------------------------------------------------

	public CmdViewer()
	{   super("DocViewer");
		
		Container content=getContentPane( );
		cmdDoc.setDocController(this);
		
		// more tools
		buttonBar.add(new JLabel("     "));
		addToolBarBtn(buttonBar,"document-save.png","Save",
				"Save the current commands.",this);
		addToolBarBtn(buttonBar,"go-next.png","RunAnalysis",
				"Run FEA or MPM Analysis.",this);
		addToolBarBtn(buttonBar,"go-last.png","CheckAnalysis",
				"Check mesh for FEA or MPM Analysis.",this);
		addToolBarBtn(buttonBar,"image-x-generic.png","ShowPartner",
				"Show associated results window (if available).",this);

				
		// top displays sections and text of selected section
		cmdField=new CommandEdit();
		cmdField.textPane.addFocusListener(this);
		
		// console on the bottom
		console=new ConsolePane();
		console.soutPane.addFocusListener(this);
		
		// full window split pane
		full=new JSplitPane(JSplitPane.VERTICAL_SPLIT,cmdField,console);
		
		// add to content
		content.add(full,BorderLayout.CENTER);
		
		// size and location
		setFramePrefs("Commands Window Width",600,"Commands Window Height",800);
		setFrameLocAndSize(this);
	}
	
		
	// load the file
	public void loadTheFile(File file,String textContents) throws Exception
	{	cmdField.setCommands(textContents);
		cmdDoc.setFile(file);
		setTitle(cmdDoc.getName());
	}
	
	// save commands to existing or new file
	public void saveCommandFile(boolean getNameFirst)
	{	File saveFile=cmdField.saveCommands(cmdDoc.getFile(),getNameFirst,this);
		if(saveFile!=null)
		{	cmdDoc.setFile(saveFile);
			setTitle(cmdDoc.getName());
			NairnFEAMPMViz.appCtrl.addWindowMenuItem(this);
		}
	}

	//----------------------------------------------------------------------------
	// handle application commands
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{	String theCmd=e.getActionCommand();
	
		if(theCmd.equals("RunAnalysis"))
		{	runNFMAnalysis(false,false);
		}
		
		else if(theCmd.equals("BgAnalysis"))
		{	runNFMAnalysis(true,false);
		}
		
		else if(theCmd.equals("CheckAnalysis"))
		{	runNFMAnalysis(false,true);
		}
		
		else if(theCmd.equals("StopAnalysis"))
		{	running=false;
		}
		
		else if(theCmd.equals("Save"))
		{	saveCommandFile(false);
		}
		
		else if(theCmd.equals("SaveAs"))
		{	saveCommandFile(true);
		}
		
		else if(theCmd.equals("ShowPartner"))
		{	if(running)
			{	Toolkit.getDefaultToolkit().beep();
				JOptionPane.showMessageDialog(this,"Cannot open results while analysis is running.");
			}
			else if(partnerCtrl!=null)
				partnerCtrl.toFront();
			else
			{	File resultsFile=console.getFile();
				if(resultsFile!=null)
					NairnFEAMPMViz.appCtrl.openNFMFileNow(resultsFile, this, false);
				else
				{	Toolkit.getDefaultToolkit().beep();
					JOptionPane.showMessageDialog(this,"This window does not have an associated results.");
				}
			}
		}
		
		else if(theCmd.equals("GoToLine"))
		{	if(currentFocus==console.soutPane)
				System.out.println("In the output console");
			else if(currentFocus==cmdField.textPane)
				System.out.println("In the command editor");
			else
				System.out.println("No focus"+currentFocus);
		}
		
		else
			super.actionPerformed(e);
		
	}
	
	//----------------------------------------------------------------------------
	// Window events
	//----------------------------------------------------------------------------
	
	public void windowOpened(WindowEvent e)
	{	double splitLoc=NFMVPrefs.prefs.getDouble(NFMVPrefs.CommandsSplitKey,
													NFMVPrefs.CommandsSplitDef);
		if(splitLoc<0. || splitLoc>1.) splitLoc=NFMVPrefs.CommandsSplitDef;
	    full.setDividerLocation(splitLoc);
		JScrollBar vertBar=cmdField.getVerticalScrollBar();
		vertBar.setValue(vertBar.getMinimum());
	}
	
	// when closed
	public void windowClosing(WindowEvent e)
	{	if(!confirmClosing()) return;
	
		double loc=(double)full.getDividerLocation()/
					(double)(full.getMaximumDividerLocation()-full.getMinimumDividerLocation());
		NFMVPrefs.prefs.putDouble(NFMVPrefs.CommandsSplitKey,loc);
		super.windowClosing(e);
	}
	
	// before close, give change to same. Return false if canceled or truie if can proceed
	public boolean confirmClosing()
	{	if(cmdField.getChanged())
		{	Toolkit.getDefaultToolkit().beep();
			String message="The XML input commands have been modified\nDo you want to save the changes?";
			int result=JOptionPane.showConfirmDialog(this, message, "Save?", JOptionPane.YES_NO_CANCEL_OPTION );
			if(result==JOptionPane.CANCEL_OPTION)
				return false;
			else if (result==JOptionPane.YES_OPTION) 
				saveCommandFile(false);
		}
		return true;
	}

	//----------------------------------------------------------------------------
	// Window events
	//----------------------------------------------------------------------------
	
	public void runNFMAnalysis(boolean doBackground,boolean checkMesh)
	{
		// only allowed if the commands have been saved
		if(cmdDoc.getFile()==null)
		{	Toolkit.getDefaultToolkit().beep();
			JOptionPane.showMessageDialog(this,"The XML input commands have to be saved to a file before running an analysis.");
			return;
		}
		
		// what is process is current running?
		if(running)
		{	Toolkit.getDefaultToolkit().beep();
			String message="An FEA or MPM process is currently running.\nDo you want stop it and start a new process?";
			int result=JOptionPane.showConfirmDialog(this, message, "Continue?", JOptionPane.OK_CANCEL_OPTION);
			if(result==JOptionPane.CANCEL_OPTION) return;
			if(runThread.isAlive())
			{	running=false;
				try
				{	runThread.join();
				}
				catch(InterruptedException ie)
				{
				}
			}
		}
		
		// save the commands
		File inFile=cmdDoc.getFile();
		
		// prepare output path
		boolean mpmAnalysis=cmdField.isMPMAnalysis();
		boolean doValidate;
		String myCmd,myDTD;
		if(mpmAnalysis)
		{	if(!console.setOutputPath(inFile,"mpm")) return;
		
			// get path to NairnFEA and NairnMPM
			myCmd=NFMVPrefs.prefs.get(NFMVPrefs.NairnMPMKey,NFMVPrefs.NairnMPMDef);
			myDTD=NFMVPrefs.prefs.get(NFMVPrefs.NairnMPMDTDKey,NFMVPrefs.NairnMPMDTDDef);
			doValidate=NFMVPrefs.prefs.getBoolean(NFMVPrefs.NairnMPMValidateKey,NFMVPrefs.NairnMPMValidateDef);
		}
		else
		{	if(!console.setOutputPath(inFile,"fea")) return;
		
			// get path to NairnFEA and NairnMPM
			myCmd=NFMVPrefs.prefs.get(NFMVPrefs.NairnFEAKey,NFMVPrefs.NairnFEADef);
			myDTD=NFMVPrefs.prefs.get(NFMVPrefs.NairnFEADTDKey,NFMVPrefs.NairnFEADTDDef);
			doValidate=NFMVPrefs.prefs.getBoolean(NFMVPrefs.NairnFEAValidateKey,NFMVPrefs.NairnFEAValidateDef);
		}
		
		// insert DTD path if validating
		if(doValidate)
		{	if(!cmdField.insertDTD(myDTD))
			{	Toolkit.getDefaultToolkit().beep();
				JOptionPane.showMessageDialog(this,"The file does not start with the required <?xml ...?> element.");
				return;
			}
		}
		
		// copy command file (if needed)
		saveCommandFile(false);
		File outFile=console.getFile();
		if(!inFile.getParent().equals(outFile.getParent()))
			inFile=cmdField.saveCopyOfCommands(new File(outFile.getParent()+"/"+inFile.getName()),this);
		if(inFile==null) return;
			
		// set commands and options
		ArrayList<String> cmds=new ArrayList<String>(20);
		
		// start login bash shell
		cmds.add(NFMVPrefs.prefs.get(NFMVPrefs.ShellKey,NFMVPrefs.ShellDef));
		cmds.add("--login");			// only needed for windows
		cmds.add("-c");
		
		// build shell command
		StringBuffer shell=new StringBuffer();
		
		// cd to file directory
		shell.append("cd ");
		String shellCD=outFile.getParent();
		if(NairnFEAMPMViz.appCtrl.isWindows())
			shellCD=shellCD.replace('\\','/');
		if(shellCD.indexOf(' ')>=0)
			shell.append("'"+shellCD+"'; ");
		else
			shell.append(shellCD+"; ");

		// executable (as shell command)
		if(NairnFEAMPMViz.appCtrl.isWindows())
			myCmd=myCmd.replace('\\','/');
		if(myCmd.indexOf(' ')>=0)
			shell.append("'"+myCmd+"'");
		else
			shell.append(myCmd);
			
		// options (-a to abort after setup, -v to validate)
		if(checkMesh) shell.append(" -a");
		if(doValidate) shell.append(" -v");
			
		String inName=inFile.getName();
		if(inName.indexOf(' ')>=0)
			shell.append(" '"+inName+"'");
		else
			shell.append(" "+inName);
		
		if(doBackground)
		{	shell.append(" >& ");
			
			// output file
			String outFileName=outFile.getName();
			if(outFileName.indexOf(' ')>=0)
				shell.append("'"+outFileName+"'");
			else
				shell.append(outFileName);
				
			shell.append(" &");
		}
		
		// add the shell command
		System.out.println(shell);
		cmds.add(shell.toString());

		// create the process and set working directory
		builder = new ProcessBuilder(cmds);
		//Map<String, String> environ = builder.environment();
		//builder.directory(myFile.getParentFile());
		builder.redirectErrorStream(true);
		
		runThread=new Thread(this);
		running=true;
		openMesh=checkMesh;
		wasSubmitted=doBackground;
		runThread.start();
	}
	
	//----------------------------------------------------------------------------
	// detachable to run the process
	//----------------------------------------------------------------------------
	
	public void run()
	{	
		try
		{	Process process = builder.start();
			InputStream is = process.getInputStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			
			console.clear();
			String line;
			while((line = br.readLine()) != null)
			{	console.appendLine(line);
				if(!running) break;
			}
			
			// if done and no errors, save and open for visualization
			if(running)
			{	int result=1;
				try
				{	process.waitFor();
					result=process.exitValue();
				}
				catch (InterruptedException ie) {}
				if(result==0)
				{	if(wasSubmitted)
					{	JOptionPane.showMessageDialog(this,"FEA or MPM job submitted"+console.processID());
					}
					else if(console.saveOutput(this))
					{	if(partnerCtrl!=null)
						{	partnerCtrl.windowClosing(null);
							partnerCtrl=null;
						}
						NairnFEAMPMViz.appCtrl.openNFMFileNow(console.getFile(),this,openMesh);
					}
				}
				else
					console.saveOutput(this);		// save in case needed
			}
			else
				console.saveOutput(this);			// save in case needed
			is.close();
			process.destroy();
		}
		catch(IOException tpe)
		{	Toolkit.getDefaultToolkit().beep();
			JOptionPane.showMessageDialog(this,tpe.toString());
		}
		running=false;
	}
	
	//----------------------------------------------------------------------------
	// Track which text field has focus (or null if neither
	//----------------------------------------------------------------------------
	
	public void focusGained(FocusEvent e)
	{	currentFocus=e.getComponent();
	}

	public void focusLost(FocusEvent e)
	{	if(!e.isTemporary())
			currentFocus=null;
	}

	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	public String getFullPath() { return cmdDoc.getFullPath(); }
	
	public boolean isCommandViewer() { return true; }
	
	public boolean isReuseable() { return !cmdField.getChanged() && cmdDoc.getFile()==null ;}

}	
