/*
 * NFMAnalysis.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 15 Dec 2012.
 * Copyright (c) 2012 RSAC Software. All rights reserved.
 */

import java.io.*;
import java.util.*;
import javax.swing.*;

import geditcom.JNFramework.*;

public class NFMAnalysis  implements Runnable
{
	protected CmdViewer doc;
	protected String cmds;
	protected ConsolePane soutConsole;
	protected boolean running=false;
	
	private ProcessBuilder builder;
	private Thread runThread;
	private int openMesh=FULL_ANALYSIS;
	private boolean wasSubmitted=false;
	
	private File tmpFile = null;
	
	public static final int INTERPRET_ONLY=0;
	public static final int RUN_CHECK_MESH=1;
	public static final int FULL_ANALYSIS=2;

	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public NFMAnalysis(CmdViewer cmdDoc)
	{	// save parent CmdViewer
		doc = cmdDoc;
	}
	
	//----------------------------------------------------------------------------
	// Interpret
	//----------------------------------------------------------------------------
	
	// Run FEA or MPM analysis
	public void runNFMAnalysis(boolean doBackground,int runType,String xmlData,ConsolePane sout,int processors)
	{
		// set the xml data
		cmds = xmlData;
		soutConsole = sout;
		
		// read command file name
		File inFile=doc.getFile();
		
		// prepare output path and check on DTD validation
		boolean mpmAnalysis=isMPMAnalysis();
		boolean doValidate;
		String myCmd,myDTD;
		if(mpmAnalysis)
		{	// get path to NairnFEA and NairnMPM
			myCmd=NFMVPrefs.prefs.get(NFMVPrefs.NairnMPMKey,NFMVPrefs.NairnMPMDef);
			myDTD=NFMVPrefs.prefs.get(NFMVPrefs.NairnMPMDTDKey,NFMVPrefs.NairnMPMDTDDef);
			doValidate=NFMVPrefs.prefs.getBoolean(NFMVPrefs.NairnMPMValidateKey,NFMVPrefs.NairnMPMValidateDef);
		}
		else
		{	// get path to NairnFEA and NairnMPM
			myCmd=NFMVPrefs.prefs.get(NFMVPrefs.NairnFEAKey,NFMVPrefs.NairnFEADef);
			myDTD=NFMVPrefs.prefs.get(NFMVPrefs.NairnFEADTDKey,NFMVPrefs.NairnFEADTDDef);
			doValidate=NFMVPrefs.prefs.getBoolean(NFMVPrefs.NairnFEAValidateKey,NFMVPrefs.NairnFEAValidateDef);
		}
		
		// insert DTD path if validating
		if(doValidate)
		{	if(!insertDTD(myDTD))
			{	JNApplication.appBeep();
				JOptionPane.showMessageDialog(doc,"The XML input commands do not start with the required <?xml ...?> element.");
				return;
			}
		}
		
		// just append XML commands
		if(runType==INTERPRET_ONLY)
		{	soutConsole.appendText(cmds);
			return;
		}
		
		// get output file
		if(mpmAnalysis)
		{	if(!soutConsole.setOutputPath(inFile,"mpm")) return;
		}
		else
		{	if(!soutConsole.setOutputPath(inFile,"fea")) return;
		}
		
		// get current output file from the output console object
		File outFile = soutConsole.getFile();
		
		// if not in the same folder, copy commands to that folder too
		tmpFile = saveCopyOfCommands(new File(outFile.getParent()+"/"+inFile.getName()));
		if(tmpFile==null) return;
			
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
		if(NairnFEAMPMViz.isWindowsOS())
		{	//shellCD=shellCD.replace('\\','/');
			shellCD=PathToCygwin(shellCD);
		}
		if(shellCD.indexOf(' ')>=0)
			shell.append("'"+shellCD+"'; ");
		else
			shell.append(shellCD+"; ");

		// executable (as shell command)
		if(NairnFEAMPMViz.isWindowsOS())
		{	//myCmd=myCmd.replace('\\','/');
			myCmd=PathToCygwin(myCmd);
		}
		if(myCmd.indexOf(' ')>=0)
			shell.append("'"+myCmd+"'");
		else
			shell.append(myCmd);
			
		// options (-a to abort after setup, -v to validate)
		if(runType==RUN_CHECK_MESH) shell.append(" -a");
		if(doValidate) shell.append(" -v");
		
		// processors
		if(processors>1)
		{	shell.append(" -np "+processors);
		}
			
		String inName = tmpFile.getName();
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
		openMesh=runType;
		wasSubmitted=doBackground;
		runThread.start();
	}
	
	// insert DTD path into commands
	public boolean insertDTD(String dtdPath)
	{
		int offset=cmds.indexOf("<!DOCTYPE"),endOffset=0;
		int docOffset=offset;
		if(offset>0)
		{	offset=cmds.indexOf("JANFEAInput",offset+9);
			if(offset>0)
			{	offset=cmds.indexOf("SYSTEM",offset+6);
			
				// if filed <!DOCTYPE JANFEAInput SYSEM located quoted path name
				if(offset>0)
				{	char startChar=' ';
					while(offset<cmds.length())
					{	startChar=cmds.charAt(offset);
						if(startChar=='\'' || startChar=='\"') break;
						if(startChar=='\n' || startChar=='\r')
						{	offset=0;
							break;
						}
						offset++;
					}
					if(offset>0)
					{	offset++;
						endOffset=offset;
						while(endOffset<cmds.length())
						{	char endChar=cmds.charAt(endOffset);
							if(endChar==startChar) break;
							if(endChar=='\n' || endChar=='\r')
							{	offset=0;
								break;
							}
							endOffset++;
						}
					}
				}
			}
		}
		dtdPath=PathToCygwin(dtdPath);
		if(offset>0)
		{	// insert path if needed
			String oldDTD=cmds.substring(offset,endOffset);
			if(!oldDTD.equals(dtdPath))
			{	cmds = cmds.substring(0,offset)+dtdPath+
							cmds.substring(endOffset,cmds.length());
			}
		}
		else if(docOffset>0)
		{	docOffset+=9;
			cmds = cmds.substring(0,docOffset)+" JANFEAInput SYSTEM '"
						+dtdPath+"'"+cmds.substring(docOffset,cmds.length());
		}
		else
		{	// insert <!DOCTYPE
			docOffset=cmds.indexOf("?>");
			if(docOffset<=0) return false;
			docOffset+=2;
			cmds = cmds.substring(0,docOffset)
						+"\n<!DOCTYPE JANFEAInput SYSTEM '"+dtdPath+"'>"
						+cmds.substring(docOffset,cmds.length());
		}
		
		return true;
	}

	// convert Windows path to cygwin
	private String PathToCygwin(String path)
	{	if(JNApplication.isWindowsOS())
		{	path=path.replace('\\','/');
			// Replace C: by /cygdrive/c (if needed)
			if(path.charAt(1)==':')
			{	char drive=path.charAt(0);
				if(drive<'a') drive+=32;
				path="/cygdrive/"+drive+path.substring(2);
			}
		}
		return path;
	}

	// save commands and return saved file (or null)
	public File saveCopyOfCommands(File origFile)
	{
		// find empty file
		String ext = JNUtilities.getExtension(origFile);
		String savePath = origFile.getPath();
		String root = savePath.substring(0, savePath.length()-ext.length());
		int num = 1;
		File saveFile = null;
		while(true)
		{	saveFile = new File(root+"_XML_"+num+ext);
			System.out.println(saveFile.getPath());
			if(!saveFile.exists()) break;
			num++;
		}
		
		// save to saveFile
		try
		{	FileWriter theFile=new FileWriter(saveFile);
			theFile.write(cmds);
			theFile.flush();
			theFile.close();
		}
		catch (Exception fe)
		{	JNApplication.appBeep();
			JOptionPane.showMessageDialog(doc,"Error writing XML input commands to temporary file: " + fe);
			return null;
		}
		
		return saveFile;
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
			
			soutConsole.clear();
			String line;
			while((line = br.readLine()) != null)
			{	soutConsole.appendLine(line);
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
					{	JOptionPane.showMessageDialog(doc,"FEA or MPM job submitted"+soutConsole.processID());
					}
					else if(soutConsole.saveOutput(doc))
					{	DocViewer newResults = doc.linkToResults();
						if(newResults!=null && openMesh==RUN_CHECK_MESH)
						{	newResults.checkMeshNow();
						}
					}
				}
				else
					soutConsole.saveOutput(doc);		// save in case needed
			}
			else
				soutConsole.saveOutput(doc);			// save in case needed
			is.close();
			process.destroy();
		}
		catch(IOException tpe)
		{	JNApplication.appBeep();
			JOptionPane.showMessageDialog(doc,tpe.toString());
		}
		// remove tmpFile
		if(tmpFile != null)
		{	try
			{	tmpFile.delete();
				tmpFile = null;
			}
			catch (Exception fe)
			{	JNApplication.appBeep();
				JOptionPane.showMessageDialog(doc,"Error deleting the temporary XML file: " + fe);
			}
		}
		running=false;
	}
	
	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	public boolean isRunning() { return running; }
	
	public void stopRunning()
	{	if(runThread.isAlive())
		{	running=false;
			try
			{	runThread.join();
			}
			catch(InterruptedException ie)
			{
			}
		}
	}

	// look for required MPMHeader element
	public boolean isMPMAnalysis()
	{	return cmds.indexOf("<MPMHeader>")>0;
	}
	

}
