/*
 * NFMAnalysis.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 15 Dec 2012.
 * Copyright (c) 2012 RSAC Software. All rights reserved.
 */

import geditcom.JNFramework.*;

import java.io.*;
import java.util.ArrayList;

import javax.swing.JOptionPane;

import com.jcraft.jsch.JSchException;


public class NFMAnalysis  implements Runnable
{
	protected CmdViewer doc;
	protected String cmds;
	protected ConsolePane soutConsole;
	protected boolean running=false;
	
	private ProcessBuilder builder = null;
	private Thread runThread;
	private int openMesh=FULL_ANALYSIS;
	private boolean wasSubmitted=false;
	
	private File tmpFile = null;
	
	String remoteFolder = null;
	String shellRemote = null;
	String outputFolder = null;
	String remoteName = null;
	
	public static final int INTERPRET_ONLY=0;
	public static final int RUN_CHECK_MESH=1;
	public static final int FULL_ANALYSIS=2;
	public static final int SCRIPT_ONLY=3;
	
	public static final int MAC_UNIX=0;
	public static final int WINDOWS_CYGWIN=1;
	public static final int WINDOWS_EXE=2;

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
	// if scriptInfo is not null, it has script params
	//		local run (0) is path to output folder
	// If null, use folder that contains the commands file or run in remote mode
	public void launchNFMAnalysis(boolean doBackground,int runType,String xmlData,ConsolePane sout,
					int processors,ArrayList<String>scriptInfo)
	{
		// set the xml data
		cmds = xmlData;
		soutConsole = sout;
		
		// read command file name
		File inFile=doc.getFile();
		System.out.println("Launching Analysis for file: "+inFile.getPath());
		
		// get shell command (if needed) and set command style
		int commandStyle = MAC_UNIX;
		String commandSep = ";";
		String pathDelim = "/";
		String bracket = "'";
		String bashPath = NFMVPrefs.prefs.get(NFMVPrefs.ShellKey,NFMVPrefs.ShellDef);
		if(NairnFEAMPMViz.isWindowsOS())
		{	if(bashPath.indexOf("$(windows)")<0)
				commandStyle = WINDOWS_CYGWIN;
			else
			{	commandStyle = WINDOWS_EXE;
				commandSep = " &";
				pathDelim = "\\";
				bracket = "\"";
			}
		}
		
		// trap background run in windows
		if(commandStyle==WINDOWS_EXE && doBackground)
		{	String msg = "Windows command line binary cannot yet be run in background.\n";
			msg = msg+"Use 'Run FEA/MPM Analysis' menu command instead.";
			JNApplication.appBeep();
			JOptionPane.showMessageDialog(doc,msg);
			return;
		}

		// prepare output path and check on DTD validation by finding
		// 		myCmd = path to executable
		//		myDTD - path to DTD file
		//		doValidate = if should validate input file
		boolean mpmAnalysis=isMPMAnalysis();
		boolean doValidate;
		String myCmd,myDTD;
		if(!NFMVPrefs.getRemoteMode())
		{	// LOCAL_EXECUTION ----------------------
			if(mpmAnalysis)
			{	// get path to NairnFEA and NairnMPM
				myCmd=NFMVPrefs.prefs.get(NFMVPrefs.NairnMPMKey,NFMVPrefs.NairnMPMDef);
				myDTD=NFMVPrefs.prefs.get(NFMVPrefs.NairnMPMDTDKey,NFMVPrefs.NairnMPMDTDDef);
				doValidate=NFMVPrefs.prefs.getBoolean(NFMVPrefs.NairnMPMValidateKey,NFMVPrefs.NairnMPMValidateDef);
				if(myCmd.indexOf("$(bundle)")>=0)
					myCmd=NairnFEAMPMViz.jarFolder+"bundle"+pathDelim+"NairnMPM.exe";
				if(myDTD.indexOf("$(bundle)")>=0)
					myDTD=NairnFEAMPMViz.jarFolder+"bundle"+pathDelim+"NairnMPM.dtd";
			}
			else
			{	// get path to NairnFEA and NairnMPM
				myCmd=NFMVPrefs.prefs.get(NFMVPrefs.NairnFEAKey,NFMVPrefs.NairnFEADef);
				myDTD=NFMVPrefs.prefs.get(NFMVPrefs.NairnFEADTDKey,NFMVPrefs.NairnFEADTDDef);
				doValidate=NFMVPrefs.prefs.getBoolean(NFMVPrefs.NairnFEAValidateKey,NFMVPrefs.NairnFEAValidateDef);
				if(myCmd.indexOf("$(bundle)")>=0)
					myCmd=NairnFEAMPMViz.jarFolder+"bundle"+pathDelim+"NairnFEA.exe";
				if(myDTD.indexOf("$(bundle)")>=0)
					myDTD=NairnFEAMPMViz.jarFolder+"bundle"+pathDelim+"NairnFEA.dtd";
			}
		}
		else
		{	// REMOTE_ACCESS ---------------------
			// find variables and also check for entered user name and password
			if(mpmAnalysis) 
			{	myCmd = NFMVPrefs.prefs.get(NFMVPrefs.RemoteMPMPathKey,NFMVPrefs.RemoteMPMPathDef);
				myDTD = NFMVPrefs.prefs.get(NFMVPrefs.RemoteMPMDTDKey,NFMVPrefs.RemoteMPMDTDDef);
				doValidate = NFMVPrefs.prefs.getBoolean(NFMVPrefs.NairnMPMValidateKey,NFMVPrefs.NairnMPMValidateDef);
			} 
			else
			{	myCmd = NFMVPrefs.prefs.get(NFMVPrefs.RemoteFEAPathKey,NFMVPrefs.RemoteFEAPathDef);
				myDTD = NFMVPrefs.prefs.get(NFMVPrefs.RemoteFEADTDKey,NFMVPrefs.RemoteFEADTDDef);
				doValidate = NFMVPrefs.prefs.getBoolean(NFMVPrefs.NairnFEAValidateKey,NFMVPrefs.NairnFEAValidateDef);
			}
			
			// make sure user name and server have been entered
			String remoteUser = NFMVPrefs.prefs.get(NFMVPrefs.RemoteUserKey,NFMVPrefs.RemoteUserDef);
			String remoteServer = NFMVPrefs.prefs.get(NFMVPrefs.RemoteServerKey,NFMVPrefs.RemoteServerDef);
			if(remoteUser.length()==0 || remoteServer.length()==0)
			{	JNApplication.appBeep();
				JOptionPane.showMessageDialog(null,"You must enter server and user names in application preferences to be able to run in remote mode.");
				return;
			}
		}
				
		// insert DTD path if validating
		if(doValidate)
		{	if(!insertDTD(myDTD,commandStyle))
			{	JNApplication.appBeep();
				JOptionPane.showMessageDialog(doc,"The XML input commands do not start with the required <?xml ...?> element.");
				return;
			}
		}
		
		// just append XML commands and exit
		if(runType==INTERPRET_ONLY)
		{	soutConsole.appendText(cmds);
			return;
		}
		
		// get output file or remote settings for output files
		File outFile = null;
		if(!NFMVPrefs.getRemoteMode())
		{	// LOCAL_EXECUTION ----------------------
			File useOutput = scriptInfo!=null ? new File(scriptInfo.get(0)) : null ;
			if(mpmAnalysis)
			{	if(!soutConsole.setOutputPath(inFile,"mpm",useOutput)) return;
			}
			else
			{	if(!soutConsole.setOutputPath(inFile,"fea",useOutput)) return;
			}
		
			// get current output file from the output console object
			outFile = soutConsole.getFile();
			System.out.println("...output: "+outFile.getPath());
			
			// write temporary file to the selected output folder
			String tmpName = outFile.getParent()+pathDelim+inFile.getName();
			tmpFile = saveCopyOfCommands(new File(tmpName));
		}
		else if(scriptInfo!=null)
		{	// REMOTE_ACCESS script mode -----------------
			
			// remote name (script code already checked to "/" so must be there
			String path = scriptInfo.get(0);
			int offset = path.lastIndexOf("/");
			remoteFolder = path.substring(0, offset);
			remoteName = path.substring(offset+1);
			if(mpmAnalysis)
			{	if(!remoteName.endsWith(".mpm")) remoteName += ".mpm";
			}
			else
			{	if(!remoteName.endsWith(".fea")) remoteName += ".fea";
			}
			sout.setRemoteFilePath(remoteFolder+"/"+remoteName);
			
			// is it unique or to be cleared? (assumed valid option)
			String remoteOption = scriptInfo.get(1);
			soutConsole.uniqueOutput = false;
			soutConsole.clearPriorContents = false;
			if(remoteOption.equals("unique"))
				soutConsole.uniqueOutput = true;
			else if(remoteOption.equals("clear"))
				soutConsole.clearPriorContents = true;
			
			// save options (assumed valid)
			outputFolder = scriptInfo.get(2);
			String outputOption = scriptInfo.get(3);
			if(outputOption.equals("download"))
				soutConsole.downloadResults = LaunchRemoteCalc.DOWNLOAD_TO_FOLDER;
			else if(outputOption.equals("nodownload"))
			{	soutConsole.downloadResults = LaunchRemoteCalc.NO_DOWNLOAD;
				outputFolder = "";
			}
			else
				soutConsole.downloadResults = LaunchRemoteCalc.OPEN_ON_SERVER;
			
			// write temporary file to the selected input file folder
			tmpFile = saveCopyOfCommands(new File(inFile.getPath()));
		}
		else
		{	// REMOTE_ACCESS --------------------------------
			
			// read current settings, then change in dialog
			outFile = soutConsole.getFile();
			
			// get remote file name
			String remotePath = soutConsole.getRemoteFilePath();
			if(remotePath==null)
			{	String outname = doc.getFile().getName();
			
				// remove extension and append one needed now (may be the same)
				int eoff = outname.lastIndexOf(".");
				if(eoff>1) outname = outname.substring(0,eoff);
				if(mpmAnalysis)
					outname += ".mpm";
				else
					outname += ".fea";
				remotePath="RemoteOutput/"+outname;
			}
			
			// get output folder from previous output or working directory
			if(outFile!=null)
			{	if(soutConsole.downloadResults == LaunchRemoteCalc.OPEN_ON_SERVER)
				{	// if opened on server, look for save path to home directory on the server
					outputFolder = soutConsole.getRemoteHomePath();
				}
				else
				{	// back up one folder and one zipped folder to original parent folder
					try
					{	outputFolder = outFile.getParentFile().getParentFile().getPath();
					}
					catch(Exception e) {}
				}
			}
			if(outputFolder == null)
				outputFolder = NFMVPrefs.prefs.get(NFMVPrefs.WorkSpaceKey,NFMVPrefs.WorkSpaceKey);
			
			// Run dialog to get items listed after the dialo
			LaunchRemoteCalc lrc = new LaunchRemoteCalc(doc,remotePath,soutConsole.uniqueOutput,
										outputFolder,soutConsole.downloadResults,soutConsole.clearPriorContents);
			if(lrc.getClickedButton()==JNDialog.CANCEL_BUTTON) return;
			
			// --------- Remote folder and name. These were separated in the dialog.
			// remoteFolder is folder (relative to home directory) and name and the file name,
			// which is last component of the enerter remote path
			remoteFolder = lrc.getRemoteFolder();
			remoteName = lrc.getRemoteFileName();
			
			// make sure extension is correct
			if(mpmAnalysis)
			{	if(!remoteName.endsWith(".mpm")) remoteName += ".mpm";
			}
			else
			{	if(!remoteName.endsWith(".fea")) remoteName += ".fea";
			}
			
			// save full path to remote file relative to home directory
			sout.setRemoteFilePath(remoteFolder+"/"+remoteName);
			
			// -------Is it unique? (i.e. Create Unique Subfolder box checked)
			soutConsole.uniqueOutput = lrc.getMakeUnique();
			
			// ------- Is it to be deleted first? (i.e. Clear Parent Folder First checked)
			soutConsole.clearPriorContents = lrc.getClearContents();
			
			// -------- Where to save the results (i.e., download zip, do not download, open on server)
			// If download save local folder in outputFolder or if open on server store home folder in outputFolder
			soutConsole.downloadResults = lrc.getDoDownload();
			if(soutConsole.downloadResults != LaunchRemoteCalc.NO_DOWNLOAD)
				outputFolder = lrc.getLocalFolder();
			else
				outputFolder = "";
			
			// write temporary file to the selected input file folder
			tmpFile = saveCopyOfCommands(new File(inFile.getPath()));
		}
		
		// if did not write to file, then exit
		if(tmpFile==null) return;
		System.out.println("...tmp commands: "+tmpFile.getPath());
			
		// set commands and options
		// bashcmds - list of commands to launch bash shell for cygwin or mac in background
		ArrayList<String> bashcmds=new ArrayList<String>(20);
		// pbcmds for process building
		ArrayList<String> pbcmds=new ArrayList<String>(20);
		
		// start login bash shell (only used in cygwin or mac in backgorund)
		if(commandStyle!=WINDOWS_EXE)
		{	bashcmds.add(NFMVPrefs.prefs.get(NFMVPrefs.ShellKey,NFMVPrefs.ShellDef));
			bashcmds.add("--login");
			bashcmds.add("-c");
		}
		
		// build shell command which will be
		// cd (parent folder); (executable)
		StringBuffer shell=new StringBuffer();
		
		// Get command to go to the parent folder directory
		if(!NFMVPrefs.getRemoteMode())
		{	// LOCAL_EXECUTION
			String shellCD = outFile.getParent();
			if(commandStyle==WINDOWS_EXE)
			{	if(shellCD.indexOf(':')==1)
				{	shell.append(shellCD.substring(0,2)+commandSep+" ");
				}
				shell.append("CD ");
				
			}
			else
			{	if(commandStyle==WINDOWS_CYGWIN)
					shellCD = PathToCygwin(shellCD);
				shell.append("cd ");
			}
			if(shellCD.indexOf(' ') >= 0)
				shell.append(bracket + shellCD + bracket);
			else
				shell.append(shellCD);
			shell.append(commandSep+" ");
		}

		// executable (Mac and Exe to process builder, all to shell)
		// pbcmds will be [(cmd),[options],(input)]
		// shell will be cd (parent) ; (cmd) (options) (input)
		if(commandStyle==WINDOWS_CYGWIN)
		{	myCmd=PathToCygwin(myCmd);
		}
		else
		{	// as direct process builder command
			pbcmds.add(myCmd);
			if(runType==RUN_CHECK_MESH) pbcmds.add("-a");
			if(doValidate) pbcmds.add("-v");
			if(processors>1)
			{	pbcmds.add("-np");
				pbcmds.add(""+processors);
			}
		}
		if(myCmd.indexOf(' ')>=0)
			shell.append(bracket+myCmd+bracket);
		else
			shell.append(myCmd);
					
		// shell options (-a to abort after setup, -v to validate, -np processors)
		if(runType==RUN_CHECK_MESH) shell.append(" -a");
		if(doValidate) shell.append(" -v");
		if(processors>1) shell.append(" -np "+processors); 
		
		// input file name to shell command
		String inName = tmpFile.getName();
		if(inName.indexOf(' ')>=0)
			shell.append(" "+bracket+inName+bracket);
		else
			shell.append(" "+inName);
		
		// and to direct file name command
		if(commandStyle!=WINDOWS_CYGWIN)
			pbcmds.add(inName);
		
		// Finish building commands then launch thread (see run() method)
		openMesh=runType;
		
		if(NFMVPrefs.getRemoteMode())
		{	// REMOTE_ACCESS ---------------------
			
			// add output file name
			if(remoteName.indexOf(' ')>=0)
				shell.append(" > '"+remoteName+"'");
			else
				shell.append(" > "+remoteName);
			
			// convert to string of commands
			shellRemote = shell.toString();

			// start thread
			builder = null;
			runThread=new Thread(this);
			running=true;
			runThread.start();
		}
		
		else
		{	// LOCAL_EXECUTION
			if(doBackground)
			{	shell.append(" >& ");
				
				// append output file
				String outFileName=outFile.getName();
				if(outFileName.indexOf(' ')>=0)
					shell.append("'"+outFileName+"'");
				else
					shell.append(outFileName);
				
				// append background command
				shell.append(" &");
			}
			
			// add shell command to the bash shell login
			bashcmds.add(shell.toString());
			
			// create the process and set working directory
			if(commandStyle!=WINDOWS_CYGWIN && !doBackground)
			{	builder = new ProcessBuilder(pbcmds);
				builder.directory(outFile.getParentFile());
				displayCommands(pbcmds);
			}
			else
			{	builder = new ProcessBuilder(bashcmds);
				displayCommands(bashcmds);
			}
			//builder.redirectErrorStream(true);
			
			runThread=new Thread(this);
			running=true;
			wasSubmitted=doBackground;
			runThread.start();
		}		
	}
	
	// display commands
	public void displayCommands(ArrayList<String> cmds)
	{	System.out.println("Launch Commands:");
		for(int i=0;i<cmds.size();i++)
			System.out.println("  "+i+": "+cmds.get(i));
	}
	
	// insert DTD path into commands
	public boolean insertDTD(String dtdPath,int commandStyle)
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
		
		// REMOTE_ACCESS - convert if not in root
		if(!NFMVPrefs.getRemoteMode() && commandStyle==WINDOWS_CYGWIN)
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
			if(!saveFile.exists()) break;
			num++;
		}
		
		// save to saveFile
		try
		{	FileWriter theFile=new FileWriter(saveFile);
			theFile.write(cmds);
			theFile.flush();
			theFile.close();
			theFile = null;
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
		if(builder!=null)
		{	// LOCAL_EXECUTION ---------------
			try
			{	Process process = builder.start();
				InputStream is = process.getInputStream();
				InputStreamReader isr = new InputStreamReader(is);
				BufferedReader br = new BufferedReader(isr);
				
				InputStream es = process.getErrorStream();
				InputStreamReader esr = new InputStreamReader(es);
				BufferedReader ebr = new BufferedReader(esr);

				soutConsole.clear();
				String line;
				String errMsg = "";
				//boolean isFinished = false;
				//boolean isErrorFinished = false;
				
				// read standard output
				while((line = br.readLine()) != null)
				{	soutConsole.appendLine(line);
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
				
				// get results code
				int result=0;
				if(running)
				{	try
					{	process.waitFor();
						result = process.exitValue();
					}
					catch(InterruptedException e)
					{	errMsg = e.getLocalizedMessage();
						result = 1;
					}
				}
				
				// save and open for visualization
				if(wasSubmitted)
				{	JOptionPane.showMessageDialog(doc,"FEA or MPM job submitted"+soutConsole.processID());
				}
				else if(soutConsole.saveOutput(doc))
				{	if(result==0)
					{	DocViewer newResults = doc.linkToResults();
						if(newResults!=null && openMesh==RUN_CHECK_MESH && !newResults.resDoc.is3D())
						{	newResults.checkMeshNow();
						}
					}
					else
					{	if(errMsg=="") errMsg = "Unknown error (perhaps Windows exe failed to launch)";
						JNApplication.appBeep();
						JOptionPane.showMessageDialog(doc,"Error: "+errMsg);
					}
				}
				
				// close all
				is.close();
				es.close();
				process.destroy();
			}
			catch(Exception tpe)
			{	JNApplication.appBeep();
				JOptionPane.showMessageDialog(doc,tpe.getLocalizedMessage());
			}
		}
		
		else
		{	// REMOVE_ACCESS ----------------------
			
			// connect to server
			RemoteConnection remoteConn = null;
			try
			{
				// this were verified as entered earlier
				String remoteUser = NFMVPrefs.prefs.get(NFMVPrefs.RemoteUserKey,NFMVPrefs.RemoteUserDef);
				String remoteServer = NFMVPrefs.prefs.get(NFMVPrefs.RemoteServerKey,NFMVPrefs.RemoteServerDef);
				
				// create session
				soutConsole.clear();
				try
				{	soutConsole.appendLine("Connecting to "+remoteServer);
					String userPass = NFMVPrefs.prefs.get(NFMVPrefs.RemoteUserPassKey,NFMVPrefs.RemoteUserPassDef);
					remoteConn = new RemoteConnection(remoteUser,userPass,remoteServer,22,soutConsole);
					remoteConn.setStrictHostKeyChecking(false);
				}
				catch(JSchException je)
				{	remoteConn = null;			// no need to disconnect
					throw new Exception("Failed to create session: "+je.getLocalizedMessage());
				}
				catch(Exception e)
				{	remoteConn = null;			// no need to disconnect
					throw new Exception("exit");
				}
				
				// exit if asked
				if(!running)
				{	remoteConn = null;
					throw new Exception("exit");
				}
				
				// connect to server
				try
				{	if(!remoteConn.connect(30000))
						throw new Exception("timeout error");
				}
				catch(Exception e)
				{	remoteConn = null;			// no need to disconnect
					throw new Exception("Failed to connect: "+e.getLocalizedMessage());
				}
				
				// exit if asked
				if(!running) throw new Exception("exit");

				// upload XML commands file
				String lastFolder=null,path = null;
				try
				{	soutConsole.appendLine("Uploading "+tmpFile.getName());
				
					// upload, return lastFolder is the folder that was used
					lastFolder = remoteConn.uploadFile(tmpFile.getPath(), remoteFolder,
															soutConsole.uniqueOutput,soutConsole.clearPriorContents);
					
					// get path to output folder (needs to update if making unique)
					if(soutConsole.uniqueOutput)
					{	if(remoteFolder.charAt(remoteFolder.length()-1)!='/')
							path = remoteFolder+"/"+lastFolder;
						else
							path = remoteFolder+lastFolder;
					}
					else
						path = remoteFolder;
				}
				catch(Exception e)
				{	throw new Exception("Failed to upload file: "+e.getLocalizedMessage());
				}
				
				// exit if asked
				if(!running) throw new Exception("exit");

				// build command to run the calculations remotely
				String command = "cd ";
				if(path.indexOf(" ")>0)
					command += "'"+path+"';";
				else
					command += path+";";
				command += shellRemote;
				soutConsole.appendLine("Running: "+command);
				
				// execute the command
				int exitStatus = 0;
				try
				{	exitStatus = remoteConn.execCommands(command, false);
					soutConsole.appendLine("         Calculations done: exit status = "+exitStatus);
					// Matt Change - errors not reported in RemoteConnection
					//if(exitStatus!=0)
					//	soutConsole.appendLine("...check system console for error message");
				}				
				catch (Exception e)
				{	throw new Exception("Remote execution failed: "+e.getLocalizedMessage());
				}
				
				// exit if asked
				if(!running) throw new Exception("exit");

				// decide what to do with final output
				String saveOutput = null;
				if(soutConsole.downloadResults==LaunchRemoteCalc.DOWNLOAD_TO_FOLDER && exitStatus==0)
				{	// zip the folder of results (cd to path and zip the folder) and then download it
					command = "cd ";
					if(path.indexOf(" ")>0)
						command += "'"+path+"';";
					else
						command += path+";";
					
					// remove previous zip and then create new one
					if(lastFolder.indexOf(" ")>0)
						command += "cd ..;rm -f '"+lastFolder+".zip';zip -r '"+lastFolder+".zip' '" +lastFolder+"'";
					else
						command += "cd ..;rm -f "+lastFolder+".zip;zip -r "+lastFolder+".zip " +lastFolder;
					soutConsole.appendLine("Zipping: "+command);
					
					// execute the zip commands
					try
					{	exitStatus = remoteConn.execCommands(command, false);
						soutConsole.appendLine("         Zip done: exit status = "+exitStatus);
						if(exitStatus!=0)
							soutConsole.appendLine("...check system console for error message");
					}				
					catch (Exception e)
					{	throw new Exception("Remote zipping failed: "+e.getLocalizedMessage());
					}
					
					// exit if asked
					if(!running) throw new Exception("exit");
					
					// download zipped file (if saved)
					if(exitStatus==0)
					{	try
						{	//soutConsole.appendLine("Downloading results to "+outputFolder);
							remoteConn.downloadExtractZip(path+".zip",outputFolder);
							saveOutput = outputFolder+File.separator+lastFolder+File.separator+remoteName;
						}
						catch (Exception e)
						{	throw new Exception("Failed to download the results: "+e.getLocalizedMessage());
						}
					}
				}
				else if(soutConsole.downloadResults==LaunchRemoteCalc.OPEN_ON_SERVER)
				{	// open on server by assuming outputFolder is path to mounted home directory on server
					// save in the console pane class
					// windows may need to parse path for file separators
					saveOutput = outputFolder+File.separator + path + File.separator + remoteName;
					soutConsole.setRemoteHomePath(outputFolder);
				}
				
				// disconnect
				try
				{	soutConsole.appendLine("Disconnecting");
					remoteConn.disconnect();
				}
				catch (Exception e)
				{	remoteConn = null;			// no need to disconnect again
					throw new Exception("Failed to disconnect: "+e.getLocalizedMessage());
				}
				
				// all done
				soutConsole.appendLine("Remote execution done");
				
				// if saved to mounted disk (saveOutput is not null), open it now
				if(saveOutput!=null)
				{	// save path to console and tell doc to open it
					soutConsole.setOutputPath(saveOutput);
					soutConsole.appendLine("Open results at "+saveOutput);
					DocViewer newResults = doc.linkToResults();
					
					// if check mesh, show mesh
					if(newResults!=null && openMesh==RUN_CHECK_MESH && !newResults.resDoc.is3D())
					{	newResults.checkMeshNow();
					}			
				}
				
			}
			catch(Exception reme)
			{	if(remoteConn!=null) remoteConn.disconnect();
				String msg = reme.getLocalizedMessage();
				if(msg.equals("exit"))
					soutConsole.appendLine("Remote execution cancelled");
				else
				{	JNApplication.appBeep();
					JOptionPane.showMessageDialog(doc,msg);
				}
			}
		}
		
		// remove tmpFile
		if(tmpFile != null)
		{	try
			{	if(!tmpFile.delete())
				{	Thread.sleep(250);
					if(!tmpFile.delete())
						System.out.println("Failed to delete temporary file '"+tmpFile.getName()+"'");
				}
			}
			catch (Exception fe)
			{	JNApplication.appBeep();
				JOptionPane.showMessageDialog(doc,"Error deleting the temporary XML file: " + fe.getLocalizedMessage());
			}
			tmpFile = null;
		}
		
		// all done
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
			{	// this waits until the thres is done
				runThread.join();
			}
			catch(InterruptedException ie)
			{	System.out.println("Failed to stop thread");
			}
		}
	}

	// look for required MPMHeader element
	public boolean isMPMAnalysis()
	{	return cmds.indexOf("<MPMHeader>")>0;
	}
	
}
