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
		
		// prepare output path and check on DTD validation
		boolean mpmAnalysis=isMPMAnalysis();
		boolean doValidate;
		String myCmd,myDTD;
		if(!NFMVPrefs.getRemoteMode())
		{	// running localling
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
		}
		else
		{	// REMOTE_ACCESS - running remotely
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
		
		// get output file or remote settings
		File outFile = null;
		if(!NFMVPrefs.getRemoteMode())
		{	File useOutput = scriptInfo!=null ? new File(scriptInfo.get(0)) : null ;
			if(mpmAnalysis)
			{	if(!soutConsole.setOutputPath(inFile,"mpm",useOutput)) return;
			}
			else
			{	if(!soutConsole.setOutputPath(inFile,"fea",useOutput)) return;
			}
		
			// get current output file from the output console object
			outFile = soutConsole.getFile();
			
			// write temporary file to the selected output folder
			tmpFile = saveCopyOfCommands(new File(outFile.getParent()+"/"+inFile.getName()));
		}
		else if(scriptInfo!=null)
		{	// REMOTE_ACCESS script mode
			
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
		{	// REMOTE_ACCESS
			outFile = soutConsole.getFile();
			
			// get remote file name
			String remotePath = soutConsole.getRemoteFilePath();
			if(remotePath==null)
			{	String outname = doc.getFile().getName();
				int eoff = outname.lastIndexOf(".");
				if(eoff>1) outname = outname.substring(0,eoff);
				if(mpmAnalysis)
					outname += ".mpm";
				else
					outname += ".fea";
				remotePath="RemoteOuput/"+outname;
			}
			
			// change extension id needed
			
			// get output folder from previous output or working directory
			if(outFile!=null)
			{	if(soutConsole.downloadResults == LaunchRemoteCalc.OPEN_ON_SERVER)
				{	outputFolder = soutConsole.getRemoteHomePath();
				}
				else
				{	try
					{	outputFolder = outFile.getParentFile().getParentFile().getPath();
					}
					catch(Exception e) {}
				}
			}
			if(outputFolder == null)
				outputFolder = NFMVPrefs.prefs.get(NFMVPrefs.WorkSpaceKey,NFMVPrefs.WorkSpaceKey);
			
			// dialog to get these
			LaunchRemoteCalc lrc = new LaunchRemoteCalc(doc,remotePath,soutConsole.uniqueOutput,
										outputFolder,soutConsole.downloadResults,soutConsole.clearPriorContents);
			if(lrc.getClickedButton()==JNDialog.CANCEL_BUTTON) return;
			
			// remote name
			remoteFolder = lrc.getRemoteFolder();
			remoteName = lrc.getRemoteFileName();
			if(mpmAnalysis)
			{	if(!remoteName.endsWith(".mpm")) remoteName += ".mpm";
			}
			else
			{	if(!remoteName.endsWith(".fea")) remoteName += ".fea";
			}
			sout.setRemoteFilePath(remoteFolder+"/"+remoteName);
			
			// is it unique?
			soutConsole.uniqueOutput = lrc.getMakeUnique();
			
			// is it to be deleted first
			soutConsole.clearPriorContents = lrc.getClearContents();
			
			// save folder
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
			
		// set commands and options
		ArrayList<String> cmds=new ArrayList<String>(20);
		
		// start login bash shell
		cmds.add(NFMVPrefs.prefs.get(NFMVPrefs.ShellKey,NFMVPrefs.ShellDef));
		cmds.add("--login");			// only needed for windows
		cmds.add("-c");
		
		// build shell command
		StringBuffer shell=new StringBuffer();
		
		// REMOTE_ACCESS - this used when not remote
		if(!NFMVPrefs.getRemoteMode())
		{	// running locally
			// cd to local file directory
			shell.append("cd ");
			String shellCD = outFile.getParent();
			if (NairnFEAMPMViz.isWindowsOS())
			{	// shellCD=shellCD.replace('\\','/');
				shellCD = PathToCygwin(shellCD);
			}
			if (shellCD.indexOf(' ') >= 0)
				shell.append("'" + shellCD + "'; ");
			else
				shell.append(shellCD + "; ");
		}

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
			shell.append(" -np "+processors); 
			
		String inName = tmpFile.getName();
		if(inName.indexOf(' ')>=0)
			shell.append(" '"+inName+"'");
		else
			shell.append(" "+inName);
		
		// start remote to local runnings
		openMesh=runType;
		
		if(NFMVPrefs.getRemoteMode())
		{	// REMOTE_ACCESS
			
			// add the shell command
			if(remoteName.indexOf(' ')>=0)
				shell.append(" > '"+remoteName+"'");
			else
				shell.append(" > "+remoteName);
			
			shellRemote = shell.toString();

			// start thread
			builder = null;
			runThread=new Thread(this);
			running=true;
			runThread.start();
		}
		
		else
		{	// running locally
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
			
			System.out.println(shell);
			cmds.add(shell.toString());
			
			// create the process and set working directory
			builder = new ProcessBuilder(cmds);
			//Map<String, String> environ = builder.environment();
			//builder.directory(myFile.getParentFile());
			builder.redirectErrorStream(true);
			
			runThread=new Thread(this);
			running=true;
			wasSubmitted=doBackground;
			runThread.start();
		}		
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
		
		// REMOTE_ACCESS - convert if not in root
		if(!NFMVPrefs.getRemoteMode())
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
		{	// local execution
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
							if(newResults!=null && openMesh==RUN_CHECK_MESH && !newResults.resDoc.is3D())
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
				JOptionPane.showMessageDialog(doc,tpe.getLocalizedMessage());
			}
		}
		
		else
		{	// Run remotely
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
					remoteConn = new RemoteConnection(remoteUser,userPass,remoteServer,22);
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

				// upload XML file
				String lastFolder=null,path = null;
				try
				{	soutConsole.appendLine("Uploading "+tmpFile.getName());
					lastFolder = remoteConn.uploadFile(tmpFile.getPath(), remoteFolder,
															soutConsole.uniqueOutput,soutConsole.clearPriorContents);
					
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

				// run the calculations remotely
				String command = "cd ";
				if(path.indexOf(" ")>0)
					command += "'"+path+"';";
				else
					command += path+";";
				command += shellRemote;
				soutConsole.appendLine("Running: "+command);
				
				int exitStatus = 0;
				try
				{	exitStatus = remoteConn.execCommands(command, false);
					soutConsole.appendLine("         Calculations done: exit status = "+exitStatus);
					if(exitStatus!=0)
						soutConsole.appendLine("...check system console for error message");
				}				
				catch (Exception e)
				{	throw new Exception("Remote execution failed: "+e.getLocalizedMessage());
				}
				
				// exit if asked
				if(!running) throw new Exception("exit");

				// create zipped file (if requested)
				String saveOutput = null;
				if(soutConsole.downloadResults==LaunchRemoteCalc.DOWNLOAD_TO_FOLDER && exitStatus==0)
				{	command = "cd ";
					if(path.indexOf(" ")>0)
						command += "'"+path+"';";
					else
						command += path+";";
					if(lastFolder.indexOf(" ")>0)
						command += "cd ..;rm -f '"+lastFolder+".zip';zip -r '"+lastFolder+".zip' '" +lastFolder+"'";
					else
						command += "cd ..;rm -f "+lastFolder+".zip;zip -r "+lastFolder+".zip " +lastFolder;
					soutConsole.appendLine("Zipping: "+command);
					
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
						{	soutConsole.appendLine("Downloading results to "+outputFolder);
							remoteConn.downloadExtractZip(path+".zip",outputFolder);
							saveOutput = outputFolder+File.separator+lastFolder+File.separator+remoteName;
						}
						catch (Exception e)
						{	throw new Exception("Failed to download the results: "+e.getLocalizedMessage());
						}
					}
				}
				else if(soutConsole.downloadResults==LaunchRemoteCalc.OPEN_ON_SERVER)
				{	// windows may need to parse path for file separators
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
				
				soutConsole.appendLine("Remote execution done");
				
				// if saved to mounted disk, open it now
				if(saveOutput!=null)
				{	soutConsole.setOutputPath(saveOutput);
					DocViewer newResults = doc.linkToResults();
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
