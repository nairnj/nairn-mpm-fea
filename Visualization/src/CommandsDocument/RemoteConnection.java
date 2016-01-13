/*
 * RemoteConnection.java
 * NairnFEAMPMViz
 * 
 * Created by Matt Viehdorfer, Nov 2014.
 * 
 * REMOTE_ACCESS - this entire class
 */


import geditcom.JNFramework.JNApplication;

import java.awt.GridLayout;
import java.io.*;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.event.AncestorEvent;
import javax.swing.event.AncestorListener;

import com.jcraft.jsch.*;

public class RemoteConnection {
	private JSch jsch;
	public Session session;
	public Boolean sessionOpen;
	private UserInfo ui;
	private ConsolePane soutConsole;

	public RemoteConnection(String user, String userPass,String host, int port, ConsolePane sout) throws Exception
	{	jsch = new JSch();
		session = jsch.getSession(user, host, port);
		generateUI();
		getUserPass(user,host,userPass);
		soutConsole = sout;
	}

	public boolean connect() throws Exception
	{	session.connect();
		return testConnection();
	}

	public boolean connect(int timeoutDuration) throws Exception
	{	session.connect(timeoutDuration);
		return testConnection();
	}

	public void disconnect() {
		if (testConnection())
			session.disconnect();
		else
			soutConsole.appendLine("Error: No open sessions to disconnect.");
	}

	public void setStrictHostKeyChecking(boolean enable) {
		try {
			String result = enable ? "yes" : "no"; // since the api requires a string...
			session.setConfig("StrictHostKeyChecking", result);
		} catch (Exception e) {
			soutConsole.appendLine("Error: "+e);
		}
	}

	public void setKnownHosts(String hostsFile) {
		try {
			jsch.setKnownHosts(hostsFile);
		} catch (Exception e) {
			soutConsole.appendLine("Error: "+e);
		}
	}

	public void writeToRemote(String fileContents, String directory, String fileName) {
		try {
			Channel channel = session.openChannel("sftp");
			channel.connect();
			ChannelSftp c = (ChannelSftp) channel;
			InputStream obj_InputStream = new ByteArrayInputStream(
					fileContents.getBytes());
			c.put(obj_InputStream, directory + fileName);
			c.exit();
			channel.disconnect();
		} catch (Exception e) {
			soutConsole.appendLine("Error: "+e);
		}
	}

	public String readFromRemote(String directory, String fileName) {
		StringBuilder fileContents = new StringBuilder();
		if (directory == null) 
			directory = "";
		
		try {
			Channel channel = session.openChannel("sftp");
			channel.connect();
			ChannelSftp c = (ChannelSftp) channel;
			if (directory.length() > 0)
				c.cd(directory);
			InputStream inputStream = c.get(fileName);
			char[] ch_Buffer = new char[0x10000];
			Reader reader = new InputStreamReader(inputStream, "UTF-8");
			int int_Line = 0;
			do {
				int_Line = reader.read(ch_Buffer, 0, ch_Buffer.length);
				if (int_Line > 0) {
					fileContents.append(ch_Buffer, 0, int_Line);
				}
			} while (int_Line >= 0);
			reader.close();
			inputStream.close();
			c.exit();
			channel.disconnect();
		} catch (Exception e) {
			soutConsole.appendLine("Error: "+e);
		}

		return fileContents.toString();
	}

	// upload local file to remote directory on the server
	// the upload will use name of local file and overwrite a file if it is there
	// return file folder name in the remote destination
	public String uploadFile(String localFilePath, String remoteDestinationDir,
					boolean makeUnique,boolean clearContents) throws Exception
	{	Channel channel = session.openChannel("sftp");
		channel.connect();
		ChannelSftp c = (ChannelSftp) channel;
		// implied "./" prefix on path and if starts with "/" it is ignored
		String[] folders = remoteDestinationDir.split("/");
		String lastFolder = "";
		boolean lastExists = true;
		for (String folder : folders) {
			if (folder.length() > 0) {
				try {
					c.cd(folder);
				} catch (SftpException e) {
					c.mkdir(folder);
					c.cd(folder);
					lastExists = false;
				}
				lastFolder = folder;
			}
		}
		
		// create a new sub folder is requiring unique
		if(makeUnique)
		{	int num=1;
			while(true)
			{	lastFolder = "out"+num;
				try
				{	c.cd(lastFolder);
					c.cd("..");
					num++;
				}
				catch (SftpException e)
				{	c.mkdir(lastFolder);
					c.cd(lastFolder);
					break;
				}
			}
		}
		else if(clearContents && lastExists && lastFolder.length()>0)
		{	// if requested, delete last folder and create it again
			c.cd("..");
			String cmd;
			if(lastFolder.indexOf(" ")>0)
				cmd = "rm -r '"+lastFolder+"'";
			else
				cmd = "rm -r "+lastFolder;
			execCommands(cmd,false);
			c.mkdir(lastFolder);
			c.cd(lastFolder);
		}

		File f1 = new File(localFilePath);
		c.put(new FileInputStream(f1), f1.getName(), ChannelSftp.OVERWRITE);

		c.exit();
		channel.disconnect();
		
		return lastFolder;
	}
	
	public void downloadExtractZip(String remoteFilePath, String localDestinationDir) throws Exception
	{	Channel channel = session.openChannel("sftp");
		channel.connect();
		ChannelSftp c = (ChannelSftp) channel;
		c.lcd(localDestinationDir);
		
		File tempFile = new File(remoteFilePath);	
		//c.get(remoteFilePath, localDestinationDir, new ProgressModel());
		SftpProgressMonitor monitor = new SystemOutProgressMonitor();
		c.get(remoteFilePath, localDestinationDir, monitor, ChannelSftp.OVERWRITE);
			
		if(localDestinationDir.charAt(localDestinationDir.length()-1)!=File.separatorChar){
			localDestinationDir += File.separator;
		}
		
		String zipFile = localDestinationDir+tempFile.getName();		
		unzipDownloadedFile(zipFile);
		removeFile(zipFile);

		c.exit();
		channel.disconnect();
	}

	public int execCommands(String commands, boolean showOutput) throws Exception
	{	
		PrintStream stderr = System.err;
		ByteArrayOutputStream allOutput = new ByteArrayOutputStream();
		PrintStream err = new PrintStream(allOutput);	
		System.setErr(err);
			
		Channel channel = session.openChannel("exec");
		((ChannelExec) channel).setCommand(commands);
		channel.setInputStream(null);
		((ChannelExec) channel).setErrStream(System.err);
		
		InputStream in = channel.getInputStream();
		channel.connect();
		printCommandOutput(in, channel, showOutput);
		int exitStatus = channel.getExitStatus();
		channel.disconnect();
		
		if (allOutput.size() > 0) {
			soutConsole.appendLine("Error executing remote command: ");
			soutConsole.appendLine(allOutput.toString());
		}
		System.setErr(stderr);
		
		return exitStatus;
	}
	
	private boolean testConnection() {
		return session.isConnected();
	}
	
	//Recursively decompress zip contents into a folder of the same name in the same dir
	private void unzipDownloadedFile(String zipFile) throws ZipException, IOException{
	    int BUFFER = 2048;
	    
	    // get zip file
	    File file = new File(zipFile);
	    ZipFile zip = new ZipFile(file);
	    
	    // unzip to directory of the same name
	    // When sip is a directory, this gets two folders
	    //String newPath = zipFile.substring(0, zipFile.length() - 4);
	    //new File(newPath).mkdir();
	    
	    // unzip to parent directory of the zip file
	    // This is assuming zip if of a directory
	    String newPath = file.getParent();
	    
	    // Process each entry
	    Enumeration<? extends ZipEntry> zipFileEntries = zip.entries();
	    while (zipFileEntries.hasMoreElements())
	    {
	        // grab a zip file entry
	        ZipEntry entry = (ZipEntry) zipFileEntries.nextElement();
	        String currentEntry = entry.getName();
	        File destFile = new File(newPath, currentEntry);
	        File destinationParent = destFile.getParentFile();

	        // create the parent directory structure if needed
	        destinationParent.mkdirs();

	        if (!entry.isDirectory())
	        {
	            BufferedInputStream is = new BufferedInputStream(zip.getInputStream(entry));
	            int currentByte;
	            // establish buffer for writing file
	            byte data[] = new byte[BUFFER];

	            // write the current file to disk
	            FileOutputStream fos = new FileOutputStream(destFile);
	            BufferedOutputStream dest = new BufferedOutputStream(fos,
	            BUFFER);

	            // read and write until last byte is encountered
	            while ((currentByte = is.read(data, 0, BUFFER)) != -1) {
	                dest.write(data, 0, currentByte);
	            }
	            dest.flush();
	            dest.close();
	            is.close();
	        }

	        if (currentEntry.endsWith(".zip"))
	        {
	            // found a zip file, try to open
	        	unzipDownloadedFile(destFile.getAbsolutePath());
	        }
	    }
	    zip.close();
	}
	
	public void openFile(String directory, String fileName) {
		File dir = new File(directory);

		if (dir.isDirectory()) {
			search(dir, fileName);
		} else {
			soutConsole.appendLine("Error: " + dir.getAbsoluteFile() + " is not a directory");
		}
	}

	private void search(File file, String fileName) {
		if (file.isDirectory()) {
			// do you have permission to read this directory?
			if (file.canRead()) {
				for (File temp : file.listFiles()) {
					if (temp.isDirectory()) {
						search(temp, fileName);
					} else {
						if (fileName.equals(temp.getName().toLowerCase())) {
							//result.add(temp.getAbsoluteFile().toString());
							JNApplication.main.openDocument(temp);
						}

					}
				}

			} else {
				soutConsole.appendLine("Error: "+ file.getAbsoluteFile() + "Permission Denied");
			}
		}
	}
	
	private void removeFile(String filePath) {
		File delFile = new File(filePath);
		try {
		    delFile.delete();
		} catch (Exception x) {
		    soutConsole.appendLine("Error deleting file: "+ filePath);
		} 
	}

	// monitor exec command and optionally dump any output
	private void printCommandOutput(InputStream in, Channel channel,boolean showOutput) throws Exception
	{
		byte[] tmp = new byte[1024];
		while(true)
		{	while(in.available() > 0)
			{	int i = in.read(tmp, 0, 1024);
				if (i < 0) break;
				if(showOutput) soutConsole.appendLine(new String(tmp, 0, i));
			}
			if (channel.isClosed())
			{	if(showOutput)
					soutConsole.appendLine("exit-status: "+ channel.getExitStatus());
				break;
			}
			try
			{	Thread.sleep(1000);
			}
			catch (Exception e) {
				soutConsole.appendLine("Error: "+e.getMessage());
			}
		}
	}

	private void generateUI() {
		ui = new MyUserInfo() {
			public void showMessage(String message) {
				JOptionPane.showMessageDialog(null, message);
			}

			public boolean promptYesNo(String message) {
				Object[] options = { "yes", "no" };
				int foo = JOptionPane.showOptionDialog(null, message,
						"Warning", JOptionPane.DEFAULT_OPTION,
						JOptionPane.WARNING_MESSAGE, null, options, options[0]);
				return foo == 0;
			}

			// If password is not given before the invocation of
			// Session#connect(),
			// implement also following methods,
			// * UserInfo#getPassword(),
			// * UserInfo#promptPassword(String message) and
			// * UIKeyboardInteractive#promptKeyboardInteractive()
		};
		session.setUserInfo(ui);
	}

	// if have password, set it, otherwise dialog to get it
	private void getUserPass(String remoteUser,String remoteHost,String userPass) throws Exception
	{	if(userPass!=null)
		{	if(userPass.length()>0)
			{	session.setPassword(userPass);
				return;
			}
		}
		
		JPanel passPanel = new JPanel(new GridLayout(2,1));
		JLabel label = new JLabel("Enter password for "+remoteUser+" on "+remoteHost+":");
		JPasswordField passwd = new JPasswordField(35);
		passwd.addAncestorListener(new RequestFocusListener());
		passPanel.add(label);
		passPanel.add(passwd);
		int option = JOptionPane.showConfirmDialog(null, passPanel,
					"Enter Password", JOptionPane.OK_CANCEL_OPTION,
					JOptionPane.PLAIN_MESSAGE);
		if(option != JOptionPane.OK_OPTION)
		{	throw new Exception("user canceled password");
		}
	
		String pass = new String(passwd.getPassword());
		session.setPassword(pass);
	}

	// User Info Class for remote auth
	public static abstract class MyUserInfo implements UserInfo,
			UIKeyboardInteractive {

		public String getPassword() {
			return null;
		}

		public boolean promptYesNo(String str) {
			return false;
		}

		public String getPassphrase() {
			return null;
		}

		public boolean promptPassphrase(String message) {
			return false;
		}

		public boolean promptPassword(String message) {
			return false;
		}

		public void showMessage(String message) {
		}

		public String[] promptKeyboardInteractive(String destination,
				String name, String instruction, String[] prompt, boolean[] echo) {
			return null;
		}
	}
	
	public class RequestFocusListener implements AncestorListener
	{
		public void ancestorAdded(AncestorEvent ae)
		{	ae.getComponent().requestFocusInWindow();
		}
		public void ancestorMoved(AncestorEvent ae)
		{
		}
		public void ancestorRemoved(AncestorEvent ae)
		{
		}
	}
	
	public class SystemOutProgressMonitor implements SftpProgressMonitor
	{
	    public SystemOutProgressMonitor() {;}	    
	    private long fileSize;
	    private long downloaded;

	    public void init(int op, java.lang.String src, java.lang.String dest, long max) 
	    {  	fileSize = max;
	    	downloaded = 0;
	        soutConsole.appendLine("Downloading: "+src+" -> "+dest+" total: "+max);
	        soutConsole.appendLine();			// This line gets replaced below
	    }

	    public boolean count(long bytes)
	    {  	downloaded += bytes;
	    	double pdone = 100.*(double)downloaded/(double)fileSize ;
	    	soutConsole.replaceLastLine("... " + String.valueOf((int)pdone) + "%" );
	        return(true);
	    }

	    public void end()
	    {   soutConsole.replaceLastLine("... 100% done");
	    }
	}
	    
}
