/*
	NairnFEAMPMViz.java
	NairnFEAMPMViz Application	

	Created by John Nairn on 3/7/08.
	Copyright 2008 RSAC Software. All rights reserved.
*/

import java.io.File;
import java.io.InputStream;
import java.net.URL;

import javax.swing.JFileChooser;
import geditcom.JNFramework.*;

public class NairnFEAMPMViz extends JNApplication
{
	public static String jarFolder = null;
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public NairnFEAMPMViz()
	{	super("NairnFEAMPMViz","Version 7.3",
				"Java application for running and visualizing NairnMPM and NairnFEA calculations.");
		NFMVPrefs.setWorkspace(chooser);
		
		// path to folder containing this jar file ending in slash
		URL jarURL = getClass().getProtectionDomain().getCodeSource().getLocation();
		File jarFile = new java.io.File(jarURL.getPath()).getAbsoluteFile();
		jarFolder = jarFile.getPath();
		int slash = jarFolder.lastIndexOf('/');
		if(slash<0) slash = jarFolder.lastIndexOf('\\');
		if(slash>=0) jarFolder = jarFolder.substring(0, slash+1);

		if(!finishLaunching(false))
			new NFMVStartUp();
	}
		
	//----------------------------------------------------------------------------
	// Override document methods
	//----------------------------------------------------------------------------
	
	// create new document
	public void newDocument(String docType)
	{	// initial text in new document
		String readData=new String();
	
		if(docType.equals("FEACmd") || docType.equals("MPMCmd"))
		{	// read default commands
			InputStream ins=null;
			if(docType.equals("MPMCmd"))
				ins=NairnFEAMPMViz.class.getResourceAsStream("Resources/MPMCommands.fmcmd");
			else
				ins=NairnFEAMPMViz.class.getResourceAsStream("Resources/FEACommands.fmcmd");
		
			// read commands
			try
			{	if(ins==null)
					throw new Exception("resource not found");
				int remaining;
				while(true)
				{	remaining=ins.available();
					if(remaining==0) break;
					byte [] buffer=new byte[remaining];
					ins.read(buffer,0,remaining);
					readData=readData+(new String(buffer));
				}
			}
			catch (Exception e)
			{	System.out.println("Error reading commands resource: " + e);
				readData="";
			}
		}
		else if(docType.equals("SCRPCmd"))
			readData="!********** Control Script **********\nScript\n";
		else
			readData="Welcome to "+appNameReadable+", "+versionReadable;

		// create viewer (or use empty one in the front) and add to application
		docType="MPMCmd";
		CmdViewer docCtrl=null;
		JNDocument frontDoc=(JNDocument)frontDocument();
		if(frontDoc==null)
			docCtrl=new CmdViewer(docType);
		else if(frontDoc.getDocType().equals(docType) && frontDoc.isEmptyDocument())
			docCtrl=(CmdViewer)frontDoc;
		else
			docCtrl=new CmdViewer(docType);
		docCtrl.loadTextFromFile(readData);
		openUntitledDocument(docCtrl);
	}
	
	public void openDocumentType(String docType,String docText,File docFile)
	{	// open that type
		JNDocument front=JNApplication.main.frontDocument();
		super.openDocumentType(docType,docText,docFile);
		if(front!=JNApplication.main.frontDocument() && front!=null)
		{	if(front.isEmptyDocument() && !front.getChanged())
				front.closeDocument(false,null);
		}
	}

	// called when open document
	public JNDocument createDocumentObject(String docType)
	{	JNDocument front=JNApplication.main.frontDocument();
		if(front!=null)
		{	if(docType.equals(front.getDocType()) && front.isEmptyDocument())
				return front;
		}
		
		if(docType.equals("FEACmd") || docType.equals("MPMCmd"))
			return new CmdViewer(docType);
		else if(docType.equals("DocViewer"))
			return new DocViewer(false);
		else
		{	// must be unknown, type to open the text
			return new CmdViewer("MPMCmd");
		}
	}
	
	public boolean tryToContinueApplication()
	{	// display small window with limited commands
		new NFMVStartUp();
		return true;
	}
	
	//----------------------------------------------------------------------------
	// Override preferences and help methods
	//----------------------------------------------------------------------------
	
	// called first time preferences window needs to open
	public JNPreferences openPreferencesWindow()
	{	return new NFMVPrefs();
	}

	public void openHelp()
	{	NFMVHelp helpWindow=NFMVHelp.customHelpWindow(this);
		helpWindow.setVisible(true);
		helpWindow.toFront();
	}

	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	public static JFileChooser GetAppChooser() { return JNApplication.chooser; }
	
	//----------------------------------------------------------------------------
	// main entry as static methods
	//----------------------------------------------------------------------------
	
	public static void main(String[] args)
	{	javax.swing.SwingUtilities.invokeLater(new Runnable()
		{	public void run()
            {	createAndShowGUI();
            }
        });
    }
	
	// Create the GUI and show it.  For thread safety,
	// this method should be invoked from the
	// event-dispatching thread.
    private static void createAndShowGUI()
	{	// initial preferences first
		NFMVPrefs.initializePrefs();
				
		// select application options
		// first is to use native LNF and second is menu bar for Mac
		setLookAndFeel(true,true);
		JNHelpWindow.setHelpResources("Resources/help.html",JNApplication.appNameReadable+" Help");
		
		// optional info strings
		JNApplication.iconResource="Resources/AboutIcon.png";
		JNApplication.copyright="Copyright 2004-2019, John A. Nairn, All Rights Reserved";
		JNApplication.author="Written and documented by John A. Nairn";
		JNApplication.webSite="http://www.cof.orst.edu/cof/wse/faculty/Nairn/";
		
		// set document types
		String[] exts1={"fmcmd","cmd"};
		JNDocument.setDocumentType("MPMCmd","MPM Commands Document",exts1);
		String[] exts2={"fmcmd","cmd"};
		JNDocument.setDocumentType("FEACmd","FEA Commands Document",exts2);
		String[] exts3={"fmcmd","cmd"};
		JNDocument.setDocumentType("SCRPCmd","Control Script Document",exts3);
		String[] exts4={"mpm","fea"};
		JNDocument.setReadOnlyDocumentType("DocViewer","MPM or FEA Results Document",exts4);
		
		// set save workspace
		NFMVPrefs.setWorkspace(JNDocument.getChooser("MPMCmd"));
		
		// create custom instance of the application class
		new NairnFEAMPMViz();
    }
    
}
