/*
	NairnFEAMPMViz.java
	NairnFEAMPMViz Application	

	Created by John Nairn on 3/7/08.
	Copyright 2008 RSAC Software. All rights reserved.
*/

import java.awt.event.ActionEvent;
import java.io.File;
import java.io.InputStream;
import java.net.URL;

import javax.swing.JFileChooser;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

import geditcom.JNFramework.*;

public class NairnFEAMPMViz extends JNApplication
{
	public static String jarFolder = null;
	//----------------------------------------------------------------------------
	// Initialize
	//----------------------------------------------------------------------------
	
	public NairnFEAMPMViz()
	{	super("NairnFEAMPMViz","Version 7.4",
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
	
	// open example files
	public void actionPerformed(ActionEvent e)
	{	String theCmd=e.getActionCommand();
	
		if(theCmd.startsWith("$SAMPLE$"))
		{	String exFile = theCmd.substring(8);
			newDocument("MPMCmd"+exFile);
		}
		else
			super.actionPerformed(e);;
	}
		

	// create new document
	public void newDocument(String docType)
	{	// initial text in new document
		String readData=new String();
	
		if(docType.startsWith("FEACmd") || docType.startsWith("MPMCmd"))
		{	// read default commands or an example file
			docType = docType.length()>6 ? docType.substring(6) : docType ;
			
			InputStream ins=null;
			if(docType.equals("MPMCmd"))
				ins=NairnFEAMPMViz.class.getResourceAsStream("Resources/MPMCommands.fmcmd");
			else if(docType.equals("FEACmd"))
				ins=NairnFEAMPMViz.class.getResourceAsStream("Resources/FEACommands.fmcmd");
			else
			{	// an example files
				String resName = "Resources/"+docType.replace(' ','_')+".fmcmd";
				ins=NairnFEAMPMViz.class.getResourceAsStream(resName);
			}
		
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
    
    // add menu of examples
    public static void addExamplesMenu(JMenuBar menuBar,String menuName)
    {
		// get menu
		JMenu fileMenu = null;
		int count = menuBar.getMenuCount();
		for(int i=0;i<count;i++)
		{	JMenu testMenu = menuBar.getMenu(i);
			if(testMenu.getText().equals(menuName))
			{	fileMenu = testMenu;
				break;
			}
		}
		if(fileMenu==null)
		{	System.out.println(menuName+" menu not found");
			return;
		}
		
		// add examples menu
		JMenu examplesMenu = new JMenu("Examples");
		fileMenu.add(examplesMenu);

		// read list of examples
		String readData=new String();
		InputStream ins=NairnFEAMPMViz.class.getResourceAsStream("Resources/examples.txt");
		
		// erro if not found = insert into menu
		if(ins==null)
		{	JMenuItem menuItem = new JMenuItem("(example files list not found)");
			examplesMenu.add(menuItem);
			return;
		}
		
		try
		{	int remaining;
			while(true)
			{	remaining=ins.available();
				if(remaining==0) break;
				byte [] buffer=new byte[remaining];
				ins.read(buffer,0,remaining);
				readData=readData+(new String(buffer));
			}
			String[] examples = readData.split("\\r\\n|\\n|\\r");
			
			if(examples.length==0)
			{	JMenuItem menuItem = new JMenuItem("(example files list is empty)");
				examplesMenu.add(menuItem);
			}
			else
			{	// add examples to submenu
				for(int i=0;i<examples.length;i++)
				{	JMenuItem menuItem = new JMenuItem(examples[i]);
					menuItem.setActionCommand("$SAMPLE$"+examples[i]);
					menuItem.addActionListener(JNApplication.main);
					examplesMenu.add(menuItem);
				}
			}
		}
		catch (Exception e)
		{	System.out.println("Error reading list of examples: " + e);
			readData="";
		}
    }
}
