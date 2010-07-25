/*******************************************************************
	NairnFEAMPMViz.java
	NairnFEAMPMViz

	Created by John Nairn on Fri Mar 05 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.util.*;
import java.awt.event.*;
import javax.swing.*;

import java.io.*;

public class NairnFEAMPMViz implements ActionListener
{
	static final long serialVersionUID=14L;
	
	//----------------------------------------------------------------------------
	// constants and variables
	//----------------------------------------------------------------------------
	
	//private static final int SPACING=6;
	private static final int MACOS=0;
	private static final int WINOS=1;
	private static final int UNIXOS=2;
	
	private static final int LOAD_EXAMPLE=0;
	private static final int WELCOME_MSG=1;
	private static final int QUIT_MSG=2;
	
	public static final String appNameReadable="NairnFEAMPMViz";
	public static final String versionReadable="Version 3.1.0";
	public static final String copyright="Copyright 2004-2010, John A. Nairn, All Rights Reserved";
	
	public static NairnFEAMPMViz appCtrl;	// global to application controller
	private ArrayList<NFMVViewer> ctrls;		// list of active controllers
	private int docCount=0;
	private AboutWindow aboutWindow=null;
	private HelpWindow helpWindow=null;
	private NFMVPrefs prefWindow=null;
	private int osType;
	public JMenu windowMenu;
	
	public JFileChooser chooser = new JFileChooser( );
  
	//----------------------------------------------------------------------------
	// constructor
	//----------------------------------------------------------------------------

	public NairnFEAMPMViz( )
	{   // global link to this controller
		appCtrl=this;
		
		// hold window menu data
		windowMenu = new JMenu("Window");
		
		// list of controllers
		ctrls=new ArrayList<NFMVViewer>(10);

		// current OS in use
		checkCurrentOS();
		
		// set file filter
		NFMVFilter filter=new NFMVFilter();
		filter.addExtension("mpm");
		filter.addExtension("fea");
		filter.setDescription("MPM or FEA Output Files");
		chooser.setFileFilter(filter);
		filter=new NFMVFilter();
		filter.addExtension("fmcmd");
		filter.addExtension("fcmd");
		filter.addExtension("mcmd");
		filter.addExtension("cmd");
		filter.setDescription("MPM or FEA Command Files");
		chooser.addChoosableFileFilter(filter);
		filter=new NFMVFilter();
		filter.addExtension("mpm");
		filter.addExtension("fea");
		filter.addExtension("fmcmd");
		filter.addExtension("fcmd");
		filter.addExtension("mcmd");
		filter.addExtension("cmd");
		chooser.addChoosableFileFilter(filter);
		
		NFMVPrefs.setWorkspace(chooser);
		
		// preferences
		initializePreferences();
		
		// open new window
		newNFMFile(true,WELCOME_MSG);
	}
	
	private void initializePreferences()
	{
		ColorPicker.setSpectrumType();
		ColorPicker.setNumberOfContours();
	}
	
	//----------------------------------------------------------------------------
	// handle application commands
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{   String theCmd=e.getActionCommand();
	
		if(theCmd.equals("OpenFile"))
		{	openNFMFile();
		}
		
		else if(theCmd.equals("NewFEA"))
		{	newNFMFile(false,LOAD_EXAMPLE);
		}
		
		else if(theCmd.equals("NewMPM"))
		{	newNFMFile(true,LOAD_EXAMPLE);
		}
		
		else if(theCmd.equals("About"))
		{	if(aboutWindow==null)
			{	aboutWindow=new AboutWindow();
				aboutWindow.setLocation(50,screenTopMargin()+10);
			}
			aboutWindow.setVisible(true);
			aboutWindow.toFront();
		}
		
		else if(theCmd.equals("Preferences"))
		{	if(prefWindow==null)
			{	prefWindow=new NFMVPrefs();
				prefWindow.setLocation(40,screenTopMargin()+30);
			}
			prefWindow.setVisible(true);
			prefWindow.toFront();
		}

		else if(theCmd.equals("Help"))
		{	if(helpWindow==null)
			{	helpWindow=new HelpWindow(Main.class.getResource("help.html"));
				helpWindow.setLocation(20,screenTopMargin()+30);
			}
			helpWindow.setVisible(true);
			helpWindow.toFront();
		}
			
		
		else if(theCmd.equals("Quit"))
		{	quit();
		}
	}
	
	// open file with chooser
	public void openNFMFile()
	{
		int result=chooser.showOpenDialog(null);
		if(result==JFileChooser.CANCEL_OPTION) return;
		openNFMFileNow(chooser.getSelectedFile(),null,false);
	}
	
	// open file now
	public void openNFMFileNow(File file,NFMVViewer partner,boolean openMesh)
	{	// exit if not file
		if(file==null) return;
		
		// if already open, show it, or use it
		NFMVViewer docCtrl=isAlreadyOpened(file.getPath());
		if(docCtrl!=null && partner==null) return;
		
		// open the file
		char [] buffer;
		try
		{	FileReader fr=new FileReader(file);
			buffer=new char [(int)file.length()];
			fr.read(buffer);
			String readData=new String(buffer);
			String fname=file.getName();
			if(!fname.endsWith(".mpm") && !fname.endsWith(".fea"))
			{	if(readData.indexOf((char)10)<0)
					readData=readData.replace((char)13,(char)10);
			}
			
			// select appropriate viewer object
			if(docCtrl==null)
			{	if(fname.endsWith(".mpm") || fname.endsWith(".fea"))
				{	// close empty one
					docCtrl=new DocViewer(openMesh);
					ctrls.add(docCtrl);
					NFMVViewer oldCtrl=emptyDocument();
					if(oldCtrl!=null && oldCtrl!=docCtrl) oldCtrl.windowClosing(null);
				}
				else
				{	// is there an unsed one?
					docCtrl=emptyDocument();
					if(docCtrl==null)
					{	docCtrl=new CmdViewer();
						ctrls.add(docCtrl);
					}
				}
			}
			docCtrl.loadTheFile(file,readData);
			docCtrl.setVisible(true);
			docCtrl.setPartner(partner);
			if(partner!=null) partner.setPartner(docCtrl);
			addWindowMenuItem(docCtrl);
		}
		catch (Exception re)
		{	if(file!=null)
			{	Toolkit.getDefaultToolkit().beep();
				JOptionPane.showMessageDialog(null,"Error loading file '"+file.getName()+"':\n   "+re.getMessage());
			}
			else
			{	Toolkit.getDefaultToolkit().beep();
				JOptionPane.showMessageDialog(null,"Error loading the file:\n   "+re.getMessage());
			}
			if(docCtrl!=null) docCtrl.windowClosing(null);
		}
		catch(OutOfMemoryError me)
		{	Toolkit.getDefaultToolkit().beep();
			JOptionPane.showMessageDialog(null,"Out of memory error: "+me.getMessage());
			if(docCtrl!=null) docCtrl.windowClosing(null);
		}
	}
	
	// new commands file
	public void newNFMFile(boolean isMPM,int loadType)
	{	// start with empty string
		String readData=new String();
		
		if(loadType==LOAD_EXAMPLE)
		{	// read commands
			InputStream ins=null;
			
			if(isMPM)
				ins=Main.class.getResourceAsStream("MPMCommands.fmcmd");
			else
				ins=Main.class.getResourceAsStream("FEACommands.fmcmd");
		
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
		else if(loadType==WELCOME_MSG)
			readData="Welcome to "+appNameReadable+", "+versionReadable;
		else
			readData="Quit to exit or open/create new set of commands";

		// create viewer
		try
		{	CmdViewer docCtrl=(CmdViewer)emptyDocument();
			if(docCtrl==null)
			{	docCtrl=new CmdViewer();
				ctrls.add(docCtrl);
			}
			docCtrl.loadTheFile(null,readData);
			docCtrl.setVisible(true);
			addWindowMenuItem(docCtrl);
		}
		catch(Exception e) {}
	}
	
    public void quit()
	{	if(quitConfirmed())
            System.exit(0);
    }
	
	public boolean quitConfirmed()
	{	// check all current controllers
		for(int i=0;i<ctrls.size();i++)
		{	if(!ctrls.get(i).confirmClosing())
				return false;
		}
		return true;
	}
	
	//----------------------------------------------------------------------------
	// methods
	//----------------------------------------------------------------------------
	
	// when window closes, remove from controller list
	public void closeViewer(NFMVViewer docCtrl)
	{	// remove partners
		int i;
		for(i=0;i<ctrls.size();i++)
			ctrls.get(i).removePartner(docCtrl);
		
		// remove from Window menu
		for(i=0;i<windowMenu.getMenuComponentCount();i++)
		{	NFMVWindowMenuItem mItem=(NFMVWindowMenuItem)(windowMenu.getMenuComponent(i));
			if(mItem.getController()==docCtrl)
			{	windowMenu.remove(mItem);
				break;
			}
		}
		
		// other controllers
		for(i=0;i<ctrls.size();i++) ctrls.get(i).removeWindowMenuItem(docCtrl);

		// remove the viewer object
		if(!ctrls.remove(docCtrl))
			System.out.println("Could not find window's viewer in the list of controllers");
		
		// was it the last one?
		if(ctrls.size()==0)
		{	Toolkit.getDefaultToolkit().beep();
			String message="Click 'OK' to quit/exit NairnFEAMPMViz or use 'Cancel' to continue with a new window.";
			int result=JOptionPane.showConfirmDialog(null, message, quitText()+"?", JOptionPane.OK_CANCEL_OPTION);
			if(result==JOptionPane.OK_OPTION) System.exit(0);
			
			// open new window
			newNFMFile(true,QUIT_MSG);
		}
	}
	
	// add window menu item
	public void addWindowMenuItem(NFMVViewer docCtrl)
	{	
		NFMVWindowMenuItem mItem=null;
		int i;
		for(i=0;i<windowMenu.getMenuComponentCount();i++)
		{	mItem=(NFMVWindowMenuItem)(windowMenu.getMenuComponent(i));
			if(mItem.getController()==docCtrl) break;
			mItem=null;
		}
		
		// create or change the menu item
		if(mItem==null)
		{	mItem=new NFMVWindowMenuItem(docCtrl,docCtrl.getTitle());
			windowMenu.add(mItem);
		}
		else
			mItem.setText(docCtrl.getTitle());
		
		// other controllers
		for(i=0;i<ctrls.size();i++) ctrls.get(i).addWindowMenuItem(docCtrl);
	}
	
	// if already opened, bring to front and return true
	private NFMVViewer isAlreadyOpened(String path)
	{	for(int i=0;i<ctrls.size();i++)
		{	if(path.equals(ctrls.get(i).getFullPath()))
			{	ctrls.get(i).toFront();
				return ctrls.get(i);
			}
		}
		return null;
	}
	
	// look for unused new document
	private NFMVViewer emptyDocument()
	{	for(int i=0;i<ctrls.size();i++)
		{	NFMVViewer docCtrl=ctrls.get(i);
			if(docCtrl.isReuseable())
				return docCtrl;
		}
		return null;
	}

	
	// Shift window rectangle to new location and make sure fits in the window
	public Rectangle getViewerBounds(Dimension d,Rectangle gcBounds)
	{
		// height
		int shift=screenTopMargin()+4+docCount*15;
		gcBounds.y+=shift;
		gcBounds.height-=(shift+screenBottomMargin()+10);		// for windows start bar
		if(gcBounds.height>d.height) gcBounds.height=d.height;
		
		// width
		shift=10+30*docCount;
		gcBounds.x+=shift;
		gcBounds.width-=(shift+10);
		if(gcBounds.width>d.width) gcBounds.width=d.width;
		
		// doc count
		docCount++;
		if(docCount>5) docCount=0;
		
		return gcBounds;
	}
	
	//----------------------------------------------------------------------------
	// OS specific stuff
	//----------------------------------------------------------------------------
	
	// save look and feel for some specific settings
	public void checkCurrentOS()
	{	String os=System.getProperty("os.name").toLowerCase();
		if(os.indexOf("mac")>=0)
			osType=MACOS;
		else if(os.indexOf("win")>=0)
			osType=WINOS;
		else
			osType=UNIXOS;
	}
	
	public boolean isWindows() { return osType==WINOS; }
	
	// bottom margin - look out for windows start bar
	public int screenBottomMargin() { return osType==WINOS ? 24 : 0 ; }
	public int menuKeyMask() { return osType==MACOS ? ActionEvent.META_MASK : ActionEvent.CTRL_MASK ; }
	public int screenTopMargin() { return osType==WINOS ? 0 : 24 ; }

	// try to open browser
	public boolean showInBrowser(String url)
	{
		Runtime rt = Runtime.getRuntime();
		try
		{	switch(osType)
			{	case MACOS:
					rt.exec( "open " + url);
					break;
				case WINOS:
					rt.exec( "rundll32 url.dll,FileProtocolHandler " + url);
					break;
				default:
					return false;
			}
		}
		catch (IOException e)
		{	return false;
        }
		return true;
	}
	
	// return quit string (which is Exit in Windows)
	public String quitText()
	{	if (osType == WINOS) return "Exit";
		return "Quit";
	}
	
	//----------------------------------------------------------------------------
	// class methods
	//----------------------------------------------------------------------------
	
	public static File CheckFileStatus(File newFile,JFrame fileWindow,String ext)
	{	// force extension
		if(ext!=null)
		{	String fileName=newFile.getPath();
			int offset = fileName.lastIndexOf(".");
			if(offset>0)
			{	String fext=fileName.substring(offset+1);
				if(!fext.equalsIgnoreCase(ext))
					newFile=new File(fileName.substring(0,offset+1)+ext);
			}
			else
				newFile=new File(fileName+"."+ext);
		}
		
		// OK if does not exist
		if(!newFile.exists()) return newFile;
		
		// confirm replacement
		Toolkit.getDefaultToolkit().beep();
		int result=JOptionPane.showConfirmDialog(fileWindow,"The selected file ("+newFile.getName()+") already exists. Do you want to replace it?",
                               "Save File",JOptionPane.YES_NO_OPTION);
		if(result==JOptionPane.YES_OPTION)
			return newFile;
		else
			return null;
	}

}

