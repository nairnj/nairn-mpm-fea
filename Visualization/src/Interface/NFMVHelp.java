/*
 * NFMVHelp.java
 * NairnFEAMPMViz Application
 * 
 * Created by John Nairn on 12/18/12
 * Copyright (c) 2008-2012 RSAC Software. All rights reserved
 */

import java.awt.event.*;
import java.net.URL;

import javax.swing.*;

import geditcom.JNFramework.*;

public class NFMVHelp extends JNHelpWindow implements ActionListener
{	
	private static final long serialVersionUID = 1L;
	
	// constructor
	public NFMVHelp(URL helpURL,String name)
	{	super(helpURL,name);
	
		// tool bar
		Class<?> baseClass=JNApplication.main.getClass();
		ImageIcon showTopic = new ImageIcon(baseClass.getResource("geditcom/JNFramework/Resources/app-help.png"));
		addToolBarIcon(showTopic,"apphelp","Show main application help information.",this);
		showTopic=new ImageIcon(baseClass.getResource("Resources/commands-editor.png"));
		addToolBarIcon(showTopic,"genhelp","Show help on scripting commands.",this);
		showTopic=new ImageIcon(baseClass.getResource("Resources/materialtypes.png"));
		addToolBarIcon(showTopic,"mathelp","Show help on material types.",this);
		showTopic=new ImageIcon(baseClass.getResource("Resources/scripticon.png"));
		addToolBarIcon(showTopic,"scrhelp","Show help on script control of calculations.",this);
		showTopic=new ImageIcon(baseClass.getResource("Resources/osupdoc-wiki.png"));
		addToolBarIcon(showTopic,"osupdocs","Open OSUPDocs wiki site in your browser.",this);
	}
	
	// pick new help topic
	public void actionPerformed(ActionEvent e)
	{   String theCmd=e.getActionCommand();

		URL helpURL = null;
		if(theCmd.equals("apphelp"))
		{	helpURL = JNApplication.main.getClass().getResource("Resources/help.html");
		}
		
		else if(theCmd.equals("genhelp"))
		{	helpURL = JNApplication.main.getClass().getResource("Resources/commands.html");
		}
		
		else if(theCmd.equals("mathelp"))
		{	helpURL = JNApplication.main.getClass().getResource("Resources/materials.html");
		}
		
		else if(theCmd.equals("scrhelp"))
		{	helpURL = JNApplication.main.getClass().getResource("Resources/scripting.html");
		}
		
		else if(theCmd.equals("osupdocs"))
		{	JNApplication.showInBrowser("http://osupdocs.forestry.oregonstate.edu/index.php/Main_Page");
			return;
		}
		
		if(helpURL==null)
			JNApplication.appBeep();
		else
			openURL(helpURL);
	}
	
	// find current help window if available or open new one if not
	public static NFMVHelp customHelpWindow(JNApplication appCtrl)
	{	NFMVHelp currentHelpWindow = (NFMVHelp)getTheHelpWindow();
		if(currentHelpWindow==null)
		{	URL helpURL = getHelpResource()!= null ? appCtrl.getClass().getResource(getHelpResource()) : null ;
			currentHelpWindow = new NFMVHelp(helpURL,"NairnFEAMPMViz Help");
		}
		return currentHelpWindow;
	}
}
