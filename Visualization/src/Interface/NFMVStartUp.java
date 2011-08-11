/*
 * NFMVStartUp.java
 * NairnFEAMPMViz Application
 * 
 * Created by John Nairn on 8/10/11
 * Copyright (c) 2008-2011 RSAC Software. All rights reserved
 */

import java.awt.*;
import java.net.URL;

import javax.swing.*;

import geditcom.JNFramework.*;

public class NFMVStartUp extends JNStartUpBase
{
	private static final long serialVersionUID = 1L;
	
	// Create the window
	public NFMVStartUp()
	{	super("NairnFEAMPMViz");
	
		// set window size
		setFramePrefs("NFMVStart Width",375,"NFMVStart Height",330+JNApplication.windowMenuBarHeight());
		
		// menus
		makeMenuBar(false);
		
		setLayout(new BorderLayout());
		
		// icon (optional) and name in large font
		JLabel iconLabel;
		String iconResource=JNApplication.iconResource;
		if(iconResource!=null)
		{	try
			{	Class<?> baseClass=JNApplication.main.getClass();
				URL iconImage=baseClass.getResource(iconResource);
				iconLabel=new JLabel("  "+JNApplication.appNameReadable,new ImageIcon(iconImage),JLabel.CENTER);
			}
			catch (Exception e)
			{	System.out.println("Icon resource '"+iconResource+"' failed to load");
				iconLabel=new JLabel("  "+JNApplication.appNameReadable,JLabel.CENTER);
			}
		}
		else
			iconLabel=new JLabel("  "+JNApplication.appNameReadable,JLabel.CENTER);
		iconLabel.setFont(new Font("Sans Serif",Font.BOLD,24));
		Dimension size=iconLabel.getPreferredSize();
		
		// a horizontal box
		Container box=Box.createHorizontalBox();
		Dimension hrigid=new Dimension(12,size.height+24);
		box.add(Box.createRigidArea(hrigid));
		box.add(iconLabel);
		box.add(Box.createRigidArea(hrigid));
		add(box,BorderLayout.NORTH);
		
		// Panel of buttons
		JPanel buttons = new JPanel();
		buttons.setLayout(new GridLayout(5,1));
		
		JButton openDoc = new JButton("Open existing document");
		openDoc.setActionCommand("openDocument");
		openDoc.addActionListener(JNApplication.main);
		openDoc.setToolTipText("Click to open a previous commands or results document");
		
		JButton newMPM = new JButton("New MPM commands document");
		newMPM.setActionCommand("newDocumentMPMCmd");
		newMPM.addActionListener(JNApplication.main);
		newMPM.setToolTipText("Click to create a new document for MPM calculations");
		
		JButton newFEA = new JButton("New FEA commands document");
		newFEA.setActionCommand("newDocumentFEACmd");
		newFEA.addActionListener(JNApplication.main);
		newFEA.setToolTipText("Click to create a new document for FEA calculations");

		JButton quitApp = new JButton("Quit NairnFEAMPMViz");
		quitApp.setActionCommand("Quit");
		quitApp.addActionListener(JNApplication.main);
		quitApp.setToolTipText("Click to exit the NairnFEAMPMViz application");
		
		JLabel seeAlso = new JLabel("(See menu commands for more options)",JLabel.CENTER);
		
		// add buttons and add to window
		buttons.add(openDoc);
		buttons.add(newMPM);
		buttons.add(newFEA);
		buttons.add(quitApp);
		buttons.add(seeAlso);
		add(buttons,BorderLayout.CENTER);

		// finish up
		finishFrameworkWindow(true,new Point(40,10+JNApplication.menuBarHeight()));
		
		// fixed size
		setSize(new Dimension(385,325+JNApplication.windowMenuBarHeight()));
	}
	
}
