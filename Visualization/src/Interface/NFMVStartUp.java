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
	{	super(JNApplication.appNameReadable);
	
		// set window size
		setFramePrefs("NFMVStart Width",375,"NFMVStart Height",330+JNApplication.windowMenuBarHeight());
		
		// menus
		makeMenuBar(false);
	    NairnFEAMPMViz.addExamplesMenu(getJMenuBar(),"File");
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
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		buttons.setLayout(gridbag);

		int top=0;
		int bottom = 3;
		int left = 2;
		int right = 2;
		c.insets=new Insets(top,left,bottom,right);			// tlbr
		c.gridx=0;
		c.gridwidth = 1;
		c.weightx = 1.0;		
		c.weighty = 1.0;
		c.fill = GridBagConstraints.HORIZONTAL;
		
		JButton openDoc = new JButton("Open existing document");
		openDoc.setActionCommand("openDocument");
		openDoc.addActionListener(JNApplication.main);
		openDoc.setToolTipText("Click to open a previous commands or results document");
		openDoc.setFocusPainted(false);
		
		JButton newMPM = new JButton("New MPM commands document");
		newMPM.setActionCommand("newDocumentMPMCmd");
		newMPM.addActionListener(JNApplication.main);
		newMPM.setToolTipText("Click to create a new document for MPM calculations");
		newMPM.setFocusPainted(false);
		
		JButton newFEA = new JButton("New FEA commands document");
		newFEA.setActionCommand("newDocumentFEACmd");
		newFEA.addActionListener(JNApplication.main);
		newFEA.setToolTipText("Click to create a new document for FEA calculations");
		newFEA.setFocusPainted(false);

		JButton quitApp = new JButton("Quit NairnFEAMPMViz");
		quitApp.setActionCommand("Quit");
		quitApp.addActionListener(JNApplication.main);
		quitApp.setToolTipText("Click to exit the NairnFEAMPMViz application");
		quitApp.setFocusPainted(false);
		
		JLabel seeAlso = new JLabel("(See menu commands for more options)",JLabel.CENTER);
		
		// add buttons and add to window
		gridbag.setConstraints(openDoc, c);
		buttons.add(openDoc);
		gridbag.setConstraints(newMPM, c);
		buttons.add(newMPM);
		gridbag.setConstraints(newFEA, c);
		buttons.add(newFEA);
		gridbag.setConstraints(quitApp, c);
		buttons.add(quitApp);
		gridbag.setConstraints(seeAlso, c);
		buttons.add(seeAlso);
		add(buttons,BorderLayout.CENTER);

		// finish up
		finishFrameworkWindow(true,new Point(40,10+JNApplication.menuBarHeight()));
		
		// fixed size
		setSize(new Dimension(435,325+JNApplication.windowMenuBarHeight()));
	}
	
}
