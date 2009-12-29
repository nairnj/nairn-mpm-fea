/*******************************************************************
	AboutWindow.java
	NairnFEAMPMViz

	Created by John Nairn on 8/17/07.
	Copyright (c) 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import javax.swing.*;
import java.net.*;
import java.awt.event.*;

public class AboutWindow extends JFrame implements ActionListener
{
	static final long serialVersionUID=1L;
	
	private static final int HMARGIN=12;
	private static final int VMARGIN=12;
	private static final int LINE_SPACING=6;
	
	public AboutWindow()
	{	super("About NairnFEAMPMViz");
	
		// icon and name in large font
		URL iconImage=Main.class.getResource("AboutIcon.png");
		JLabel iconLabel=new JLabel("  "+NairnFEAMPMViz.appNameReadable,new ImageIcon(iconImage),JLabel.CENTER);
		iconLabel.setFont(new Font("Sans Serif",Font.BOLD,24));
		iconLabel.getPreferredSize();
		Dimension size=iconLabel.getPreferredSize();
		
		// a horiztonal box
		Container box=Box.createHorizontalBox();
		Dimension hrigid=new Dimension(HMARGIN,size.height+2*VMARGIN);
		box.add(Box.createRigidArea(hrigid));
		box.add(iconLabel);
		box.add(Box.createRigidArea(hrigid));
		add(box,BorderLayout.CENTER);
		
		// text and links in the bottom
		Font small=new Font("Sans Serif",Font.PLAIN,12);
		JLabel type=new JLabel("   Java Application for NairnFEA and NairnMPM Visualization");
		type.setFont(small);
		type.setAlignmentX(Component.LEFT_ALIGNMENT);
		
		JLabel version=new JLabel("   "+NairnFEAMPMViz.versionReadable);
		version.setFont(small);
		
		JLabel developer=new JLabel("   Written and documented by John A. Nairn");
		developer.setFont(small);
		
		JButton homePage=addURLButton("http://woodscience.oregonstate.edu/faculty/nairn",small);
		
		JButton webSite=addURLButton("http://oregonstate.edu/~nairnj/NairnFEAMPM",small);
		
		JLabel copy=new JLabel("          "+NairnFEAMPMViz.copyright);
		copy.setFont(new Font("Sans Serif",Font.PLAIN,10));
		
		// vertical box
		Container vbox=Box.createVerticalBox();
		vbox.add(Box.createVerticalStrut(LINE_SPACING));
		vbox.add(type);
		vbox.add(Box.createVerticalStrut(LINE_SPACING));
		vbox.add(version);
		vbox.add(Box.createVerticalStrut(LINE_SPACING));
		vbox.add(developer);
		vbox.add(Box.createVerticalStrut(LINE_SPACING));
		vbox.add(homePage);
		vbox.add(Box.createVerticalStrut(LINE_SPACING));
		vbox.add(webSite);
		vbox.add(Box.createVerticalStrut(LINE_SPACING));
		vbox.add(copy);
		vbox.add(Box.createVerticalStrut(VMARGIN));
		add(vbox,BorderLayout.SOUTH);
		
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		pack();
		setResizable(false);
	}
	
	// create button to link to a URL
	private JButton addURLButton(String linkText,Font linkFont)
	{
		JButton urlLink=new JButton("        "+linkText);
		urlLink.setFont(linkFont);
		urlLink.setForeground(Color.blue);
		urlLink.setContentAreaFilled(false);
		urlLink.setFocusPainted(false);
		urlLink.setBorderPainted(false);
		urlLink.setActionCommand(linkText);
		urlLink.addActionListener(this);
		return urlLink;
	}

	//----------------------------------------------------------------------------
	// try to open browser
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{   String theCmd=e.getActionCommand();
		NairnFEAMPMViz.appCtrl.showInBrowser(theCmd);
	}
}
