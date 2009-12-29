/*******************************************************************
	HelpWindow.java
	NairnFEAMPMViz

	Created by John Nairn on 8/20/07.
	Copyright (c) 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import javax.swing.*;
import javax.swing.event.*;
import java.net.*;

public class HelpWindow extends JFrame
{
	static final long serialVersionUID=6L;
	
	protected JEditorPane editorPane;
	
	HelpWindow(URL helpURL)
	{	super("NairnFEAMPMViz Help");
	
		Container content=getContentPane();
		content.setLayout(new BorderLayout());
		
		editorPane=new JEditorPane();
		editorPane.setEditable(false);
		JScrollPane scroller=new JScrollPane(editorPane);
		content.add(scroller,BorderLayout.CENTER);
		
		openURL(helpURL);
		
		editorPane.addHyperlinkListener(new LinkActivator());
		
		setSize(600,600);
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
	}
	
	// load URL into the window
	protected void openURL(URL helpURL)
	{
		try
		{	editorPane.setPage(helpURL);
		}
		catch(Exception e)
		{	System.out.println("Could not open "+helpURL+" : "+e.getMessage());
		}
	}
	
	class LinkActivator implements HyperlinkListener
	{
		public void hyperlinkUpdate(HyperlinkEvent he)
		{	HyperlinkEvent.EventType type=he.getEventType();
			if(type==HyperlinkEvent.EventType.ACTIVATED)
				openURL(he.getURL());
		}
	}
}
