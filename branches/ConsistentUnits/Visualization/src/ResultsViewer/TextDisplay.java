/*******************************************************************
	TextDisplay.java
	NairnFEAMPMViz

	Created by John Nairn on Sat Mar 06 2004.
	Copyright (c) 2004 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.io.*;
import javax.swing.*;
import javax.swing.event.*;
import java.util.*;
import java.awt.event.*;

import geditcom.JNFramework.*;

public class TextDisplay extends JSplitPane
{
	static final long serialVersionUID=20L;
	
	//-----------------------------------------------------------------
	// Class variables and constants
	//-----------------------------------------------------------------
	
	// members
	public JNTextArea textPane = new JNTextArea( );
	public JTable sectionList;
	private ResultsDocument resDoc;
	private JScrollPane right;
	private boolean newSelection=false;
  
	//-----------------------------------------------------------------
	// Initialize - a split pane with table of sections on the
	// left and scrollable text content on the right
	//-----------------------------------------------------------------
	
	public TextDisplay(ResultsDocument gResDoc)
	{   super(JSplitPane.HORIZONTAL_SPLIT);
		
		// left side is table with sections
		resDoc=gResDoc;
		sectionList=new JTable(resDoc);
		sectionList.setFont(new Font("Monospaced",Font.PLAIN,12));
		sectionList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		ListSelectionModel rowSM=sectionList.getSelectionModel();
		rowSM.addListSelectionListener(new ListSelectionListener()
		{   public void valueChanged(ListSelectionEvent e)
			{   int index;
			    if(e.getValueIsAdjusting()) return;
				ListSelectionModel lsm=(ListSelectionModel)e.getSource();
				index=lsm.getMinSelectionIndex();
				textPane.setText(resDoc.section(index));
				newSelection=true;
				textPane.goTopGetFocus();
			}
		});
		sectionList.setShowGrid(false);
		Component left=new JScrollPane(sectionList);
		setLeftComponent(left);
		
		// right side is text view
		textPane.setEditable(false);
		textPane.setFont(new Font("Monospaced",Font.PLAIN,12));
		right=new JScrollPane(textPane);
		setRightComponent(right);
		
		// listener to scroll new selection to the top
		textPane.addComponentListener(new ComponentListener()
		{	public void componentHidden(ComponentEvent e) {  }
			public void componentMoved(ComponentEvent e) {  }
			public void componentResized(ComponentEvent e)
			{	if(newSelection)
				{	JScrollBar vertBar=right.getVerticalScrollBar();
					vertBar.setValue(vertBar.getMinimum());
					newSelection=false;
				}
			}
			public void componentShown(ComponentEvent e) {  }
		});
}

	//--------------------------------------------------------------------
	// Opening new file starts here. This class decodes generic sections
	// and then calls ResultsDocument to decode specific data
	//--------------------------------------------------------------------
	
	// decode into sections marked by "***** #."
	public void readFile(File file,String fileText) throws Exception
	{
		if(FindGenericSections(fileText)<5)
			throw new Exception("The file does not appear to be output from NairnFEA or NairnMPM.\n   (Note: you should not open MPM binary archive files directly)");
		
		// if nodes and elements, etc. in files, load them now
		loadSubFile(file,"NODAL POINT COORDINATES",
					" Node         x               y               z\n",
					"--------------------------------------------------------\n");
					
		loadSubFile(file,"ELEMENT DEFINITIONS",
					"  No. ID    Nd1   Nd2   Nd3   Nd4   Nd5   Nd6   Nd7   Nd8\n",
					"---------------------------------------------------------------------------------\n");
					
		loadSubFile(file,"NODAL POINTS WITH FIXED DISPLACEMENTS",
					" Node DOF ID  Vel (mm/sec)    Arg (ms/ms^-1)  Angle  Function\n",
					"----------------------------------------------------------------\n");
					
		// decode data needed for plotting
		try
		{	resDoc.DecodeFileSections(file);
		}
		catch (Exception e)
		{	String emsg = e.getMessage();
			if(emsg == null)
				emsg = "Scanner error probably due to corrupted or misformatted data.";
			throw new Exception("Could not find all the required data:\n   " + emsg);
		}
		
		// load first section and exit
		textPane.setText(resDoc.section(0));
	}

	// decode file into sections marked by '*****'
	public int FindGenericSections(String mpm)
	{
		// clear previous one
		resDoc.clear(true);
		
		String section,title;
		int titleEnd,sectionNum=1,crEnd;
		Scanner s=new Scanner(mpm);
		s.useDelimiter("\\*\\*\\*\\*\\*\\s+\\d+\\.\\s+");
		s.useLocale(Locale.US);
		if(s.hasNext()) resDoc.add("TITLE",s.next());
		while(s.hasNext())
		{	section=s.next();
			titleEnd=section.indexOf("\n");
			crEnd=section.indexOf("\r");
			if(titleEnd<0)
				titleEnd=crEnd;
			else if(crEnd>0 && crEnd<titleEnd)
				titleEnd=crEnd;
			if(titleEnd>0)
				title=section.substring(0,titleEnd);
			else
				title="Unknown section";
			resDoc.add(title,sectionNum+". "+section);
			sectionNum++;
		}
		s.close();
		if(sectionNum<=1) return 0;
		
		// redisplay the table
		sectionList.revalidate();
		sectionList.repaint();
		sectionList.setRowSelectionInterval(0,0);
		
		return sectionNum;
	}
	
	// if there is subfile with nodes or elements, load now into the section
	public void loadSubFile(File file,String title,String head1,String head2) throws Exception
	{
		String fdata=resDoc.section(title);
		int index=fdata.indexOf("File:");
		if(index>0)
		{	// get full title
			Scanner s=new Scanner(fdata);
			s.useDelimiter("\\r\\n|\\n|\\r");
			title = s.next();
			// remove the number Assum #. or ##.
			if(title.charAt(1)=='.')
				title = title.substring(3);
			else
				title = title.substring(4);
			try
			{	int endFile=fdata.indexOf("\n",index);
				int crEnd=fdata.indexOf("\r",index);
				if(endFile<0)
					endFile=crEnd;
				else if(crEnd>0 && crEnd<endFile)
					endFile=crEnd;
				if(endFile>0)
				{	File subFile=new File(file.getParentFile(),fdata.substring(index+6,endFile));
					FileReader fr=new FileReader(subFile);
					char [] buffer=new char [(int)subFile.length()];
					fr.read(buffer);
					fr.close();
					String newSection=fdata.substring(0,index).concat(head1).concat(head2);
					resDoc.replace(title,newSection.concat(new String(buffer)));
				}
			}
			catch (Exception e)
			{	throw new Exception("Could not load sub file:\n   " + e.getMessage());
			}
		}
	}
	
}

