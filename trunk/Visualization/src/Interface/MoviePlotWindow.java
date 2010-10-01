/*******************************************************************
	MoviePlotWindow.java
	NairnFEAMPMViz

	Created by John Nairn on 2/27/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;

import javax.swing.*;

import java.text.*;
import java.awt.event.*;
import java.util.*;
import java.awt.image.*;
import java.io.*;

import javax.imageio.*;
import javax.imageio.stream.ImageOutputStream;
import javax.imageio.event.*;

public class MoviePlotWindow extends NFMVFrame implements Runnable, IIOWriteProgressListener, JNNotificationListener
{
	static final long serialVersionUID=12L;
	
	// constants
	public static final int FRAME_RATE=30;
	
	// variables
	protected ResultsDocument resDoc;
	protected DocViewer docView;
	protected MeshPlotView plotView;
	protected MovieControls movieControls;
	protected boolean runFlag,imageStillWriting;
	protected int movieComponent;
	public int plotType;
	protected String frameFileRoot=null;
	protected MeshPlotScroll plotScroll;
	protected JCheckBoxMenuItem checkedZoomItem;
	
	// command file chooser
	private JFileChooser chooser=new JFileChooser();
	
	// initialize
	public MoviePlotWindow(ResultsDocument gResDoc,DocViewer gDocView)
	{	super(gResDoc.name);
		resDoc=gResDoc;
		docView=gDocView;
		setFramePrefs("Mesh Window Width",800,"Mesh Window Height",600);
		
		// size and location
		Dimension d=setFrameLocAndSize(this);
		
		// add plot view (entire window for now)
		Container content=getContentPane();
		plotView=new MeshPlotView(resDoc);
		//content.add(plotView,BorderLayout.CENTER);
		
		// scroll pane
		plotScroll=new MeshPlotScroll(plotView);
		content.add(plotScroll,BorderLayout.CENTER);
		
		// add controls
		movieControls=new MovieControls(d.width,resDoc,this,gDocView);
		content.add(movieControls,BorderLayout.SOUTH);
		
		NFMVPrefs.setWorkspace(chooser);
		
		// notifications
		JNNotificationCenter.getInstance().addNameAndObjectForTarget("TimeSliderChanged",gDocView,this);
		JNNotificationCenter.getInstance().addNameAndObjectForTarget("PlotQuantityChanged",gDocView,this);
	}
	
	// dispose - remove the notification
	public void dispose()
	{	JNNotificationCenter.getInstance().removeAllForTarget(this);
		super.dispose();
	}

	// make menu bar on launch
	protected void makeMenuBar()
	{	//Create the menu bar
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		
		// File menu
		JMenu menu = new JMenu("File");
		menuBar.add(menu);
		makeMenuItem(menu,"Open...","OpenFile",KeyEvent.VK_O,NairnFEAMPMViz.appCtrl);
		makeMenuItem(menu,"Export Graphics...","ExportFrame",KeyEvent.VK_E,this);
		makeMenuItem(menu,"Export Movie Frames...","ExportMovie",KeyEvent.VK_L,this);
		menu.addSeparator();
		makeMenuItem(menu,"Close","Close",KeyEvent.VK_W,this);
		
		// Window menu
		menu = new JMenu("Window");
		menuBar.add(menu);
		makeMenuItem(menu,"Analysis Results","ShowResults",KeyEvent.VK_D,this);
		
		menu = new JMenu("Zoom");
		menuBar.add(menu);
		checkedZoomItem=makeCheckBoxMenuItem(menu,"100%","ZoomPlot",KeyEvent.VK_1,this);
		makeCheckBoxMenuItem(menu,"150%","ZoomPlot",KeyEvent.VK_2,this);
		makeCheckBoxMenuItem(menu,"200%","ZoomPlot",KeyEvent.VK_3,this);
		makeCheckBoxMenuItem(menu,"300%","ZoomPlot",KeyEvent.VK_4,this);
		makeCheckBoxMenuItem(menu,"400%","ZoomPlot",KeyEvent.VK_5,this);
		makeCheckBoxMenuItem(menu,"500%","ZoomPlot",KeyEvent.VK_6,this);
		makeCheckBoxMenuItem(menu,"700%","ZoomPlot",KeyEvent.VK_7,this);
		makeCheckBoxMenuItem(menu,"1000%","ZoomPlot",KeyEvent.VK_8,this);
		checkedZoomItem.setSelected(true);
}
	
	//----------------------------------------------------------------------------
	// detachable thread and methods for controlling movies
	//----------------------------------------------------------------------------
	
	public void run()
	{
		Date theTime,endTime;
		long pause;
		
		movieControls.setPlaying(true);
		
		try
		{	while(movieControls.incrementArchiveIndex() && runFlag)
			{	theTime=new Date();
			
				try
				{	loadPlotData();
				}
				catch(Exception e)
				{	dispose();
					Toolkit.getDefaultToolkit().beep();
					JOptionPane.showMessageDialog(null,"Error loading plot data: "+e.getMessage());
					break;
				}
				catch(OutOfMemoryError me)
				{	dispose();
					Toolkit.getDefaultToolkit().beep();
					JOptionPane.showMessageDialog(null,"Out of memory error: "+me.getMessage());
					break;
				}
				plotView.repainting=true;
				plotView.repaint();
				
				endTime=new Date();
				pause=FRAME_RATE-endTime.getTime()+theTime.getTime();
				if(pause>0) Thread.sleep(pause);
				while(plotView.repainting)
				{	Thread.sleep(0);
				}
				
				// on movie export, save the frame
				if(frameFileRoot!=null)
				{	File saveFile=new File(frameFileRoot+String.format("%04d.png",movieControls.getArchiveIndex()));
					if(!exportPlotFrame(saveFile)) break;
				}
			}
		}
		catch(InterruptedException e) {}
		
		// reset flags and buttons before done
		movieControls.setPlaying(false);
		runFlag=false;
		frameFileRoot=null;
	}
	
	// change to new time unless in a movie
	public void changeArchiveIndex(int newIndex)
	{	if(runFlag) return;
		try
		{	loadPlotData();
			plotView.repaint();
		}
		catch(Exception e)
		{	dispose();
			Toolkit.getDefaultToolkit().beep();
			JOptionPane.showMessageDialog(null,"Error loading plot data: "+e.getMessage());
		}
		catch(OutOfMemoryError me)
		{	dispose();
			Toolkit.getDefaultToolkit().beep();
			JOptionPane.showMessageDialog(null,"Out of memory error: "+me.getMessage());
		}
	}
	
	// start plot over at new initial index and possibly new component and options
	public void beginNewIndexNewComponent(int newIndex,int newComponent)
	{	movieComponent=newComponent;
		ControlPanel controls=resDoc.docCtrl.controls;
		plotView.setOptions(controls.getOptions(),movieComponent,plotType,controls.getParticleSize());
		
		boolean sameIndex=!movieControls.setArchiveIndex(newIndex);
		if(sameIndex || resDoc.isFEAAnalysis())
		{	// means the index is same OR an FEA plot
			// This will load data and replot (same index reload needed in case other settings changed)
			changeArchiveIndex(newIndex);
		}
		else
		{	// Means MPM plot where index has changed
			// The set archive index above will have already loaded the new plot data
			// but still need to repaint
			plotView.repaint();
		}
	}
	
	// load everything needed to plot or replat data - must override
	public void loadPlotData() throws Exception { }
	
	//----------------------------------------------------------------------------
	// handle movie commands
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{   String theCmd=e.getActionCommand();
	
		if(theCmd.equals("Play"))
		{	if(runFlag)
			{   runFlag=false;
				movieControls.setPlaying(false);
			}
			else
			{	// movie thread
				runFlag=true;
				Thread movieThread=new Thread(this);
				movieThread.start();
			}
		}
		
		else if(theCmd.equals("Rewind"))
			movieControls.setArchiveIndex(0);
			
		else if(theCmd.equals("ShowResults"))
		{	docView.toFront();
		}
		
		else if(theCmd.equals("ExportFrame"))
		{	// get file if need it
			int result = chooser.showSaveDialog(this);
			if(result != JFileChooser.APPROVE_OPTION) return;
			File saveFile=NairnFEAMPMViz.CheckFileStatus(chooser.getSelectedFile(),this,"png");
			if(saveFile==null) return;
			
			// save to png file
			exportPlotFrame(saveFile);
		}
		
		else if(theCmd.equals("ExportMovie"))
		{	if(runFlag || !resDoc.isMPMAnalysis())
			{	// can't do frames if movie is already running
				Toolkit.getDefaultToolkit().beep();
				return;
			}
		
			// get base name without .png extension and not ending in . either
			int result = chooser.showSaveDialog(this);
			if(result != JFileChooser.APPROVE_OPTION) return;
			frameFileRoot=chooser.getSelectedFile().getPath();
			int offset = frameFileRoot.lastIndexOf(".");
			if(offset>0)
			{	String fext=frameFileRoot.substring(offset+1);
				if(fext.equalsIgnoreCase("png") || offset+1==frameFileRoot.length())
					frameFileRoot=frameFileRoot.substring(0,offset);
			}
			
			// export current frame to start
			File saveFile=new File(frameFileRoot+String.format("%04d.png",movieControls.getArchiveIndex()));
			if(!exportPlotFrame(saveFile)) return;
			
			// movie thread
			runFlag=true;
			Thread movieThread=new Thread(this);
			movieThread.start();
		}
		
		else if(theCmd.equals("ZoomPlot"))
		{	checkedZoomItem.setSelected(false);
			checkedZoomItem=(JCheckBoxMenuItem)(e.getSource());
			checkedZoomItem.setSelected(true);
			String theZoomCmd=checkedZoomItem.getText();
			double scale=1.0;
			if(theZoomCmd.equals("150%"))
				scale=1.5;
			else if(theZoomCmd.equals("200%"))
				scale=2.0;	
			else if(theZoomCmd.equals("300%"))
				scale=3.0;	
			else if(theZoomCmd.equals("400%"))
				scale=4.0;	
			else if(theZoomCmd.equals("500%"))
				scale=7.0;	
			else if(theZoomCmd.equals("700%"))
				scale=7.0;	
			else if(theZoomCmd.equals("1000%"))
				scale=10.0;
			plotScroll.setScale(scale);
		}
		
		else
			super.actionPerformed(e);
	}
	
	public void receiveNotification(JNNotificationObject obj)
	{	if(obj.getName().equals("TimeSliderChanged"))
		{	JSlider select=(JSlider)obj.getUserInfo();
			if(select==movieControls.selectTime.select)
			{	// if slider in movie controller changes, load data and synch with doc view
				changeArchiveIndex(select.getValue());
				docView.controls.setArchiveIndex(select.getValue());
			}
			else
			{	// if slider in doc view changed, change movie controller
				// if it causes a changed, this will trigger notification to call again and
				// load the data with above code
				movieControls.setArchiveIndex(select.getValue());
			}
		}
	
		else if(obj.getName().equals("PlotQuantityChanged"))
		{	JComboBox qmenu=(JComboBox)obj.getUserInfo();
			int ctrlPlotType=docView.controls.getPlotType();
			if(qmenu==movieControls.pquant)
			{	// changed in plot window, synch with results window
				if(ctrlPlotType==plotType)
				{	int newIndex=qmenu.getSelectedIndex();
					JComboBox plotQuant=docView.controls.getQuantityMenu();
					if(plotQuant.getSelectedIndex()!=newIndex)
						plotQuant.setSelectedIndex(newIndex);
				}
				if(!movieControls.disableStartPlot)
					docView.startNewPlot(plotType);
			}
			else if(qmenu==movieControls.pcmpnt)
			{	// changed in plot window, synch with results window
				if(ctrlPlotType==plotType)
				{	int newIndex=qmenu.getSelectedIndex();
					JComboBox plotCmpnt=docView.controls.getComponentMenu();
					if(plotCmpnt.getSelectedIndex()!=newIndex)
						plotCmpnt.setSelectedIndex(newIndex);
				}
				if(!movieControls.disableStartPlot)
					docView.startNewPlot(plotType);
			}
			else
			{	// changed in controll panel - synch and replot
				movieControls.syncPlotQuantityMenus();
				docView.startNewPlot(plotType);
			}
		}
	}

	// export current frame to saveFile and return true of false if done OK
	public boolean exportPlotFrame(File saveFile)
	{	try
		{	// create my own ImageWriter
			Iterator<ImageWriter> writers=ImageIO.getImageWritersByFormatName("png");
			ImageWriter writer=writers.next();
			ImageOutputStream ios=ImageIO.createImageOutputStream(saveFile);
			writer.setOutput(ios);
			BufferedImage frameCopy=plotView.frameImage();
			writer.addIIOWriteProgressListener(this);
			imageStillWriting=true;
			writer.write(frameCopy);
			while(imageStillWriting)
			{	Thread.sleep(0);
			}
			frameCopy=null;
		
			// use and arbitrary ImageWriter
			//BufferedImage frameCopy=plotView.frameImage();
			//ImageIO.write(frameCopy,"png",saveFile);
		}
		catch(Exception fe)
		{	Toolkit.getDefaultToolkit().beep();
			JOptionPane.showMessageDialog(this,"Error exporting graphic image: " + fe);
			return false;
		}
		return true;
	}
	
	//----------------------------------------------------------------------------
	// IIOWriteProgressListener interface to verify image is done writing
	//----------------------------------------------------------------------------

	public void imageComplete(ImageWriter source) { imageStillWriting=false; }
	public void imageProgress(ImageWriter source, float percentageDone) {}
	public void imageStarted(ImageWriter source, int imageIndex) {}
	public void thumbnailComplete(ImageWriter source) {}
	public void thumbnailProgress(ImageWriter source, float percentageDone) {}
	public void thumbnailStarted(ImageWriter source, int imageIndex, int thumbnailIndex) {}
	public void writeAborted(ImageWriter source) {}
	
	//----------------------------------------------------------------------------
	// Window events
	//----------------------------------------------------------------------------
	
	public void windowOpened(WindowEvent e)
	{	int newComponent=docView.controls.getPlotComponent();
		beginNewIndexNewComponent(docView.controls.getArchiveIndex(),newComponent);
		plotView.setFirstLoad(true);
	}
	public void	windowClosed(WindowEvent e)
	{	docView.closePlotFrame(this);
	}
	
	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------
	
	public MeshPlotView getPlotView() { return plotView; }
	public boolean isMovieRunning() { return runFlag; }
	public int getPlotType() { return plotType; }
	
	//----------------------------------------------------------------------------
	// static methods
	//----------------------------------------------------------------------------
	
	// nice format for any double
	public static String formatDouble(double d) { return formatDouble(d,false); }
	
	// nice format for any double
	public static String formatDouble(double d,boolean noFraction)
	{	StringBuffer ddec;
		int dlog;
		double da;
		NumberFormat nf=NumberFormat.getInstance();
		
		// display number in good format
		da=d;
		if(d<0) da=-d;
		if(da>1e9)
		{	dlog=(int)(Math.log(da)/Math.log(10.));
			d*=Math.pow(10,-dlog);
			nf.setMaximumFractionDigits(5);
			nf.setMinimumFractionDigits(1);
			ddec=new StringBuffer(nf.format(d));
			ddec.append("e+");
			nf.setMaximumFractionDigits(0);
			ddec.append(nf.format(dlog));
			return ddec.toString();
		}
		else if(da>1.e6)
		{	nf.setMaximumFractionDigits(0);
		}
		else if(da>1.e3)
		{	nf.setMaximumFractionDigits(2);
			if(noFraction) nf.setMaximumFractionDigits(0);
		}
		else if(da>1.)
		{	nf.setMaximumFractionDigits(4);
			if(noFraction) nf.setMaximumFractionDigits(0);
		}
		else if(da>1.e-2)
		{	nf.setMaximumFractionDigits(6);
		}
		else if(da==0.)
		{	nf.setMaximumFractionDigits(0);
		}
		else
		{	if(noFraction) return "0";
			dlog=(int)(Math.log(da)/Math.log(10.));
			dlog--;
			d*=Math.pow(10,-dlog);
			nf.setMaximumFractionDigits(5);
			nf.setMinimumFractionDigits(1);
			ddec=new StringBuffer(nf.format(d));
			ddec.append("e");
			nf.setMaximumFractionDigits(0);
			ddec.append(nf.format(dlog));
			return ddec.toString();
		}
		return nf.format(d);
	}
}
