/*
 * MoviePlotWindow.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 2/27/07.
 * Copyright 2007 RSAC Software. All rights reserved.
 */

import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.*;
import javax.imageio.*;
import javax.imageio.event.IIOWriteProgressListener;
import javax.imageio.stream.ImageOutputStream;
import javax.swing.*;

import geditcom.JNFramework.*;

public class MoviePlotWindow extends JPanel implements  Runnable,
				IIOWriteProgressListener, JNNotificationListener, ActionListener
{
	private static final long serialVersionUID = 1L;

	// constants
	public static final int FRAME_RATE=30;

	// variables
	protected ResultsDocument resDoc;
	protected DocViewer docView;
	protected MeshPlotView plotView;
	protected MovieControls movieControls;
	protected int movieComponent;
	protected boolean runFlag,imageStillWriting;
	protected MeshPlotScroll plotScroll;
	protected String frameFileRoot=null;
	double currentScale = 1.0;

	// command file chooser
	private JFileChooser chooser=new JFileChooser();
	
	private static Cursor zoomIn = null;
	private static Cursor zoomOut = null;
	private static Cursor zoomNo = null;

	// initialize
	public MoviePlotWindow(ResultsDocument gResDoc,DocViewer gDocView)
	{	super(new BorderLayout());
		resDoc=gResDoc;
		docView=gDocView;
		
		// add plot view (will be to center)
		plotView=new MeshPlotView(resDoc);
		
		// listen to clicks in the plot view
		plotView.addMouseListener(new MouseAdapter()
		{	public void mouseClicked(MouseEvent e)
		    {	if(e.isMetaDown())
		    	{	int newItem = 0;
		    		if(e.isAltDown())
		    		{	// reduce 1, 1.5, 2,3,4,5,7,10
		    			if(currentScale<1.25)
		    			{	JNApplication.appBeep();
		    				return;
		    			}
		    			else if(currentScale<2.5)
		    				currentScale -= 0.5;
		    			else if(currentScale<6.)
		    				currentScale -= 1.0;
		    			else if(currentScale<8.0)
		    				currentScale = 5.;
		    			else
		    				currentScale = 7.;
		    		}
		    		else
		    		{	// zooom 1, 1.5, 2,3,4,5,7,10
		    			if(currentScale>8.)
		    			{	JNApplication.appBeep();
		    				return;
		    			}
		    			else if(currentScale>6.)
		    				currentScale = 10.;
		    			else if(currentScale>4.5)
		    				currentScale = 7.;
		    			else if(currentScale>1.75)
		    				currentScale += 1.;
		    			else
		    				currentScale +=0.5;
		    		}
		    	
		    		// get new checked scale item
		    		docView.checkedZoomItem.setSelected(false);
		    		if(currentScale<2.5)
		    			newItem = (int)(2*currentScale-1.99);
		    		else if(currentScale<6.)
		    			newItem = (int)(currentScale+0.1);
		    		else if(currentScale<8.)
		    			newItem = 6;
		    		else
		    			newItem = 7;
		    		docView.checkedZoomItem = (JCheckBoxMenuItem)docView.zoomMenu.getMenuComponent(newItem);
		    		docView.checkedZoomItem.setSelected(true);
		    	
		    		// change scale
		    		plotScroll.setScale(currentScale,e.getPoint());
		    	}
		    }
		});
		
		// listen to clicks in the plot view
		plotView.addMouseMotionListener(new MouseMotionListener()
		{	public void mouseMoved(MouseEvent e)
			{	if(e.isMetaDown())
				{	if(e.isAltDown())
					{	if(currentScale<1.25)
							plotView.setCursor(zoomNo);
						else
							plotView.setCursor(zoomOut);
					}
					else
					{	if(currentScale>8.)
							plotView.setCursor(zoomNo);
						else
							plotView.setCursor(zoomIn);
					}
				}
				else
				{	plotView.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
				}
		    }
		
			public void mouseDragged(MouseEvent e)
			{
			}
		});
		
		// scroll pane
		plotScroll=new MeshPlotScroll(plotView);
		add(plotScroll,BorderLayout.CENTER);
		
		// add controls
		movieControls=new MovieControls(500,resDoc,this,gDocView);
		add(movieControls,BorderLayout.SOUTH);
		
		NFMVPrefs.setWorkspace(chooser);
		
		// notifications
		JNNotificationCenter.getInstance().addNameAndObjectForTarget("PlotQuantityChanged",gDocView,this);
		JNNotificationCenter.getInstance().addNameAndObjectForTarget("ParticleSizeChanged",gDocView,this);
		JNNotificationCenter.getInstance().addNameAndObjectForTarget("MaxElongChanged",gDocView,this);
		
		// load cursors
		if(zoomIn==null)
		{	Toolkit toolkit = Toolkit.getDefaultToolkit();
			Class<?> baseClass=JNApplication.main.getClass();
			ImageIcon cIcon = new ImageIcon(baseClass.getResource("Resources/ZoomMinus.png"));
			zoomOut = toolkit.createCustomCursor(cIcon.getImage(),new Point(8,8),"zoomOut");
			cIcon = new ImageIcon(baseClass.getResource("Resources/ZoomPlus.png"));
			zoomIn = toolkit.createCustomCursor(cIcon.getImage(),new Point(8,8),"zoomIn");
			cIcon = new ImageIcon(baseClass.getResource("Resources/ZoomNo.png"));
			zoomNo = toolkit.createCustomCursor(cIcon.getImage(),new Point(8,8),"zoomNo");
		}
		
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
		{	while(docView.controls.incrementArchiveIndex() && runFlag)
			{	theTime=new Date();
			
				try
				{	loadPlotData();
				}
				catch(Exception e)
				{	Toolkit.getDefaultToolkit().beep();
					JNUtilities.showMessage(null,"Error loading plot data: "+e.getMessage());
					break;
				}
				catch(OutOfMemoryError me)
				{	Toolkit.getDefaultToolkit().beep();
					JNUtilities.showMessage(null,"Out of memory error: "+me.getMessage());
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
				{	File saveFile=new File(frameFileRoot+String.format("%04d.png",docView.controls.getArchiveIndex()));
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
	
	// stop movie if it is running
	protected void stopMovie()
	{	if(!runFlag) return;
		runFlag=false;
		movieControls.setPlaying(false);
	}
	
	// change to new time unless in a movie
	public void changeArchiveIndex(int newIndex)
	{	if(runFlag) return;
		try
		{	loadPlotData();
			plotView.repaint();
		}
		catch(Exception e)
		{	Toolkit.getDefaultToolkit().beep();
			JNUtilities.showMessage(null,"Error loading plot data: "+e.getMessage());
		}
		catch(OutOfMemoryError me)
		{	Toolkit.getDefaultToolkit().beep();
			JNUtilities.showMessage(null,"Out of memory error: "+me.getMessage());
		}
	}
	
	// start plot over at new initial index and possibly new component and options
	public void beginNewIndexNewComponent(int newIndex,int newComponent)
	{	movieComponent=newComponent;
		ControlPanel controls=docView.controls;
		plotView.setOptions(controls.getOptions(),movieComponent,getPlotType(),controls.getParticleSize());
		
		boolean sameIndex=!controls.setArchiveIndex(newIndex);
		if(sameIndex || resDoc.isFEAAnalysis())
		{	// means the index is same OR an FEA plot
			// This will load data and replot (same index reload needed in case other settings changed)
			stopMovie();
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
			else if(!docView.didInternalPlot)
			{	Toolkit.getDefaultToolkit().beep();
			}
			else
			{	// movie thread
				runFlag=true;
				Thread movieThread=new Thread(this);
				movieThread.start();
			}
		}
		
		else if(theCmd.equals("Rewind"))
			docView.controls.setArchiveIndex(0);
			
		else if(theCmd.equals("ExportFrame"))
		{	// get file if need it
			int result = chooser.showSaveDialog(this);
			if(result != JFileChooser.APPROVE_OPTION) return;
			File saveFile=chooser.getSelectedFile();
			String ext = JNUtilities.getExtension(saveFile);
			if(!ext.equalsIgnoreCase(".png"))
				saveFile = new File(saveFile.getParent(),saveFile.getName()+".png");
			saveFile = JNUtilities.CheckFileStatus(docView, saveFile);
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
			File saveFile=new File(frameFileRoot+String.format("%04d.png",docView.controls.getArchiveIndex()));
			if(!exportPlotFrame(saveFile)) return;
			
			// movie thread
			runFlag=true;
			Thread movieThread=new Thread(this);
			movieThread.start();
		}
		
		else if(theCmd.equals("ZoomPlot"))
		{	docView.checkedZoomItem.setSelected(false);
			docView.checkedZoomItem=(JCheckBoxMenuItem)(e.getSource());
			docView.checkedZoomItem.setSelected(true);
			String theZoomCmd=docView.checkedZoomItem.getText();
			if(theZoomCmd.equals("150%"))
				currentScale=1.5;
			else if(theZoomCmd.equals("200%"))
				currentScale=2.0;	
			else if(theZoomCmd.equals("300%"))
				currentScale=3.0;	
			else if(theZoomCmd.equals("400%"))
				currentScale=4.0;	
			else if(theZoomCmd.equals("500%"))
				currentScale=5.0;	
			else if(theZoomCmd.equals("700%"))
				currentScale=7.0;	
			else if(theZoomCmd.equals("1000%"))
				currentScale=10.0;
			else
				currentScale = 1.0;
			plotScroll.setScale(currentScale,null);
		}
		
		else if(theCmd.equals("ShowPrevious"))
		{	if(!docView.controls.decrementArchiveIndex())
				JNApplication.appBeep();
		}
		
		else if(theCmd.equals("ShowNext"))
		{	if(!docView.controls.incrementArchiveIndex())
				JNApplication.appBeep();
		}
	}
	
	public void receiveNotification(JNNotificationObject obj)
	{	if(obj.getName().equals("PlotQuantityChanged"))
		{	// first two when choose new quantity or component in a movie plot window
			JComboBox<?> qmenu = (JComboBox<?>)obj.getUserInfo();
			if(qmenu==movieControls.pquant)
			{	// change control panel to match type and index
				docView.controls.changePlotType(getPlotType());
				int newIndex = qmenu.getSelectedIndex();
				JComboBox<PlotMenuItem> plotQuant=docView.controls.getQuantityMenu();
				if(plotQuant.getSelectedIndex()!=newIndex)
					plotQuant.setSelectedIndex(newIndex);
			
				// adjust component menus without reploting
				boolean oldDisable = movieControls.disableStartPlot;
				movieControls.disableStartPlot = true;
				movieControls.setComponentMenu();
				movieControls.disableStartPlot = oldDisable;
			
				// plot if not disabled
				if(!movieControls.disableStartPlot)
					docView.startNewPlot(getPlotType());
			}
			else if(qmenu==movieControls.pcmpnt)
			{	// changed in plot window, synch with results window
				docView.controls.changePlotType(getPlotType());
				int newIndex=movieControls.pcmpnt.getSelectedIndex();
				JComboBox<String> plotCmpnt=docView.controls.getComponentMenu();
				if(plotCmpnt.getSelectedIndex()!=newIndex)
					plotCmpnt.setSelectedIndex(newIndex);
				if(!movieControls.disableStartPlot)
					docView.startNewPlot(getPlotType());
			}
			else
			{	// changed in control panel - synch and replot (might never be used)
				movieControls.syncPlotQuantityMenus();
				docView.startNewPlot(getPlotType());
			}
		}
	
		else if(obj.getName().equals("ParticleSizeChanged"))
		{	if((JSlider)obj.getUserInfo()==movieControls.mpmParticleSize)
			{	int newSize = movieControls.mpmParticleSize.getValue();
				JSlider masterSize = docView.controls.getParticleSizeSlider();
				masterSize.setValue(newSize);
				if(!movieControls.disableStartPlot)
					docView.startNewPlot(getPlotType());
			}
			else
			{	// changed in control panel - synch and replot (might never be used)
				movieControls.syncPlotQuantityMenus();
				docView.startNewPlot(getPlotType());
			}
		}
	
		else if(obj.getName().equals("MaxElongChanged"))
		{	if(!movieControls.disableStartPlot)
				docView.startNewPlot(getPlotType());
		}
	}
	
	// just picked new quantity in the menu the movie frame
	protected void plotQuantityChanged(JComboBox<PlotMenuItem> qmenu)
	{	System.out.println("new plot index="+qmenu.getSelectedIndex());
	
		// change control panel to match type and index
		docView.controls.changePlotType(getPlotType());
		int newIndex = qmenu.getSelectedIndex();
		JComboBox<PlotMenuItem> plotQuant=docView.controls.getQuantityMenu();
		if(plotQuant.getSelectedIndex()!=newIndex)
			plotQuant.setSelectedIndex(newIndex);
		
		// adjust component menus without reploting
		boolean oldDisable = movieControls.disableStartPlot;
		movieControls.disableStartPlot = true;
		movieControls.setComponentMenu();
		movieControls.disableStartPlot = oldDisable;
		
		System.out.println("plot="+newIndex+", quant="+plotQuant.getSelectedIndex());
		// plot if not disabled
		if(!movieControls.disableStartPlot)
			docView.startNewPlot(getPlotType());
		
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
			JNUtilities.showMessage(docView,"Error exporting graphic image: " + fe);
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
	// accessors
	//----------------------------------------------------------------------------
	
	public MeshPlotView getPlotView() { return plotView; }
	public boolean isMovieRunning() { return runFlag; }
	public int getPlotType() { return LoadArchive.PARTICLE_PLOT; }
	
}
