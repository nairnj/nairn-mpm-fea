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

public class MoviePlotWindow extends JNChildWindow implements  Runnable, IIOWriteProgressListener
{
	private static final long serialVersionUID = 1L;

	// constants
	public static final int FRAME_RATE=30;

	// variables
	protected ResultsDocument resDoc;
	protected MeshPlotView plotView;
	protected MovieControls movieControls;
	protected int movieComponent;
	protected boolean runFlag,imageStillWriting;
	protected MeshPlotScroll plotScroll;
	protected String frameFileRoot=null;
	protected JMenu zoomMenu;
	protected JCheckBoxMenuItem checkedZoomItem;
	double currentScale;

	// command file chooser
	private JFileChooser chooser=new JFileChooser();
	
	private static Cursor zoomIn = null;
	private static Cursor zoomOut = null;
	private static Cursor zoomNo = null;

	// initialize
	public MoviePlotWindow(ResultsDocument gResDoc,DocViewer gDocView)
	{	super(gDocView,gResDoc.name);
		setChildType("movieFrame");
		resDoc=gResDoc;
		setFramePrefs("Mesh Window Width",800,"Mesh Window Height",600);
		
		// size and location
		finishFrameworkWindow(false);
		Dimension d=getPreferredSize();
		
		makeMenuBar();
		
		// add plot view (entire window for now)
		Container content=getContentPane();
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
		    		checkedZoomItem.setSelected(false);
		    		if(currentScale<2.5)
		    			newItem = (int)(2*currentScale-1.99);
		    		else if(currentScale<6.)
		    			newItem = (int)(currentScale+0.1);
		    		else if(currentScale<8.)
		    			newItem = 6;
		    		else
		    			newItem = 7;
		    		checkedZoomItem = (JCheckBoxMenuItem)zoomMenu.getMenuComponent(newItem);
		    		checkedZoomItem.setSelected(true);
		    	
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
		content.add(plotScroll,BorderLayout.CENTER);
		
		// add controls
		movieControls=new MovieControls(d.width,resDoc,this,gDocView);
		content.add(movieControls,BorderLayout.SOUTH);
		
		NFMVPrefs.setWorkspace(chooser);
		
		// notifications
		JNNotificationCenter.getInstance().addNameAndObjectForTarget("TimeSliderChanged",gDocView,this);
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
		
		finishFrameworkWindow(false);
	}
	
	// make menu bar on launch
	protected void makeMenuBar()
	{	//Create the menu bar
		JMenuBar menuBar = new JMenuBar();
		
		// File menu
		JMenu menu = new JMenu("File");
		menuBar.add(menu);
		JNDocument.makeMenuItem(menu,"Export Graphics...","ExportFrame",this,KeyEvent.VK_G,true);
		JNDocument.makeMenuItem(menu,"Export Movie Frames...","ExportMovie",this,KeyEvent.VK_M,true);
		menu.addSeparator();
		JNDocument.makeMenuItem(menu,"Close","closeWindow",this,KeyEvent.VK_W);
		
		// Zoom menu
		zoomMenu = new JMenu("Zoom");
		menuBar.add(zoomMenu);
		checkedZoomItem=makeCheckBoxMenuItem(zoomMenu,"100%","ZoomPlot",this,KeyEvent.VK_1);
		makeCheckBoxMenuItem(zoomMenu,"150%","ZoomPlot",this,KeyEvent.VK_2);
		makeCheckBoxMenuItem(zoomMenu,"200%","ZoomPlot",this,KeyEvent.VK_3);
		makeCheckBoxMenuItem(zoomMenu,"300%","ZoomPlot",this,KeyEvent.VK_4);
		makeCheckBoxMenuItem(zoomMenu,"400%","ZoomPlot",this,KeyEvent.VK_5);
		makeCheckBoxMenuItem(zoomMenu,"500%","ZoomPlot",this,KeyEvent.VK_6);
		makeCheckBoxMenuItem(zoomMenu,"700%","ZoomPlot",this,KeyEvent.VK_7);
		makeCheckBoxMenuItem(zoomMenu,"1000%","ZoomPlot",this,KeyEvent.VK_8);
		checkedZoomItem.setSelected(true);
		currentScale = 1.0;
		zoomMenu.addSeparator();
		JNDocument.makeMenuItem(zoomMenu,"Previous Frame","ShowPrevious",this,KeyEvent.VK_LEFT);
		JNDocument.makeMenuItem(zoomMenu,"Next Frame","ShowNext",this,KeyEvent.VK_RIGHT);
		
		// Window
		menu = new JMenu("Window");
		menuBar.add(menu);
		if(JNApplication.isMacLNF())
		{	menu.add(JNApplication.main.getOpenHelpAction());
		}
		JNDocument.makeMenuItem(menu,"Analysis Results","ShowResults",this,KeyEvent.VK_D);
		menu.addSeparator();
		setWindowMenu(menu);
		
		// add the menu
		setJMenuBar(menuBar);
	}
	
	// create check box item for the Zoom menu
	protected JCheckBoxMenuItem makeCheckBoxMenuItem(JMenu menu,String menuTitle,String menuAction,ActionListener target,int mKey)
	{	JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem(menuTitle);
		menuItem.setActionCommand(menuAction);
		menuItem.addActionListener(target);
		if(mKey!=0)
			menuItem.setAccelerator(KeyStroke.getKeyStroke(mKey,JNApplication.menuKeyMask()));
		menu.add(menuItem);
		return menuItem;
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
					JNUtilities.showMessage(null,"Error loading plot data: "+e.getMessage());
					break;
				}
				catch(OutOfMemoryError me)
				{	dispose();
					Toolkit.getDefaultToolkit().beep();
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
			JNUtilities.showMessage(null,"Error loading plot data: "+e.getMessage());
		}
		catch(OutOfMemoryError me)
		{	dispose();
			Toolkit.getDefaultToolkit().beep();
			JNUtilities.showMessage(null,"Out of memory error: "+me.getMessage());
		}
	}
	
	// start plot over at new initial index and possibly new component and options
	public void beginNewIndexNewComponent(int newIndex,int newComponent)
	{	movieComponent=newComponent;
		ControlPanel controls=resDoc.docCtrl.controls;
		plotView.setOptions(controls.getOptions(),movieComponent,getPlotType(),controls.getParticleSize());
		
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
		{	document.toFront();
		}
		
		else if(theCmd.equals("ExportFrame"))
		{	// get file if need it
			int result = chooser.showSaveDialog(this);
			if(result != JFileChooser.APPROVE_OPTION) return;
			ArrayList<String> exts=new ArrayList<String>(1);
			exts.add("png");
			File saveFile=JNUtilities.CheckFileExtension(this,chooser.getSelectedFile(),exts,"Exported movie frame file");
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
		{	if(!movieControls.decrementArchiveIndex())
				JNApplication.appBeep();
		}
		
		else if(theCmd.equals("ShowNext"))
		{	if(!movieControls.incrementArchiveIndex())
				JNApplication.appBeep();
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
				((DocViewer)document).controls.setArchiveIndex(select.getValue());
			}
			else
			{	// if slider in doc view changed, change movie controller
				// if it causes a changed, this will trigger notification to call again and
				// load the data with above code
				movieControls.setArchiveIndex(select.getValue());
			}
		}
	
		else if(obj.getName().equals("PlotQuantityChanged"))
		{	// called when choose new quantity on component in a move plot window
			JComboBox<?> qmenu = (JComboBox<?>)obj.getUserInfo();
			if(qmenu==movieControls.pquant)
			{	// changed in plot window, synch with results window
				((DocViewer)document).controls.changePlotType(getPlotType());
				int newIndex=movieControls.pquant.getSelectedIndex();
				JComboBox<PlotMenuItem> plotQuant=((DocViewer)document).controls.getQuantityMenu();
				if(plotQuant.getSelectedIndex()!=newIndex)
					plotQuant.setSelectedIndex(newIndex);
				if(!movieControls.disableStartPlot)
					((DocViewer)document).startNewPlot(getPlotType());
			}
			else if(qmenu==movieControls.pcmpnt)
			{	// changed in plot window, synch with results window
				((DocViewer)document).controls.changePlotType(getPlotType());
				int newIndex=movieControls.pcmpnt.getSelectedIndex();
				JComboBox<String> plotCmpnt=((DocViewer)document).controls.getComponentMenu();
				if(plotCmpnt.getSelectedIndex()!=newIndex)
					plotCmpnt.setSelectedIndex(newIndex);
				if(!movieControls.disableStartPlot)
					((DocViewer)document).startNewPlot(getPlotType());
			}
			else
			{	// changed in control panel - synch and replot (might never be used)
				movieControls.syncPlotQuantityMenus();
				((DocViewer)document).startNewPlot(getPlotType());
			}
		}
	
		else if(obj.getName().equals("ParticleSizeChanged"))
		{	if((JSlider)obj.getUserInfo()==movieControls.mpmParticleSize)
			{	int newSize = movieControls.mpmParticleSize.getValue();
				JSlider masterSize = ((DocViewer)document).controls.getParticleSizeSlider();
				masterSize.setValue(newSize);
				if(!movieControls.disableStartPlot)
					((DocViewer)document).startNewPlot(getPlotType());
			}
			else
			{	// changed in control panel - synch and replot (might never be used)
				movieControls.syncPlotQuantityMenus();
				((DocViewer)document).startNewPlot(getPlotType());
			}
		}
	
		else if(obj.getName().equals("MaxElongChanged"))
		{	if(!movieControls.disableStartPlot)
				((DocViewer)document).startNewPlot(getPlotType());
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
			JNUtilities.showMessage(this,"Error exporting graphic image: " + fe);
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
	{	int newComponent=((DocViewer)document).controls.getPlotComponent(-1,false,null);
		beginNewIndexNewComponent(((DocViewer)document).controls.getArchiveIndex(),newComponent);
		plotView.setFirstLoad(true);
		super.windowOpened(e);
	}
	
	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------
	
	public MeshPlotView getPlotView() { return plotView; }
	public boolean isMovieRunning() { return runFlag; }
	public int getPlotType() { return LoadArchive.PARTICLE_PLOT; }
	
}
