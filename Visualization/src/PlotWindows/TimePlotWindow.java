/*
 * TimePlotWindow.java
 * NairnFEAMPMViz
 *
 * Created by John Nairn on 8/3/07.
 * Copyright 2007 RSAC Software. All rights reserved.
 */

import geditcom.JNFramework.JNUtilities;
import java.util.*;

import javax.swing.SwingUtilities;

import java.io.*;
import java.nio.*;
import java.awt.Toolkit;

public class TimePlotWindow extends TwoDPlotWindow implements Runnable
{
	static final long serialVersionUID=22L;
	
	// constants
	static final int DATA_PLOTTED=0;
	static final int FILE_ERROR=-1;
	
	// for plot thread calculations
	private int component;
	private ControlPanel controls;
	private ISDictType settings;

	// initialize
	public TimePlotWindow(DocViewer parent)
	{	super(parent);
		setChildType("timeFrame");
	}
	
	// add another plot
	public void addPlot(ControlPanel plotControls,ISDictType plotSettings) throws Exception
	{
		// the component to plot
		settings=plotSettings;
		if(settings!=null)
			settings.gcis_setObjectforKey("timeplot","plottype");
		controls=plotControls;
		component=controls.getPlotComponent(LoadArchive.TIME_PLOT,true,settings);
		if(component<0)
			throw new Exception("The result to be plotted was not recognized");
		component=controls.adjustComponent(component,settings);
		ResultsDocument resDoc=((DocViewer)document).resDoc;
		
		// plot file of global results (never in scripted plot)
		if(component==PlotQuantity.MPMGLOBALRESULTS)
		{	if(resDoc.globalArchive==null)
				throw new Exception("Global results file not found ");
				
			try
			{	FileReader fr=new FileReader(resDoc.globalArchive);
				char [] buffer=new char [(int)resDoc.globalArchive.length()];
				fr.read(buffer);
				fr.close();
				String globalTable = new String(buffer);
				plot2DView.readTable(globalTable);
				plot2DView.setXTitle("Time ("+resDoc.units.timeUnits()+")");
				plot2DView.setYTitle("Global Quantity");
				SwingUtilities.invokeLater(new Runnable()
				{	public void run()
					{	setVisible(true);
						toFront();
					}
				});
			}
			catch (Exception e)
			{	throw new Exception("Could not load global results file:\n   " + e.getMessage());
			}
			
			return;
		}
		
		// plot file of global results (never in scripted file)
		else if(component==PlotQuantity.IMPORTANDPLOTFILE)
		{	importAndPlotFile();			
			return;
		}
		
		// detach thread to gather plot information
		Thread plot2DThread=new Thread(this);
		plot2DThread.start();
	}
	
	public boolean isPlotting()
	{	return controls.isPlotting();
	}
	
	//----------------------------------------------------------------------------
	// detachable thread for loading time plot data
	//----------------------------------------------------------------------------
	
	public void run()
	{
		controls.enableProgress(((DocViewer)document).resDoc.mpmArchives.size());
		
		try
		{	switch(component)
			{	case PlotQuantity.MPMTOTSTRENERGY:
				case PlotQuantity.MPMTOTKINENERGY:
				case PlotQuantity.MPMTOTENERGY:
				case PlotQuantity.MPMTOTWORKENERGY:
				case PlotQuantity.MPMTOTPLASTICENERGY:
				case PlotQuantity.MPMTOTHEATENERGY:
				case PlotQuantity.MPMTOTELEMENTCROSSINGS:
					plotTotalEnergy(controls.getAveraging(settings));
					break;
				
				case PlotQuantity.MPMJ1:
				case PlotQuantity.MPMJ2:
				case PlotQuantity.MPMKI:
				case PlotQuantity.MPMKII:
				case PlotQuantity.MPMMODEIFB:
				case PlotQuantity.MPMMODEIIFB:
				case PlotQuantity.MPMCZLENGTH:
				case PlotQuantity.MPMLENGTH:
				case PlotQuantity.MPMDEBONDLENGTH:
				case PlotQuantity.MPMNORMALCTOD:
				case PlotQuantity.MPMSHEARCTOD:
				case PlotQuantity.MPMDEBONDNCTOD:
				case PlotQuantity.MPMDEBONDSCTOD:
					plotCrackData();
					break;
				
				default:
					plotParticleData();
					break;
			}
		}
		catch(Exception pe)
		{	Toolkit.getDefaultToolkit().beep();
			JNUtilities.showMessage(this,"Time plot error: "+pe.getMessage());
			if(plot2DView.getNumberOfPlots()==0) dispose();
			controls.disableProgress();
			return;
		}
		
		// axis labels
		if(plot2DView.getNumberOfPlots()<2)
		{	ResultsDocument resDoc=((DocViewer)document).resDoc;
			plot2DView.setXTitle("Time ("+resDoc.units.timeUnits()+")");
			plot2DView.setYTitle(PlotQuantity.plotLabel(component,resDoc));
		}
		controls.disableProgress();
		SwingUtilities.invokeLater(new Runnable()
		{	public void run()
			{	setVisible(true);
				toFront();
			}
		});
	}
	
	// plot particle data
	private void plotParticleData() throws Exception
	{
		// settings
		ResultsDocument resDoc=((DocViewer)document).resDoc;
		int ptNum=controls.getParticleNumber(settings);
		
		// trap averaged quantity
		if(ptNum<=0)
		{	plotTotalEnergy(controls.getAveraging(settings));
			return;
		}
		
		int i,npts=resDoc.mpmArchives.size();
		double angle=0.;
		
		// array for results
		ArrayList<Double> x=new ArrayList<Double>(npts);
		ArrayList<Double> y=new ArrayList<Double>(npts);
		
		// variables while decoding
		byte[] version=new byte[4];
		ByteBuffer bb;
		MaterialPoint mpm=new MaterialPoint(ptNum);
		
		// format
		char[] mpmOrder=new char[ReadArchive.ARCH_MAXMPMITEMS];
		resDoc.archFormat.getChars(0,ReadArchive.ARCH_MAXMPMITEMS,mpmOrder,0);
		
		//---------------------------------------------------------
		// loop over all archives
		for(i=0;i<npts;i++)
		{	// adjust progress bar
			controls.setProgress(i+1);
		
			// open file (try to continue on errors
			try
			{	bb=resDoc.openSelectedArchive(i);
			}
			catch(Exception bbe)
			{	continue;
			}
			
			// check version
			bb.get(version);
			int headerLength=4;
			int vernum=version[3]-'0';
			if(vernum>=4)
				headerLength=64;
			else if(vernum!=3)
				throw new Exception("Archive file is too old for this tool");
			
			// number of records (some may be cracks)
			int nummpms=(int)((bb.remaining()+4-headerLength)/resDoc.recSize);
			if(ptNum>nummpms) continue;
			
			// read record in the file and convert to MaterialPoint class
			bb.position(headerLength+(ptNum-1)*resDoc.recSize);
			mpm.readRecord(bb,mpmOrder,resDoc.units,resDoc.is3D());
			
			// find particle property and add to plot
			x.add(new Double(resDoc.mpmArchives.get(i).getTime()));
			y.add(new Double(mpm.getForPlot(component,angle,resDoc)));
		}
		
		if(x.size()==0)
			throw new Exception("No data found for that plot quantity");

		Hashtable<String,String> props = new Hashtable<String,String>();
		props.put("object.color",plot2DView.selectPlotColor());
		props.put("array.name",PlotQuantity.plotName(component,resDoc)+" (pt "+ptNum+")");
		plot2DView.plotData(x,y,props);
	}
	
	// plot particle data
	private void plotTotalEnergy(boolean average) throws Exception
	{
		// settings
		ResultsDocument resDoc=((DocViewer)document).resDoc;
		double total,angle=0.;
		int p,i,npts=resDoc.mpmArchives.size();
		int matNumOption=-controls.getParticleNumber(settings);		// 0 for all or (material #) for one material
		
		// array for results
		ArrayList<Double> x=new ArrayList<Double>(npts);
		ArrayList<Double> y=new ArrayList<Double>(npts);
		
		// variables while decoding
		byte[] version=new byte[4];
		ByteBuffer bb;
		MaterialPoint mpm=new MaterialPoint(1);
		
		// format
		char[] mpmOrder=new char[ReadArchive.ARCH_MAXMPMITEMS];
		resDoc.archFormat.getChars(0,ReadArchive.ARCH_MAXMPMITEMS,mpmOrder,0);
		
		//---------------------------------------------------------
		// loop over all archives
		for(i=0;i<npts;i++)
		{	// adjust progress bar
			controls.setProgress(i+1);
		
			// open file (try to continue on errors
			try
			{	bb=resDoc.openSelectedArchive(i);
			}
			catch(Exception bbe)
			{	continue;
			}
			
			// check version
			bb.get(version);
			int headerLength=4;
			int vernum=version[3]-'0';
			if(vernum>=4)
				headerLength=64;
			else if(vernum!=3)
				throw new Exception("Archive file is too old for this tool");
			int nummpms=(int)((bb.remaining()+4-headerLength)/resDoc.recSize);
			
			// read each point
			total=0.;
			int numadded=0;
			double totalVol = 0.;
			for(p=0;p<nummpms;p++)
			{	// read record in the file and convert to MaterialPoint class
				bb.position(headerLength+p*resDoc.recSize);
				mpm.readRecord(bb,mpmOrder,resDoc.units,resDoc.is3D());
				mpm.setNum(p+1);
				if(mpm.inReservoir()) continue;
			
				// exit if crack or continue if do not want this material
				if(mpm.material<0) break;
				if(matNumOption!=0)
				{	if(mpm.material!=matNumOption)
						continue;
				}
				
				// get total or get volume-weighted average
				double pValue = mpm.getForPlot(component,angle,resDoc);
				if(!average)
				{	total += pValue;
				}
				else
				{	double Vp = mpm.getDeformedVolume(resDoc);
					total += Vp*pValue;
					totalVol += Vp;
				}
				
				numadded++;
			}
			
			// convert to average by dividing by total volume
			if(average && numadded>0) total/=totalVol;
			
			// add total calculation to the  plot
			x.add(new Double(resDoc.mpmArchives.get(i).getTime()));
			y.add(new Double(total));
		}
		
		if(x.size()==0)
			throw new Exception("No data found for that plot quantity");

		Hashtable<String,String> props = new Hashtable<String,String>();
		props.put("object.color",plot2DView.selectPlotColor());
		String extra="";
		if(matNumOption!=0) extra=" (material "+matNumOption+")";
		props.put("array.name",PlotQuantity.plotName(component,resDoc)+extra);
		plot2DView.plotData(x,y,props);
	}
	
	// plot particle data
	private void plotCrackData() throws Exception
	{
		// settings
		ResultsDocument resDoc=((DocViewer)document).resDoc;
		int npts=resDoc.mpmArchives.size();
		
		// read crack number and tip
		int crackNum=controls.getCrackNumber(settings);
		int tipNum=controls.getCrackTip(settings);
		
		// array for results
		ArrayList<Double> x=new ArrayList<Double>(npts);
		ArrayList<Double> y=new ArrayList<Double>(npts);
		
		// get data
		resDoc.getTimeCrackData(controls,component,crackNum,tipNum,x,y);
		
		// plot the data
		Hashtable<String,String> props = new Hashtable<String,String>();
		props.put("object.color",plot2DView.selectPlotColor());
		String extra;
		if(component==PlotQuantity.MPMLENGTH || component==PlotQuantity.MPMDEBONDLENGTH ||
				component==PlotQuantity.MPMMODEIFB || component==PlotQuantity.MPMMODEIIFB ||
				component==PlotQuantity.MPMCZLENGTH )
			extra=" (crack "+crackNum+")";
		else if(tipNum==CrackSelector.CRACK_START)
			extra=" (start of crack "+crackNum+")";
		else
			extra=" (end of crack "+crackNum+")";
		props.put("array.name",PlotQuantity.plotName(component,resDoc)+extra);
		plot2DView.plotData(x,y,props);
	}
	
}
