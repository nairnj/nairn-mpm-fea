/*
 * TimePlotWindow.java
 * NairnFEAMPMViz
 *
 * Created by John Nairn on 8/3/07.
 * Copyright 2007 RSAC Software. All rights reserved.
 */

import javax.swing.*;
import java.util.*;
import java.io.*;
import java.nio.*;
import java.awt.Toolkit;
import java.awt.geom.*;

public class TimePlotWindow extends TwoDPlotWindow implements Runnable
{
	static final long serialVersionUID=22L;
	
	// constants
	static final int DATA_PLOTTED=0;
	static final int FILE_ERROR=-1;
	
	// for plot thread calculations
	private int component;
	private ControlPanel controls;

	// initialize
	public TimePlotWindow(DocViewer parent)
	{	super(parent);
		setChildType("timeFrame");
	}
	
	// add another plot
	public void addPlot(ControlPanel plotControls) throws Exception
	{
		// the component to plot
		controls=plotControls;
		component=controls.adjustComponent(controls.getPlotComponent());
		ResultsDocument resDoc=((DocViewer)document).resDoc;
		
		// plot file of global results
		if(component==PlotQuantity.MPMGLOBALRESULTS)
		{	if(resDoc.globalArchive==null)
				throw new Exception("Global results file not found ");
				
			try
			{	FileReader fr=new FileReader(resDoc.globalArchive);
				char [] buffer=new char [(int)resDoc.globalArchive.length()];
				fr.read(buffer);
				plot2DView.readTable(new String(buffer));
				plot2DView.setXTitle("Time ("+resDoc.timeU+")");
				plot2DView.setYTitle("Global Quantity");
				setVisible(true);
				toFront();
			}
			catch (Exception e)
			{	throw new Exception("Could not load global results file:\n   " + e.getMessage());
			}
			
			return;
		}
		
		// detach thread to gather plot information
		Thread plot2DThread=new Thread(this);
		plot2DThread.start();
	}
	
	//----------------------------------------------------------------------------
	// detachable thread for loading time plot data
	//----------------------------------------------------------------------------
	
	public void run()
	{
		controls.enableProgress(((DocViewer)document).resDoc.archives.size());
		
		try
		{	switch(component)
			{	case PlotQuantity.MPMTOTSTRENERGY:
				case PlotQuantity.MPMTOTKINENERGY:
				case PlotQuantity.MPMTOTENERGY:
				case PlotQuantity.MPMTOTEXTWORK:
				case PlotQuantity.MPMTOTPOTENERGY:
				case PlotQuantity.MPMTOTPLASTICENERGY:
				case PlotQuantity.MPMTOTTHERMALENERGY:
				case PlotQuantity.MPMTOTELEMENTCROSSINGS:
					plotTotalEnergy(false);
					break;
				
				case PlotQuantity.MPMJ1:
				case PlotQuantity.MPMJ2:
				case PlotQuantity.MPMKI:
				case PlotQuantity.MPMKII:
				case PlotQuantity.MPMCRACKRELEASE:
				case PlotQuantity.MPMCRACKABSORB:
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
			JOptionPane.showMessageDialog(this,"Time plot error: "+pe.getMessage());
			if(plot2DView.getNumberOfPlots()==0) dispose();
			controls.disableProgress();
			return;
		}
		
		// axis labels
		if(plot2DView.getNumberOfPlots()<2)
		{	ResultsDocument resDoc=((DocViewer)document).resDoc;
			plot2DView.setXTitle("Time ("+resDoc.timeU+")");
			plot2DView.setYTitle(PlotQuantity.plotLabel(component,resDoc.distU,resDoc.timeU));
		}
		controls.disableProgress();
		setVisible(true);
		toFront();
	}
	
	// plot particle data
	private void plotParticleData() throws Exception
	{
		// settings
		ResultsDocument resDoc=((DocViewer)document).resDoc;
		int ptNum=controls.getParticleNumber();
		
		// trap averaged quantity
		if(ptNum<=0)
		{	plotTotalEnergy(true);
			return;
		}
		
		int i,npts=resDoc.archives.size();
		double angle=0.;
		
		// array for results
		ArrayList<Double> x=new ArrayList<Double>(npts);
		ArrayList<Double> y=new ArrayList<Double>(npts);
		
		// variables while decoding
		byte[] version=new byte[4];
		ByteBuffer bb;
		MaterialPoint mpm=new MaterialPoint(ptNum+1);
		
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
			mpm.readRecord(bb,mpmOrder,resDoc.lengthScale,resDoc.timeScale,resDoc.is3D());
			
			// find particle property and add to plot
			x.add(new Double(resDoc.archiveTimes.get(i)));
			y.add(new Double(mpm.getForPlot(component,angle,resDoc)));
		}
		
		if(x.size()==0)
			throw new Exception("No data found for that plot quantity");

		Hashtable<String,String> props = new Hashtable<String,String>();
		props.put("object.color",plot2DView.selectPlotColor());
		props.put("array.name",PlotQuantity.plotName(component)+" (pt "+ptNum+")");
		plot2DView.plotData(x,y,props);
	}
	
	// plot particle data
	private void plotTotalEnergy(boolean average) throws Exception
	{
		// settings
		ResultsDocument resDoc=((DocViewer)document).resDoc;
		double total,angle=0.;
		int p,i,npts=resDoc.archives.size();
		int matNumOption=-controls.getParticleNumber();		// 0 for all or (material #) for one material
		
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
			for(p=0;p<nummpms;p++)
			{	// read record in the file and convert to MaterialPoint class
				bb.position(headerLength+p*resDoc.recSize);
				mpm.readRecord(bb,mpmOrder,resDoc.lengthScale,resDoc.timeScale,resDoc.is3D());
			
				// exit if crack or continue if do not want this material
				if(mpm.material<0) break;
				if(matNumOption!=0)
				{	if(mpm.material!=matNumOption)
						continue;
				}
				
				// find particle property and add to total
				total+=mpm.getForPlot(component,angle,resDoc);
			}
			
			if(average) total/=(double)nummpms;
			
			// add total calculation to the  plot
			x.add(new Double(resDoc.archiveTimes.get(i)));
			y.add(new Double(total));
		}
		
		if(x.size()==0)
			throw new Exception("No data found for that plot quantity");

		Hashtable<String,String> props = new Hashtable<String,String>();
		props.put("object.color",plot2DView.selectPlotColor());
		String extra="";
		if(matNumOption!=0) extra=" (material "+matNumOption+")";
		props.put("array.name",PlotQuantity.plotName(component)+extra);
		plot2DView.plotData(x,y,props);
	}
	
	// plot particle data
	private void plotCrackData() throws Exception
	{
		// settings
		ResultsDocument resDoc=((DocViewer)document).resDoc;
		int i,npts=resDoc.archives.size();
		int crackNum=controls.getCrackNumber();
		int tipNum=controls.getCrackTip();
		
		// array for results
		ArrayList<Double> x=new ArrayList<Double>(npts);
		ArrayList<Double> y=new ArrayList<Double>(npts);
		ArrayList<Integer> crackEnds=new ArrayList<Integer>(20);

		// variables while decoding
		byte[] version=new byte[4];
		ByteBuffer bb;
		CrackSegment seg=new CrackSegment();
		int lastoffset,offset,tipOffset,endOffset;
		boolean foundTip;
		short matnum;
		double yvalue;
		Point2D.Double pt,lastPt,cod;
		
		// format
		char[] crackOrder=new char[ReadArchive.ARCH_MAXCRACKITEMS];
		resDoc.crackFormat.getChars(0,ReadArchive.ARCH_MAXCRACKITEMS,crackOrder,0);
		
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
			
			// find the first crack (matnum>0), remember ends
			offset=headerLength+(nummpms-1)*resDoc.recSize+ReadArchive.sizeofInt+ReadArchive.sizeofDouble;
			lastoffset=offset;
			crackEnds.clear();
			crackEnds.add(new Integer(lastoffset));
			while(offset>0)
			{	bb.position(offset);
				matnum=bb.getShort();
				if(matnum>0) break;
				offset-=resDoc.recSize;
				if(matnum==-1) crackEnds.add(new Integer(offset));
			}
			
			// find start or end of desired crack (if it exists)
			if(crackNum>crackEnds.size()-1) continue;
			if(tipNum==CrackSelector.CRACK_START || component==PlotQuantity.MPMLENGTH
						|| component==PlotQuantity.MPMDEBONDLENGTH)
			{	Integer offObj=crackEnds.get(crackEnds.size()-crackNum);
				offset=offObj.intValue()+resDoc.recSize;
				offObj=crackEnds.get(crackEnds.size()-crackNum-1);
				endOffset=offObj.intValue();
			}
			else
			{	Integer offObj=crackEnds.get(crackEnds.size()-crackNum-1);
				offset=offObj.intValue();
				offObj=crackEnds.get(crackEnds.size()-crackNum);
				endOffset=offObj.intValue()+resDoc.recSize;
			}
			offset-=(ReadArchive.sizeofInt+ReadArchive.sizeofDouble);
			endOffset-=(ReadArchive.sizeofInt+ReadArchive.sizeofDouble);
			
			// read segment
			bb.position(offset);
			seg.readRecord(bb,crackOrder,resDoc.lengthScale,resDoc.timeScale);
			
			// crack tip properties
			switch(component)
			{   case PlotQuantity.MPMJ1:
					yvalue=seg.J1;
					break;
					
				case PlotQuantity.MPMJ2:
					yvalue=seg.J2;
					break;
					
				case PlotQuantity.MPMKI:
					yvalue=seg.KI;
					break;

				case PlotQuantity.MPMKII:
					yvalue=seg.KII;
					break;

				case PlotQuantity.MPMLENGTH:
				case PlotQuantity.MPMDEBONDLENGTH:
					yvalue=0.;
					double bonded=0.;
					lastPt=seg.getMedianPosition();
					while(true)
					{   offset+=resDoc.recSize;
						if(offset>lastoffset) break;
						bb.position(offset);
						seg.readRecord(bb,crackOrder,resDoc.lengthScale,resDoc.timeScale);
						if(seg.startFlag==-1) break;
						pt=seg.getMedianPosition();
						double segLength=Math.sqrt((pt.x-lastPt.x)*(pt.x-lastPt.x) +
										(pt.y-lastPt.y)*(pt.y-lastPt.y));
						yvalue+=segLength;
						if(seg.tractionMaterial>0) bonded+=segLength;
						lastPt=pt;
					}
					if(component==PlotQuantity.MPMDEBONDLENGTH) yvalue-=bonded;
					break;
				
				case PlotQuantity.MPMCRACKRELEASE:
					yvalue=seg.release;
					break;
				
				case PlotQuantity.MPMCRACKABSORB:
					yvalue=seg.absorb;
					break;
				
				case PlotQuantity.MPMDEBONDNCTOD:
				case PlotQuantity.MPMDEBONDSCTOD:
					tipOffset=offset;
					foundTip=true;
					pt=seg.getPt();
					int tlCount=0;
					lastPt=new Point2D.Double(0.,0.);
					cod=new Point2D.Double(0.,0.);
					// scan to end of debond zone from this tip
					while(seg.tractionMaterial>0)
					{	lastPt=pt;
						cod=seg.getCOD();
						pt=seg.getPt();
						tlCount++;
						if(tipNum==CrackSelector.CRACK_START)
						{	tipOffset+=resDoc.recSize;
							if(tipOffset>endOffset)
							{	foundTip=false;
								break;
							}
						}
						else
						{	tipOffset-=resDoc.recSize;
							if(tipOffset<endOffset)
							{	foundTip=false;
								break;
							}
						}
						bb.position(tipOffset);
						seg.readRecord(bb,crackOrder,resDoc.lengthScale,resDoc.timeScale);
					}
					
					// not found in the crack (entire crack is traction law)
					if(!foundTip)
					{	yvalue=0.;
						break;
					}
					
					// if tlCount==0, then no traction at this crack tip, so just fall through and use regular crack tip cod
					//   otherwise do calculations
					//	Here pt and cod are at the debond tip. lastPt is at previous traction law or at tip if tlCount==1
					if(tlCount>0)
					{	pt=seg.getPt();
						double dx=lastPt.x-pt.x;
						double dy=lastPt.y-pt.y;
						double norm=Math.sqrt(dx*dx+dy*dy);
					
						if(component==PlotQuantity.MPMDEBONDNCTOD)
							yvalue=(-cod.x*dy + cod.y*dx)/norm;
						else
							yvalue=(cod.x*dx + cod.y*dy)/norm;
						break;
					}
				
				case PlotQuantity.MPMNORMALCTOD:
				case PlotQuantity.MPMSHEARCTOD:
					cod=seg.getCOD();
					pt=seg.getMedianPosition();

					// read previous segment
					offset = tipNum==CrackSelector.CRACK_START ? offset+resDoc.recSize : offset-resDoc.recSize;
					bb.position(offset);
					seg.readRecord(bb,crackOrder,resDoc.lengthScale,resDoc.timeScale);
					lastPt=seg.getMedianPosition();
					double dx=pt.x-lastPt.x;
					double dy=pt.y-lastPt.y;
					double norm=Math.sqrt(dx*dx+dy*dy);
					
					if(component==PlotQuantity.MPMNORMALCTOD || component==PlotQuantity.MPMDEBONDNCTOD)
						yvalue=(-cod.x*dy + cod.y*dx)/norm;
					else
						yvalue=(cod.x*dx + cod.y*dy)/norm;
					break;
				
				default:
					yvalue=0.;
					break;
			}
			
			// add crack value to the  plot
			x.add(new Double(resDoc.archiveTimes.get(i)));
			y.add(new Double(yvalue));
		}
		
		if(x.size()==0)
			throw new Exception("No data found for that plot quantity");

		Hashtable<String,String> props = new Hashtable<String,String>();
		props.put("object.color",plot2DView.selectPlotColor());
		String extra;
		if(component==PlotQuantity.MPMLENGTH || component==PlotQuantity.MPMDEBONDLENGTH)
			extra=" (crack "+crackNum+")";
		else if(tipNum==CrackSelector.CRACK_START)
			extra=" (start of crack "+crackNum+")";
		else
			extra=" (end of crack "+crackNum+")";
		props.put("array.name",PlotQuantity.plotName(component)+extra);
		plot2DView.plotData(x,y,props);
	}
}
