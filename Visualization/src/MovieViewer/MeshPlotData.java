/*
 * MeshPlotData.java
 * NairnFEAMPMViz
 * 
 * Created by John Nairn on 4/8/07.
 * Copyright (c) 2007 RSAC Software. All rights reserved.
 */

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import java.awt.geom.*;
import geditcom.JNFramework.*;

public class MeshPlotData  extends PlotControl implements MouseMotionListener, JNNotificationListener
{
	static final long serialVersionUID=8L;
	
	public JLabel xlabel=new JLabel("x:");
	public JLabel ylabel=new JLabel("y:");
	public JLabel zlabel=new JLabel("z:");
	public JLabel nodelabel=new JLabel("nd:");
	public JLabel elemlabel=new JLabel("el:");
	public JLabel mplabel=new JLabel("mp:");
	public String zunits="";
	
	private MeshPlotView plotView;
	private ResultsDocument resDoc;
	
	private int currentElem=0;
	
	// initialize
	MeshPlotData(MeshPlotView gPlotView,ResultsDocument gResDoc)
	{   super(ControlPanel.WIDTH,48,gResDoc.docCtrl);
		setLayout(new GridLayout(3,2));

		Font labelFont=new Font("sanserif",Font.PLAIN,10);
		xlabel.setFont(labelFont);
		add(xlabel);
		elemlabel.setFont(labelFont);
		add(elemlabel);
		ylabel.setFont(labelFont);
		add(ylabel);
		nodelabel.setFont(labelFont);
		add(nodelabel);
		zlabel.setFont(labelFont);
		add(zlabel);
		mplabel.setFont(labelFont);
		add(mplabel);
		
		setBorder(BorderFactory.createEmptyBorder());		// of the PlotControl JPanel
		
		plotView=gPlotView;
		plotView.addMouseMotionListener(this);
		
		resDoc=gResDoc;
		
		// notifications
		JNNotificationCenter.getInstance().addNameAndObjectForTarget("PlotUnitsChanged",gResDoc.docCtrl,this);
	}
	
	// mouse motion events
	public void mouseDragged(MouseEvent e) {}
	
	public void mouseMoved(MouseEvent e)
	{	// skip until something is loaded
		if(!plotView.getFirstLoad()) return;
		
		int i;
		Point2D.Double pt=plotView.getCoords(e.getPoint());
		
		// coordinates
		xlabel.setText("x: "+JNUtilities.formatDouble(pt.x)+" "+resDoc.units.lengthUnits());
		ylabel.setText("y: "+JNUtilities.formatDouble(pt.y)+" "+resDoc.units.lengthUnits());
		
		// can't check while movie is running
		if(resDoc.docCtrl.getMovieFrame().isMovieRunning())
		{	elemlabel.setText("el: ...");
			nodelabel.setText("nd: ...");
			mplabel.setText("mp: ...");
			zlabel.setText("z: ...");
			return;
		}
		
		// element number
		if(plotView.inDisplaced() & resDoc.isFEAAnalysis())
		{	if(currentElem>0)
			{	// check previous one first
				if(!(resDoc.elements.get(currentElem-1)).PtInDispElement(pt))
					currentElem=0;
			}
			if(currentElem<=0)
			{	for(i=0;i<resDoc.elements.size();i++)
				{	if((resDoc.elements.get(i)).PtInDispElement(pt))
					{	currentElem=i+1;
						break;
					}
				}
			}
		}
		else
		{	if(currentElem>0)
			{	// check previous one first
				if(!(resDoc.elements.get(currentElem-1)).PtInElement(pt))
					currentElem=0;
			}
			if(currentElem<=0)
			{	for(i=0;i<resDoc.elements.size();i++)
				{	if((resDoc.elements.get(i)).PtInElement(pt))
					{	currentElem=i+1;
						break;
					}
				}
			}
		}
		
		int nearest;
		if(currentElem>0)
		{	elemlabel.setText("el: "+currentElem);
			nearest=(resDoc.elements.get(currentElem-1)).NearestNode(pt,resDoc,plotView.inDisplaced());
			nodelabel.setText("nd: "+nearest);
			
			// find nearest material point
			if(plotView.getPlotType()==LoadArchive.PARTICLE_PLOT)
			{	MaterialPoint mpt=resDoc.mpmPoints.get(0);
				double distance=pt.distanceSq(mpt.x,mpt.y);
				MaterialPoint nearmp=mpt;
				double approach;
				for(i=1;i<resDoc.mpmPoints.size();i++)
				{	mpt=resDoc.mpmPoints.get(i);
					approach=pt.distanceSq(mpt.x,mpt.y);
					if(approach<distance)
					{   distance=approach;
						nearmp=mpt;
					}
				}
				
				mplabel.setText("mp: "+nearmp.num);
				zlabel.setText("z: "+JNUtilities.formatDouble(nearmp.getPlotValue())+zunits);
			}
			else if(plotView.getPlotComponent()!=PlotQuantity.MESHONLY)
			{	mplabel.setText("mp: ---");
				ElementBase cElem=resDoc.elements.get(currentElem-1);
				zlabel.setText("z: "+JNUtilities.formatDouble(cElem.getValueAt(pt))+zunits);
			}
			else
			{   zlabel.setText("z: ---");
				mplabel.setText("mp: ---");
			}
		}
		else
		{	elemlabel.setText("el: ---");
			nodelabel.setText("nd: ---");
			mplabel.setText("mp: ---");
			zlabel.setText("z: ---");
		}
			
	}

	public void receiveNotification(JNNotificationObject obj)
	{	if(obj.getName().equals("PlotUnitsChanged"))
		{	setZunits();
		}
	}
	
	public void setZunits()
	{	String newUnits = PlotQuantity.plotUnits(resDoc.docCtrl.controls.getPlotComponent(-1),resDoc.units);
		if(newUnits.length()>0)
			zunits = new String(" "+newUnits);
		else
			zunits = "";
	}
}
