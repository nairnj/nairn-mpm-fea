/*******************************************************************
	NodalPoint.java
	NairnFEAMPMViz

	Created by John Nairn on 2/27/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

public class NodalPoint
{
	// variables
	private int num,numpts;
	public double x,y,z;
	public double dispx,dispy;
	private double plotValue,weight;
	public double sigxx,sigyy,sigzz,sigxy;		// FEA nodal stresses
	
	// initialize
	NodalPoint(int nodeNum,double xpt,double ypt)
	{   super();
		num=nodeNum;
		x=xpt;
		y=ypt;
		dispx=0.;
		dispy=0.;
		sigxx=sigyy=sigzz=sigxy=0.;
	}
	
	// initialize
	NodalPoint(int nodeNum,double xpt,double ypt,double zpt)
	{   super();
		num=nodeNum;
		x=xpt;
		y=ypt;
		z=zpt;
		dispx=0.;
		dispy=0.;
		sigxx=sigyy=sigzz=sigxy=0.;
		// dispz, sigxz,sigyz not used yet
	}

	// fill node with circle
	public void stroke(MeshPlotView pv,ResultsDocument doc,boolean displaced)
	{
		NodalPoint nd=doc.nodes.get(num-1);
		if(displaced)
			pv.moveTo(nd.x+nd.dispx,nd.y+nd.dispy);
		else
			pv.moveTo(nd.x,nd.y);
		pv.fillOval(0.1*doc.cellMinSide);
	}
	
	// draw the number
	public void number(MeshPlotView pv,ResultsDocument doc,boolean displaced)
	{
		NodalPoint nd=doc.nodes.get(num-1);
		if(displaced)
			pv.moveTo(nd.x+nd.dispx,nd.y+nd.dispy);
		else
			pv.moveTo(nd.x,nd.y);
		pv.nudge(2.,2.);
		pv.drawString(String.format("%d",num),MeshPlotView.JUST_DRAW);
	}
	
	// zero entries
	public void zeroPlotValue()
	{	plotValue=0.;
		weight=0.;
		numpts=0;
	}

	// increment nodal value
	public void addPlotValue(double theValue,double wt)
	{	plotValue+=theValue*wt;
		weight+=wt;
		numpts++;
	}

	// divide by weight to get final value
	public void normPlotValue()
	{	if(numpts>0)
			plotValue/=weight;
	}
	
	public double getPlotValue() { return plotValue; }
	public int getNum() { return num; }

}
