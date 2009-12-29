/*******************************************************************
	MeshPlotWindow.java
	NairnFEAMPMViz

	Created by John Nairn on 8/21/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

public class MeshPlotWindow extends MoviePlotWindow
{
	static final long serialVersionUID=10L;

	// initialize
	public MeshPlotWindow(ResultsDocument gResDoc,DocViewer gDocView)
	{	super(gResDoc,gDocView);
		plotType=MeshPlotView.MPMMESH_PLOTS;
	}

	// load everything needed to plot or replot data
	public void loadPlotData() throws Exception
	{	resDoc.readSelectedArchive(movieControls.getArchiveIndex());
		ElementBase.loadPlotData(movieComponent,resDoc,plotView.getPlotType(),false);
	}
}
