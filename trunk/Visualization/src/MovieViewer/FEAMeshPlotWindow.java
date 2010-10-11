/*
 * FEAMeshPlotWindow.java
 * NairnFEAMPMViz
 *
 * Created by John Nairn on 9/6/07.
 * Copyright 2007 RSAC Software. All rights reserved.
 */

public class FEAMeshPlotWindow extends MoviePlotWindow
{
	static final long serialVersionUID=5L;

	// initialize
	public FEAMeshPlotWindow(ResultsDocument gResDoc,DocViewer gDocView)
	{	super(gResDoc,gDocView);
	}

	// load everything needed to plot or replot data
	public void loadPlotData() throws Exception
	{	ElementBase.loadPlotData(movieComponent,resDoc,plotView.getPlotType(),plotView.inDisplaced());
	}
	
	// the plot type
	public int getPlotType() { return MeshPlotView.FEAMESH_PLOTS; }

}
