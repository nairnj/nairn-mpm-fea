/*******************************************************************
	MPMParticlePlotWindow.java
	NairnFEAMPMViz

	Created by John Nairn on 8/21/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

public class MPMParticlePlotWindow extends MoviePlotWindow
{
	static final long serialVersionUID=13L;
	
	// initialize
	public MPMParticlePlotWindow(ResultsDocument gResDoc,DocViewer gDocView)
	{	super(gResDoc,gDocView);
	}
	
	// load everything needed to plot or replat data
	public void loadPlotData() throws Exception
	{	resDoc.readSelectedArchive(movieControls.getArchiveIndex());
		MaterialPoint.loadPlotData(movieComponent,resDoc);
	}
	
}
