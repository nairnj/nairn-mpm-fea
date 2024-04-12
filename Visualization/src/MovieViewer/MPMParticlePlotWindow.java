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
	
	// load everything needed to plot or replot data
	public void loadPlotData() throws Exception
	{	int archIndex = docView.controls.getArchiveIndex();
		resDoc.readSelectedArchive(archIndex);
		plotView.dataTime = resDoc.mpmArchives.get(archIndex).getTime();
		MaterialPoint.loadPlotData(movieComponent,resDoc);
	}
	
}
