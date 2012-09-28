/*******************************************************************
	RigidMaterial.java
	NairnFEAMPMViz

	Created by John Nairn on 9/6/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

public class RigidMaterial extends MaterialBase
{

	// initialize
	RigidMaterial(String matName)
	{	super(matName,MaterialBase.RIGIDMATERIAL);
	}
	
	// this is the one rigid material
	public boolean isRigid() { return true; }

}
