/*******************************************************************
	IsotropicMat.java
	NairnFEAMPMViz

	Created by John Nairn on 9/6/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.util.*;

public class IsotropicMat extends MaterialBase
{
	protected double E;
	protected double nu;
	protected double a;

	// initialize
	IsotropicMat(String matName)
	{	super(matName,MaterialBase.ISOTROPIC);
	}

	// decode data
	public void decodeData(Scanner s)
	{	// E, nu, and G (not needed)
		Scanner sline=new Scanner(s.next());
		sline.useLocale(Locale.US);
		sline.next();			// E
		sline.next();			// =
		E=sline.nextDouble();
		sline.next();			// v
		sline.next();			// =
		nu=sline.nextDouble();
		
		// a
		sline=new Scanner(s.next());
		sline.useLocale(Locale.US);
		sline.next();			// a
		sline.next();			// =
		a=sline.nextDouble();
		
		super.decodeData(s);
	}
}
