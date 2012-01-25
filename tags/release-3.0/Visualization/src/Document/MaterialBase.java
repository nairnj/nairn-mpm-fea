/*******************************************************************
	MaterialBase.java
	NairnFEAMPMViz

	Created by John Nairn on 9/6/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.util.*;

public class MaterialBase
{
	public static final int ISOTROPIC=1;
	public static final int TRANSISO1=2;
	public static final int TRANSISO2=3;
	public static final int ORTHO=4;
	public static final int INTERFACEPARAMS=5;
	public static final int DUGDALE=6;
	public static final int VISCOELASTIC=7;
	public static final int MOONEYRIVLIN=8;
	public static final int VONMISES=9;
	public static final int BISTABLEISO=10;
	public static final int RIGIDMATERIAL=11;
	public static final int COHESIVEZONEMATERIAL=12;
	public static final int LINEARTRACTIONMATERIAL=13;
	public static final int CUBICTRACTIONMATERIAL=14;
	public static final int HILLPLASTIC=15;
	public static final int LASTMATERIALTYPE=16;
	public static final int UNKNOWNMATERIAL=17;
	
	protected int type;
	protected String name;
	protected double rho;
	
	// initialize
	MaterialBase(String matName,int matType)
	{	name=matName;
		type=matType;
	}
	
	// decode data
	public void decodeData(Scanner s)
	{	// scan to end
		while(s.hasNext())
		{	Scanner sline=new Scanner(s.next());
			sline.useLocale(Locale.US);
			if(!sline.hasNext()) break;
			if(sline.next().equals("rho="))
				rho=sline.nextDouble();
		}
	}
	
	// most are not rigid
	public boolean isRigid() { return false; }
}
