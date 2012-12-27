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
	public static final int JOHNSONCOOK=16;
	public static final int MGSCGLMATERIAL=17;
	public static final int SLMATERIAL=18;
	public static final int WOODMATERIAL=19;
	public static final int TRILINEARTRACTIONLAW=20;
	public static final int HEANISOTROPIC=21;
	public static final int IDEALGAS=22;
	public static final int COUPLESAWTOOTHTRACTIONLAW=23;
	public static final int HEISOTROPIC=24;
	public static final int LASTMATERIALTYPE=25;
	public static final int UNKNOWNMATERIAL=26;
	
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
	
	// true it this material archives plastic strain that
	// should be added to elastic strain when finding deformation
	// gradient. It is onlu low-strain plasticity materials
	public boolean hasPlasticStrainForGradient()
	{ 	switch(type)
		{	case DUGDALE:
			case VONMISES:
			case HILLPLASTIC:
			case JOHNSONCOOK:
			case MGSCGLMATERIAL:
			case SLMATERIAL:
			case WOODMATERIAL:
				return true;
			default:
				break;
		}
		return false;
	}

}
