/*******************************************************************
	ColorPicker.java
	NairnFEAMPMViz

	Created by John Nairn on 2/28/07.
	Copyright 2007 RSAC Software. All rights reserved.
	
	No instance - all static variables and methods to hande
	color spectrum of plots
*******************************************************************/

import java.awt.*;

public class ColorPicker
{
	// variables and constants
	public static final int BLUE_TO_CYAN=0;
	public static final int BLUE_TO_RED=1;
	public static final int PURPLE_TO_ORANGE=2;
	public static final int GRAY_SCALE=3;
	public static final int COOL_DIVERGING=4;
	
	private static int spectrum=BLUE_TO_RED;
	public static int numberContours=1;
	private static boolean inverted=false;
	private static float brightness=1.0f;

	public static Color PickRainbow(double fraction)
	{	// make sure between 0 and 1
		fraction=Math.min(Math.max(0.,fraction),1.);
		if(numberContours>1)
		{	int trunc=(int)(numberContours*fraction);
			fraction=Math.min((double)trunc/(double)(numberContours-1),1.);
		}
		if(inverted) fraction=1.-fraction;
		
		// get the color
		Color theColor;
		int rgb;
		float hue;
		switch(spectrum)
		{	case GRAY_SCALE:
				theColor=new Color((float)fraction,(float)fraction,(float)fraction);
				break;
			case BLUE_TO_CYAN:
				hue=(float)(0.67+0.83*fraction);
				if(hue>1.0) hue-=1.0;
				rgb=Color.HSBtoRGB(hue,(float)1.0,brightness);
				theColor=new Color(rgb);
				break;
			case BLUE_TO_RED:
				hue=(float)(0.67*(1.-fraction));
				rgb=Color.HSBtoRGB(hue,(float)1.0,brightness);
				theColor=new Color(rgb);
				break;
			case PURPLE_TO_ORANGE:
				hue=(float)(0.777-0.702*fraction);
				rgb=Color.HSBtoRGB(hue,(float)1.0,brightness);
				theColor=new Color(rgb);
				break;
			case COOL_DIVERGING:
				// diverging colors from RGB = (35,61,181) to (172,36,32) through (223,223,223)
				//                           = (0.13672,0.23828,0.70703) to (0.67188,0.14063,0.125) through (0.87109 x 3)
				float red,green,blue;
				if(fraction<0.5)
				{	red=(float)(0.13672 + 1.4688*fraction);
					green=(float)(0.23828 + 1.2656*fraction);
					blue=(float)(0.70703 + 0.32813*fraction);
				}
				else
				{	red=(float)(1.0703 - 0.39843*fraction);
					green=(float)(1.6016 - 1.4609*fraction);
					blue=(float)(1.6172 - 1.4922*fraction);
				}
				theColor=new Color(red,green,blue);
				break;
			default:
				theColor=new Color((float)fraction,(float)fraction,(float)fraction);
				break;
		}
		return theColor;
	}
	
	// get/set current spectrum type
	public static int getSpectrumType() { return spectrum; }
	public static void setSpectrumType()
	{	int newType=NFMVPrefs.prefs.getInt(NFMVPrefs.SpectrumKey,NFMVPrefs.SpectrumDef);
		if(newType>=BLUE_TO_CYAN && newType<=COOL_DIVERGING)
			spectrum=newType;
	}
	
	// get/set number of contours
	public static int getNumberOfContours() { return numberContours; }
	public static void setNumberOfContours()
	{	int newNumber=NFMVPrefs.prefs.getInt(NFMVPrefs.NumContoursKey,NFMVPrefs.NumContoursDef);
		if(newNumber>=1)
			numberContours=newNumber;
	}
}
