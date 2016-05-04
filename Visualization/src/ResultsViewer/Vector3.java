
/*
 * Vector3
 * NairnFEAMPMViz
 * 
 * Created 4/11/2016, John A. Nairn
 * Copyright 2016, RSAC Software, All Rights Reserved
 */

public class Vector3
{
	public double x,y,z;
	
	// empty and all zeroed
	Vector3()
	{	x = y = z = 0.;
	}
	
	// set components of vector
	Vector3(double vx,double vy,double vz)
	{	x = vx;
		y = vy;
		z = vz;
	}
	

}
