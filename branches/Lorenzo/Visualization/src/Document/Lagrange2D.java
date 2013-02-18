import java.awt.geom.Point2D;

/*
 * Lagrange2D.java
 * NairnFEAMPMViz
 *
 * Created by John Nairn on 4/18/11.
 * Copyright 2011 RSAC Software. All rights reserved.
 */


public class Lagrange2D extends EightNodeIsoparam
{
	// initialize
	Lagrange2D(int elemNum,int [] nds)
	{	super(elemNum,nds);
	}
	
	//--------------------------------------------------------------
	// Shape Function and (xi,eta) Methods needing overrides
	//--------------------------------------------------------------

	// shape function only
	public void getShapeFunction(double[] sfxn,Point2D.Double xiEta,NodalPoint[] ndpt)
	{
		sfxn[0]=0.25 * (xiEta.x-1.)*xiEta.x * (xiEta.y-1.)*xiEta.y;		// (-1,-1)
		sfxn[1]=0.25 * (xiEta.x+1.)*xiEta.x * (xiEta.y-1.)*xiEta.y;		// (1,-1)
		sfxn[2]=0.25 * (xiEta.x+1.)*xiEta.x * (xiEta.y+1.)*xiEta.y;		// (1,1)
		sfxn[3]=0.25 * (xiEta.x-1.)*xiEta.x * (xiEta.y+1.)*xiEta.y;		// (-1,1)
		sfxn[4]=0.5 * (1.-xiEta.x*xiEta.x) * (xiEta.y-1.)*xiEta.y;		// (0,-1)
		sfxn[5]=0.5 * (xiEta.x+1.)*xiEta.x * (1.-xiEta.y*xiEta.y);		// (1,0)
		sfxn[6]=0.5 * (1.-xiEta.x*xiEta.x) * (xiEta.y+1.)*xiEta.y;		// (0,1)
		sfxn[7]=0.5 * (xiEta.x-1.)*xiEta.x * (1.-xiEta.y*xiEta.y);		// (-1,0)
		sfxn[8]=(1.-xiEta.x*xiEta.x) * (1.-xiEta.y*xiEta.y);			// (0,0)
	}
	
	// general shape function BMatrix evaluation
	public void getShapeBMatrix(Point2D.Double xiEta,double[] xiDeriv,double[] etaDeriv,double[] asbe,NodalPoint[] eNodes,boolean derivs)
	{
		double sfxn[]=new double[9];
		sfxn[0]=0.25 * (xiEta.x-1.)*xiEta.x * (xiEta.y-1.)*xiEta.y;		// (-1,-1)
		sfxn[1]=0.25 * (xiEta.x+1.)*xiEta.x * (xiEta.y-1.)*xiEta.y;		// (1,-1)
		sfxn[2]=0.25 * (xiEta.x+1.)*xiEta.x * (xiEta.y+1.)*xiEta.y;		// (1,1)
		sfxn[3]=0.25 * (xiEta.x-1.)*xiEta.x * (xiEta.y+1.)*xiEta.y;		// (-1,1)
		sfxn[4]=0.5 * (1.-xiEta.x*xiEta.x) * (xiEta.y-1.)*xiEta.y;		// (0,-1)
		sfxn[5]=0.5 * (xiEta.x+1.)*xiEta.x * (1.-xiEta.y*xiEta.y);		// (1,0)
		sfxn[6]=0.5 * (1.-xiEta.x*xiEta.x) * (xiEta.y+1.)*xiEta.y;		// (0,1)
		sfxn[7]=0.5 * (xiEta.x-1.)*xiEta.x * (1.-xiEta.y*xiEta.y);		// (-1,0)
		sfxn[8]=(1.-xiEta.x*xiEta.x) * (1.-xiEta.y*xiEta.y);			// (0,0)

		xiDeriv[0]=(2.*xiEta.x-1.) * (xiEta.y-1.)*xiEta.y/4.;		// (-1,-1)
		xiDeriv[1]=(2.*xiEta.x+1.) * (xiEta.y-1.)*xiEta.y/4.;		// (1,-1)
		xiDeriv[2]=(2.*xiEta.x+1.) * (xiEta.y+1.)*xiEta.y/4.;		// (1,1)
		xiDeriv[3]=(2.*xiEta.x-1.) * (xiEta.y+1.)*xiEta.y/4.;		// (-1,1)
		xiDeriv[4]=(-2.*xiEta.x) * (xiEta.y-1.)*xiEta.y/2.;		// (0,-1)
		xiDeriv[5]=(2.*xiEta.x+1.) * (1.-xiEta.y*xiEta.y)/2.;		// (1,0)
		xiDeriv[6]=(-2.*xiEta.x) * (xiEta.y+1.)*xiEta.y/2.;		// (0,1)
		xiDeriv[7]=(2.*xiEta.x-1.) * (1.-xiEta.y*xiEta.y)/2.;		// (-1,0)
		xiDeriv[8]=(-2.*xiEta.x) * (1.-xiEta.y*xiEta.y);			// (0,0)
	
		etaDeriv[0]=(xiEta.x-1.)*xiEta.x * (2.*xiEta.y-1.)/4.;	// (-1,-1)
		etaDeriv[1]=(xiEta.x+1.)*xiEta.x * (2.*xiEta.y-1.)/4.;	// (1,-1)
		etaDeriv[2]=(xiEta.x+1.)*xiEta.x * (2.*xiEta.y+1.)/4.;	// (1,1)
		etaDeriv[3]=(xiEta.x-1.)*xiEta.x * (2.*xiEta.y+1.)/4.;	// (-1,1)
		etaDeriv[4]=(1.-xiEta.x*xiEta.x) * (2.*xiEta.y-1.)/2.;	// (0,-1)
		etaDeriv[5]=(xiEta.x+1.)*xiEta.x * (-2.*xiEta.y)/2.;		// (1,0)
		etaDeriv[6]=(1.-xiEta.x*xiEta.x) * (2.*xiEta.y+1.)/2.;	// (0,1)
		etaDeriv[7]=(xiEta.x-1.)*xiEta.x * (-2.*xiEta.y)/2.;		// (-1,0)
		etaDeriv[8]=(1.-xiEta.x*xiEta.x) * (-2.*xiEta.y);			// (0,0)
		
		if(derivs) return;
		
		// Get B matrix elements in xiDeriv[] and etaDeriv[] and Ni/r in sfxn[]
		// Find Jacobian Matrix
		int i;
		double temp1;
		double[][] jac={{0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.}};
		jac[1][1]=0.;
		jac[1][2]=0.;
		jac[2][1]=0.;
		jac[2][2]=0.;
		for(i=0;i<9;i++)
		{	jac[1][1]+=xiDeriv[i]*eNodes[i].x;
			jac[1][2]+=xiDeriv[i]*eNodes[i].y;
			jac[2][1]+=etaDeriv[i]*eNodes[i].x;
			jac[2][2]+=etaDeriv[i]*eNodes[i].y;
		}
		double detjac=jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1];
		temp1=jac[1][1]/detjac;
		jac[1][1]=jac[2][2]/detjac;
		jac[2][2]=temp1;
		jac[1][2]=-jac[1][2]/detjac;
		jac[2][1]=-jac[2][1]/detjac;
		
		// For radial position for axisymmetric analyses
		double asr=0.;
		for(i=0;i<9;i++)
			asr+=sfxn[i]*eNodes[i].x;
			
		/* Load J(-1)(1,1) dNi/dxi + J(-1)(1,2) dNi/deta into xiDeriv[i]
				J(-1)(2,1) dNi/dxi + J(-1)(2,2) dNi/deta into etaDeriv[i]
				Ni/r into asbe[i] */
		double temp2;
		for(i=0;i<9;i++)
		{	temp1=xiDeriv[i];
			temp2=etaDeriv[i];
			xiDeriv[i]=jac[1][1]*temp1 + jac[1][2]*temp2;
			etaDeriv[i]=jac[2][1]*temp1 + jac[2][2]*temp2;
			if(asr!=0.)
				asbe[i]=sfxn[i]/asr;
			else
				asbe[i]=0.;
		}
	}

	//--------------------------------------------------------------
	// Accessors probably needing overrides
	//--------------------------------------------------------------

	// nodes, sides, and mid side nodes (override as needed)
	public int getNumberNodes() { return 9; }
}
