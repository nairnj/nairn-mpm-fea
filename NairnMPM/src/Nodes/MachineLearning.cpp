/********************************************************************************
	MachineLearning.cpp for MPM only
 	Part of CrackVelocityFieldMulti class
	nairn-mpm-fea

	Created by John Nairn on Mar 20, 2018.
	Copyright (c) 2018 John A. Nairn, All rights reserved.
 
	Linear, Logistic, and SVN methods for contact normals
********************************************************************************/
#if defined ( _MSC_VER) || defined (__APPLE__) 
#include "stdafx.h"
#endif

#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/MaterialContactNode.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Exceptions/MPMWarnings.hpp"


#include <iostream>
using namespace std;

// To print convergence progress
//#define SHOW_CONVERGENCE

// To add penalty to normal components (only 2D for testing)
#define PENALIZE_NORMALS

// To use IRLS instead of NLLS method (2D only)
//#define USE_IRLS

// To weight by inverse number of particles (only 2D for testing)
//#define WEIGHT_BY_NUMBER

// Penalties for offset in linear regression
double offsetLambdaLinReg = 0.1;

// Penalties for normals (optional) and offset in logistic regression
// Also scale penalties after linear regression step for use in logistic regression steps
#ifdef PENALIZE_NORMALS
// normals and offset (normals in in for testing)
double normLambdaLogReg = 1.e-7;
double offsetLambdaLogReg = 0.;

// Scaling to change penalty between linear step and logistic steps
double linearPenaltyScaling = 1.;
#else
// normal components are not penalized (my notes suggest 5e-4 when only penalize offset)
double offsetLambdaLogReg = 5.e-4;

// Scaling to change penalty between linear step and logistic steps
double linearPenaltyScaling = 0.25;		// my notes suggest 0.25, but zero might be better
#endif
double linearBetaScaling = 2.;

// Logisitic regression convergence conditions
double tolerance = 1.e-5;
int maxIter = 15;

// Linear regression method
// if b<0 then treat all other materials is matrerial b
Vector CrackVelocityFieldMulti::LinearRegressionNormal(MaterialContactNode *mcn,int vfld,int a,int b,
													   bool &hasDeln,Vector *delMats)
{
	Vector norm,xp;
	double mi;
	
	// needed to get nodal point location
	NodalPoint *ndptr = mcn->GetTheNode();
	
	// get linked particles
	vector< int > list = mcn->ParticleLists(vfld);
	long linkedParticles = (long)list.size();
	
	if(fmobj->IsThreeD())
	{	// 3D linear regression matrix
		double mt[10];
		for(long i=0;i<10;i++) mt[i] = 0.;
		double rhs[4];
		for(long i=0;i<4;i++) rhs[i] = 0.;
		double lambda[4];
		for(long i=0;i<3;i++) lambda[i] = 0.;
		lambda[3] = offsetLambdaLinReg;
		
		for(long i=0;i<linkedParticles;i++)
		{	MPMBase *mpmptr = mpm[list[i]];
			const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
			int matfld = matID->GetField();
			if(matfld==a)
				mi = -1.;
			else if(matfld==b || b<0)
				mi = 1.;
			else
				continue;
			xp = mpmptr->pos;
			xp.x -= ndptr->x;
			xp.y -= ndptr->y;
			xp.z -= ndptr->z;
			
			// unweighted symmetric matrix
			mt[0] += xp.x*xp.x;
			mt[1] += xp.x*xp.y;
			mt[2] += xp.x*xp.z;
			mt[3] += xp.x;
			mt[4] += xp.y*xp.y;
			mt[5] += xp.y*xp.z;
			mt[6] += xp.y;
			mt[7] += xp.z*xp.z;
			mt[8] += xp.z;
			mt[9] += 1.;
			
			// unweighted right hand size
			rhs[0] += mi*xp.x;
			rhs[1] += mi*xp.y;
			rhs[2] += mi*xp.z;
			rhs[3] += mi;
		}
		
		// solve - normal direction in x and y components
		Matrix4 m(mt[0]+lambda[0], mt[1], 			mt[2], 			 mt[3],
								   mt[4]+lambda[1], mt[5], 			 mt[6],
													mt[7]+lambda[2], mt[8],
									   								 mt[9]+lambda[3]);
		Matrix4 minv = m.Inverse();
		double p[4];
		minv.Times(rhs,p);
		norm = MakeVector(p[0],p[1],p[2]);
		
	}
	else
	{	// 2D linear regression matrix
		double mt[6];
		for(long i=0;i<6;i++) mt[i] = 0.;
		Vector rhs;
		ZeroVector(&rhs);
		
		// penalty
		Vector lambda;
		lambda.x = 0.;
		lambda.y = 0.;
		lambda.z = offsetLambdaLinReg;

		for(long i=0;i<linkedParticles;i++)
		{	MPMBase *mpmptr = mpm[list[i]];
			const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
			int matfld = matID->GetField();
			if(matfld==a)
				mi = -1.;
			else if(matfld==b || b<0)
				mi = 1.;
			else
				continue;
			xp = mpmptr->pos;
			xp.x -= ndptr->x;
			xp.y -= ndptr->y;
			
			// unweighted symmetric matrix
			mt[0] += xp.x*xp.x;
			mt[1] += xp.x*xp.y;
			mt[2] += xp.x;
			mt[3] += xp.y*xp.y;
			mt[4] += xp.y;
			mt[5] += 1.;
			
			// unweighted right hand size
			rhs.x += mi*xp.x;
			rhs.y += mi*xp.y;
			rhs.z += mi;
		}
		
		// solve - normal direction in x and y components
		Matrix3 m(mt[0]+lambda.x, mt[1], 		  mt[2],
						 		  mt[3]+lambda.y, mt[4],
												  mt[5]+lambda.z);
		Matrix3 minv = m.Inverse();
		norm = minv.Times(&rhs);
		norm.z = 0.;
	}
	
	// normalize
	double magi=sqrt(DotVectors(&norm,&norm));
	ScaleVector(&norm,1./magi);
	
	// find separation
	double supportForA = 1.e30;
	double supportForB = 1.e30;
	for(long i=0;i<linkedParticles;i++)
	{	MPMBase *mpmptr = mpm[list[i]];
		const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
		int matfld = matID->GetField();
		Vector xp = mpmptr->pos;
		if(matfld==a)
		{	// find minimum of t (t>0 on material a side of the line) = min(-Xp.n-Rp) = -max(Xp.n+Rp)
			double t = -DotVectors(&norm,&xp)-mpmptr->GetDeformedRadius(&norm);
			if(t<supportForA) supportForA = t;
		}
		else if(matfld==b || b<0)
		{	// find minimum of -t (t<0 on material b side of the line) = min(Xj.n - Rj)
			double t = DotVectors(&norm,&xp)-mpmptr->GetDeformedRadius(&norm) ;
			if(t<supportForB) supportForB = t;
		}
	}
	
	// compared to contact paper
	//    d((xa-xi).n) = -supportForA-xi.n
	//    d((xb-xi).n) =  supportForB-xi.n
	//    dLR = d((xb-xi).n) - d((xa-xi).n) = supportForB+suppportForA
	// average distance to the node is
	//    dist = 0.5*(|d((xa-xi).n)|+|d((xb-xi).n)|)
	
	// add distances accounting for particle deformation
	delMats->x = supportForA + supportForB;
	hasDeln = true;
	
	// get signed distances to node
	double xiDotn = norm.x*ndptr->x + norm.y*ndptr->y + norm.z*ndptr->z;
	delMats->y = -supportForA - xiDotn;
	delMats->z = supportForB - xiDotn;
	
	return norm;
}

// Logistic regression method
// if b<0 then treat all other materials is matrerial b
Vector CrackVelocityFieldMulti::LogisticRegressionNormal(MaterialContactNode *mcn,int vfld,int a,int b,
														 bool &hasDeln,Vector *delMats)
{
	Vector norm,xp,linNorm;
	double deltaOneMinusCos = 1.;
	int iter=0;

	// needed to get nodal point location
	NodalPoint *ndptr = mcn->GetTheNode();

#ifdef PENALIZE_NORMALS
	// get cell^2
	double cell2 = mpmgrid.GetMinCellDimension();
	cell2 *= cell2;
#endif
	
	// get linked particles
	vector< int > list = mcn->ParticleLists(vfld);
	long linkedParticles = (long)list.size();
	int *mis = new int[linkedParticles];
	
	// 3D or 2D?
	if(fmobj->IsThreeD())
	{	// start logisitic with linear regression and get particle values (mi)
		
		// zero matrix and right hand side
		double mt[10],rhs[4],lambda[4];;
		for(long i=0;i<10;i++) mt[i] = 0.;
#ifdef PENALIZE_NORMALS
		double normPenalty = normLambdaLogReg*cell2;
		for(long i=0;i<3;i++)
		{	rhs[i] = 0.;
			lambda[i] = normPenalty;
		}
#else
		for(long i=0;i<3;i++)
		{	rhs[i] = 0.;
			lambda[i] = 0.;
		}
#endif
		rhs[3] = 0.;
		lambda[3] = offsetLambdaLogReg;

		for(long i=0;i<linkedParticles;i++)
		{	MPMBase *mpmptr = mpm[list[i]];
			const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
			int matfld = matID->GetField();
			if(matfld==a)
			{	mis[i] = -1;
			}
			else if(matfld==b || b<0)
			{	mis[i] = 1;
			}
			else
			{	mis[i] = 0;
				continue;
			}
			xp = mpmptr->pos;
			xp.x -= ndptr->x;
			xp.y -= ndptr->y;
			xp.z -= ndptr->z;
			
			// unweighted symmetric matrix
			mt[0] += xp.x*xp.x;
			mt[1] += xp.x*xp.y;
			mt[2] += xp.x*xp.z;
			mt[3] += xp.x;
			mt[4] += xp.y*xp.y;
			mt[5] += xp.y*xp.z;
			mt[6] += xp.y;
			mt[7] += xp.z*xp.z;
			mt[8] += xp.z;
			mt[9] += 1.;
			
			// unweighted right hand size
			double mi = (double)mis[i];
			rhs[0] += mi*xp.x;
			rhs[1] += mi*xp.y;
			rhs[2] += mi*xp.z;
			rhs[3] += mi;
		}
		
		// solve - normal direction in x and y components
		Matrix4 m(mt[0]+lambda[0], mt[1], 			mt[2], 			 mt[3],
								   mt[4]+lambda[1], mt[5], 			 mt[6],
												    mt[7]+lambda[2], mt[8],
																	 mt[9]+lambda[3]);
		Matrix4 minv = m.Inverse();
		double beta[4];
		minv.Times(rhs,beta);
		
		// rescale beta by 2 and lambda for improved initial logistic guess
		for(long i=0;i<4;i++)
		{	beta[i] *= linearBetaScaling;
			lambda[i] *= linearPenaltyScaling;
		}
		
		// initial norm from linear regression
		Vector prevNorm = MakeVector(beta[0],beta[1],beta[2]);
		double nmag = sqrt(DotVectors(&prevNorm,&prevNorm));
		ScaleVector(&prevNorm,1./nmag);
		linNorm = prevNorm;			// save for fallback
		
#ifdef SHOW_CONVERGENCE
		cout << "# Node " << ndptr->num << " particles " << linkedParticles << endl;
		PrintVector("# LR ",&prevNorm);
		cout << " (" << beta[0] << "," << beta[1] << "," << beta[2] << "," << beta[3] << ")" << endl;
#endif
		
		// logistic analysis nonlinear least squares iteration until done
		while(deltaOneMinusCos>tolerance && iter<maxIter)
		{	// zero matrix and vector elements
			for(long i=0;i<10;i++) mt[i] = 0.;
			for(long i=0;i<4;i++) rhs[i] = -lambda[i]*beta[i];
			
			// build J^TWJ and J^TW(M-f) (for now, w_i = 1)
			for(long i=0;i<linkedParticles;i++)
			{	if(mis[i]==0) continue;
				xp = mpm[list[i]]->pos;
				xp.x -= ndptr->x;
				xp.y -= ndptr->y;
				xp.z -= ndptr->z;
				
				double dotp = xp.x*beta[0] + xp.y*beta[1] + xp.z*beta[2] + beta[3];
				double arg = 1.+exp(-dotp);
				double phi = 2.*(arg-1.)/(arg*arg);
				double phi2 = phi*phi;
				
				// unweighted symmetric matrix
				mt[0] += phi2*xp.x*xp.x;
				mt[1] += phi2*xp.x*xp.y;
				mt[2] += phi2*xp.x*xp.z;
				mt[3] += phi2*xp.x;
				mt[4] += phi2*xp.y*xp.y;
				mt[5] += phi2*xp.y*xp.z;
				mt[6] += phi2*xp.y;
				mt[7] += phi2*xp.z*xp.z;
				mt[8] += phi2*xp.z;
				mt[9] += phi2;
				
				// unweighted right hand side
				double mifi = ((double)mis[i] - 2./arg + 1.)*phi;
				rhs[0] += mifi*xp.x;
				rhs[1] += mifi*xp.y;
				rhs[2] += mifi*xp.z;
				rhs[3] += mifi;
			}
			
			// Matrix (J^TWJ + Lambda)
			Matrix4 m(mt[0]+lambda[0], mt[1], 			mt[2], 			 mt[3],
					  				   mt[4]+lambda[1], mt[5], 			 mt[6],
					  									mt[7]+lambda[2], mt[8],
					  													 mt[9]+lambda[3]);
			// next increment added to beta
			Matrix4 minv = m.Inverse();
			double dbeta[4];
			minv.Times(rhs,dbeta);
			for(long i=0;i<4;i++) beta[i] += dbeta[i];
			
			// find normal normal and its angle with previous normal
			norm = MakeVector(beta[0],beta[1],beta[2]);
			nmag = sqrt(DotVectors(&norm,&norm));
			ScaleVector(&norm,1./nmag);
			deltaOneMinusCos = 1. - DotVectors(&norm,&prevNorm);
			
#ifdef SHOW_CONVERGENCE
			cout << "# " << iter;
			PrintVector(": ",&norm);
			cout << " (" << beta[0] << "," << beta[1] << "," << beta[2] << "," << beta[3] << ") ";
			cout << deltaOneMinusCos << endl;
#endif
			
			// on to next iteration
			iter++;
			prevNorm = norm;
		}
	}
	else
	{	// start logistic with linear regression and get particle values (mi)
		
		// zero matrix and right hand side
		double mt[6];
		for(long i=0;i<6;i++) mt[i] = 0.;
		Vector rhs = MakeVector(0.,0.,0.);
		
		// penalty
#ifdef PENALIZE_NORMALS
		double normPenalty = normLambdaLogReg*cell2;
		Vector lambda = MakeVector(normPenalty,normPenalty,offsetLambdaLogReg);
#else
		Vector lambda = MakeVector(0.,0.,offsetLambdaLogReg);
#endif

		// build matrix and vector elements
#ifdef WEIGHT_BY_NUMBER
		int numa=0,numb=0;
#endif
		for(long i=0;i<linkedParticles;i++)
		{	MPMBase *mpmptr = mpm[list[i]];
			const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
			int matfld = matID->GetField();
			if(matfld==a)
			{	mis[i] = -1;
#ifdef WEIGHT_BY_NUMBER
				numa++;
#endif
			}
			else if(matfld==b || b<0)
			{	mis[i] = 1;
#ifdef WEIGHT_BY_NUMBER
				numb++;
#endif
			}
			else
			{	mis[i] = 0;
				continue;
			}
			xp = mpmptr->pos;
			xp.x -= ndptr->x;
			xp.y -= ndptr->y;
			xp.z = 1.;
		
			// unweighted symmetric matrix
			mt[0] += xp.x*xp.x;
			mt[1] += xp.x*xp.y;
			mt[2] += xp.x;
			mt[3] += xp.y*xp.y;
			mt[4] += xp.y;
			mt[5] += 1.;
			
			// unweighted right hand size
			AddScaledVector(&rhs,&xp,(double)mis[i]);
		}
		
#ifdef WEIGHT_BY_NUMBER
		// weighting for subsequent logistic steps for material a; material b is weighted 1
		double weighta = (double)numb/(double)numa;
#endif
		// solve - normal direction in x and y components
		Matrix3 m(mt[0]+lambda.x, mt[1], 		  mt[2],
								  mt[3]+lambda.y, mt[4],
							    				  mt[5]+lambda.z);
		Matrix3 minv = m.Inverse();
		Vector beta = minv.Times(&rhs);

		// rescale beta by 2 and lambda for improved initial logistic guess
		ScaleVector(&beta,linearBetaScaling);
		ScaleVector(&lambda,linearPenaltyScaling);
		
		// initial norm from linear regression
		Vector prevNorm = beta;
		double nmag = sqrt(beta.x*beta.x+beta.y*beta.y);
		prevNorm.x /= nmag;
		prevNorm.y /= nmag;
		linNorm = prevNorm;			// save for fallback
		linNorm.z = 0.;
		
#ifdef SHOW_CONVERGENCE
		cout << "# Node " << ndptr->num << " particles " << linkedParticles << endl;
		cout << "# LR L=" << m;
		PrintVector(", R=",&rhs);
		cout << endl;
		cout << "# LR " << prevNorm.x << "," << prevNorm.y;
		PrintVector(", ",&beta);
		cout << endl;
#endif

		// logistic analysis nonlinear least squares iteration until done
		while(deltaOneMinusCos>tolerance && iter<maxIter)
		{	// zero matrix and vector elements
			for(long i=0;i<6;i++) mt[i] = 0.;
		
			// penalty term
			rhs.x = -lambda.x*beta.x;
			rhs.y = -lambda.y*beta.y;
			rhs.z = -lambda.z*beta.z;
			
			// build J^TWJ and J^TW(M-f) (for now, w_i = 1)
			for(long i=0;i<linkedParticles;i++)
			{	if(mis[i]==0) continue;
				xp = mpm[list[i]]->pos;
				xp.x -= ndptr->x;
				xp.y -= ndptr->y;
				xp.z = 1.;
				
				double arg = 1.+exp(-DotVectors(&xp,&beta));
#ifdef USE_IRLS
				double phi2 = (arg-1.)/(arg*arg);
#else
				double phi = 2.*(arg-1.)/(arg*arg);
#ifdef WEIGHT_BY_NUMBER
				// weighting for material a, material b is weighted 1
				double phi2 = mis[i]<0. ? weighta*phi*phi : phi*phi ;
#else
				double phi2 = phi*phi;
#endif
#endif
				
				// unweighted symmetric matrix
				mt[0] += phi2*xp.x*xp.x;
				mt[1] += phi2*xp.x*xp.y;
				mt[2] += phi2*xp.x;
				mt[3] += phi2*xp.y*xp.y;
				mt[4] += phi2*xp.y;
				mt[5] += phi2;
				
				// unweighted right hand side
#ifdef USE_IRLS
				double mifi = mis[i]>0. ? 1. - 1./arg : - 1./arg ;
#else
#ifdef WEIGHT_BY_NUMBER
				double mifi = mis[i]<0. ? weighta*((double)mis[i] - 2./arg + 1.)*phi : ((double)mis[i] - 2./arg + 1.)*phi ;
#else
				double mifi = ((double)mis[i] - 2./arg + 1.)*phi;
#endif
#endif
				rhs.x += mifi*xp.x;
				rhs.y += mifi*xp.y;
				rhs.z += mifi;
			}
			
			// Matrix (J^TWJ + Lambda)
			Matrix3 m(mt[0]+lambda.x, mt[1], 		  mt[2],
							 		  mt[3]+lambda.y, mt[4],
								    				  mt[5]+lambda.z);
			
			// next increment added to beta
			Matrix3 minv = m.Inverse();
			Vector dbeta = minv.Times(&rhs);
			AddVector(&beta,&dbeta);
			
			// find normal normal and its angle with previous normal
			norm = beta;
			nmag = sqrt(norm.x*norm.x+norm.y*norm.y);
			norm.x /= nmag;
			norm.y /= nmag;
			deltaOneMinusCos = 1. - norm.x*prevNorm.x - norm.y*prevNorm.y;

#ifdef SHOW_CONVERGENCE
			cout << "# " << iter << ":" << norm.x << "," << norm.y;
			PrintVector(", ",&beta);
			cout << ", " << deltaOneMinusCos << endl;
#endif
			
			// on to next iteration
			iter++;
			prevNorm = norm;
		}
		
		// final result is in norm and is normalized, but don't need offset
		norm.z = 0.;
	}
	
	// if failed, use linear one (checking one component should be enough)
	if(norm.x!=norm.x)
		norm = linNorm;
	
	// find separation from distance to arbirtary plane defined by norm
	double supportForA = 1.e30;
	double supportForB = 1.e30;
	for(long i=0;i<linkedParticles;i++)
	{	if(mis[i]==0) continue;
		xp = mpm[list[i]]->pos;
		if(mis[i]<0)
		{	// find minimum of t (t>0 on material a side of the line) = min(-Xp.n-Rp) = -max(Xp.n+Rp)
			double t = -DotVectors(&norm,&xp)-mpm[list[i]]->GetDeformedRadius(&norm) ;
			if(t<supportForA) supportForA = t;
		}
		else
		{	// find minimum of -t (t<0 on material b side of the line) = min(Xj.n - Rj)
			double t = DotVectors(&norm,&xp)-mpm[list[i]]->GetDeformedRadius(&norm) ;
			if(t<supportForB) supportForB = t;
		}
	}
	
	// compared to contact paper
	//    d((xa-xi).n) = -supportForA-xi.n
	//    d((xb-xi).n) =  supportForB-xi.n
	//    dLR = d((xb-xi).n) - d((xa-xi).n) = supportForB+suppportForA
	// average distance to the node is
	//    dist = 0.5*(|d((xa-xi).n)|+|d((xb-xi).n)|)
	
	// add distances accounting for particle deformation
	delMats->x = supportForA + supportForB;
	hasDeln = true;
	
	// if needed get distances to node
	double xiDotn = norm.x*ndptr->x + norm.y*ndptr->y + norm.z*ndptr->z;
	delMats->y = -supportForA - xiDotn;
	delMats->z = supportForB - xiDotn;

#ifdef SHOW_CONVERGENCE
	cout << "# deln = " << delMats->x << endl;
#endif
	
	// check for lack of convergence, but probably only care if in contact too
	if(iter>=maxIter && delMats->x<0.)
	{	if(warnings.Issue(MeshInfo::warnLRConvergence,-1)==GAVE_WARNING)
		{
#pragma omp critical (output)
            {	cout << "# Exceeded " << maxIter << " iterations on node seeing " << linkedParticles << " particles" << endl;
                // uncomment lines next sectino to print particle positions
				//cout << "# ";
				int numa=0,numb=0;
				for(long i=0;i<linkedParticles;i++)
				{	//if(i!=0 && i%4==0) cout << "\n# ";
					MPMBase *mpmptr = mpm[list[i]];
					const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
					int matfld = matID->GetField();
					//Vector xp = mpmptr->pos;
					if(matfld==a)
					{	//PrintVector("a: ",&xp);
						numa++;
					}
					else if(matfld==b || b<0)
					{	//PrintVector("b: ",&xp);
						numb++;
					}
				}
				//cout << endl;
				cout << "# " << numa << " of type a and " << numb;
				if(b<0)
					cout << " of other materials lumped (may only be one)" << endl;
				else
					cout << " of type b" << endl;
				PrintVector("# Final n=",&norm);
				PrintVector(", linreg n=",&linNorm);
				cout << endl;
				ndptr->Describe(false);
				cout << "#     in crack field " << fieldNum << endl;
			}
		}
	}
	
	// return normalized result
	delete [] mis;
	return norm;
}

// Get separation given normal vector (in norm) and material contact node info
void CrackVelocityFieldMulti::FindSepFromNormalAndPointCloud(Vector *norm,MaterialContactNode *mcn,int vfld,int a,int b,Vector *delMats)
{
	// get linked particles
	vector< int > list = mcn->ParticleLists(vfld);
	long linkedParticles = (long)list.size();

	// find separation
	double supportForA = 1.e30;
	double supportForB = 1.e30;
	for(long i=0;i<linkedParticles;i++)
	{	MPMBase *mpmptr = mpm[list[i]];
		const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
		int matfld = matID->GetField();
		Vector xp = mpmptr->pos;
		if(matfld==a)
		{	// find minimum of t (t>0 on material a side of the line) = min(-Xp.n-Rp) = -max(Xp.n+Rp)
			double t = -DotVectors(norm,&xp)-mpmptr->GetDeformedRadius(norm);
			if(t<supportForA) supportForA = t;
		}
		else if(matfld==b || b<0)
		{	// find minimum of -t (t<0 on material b side of the line) = min(Xj.n - Rj)
			double t = DotVectors(norm,&xp)-mpmptr->GetDeformedRadius(norm) ;
			if(t<supportForB) supportForB = t;
		}
	}
	
	// compared to contact paper
	//    d((xa-xi).n) = -supportForA-xi.n
	//    d((xb-xi).n) =  supportForB-xi.n
	//    dLR = d((xb-xi).n) - d((xa-xi).n) = supportForB+suppportForA
	// average distance to the node is
	//    dist = 0.5*(|d((xa-xi).n)|+|d((xb-xi).n)|)
	
	// separation for particle deformation
	delMats->x = supportForA + supportForB;
	
	// get signed distances to the node
	NodalPoint *ndptr = mcn->GetTheNode();
	double xiDotn = norm->x*ndptr->x - norm->y*ndptr->y - norm->z*ndptr->z;
	delMats->y = -supportForA - xiDotn;
	delMats->z = supportForB - xiDotn;
}
