//
//  ElasticMatrixCamClay.cpp
//  NairnMPM
//
//  Created by Raydel Lorenzo on 6/29/12.
//  Copyright (c) 2012 __Geotecnia-UNB__. All rights reserved.
//
// Elastic Matrix for Modified Cam Clay Model. get(line, column)

#include <iostream>

using namespace std;

class matrix
{
public:
	matrix()
	{
	}  
	matrix (double kapa, double pm, double poisson, double e0, double rho)                          
	{
        double Ecc = (3. * (1. - 2. * poisson) * (1. + e0)) * pm / kapa;                    //modulus variable with pm
        
        double aux = Ecc * (1. - poisson) / ((1. + poisson) * (1. - 2. * poisson));
    
        
        matelas[0][0] = matelas[1][1] = matelas[2][2] = aux;
        matelas[0][1] = matelas[0][2] = matelas[1][0] = matelas[1][2] = matelas[2][0] = matelas[2][1] = aux * (poisson / (1. - poisson));
        matelas[0][3] = matelas[0][4] = matelas[0][5] = matelas[3][0] = matelas[3][1] = matelas[3][2] = matelas[4][0] = matelas[4][1] =matelas[4][2] = matelas[4][3] = matelas[5][0] = matelas[5][1] = matelas[5][2] = matelas[5][3] = matelas[5][4] = matelas[1][3] = matelas[1][4] = matelas[1][5] = matelas[2][3]= 0;
        matelas[2][4] = matelas[2][5] = matelas[3][4] = matelas[3][5] = matelas [4][5] = 0;
        matelas[3][3] = matelas[4][4] = matelas[5][5] = aux * (1. - 2. * poisson) / (2. * (1. - poisson)); 
	}
    
	double get(int l, int c)
	{
		return matelas[l][c];
	}
    
private:
	int matelas[6][6];
};