/********************************************************************************
    MatRegionFEA.cpp - extra code for FEAReadHandler.cpp to generate elements
        from shape commands mapped to a mesh of elements
    NairnFEA

    Created by John Nairn on 4/20/12.
    Copyright (c) 2012 RSAC Software. All rights reserved.
********************************************************************************/

#include "Read_FEA/FEAReadHandler.hpp"

//-----------------------------------------------------------
// Check for material region command
//-----------------------------------------------------------
short FEAReadHandler::MatRegionInput(char *xName,const Attributes& attrs)
{
    // Begin material shape regions
    if(strcmp(xName,"MatRegion")==0)
    {   ValidateCommand(xName,MUST_BE_NO_BLOCK,ANY_DIM);
        block=MATREGIONBLOCK;
    }
	
    else
        return FALSE;
    
    // was handled
    return TRUE;

}

//-----------------------------------------------------------
// Subroutine to finish up mat region commands
//-----------------------------------------------------------
short FEAReadHandler::EndMatRegionInput(char *xName,int exitBlock)
{
    // finish the matregion block
    if(strcmp(xName,"MatRegion")==0)
    {   block=exitBlock;
    }
    
	// not a MatRegion block element
	else
		return FALSE;
    
	return TRUE;
}
