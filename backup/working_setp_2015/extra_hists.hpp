/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef MISSING_MASS_HISTS_H_GUARD
#define MISSING_MASS_HISTS_H_GUARD
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

//
//ogram declarations, fills, and write
//
//

TH1D *Missing_Mass = new TH1D("Missing_Mass", "Missing Mass", 500,0.0,5.0);

void Fill_Missing_Mass(double mass){
	Missing_Mass->Fill(mass);
}


void Write_found_hists(){
	Missing_Mass->Write();
}

#endif
