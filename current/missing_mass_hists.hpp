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
int bins_MM = 1000;
double MM_min = 0.0;
double MM_max = 3.0;
TH1D *Missing_Mass = new TH1D("Missing_Mass", "Missing Mass", bins_MM, MM_min, MM_max);
TH1D *Missing_Mass_e_proton_pi_only_found = new TH1D("Missing_Mass_e_proton_pi_only_found", "Missing_Mass_e_proton_pi_only_found", bins_MM, MM_min, MM_max);

void Fill_Missing_Mass(double miss_mass){
	Missing_Mass->Fill(miss_mass);
}

void Fill_Missing_Mass_P_PI(double miss_mass){
	Missing_Mass_e_proton_pi_only_found->Fill(miss_mass);
}


void Write_Missing_Mass(){
	Missing_Mass->SetXTitle("Mass (GeV)");
	Missing_Mass->Write();

	Missing_Mass_e_proton_pi_only_found->SetXTitle("Mass (GeV)");
	Missing_Mass_e_proton_pi_only_found->Write();
}

#endif
