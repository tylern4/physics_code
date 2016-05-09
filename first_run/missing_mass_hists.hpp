/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef MISSING_MASS_HISTS_H_GUARD
#define MISSING_MASS_HISTS_H_GUARD
#include "main.h"
//
//histogram declarations, fills, and write
//
//
int bins_MM = 200;
double MM_min = 0.0;
double MM_max = 3.0;
TH1D *Missing_Mass = new TH1D("Missing_Mass", "Missing Mass", bins_MM, MM_min, MM_max);
TH1D *Mass = new TH1D("Mass", "Mass", 600, 0, 6);
TH1D *Missing_Mass_square = new TH1D("Missing_Mass_square", "Missing Mass square", bins_MM, MM_min, Square(MM_max));

void Fill_Missing_Mass(double miss_mass){
	Missing_Mass->Fill(miss_mass);
}

void Fill_Mass(double mass){
	Mass->Fill(mass);
}

void Fill_Missing_Mass_square(double miss_mass_2){
	Missing_Mass_square->Fill(miss_mass_2);
}

void Write_Missing_Mass(){
	Missing_Mass->SetXTitle("Mass (GeV)");
	Missing_Mass->Write();

	Mass->SetXTitle("Mass (GeV)");
	Mass->Write();

	Missing_Mass_square->SetXTitle("Mass (GeV)");
	Missing_Mass_square->Write();
}

#endif
