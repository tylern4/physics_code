/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef HISTO2_H_GUARD
#define HISTO2_H_GUARD
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

//
//ogram declarations, fills, and write
//
//
int bins_2 = 500;
float w_min_2 = 0;
float w_max_2 = 3.25;
float q2_min_2 = 0;
float q2_max_2 = 10;

TH2D *WvsQ2_e_proton_found = new TH2D("WvsQ2_e_proton_found","W vs Q^{2} p", bins_2, w_min_2, w_max_2, bins_2, q2_min_2, q2_max_2);
TH1D *W_e_proton_found = new TH1D("W_e_proton_found","W p",bins_2,  w_min_2, w_max_2);
TH1D *Q2_e_proton_found = new TH1D("Q2_e_proton_found","Q^{2} p",bins_2, q2_min_2, q2_max_2);

TH2D *WvsQ2_e_pi_found = new TH2D("WvsQ2_e_pi_found","W vs Q^{2} #pi^{+}", bins_2, w_min_2, w_max_2, bins_2, q2_min_2, q2_max_2);
TH1D *W_e_pi_found = new TH1D("W_e_pi_found","W #pi^{+}",bins_2,  w_min_2, w_max_2);
TH1D *Q2_e_pi_found = new TH1D("Q2_e_pi_found","Q^{2} #pi^{+}",bins_2, q2_min_2, q2_max_2);

int bins_pvb = 500;
float p_min = 0;
float p_max = 2.5;
float b_min = 0.1;
float b_max = 1.2;

TH2D *MomVsBeta_e_proton_found = new TH2D("MomVsBeta_e_proton_found","Momentum Vs #beta p", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
TH2D *MomVsBeta_e_pi_found = new TH2D("MomVsBeta_e_pi_found","Momentum Vs #beta #pi^{+}", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);

TH1D *Missing_Mass = new TH1D("Missing_Mass", "Missing Mass", 500,-5.0,5.0);

void Fill_Missing_Mass(double mass){
	Missing_Mass->Fill(mass);
}


void Fill_e_proton_found(double W, double Q2,double p,double beta){
	WvsQ2_e_proton_found->Fill(W,Q2);
	W_e_proton_found->Fill(W);
	MomVsBeta_e_proton_found->Fill(p,beta);
	Q2_e_proton_found->Fill(Q2);
}

void Fill_e_pi_found(double W, double Q2,double p,double beta){
	WvsQ2_e_pi_found->Fill(W,Q2);
	W_e_pi_found->Fill(W);
	MomVsBeta_e_pi_found->Fill(p,beta);
	Q2_e_pi_found->Fill(Q2);
}

void Write_found_hists(){

	WvsQ2_e_proton_found->SetXTitle("W (GeV)");
	WvsQ2_e_proton_found->SetYTitle("Q^{2} (GeV^{2})");
	W_e_proton_found->SetXTitle("W (GeV)");
	Q2_e_proton_found->SetXTitle("Q^{2} (GeV^{2})");
	MomVsBeta_e_proton_found->SetXTitle("Momentum (GeV)");
	MomVsBeta_e_proton_found->SetYTitle("#beta");

	MomVsBeta_e_proton_found->Write();

	WvsQ2_e_pi_found->SetXTitle("W (GeV)");
	WvsQ2_e_pi_found->SetYTitle("Q^{2} (GeV^{2})");
	W_e_pi_found->SetXTitle("W (GeV)");
	Q2_e_pi_found->SetXTitle("Q^{2} (GeV^{2})");
	MomVsBeta_e_pi_found->SetXTitle("Momentum (GeV)");
	MomVsBeta_e_pi_found->SetYTitle("#beta");

	MomVsBeta_e_pi_found->Write();

}

#endif
