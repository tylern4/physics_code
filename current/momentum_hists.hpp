/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef MOMENTUM_HISTS_H_GUARD
#define MOMENTUM_HISTS_H_GUARD
#include "main.h"
//
//Histogram declarations, fills, and write
//
//
int bins_pvb = 500;
float p_min = 0;
float p_max = 2.5;
float b_min = 0.1;
float b_max = 1.2;

TH2D *MomVsBeta_hist = new TH2D("MomVsBeta","Momentum Vs #beta", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
TH2D *MomVsBeta_hist_pos = new TH2D("MomVsBeta_pos","Momentum Vs #beta Positive", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
TH2D *MomVsBeta_hist_neg = new TH2D("MomVsBeta_neg","Momentum Vs #beta Negative", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
TH1D *Mom = new TH1D("Momentum","Momentum",bins,0,2.0);
TH1D *Energy_hist = new TH1D("Energy_hist","Energy_hist",bins,0.0,2.5);

TH2D *MomVsBeta_e_proton_found = new TH2D("MomVsBeta_e_proton_found","Momentum Vs #beta p", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
TH2D *MomVsBeta_e_pi_found = new TH2D("MomVsBeta_e_pi_found","Momentum Vs #beta #pi^{+}", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);

void MomVsBeta_Fill_pos(double P, double Beta){
	MomVsBeta_hist_pos->Fill(P,Beta);
}

void MomVsBeta_Fill_neg(double P, double Beta){
	MomVsBeta_hist_neg->Fill(P,Beta);
}

void Fill_e_proton_found_P(double p,double beta){
	MomVsBeta_e_proton_found->Fill(p,beta);
}

void Fill_e_pi_found_P(double p,double beta){
	MomVsBeta_e_pi_found->Fill(p,beta);
}
void MomVsBeta_Fill(double Energy, double P, double Beta){
	Energy_hist->Fill(Energy);
	MomVsBeta_hist->Fill(P,Beta);
	Mom->Fill(P);
}
void MomVsBeta_Write(){
	MomVsBeta_hist->SetXTitle("Momentum (GeV)");
	MomVsBeta_hist->SetYTitle("#beta");
	MomVsBeta_hist_pos->SetXTitle("Momentum (GeV)");
	MomVsBeta_hist_pos->SetYTitle("#beta");
	MomVsBeta_hist_neg->SetXTitle("Momentum (GeV)");
	MomVsBeta_hist_neg->SetYTitle("#beta");
	Mom->SetXTitle("Momentum (GeV)");

	MomVsBeta_e_proton_found->SetXTitle("Momentum (GeV)");
	MomVsBeta_e_proton_found->SetYTitle("#beta");
	MomVsBeta_e_proton_found->Write();

	MomVsBeta_e_pi_found->SetXTitle("Momentum (GeV)");
	MomVsBeta_e_pi_found->SetYTitle("#beta");
	MomVsBeta_e_pi_found->Write();

	Energy_hist->Write();
	MomVsBeta_hist->Write();
	MomVsBeta_hist_pos->Write();
	MomVsBeta_hist_neg->Write();
	Mom->Write();
}

#endif
