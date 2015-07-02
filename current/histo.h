/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef HISTO_H_GUARD
#define HISTO_H_GUARD
//#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "histo2.h"

//
//Histogram declarations, fills, and write
//
//
int bins = 500;
float w_min = 0;
float w_max = 3.25;
float q2_min = 0;
float q2_max = 10;

TH2D *WvsQ2_hist = new TH2D("WvsQ2_hist","W vs Q^{2}", bins, w_min, w_max, bins, q2_min, q2_max);
TH1D *W_hist = new TH1D("W","W",bins,  w_min, w_max);
TH1D *Q2_hist = new TH1D("Q2","Q^{2}",bins, q2_min, q2_max);

TH1D *E_prime_hist = new TH1D("E_prime","Scattered Electron Energy",bins,0.0,2.0);
TH2D *Q2_vs_xb = new TH2D("Q2_vs_xb","Q^{2} vs x_{b}",bins,0.1,0.6,bins,1.0,3.5);

TH2D *MomVsBeta_hist = new TH2D("MomVsBeta","Momentum Vs #beta", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
TH2D *MomVsBeta_hist_pos = new TH2D("MomVsBeta_pos","Momentum Vs #beta Positive", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
TH2D *MomVsBeta_hist_neg = new TH2D("MomVsBeta_neg","Momentum Vs #beta Negative", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
TH1D *Mom = new TH1D("Momentum","Momentum",bins,0,2.0);
TH1D *Energy_hist = new TH1D("Energy_hist","Energy_hist",bins,0.0,2.5);

void WvsQ2_Fill(double E_prime, double W, double Q2, double xb){
	E_prime_hist->Fill(E_prime);
	WvsQ2_hist->Fill(W,Q2);
	W_hist->Fill(W);
	Q2_hist->Fill(Q2);
	Q2_vs_xb->Fill(xb,Q2);
}

void WvsQ2_Write(){
	WvsQ2_hist->SetXTitle("W (GeV)");
	WvsQ2_hist->SetYTitle("Q^{2} (GeV^{2})");

	W_hist->SetXTitle("W (GeV)");
	Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
	E_prime_hist->SetXTitle("Energy (GeV)");

	Q2_vs_xb->SetXTitle("x_{b}");
	Q2_vs_xb->SetYTitle("Q^{2}");

	E_prime_hist->Write();
	WvsQ2_hist->Write();
	W_hist->Write();
	Q2_hist->Write();
	Q2_vs_xb->Write();

}

void MomVsBeta_Fill_pos(double P, double Beta){
	MomVsBeta_hist_pos->Fill(P,Beta);
}

void MomVsBeta_Fill_neg(double P, double Beta){
	MomVsBeta_hist_neg->Fill(P,Beta);
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

	Energy_hist->Write();
	MomVsBeta_hist->Write();
	MomVsBeta_hist_pos->Write();
	MomVsBeta_hist_neg->Write();
	Mom->Write();
}

#endif
