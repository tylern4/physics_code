/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef WVSQ2_HISTS_H_GUARD
#define WVSQ2_HISTS_H_GUARD
//#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

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

TH1D *E_prime_hist = new TH1D("E_prime","Scattered Electron Energy",bins,0.0,5.0);

TH2D *Q2_vs_xb = new TH2D("Q2_vs_xb","Q^{2} vs x_{b}",bins,0.1,0.6,bins,1.0,3.5);

TH2D *WvsQ2_proton = new TH2D("WvsQ2_proton","W vs Q^{2} p", bins, w_min, w_max, bins, q2_min, q2_max);
TH1D *W_proton = new TH1D("W_proton","W p",bins,  w_min, w_max);
TH1D *Q2_proton = new TH1D("Q2_proton","Q^{2} p",bins, q2_min, q2_max);

TH2D *WvsQ2_e_proton_only_found = new TH2D("WvsQ2_e_proton_only_found","W vs Q^{2} p only", bins, w_min, w_max, bins, q2_min, q2_max);
TH1D *W_e_proton_only_found = new TH1D("W_e_proton_only_found","W p only",bins,  w_min, w_max);
TH1D *Q2_e_proton_only_found = new TH1D("Q2_e_proton_only_found","Q^{2} p only",bins, q2_min, q2_max);

TH2D *WvsQ2_pion = new TH2D("WvsQ2_pion","W vs Q^{2} p #pi^{+} only", bins, w_min, w_max, bins, q2_min, q2_max);
TH1D *W_pion = new TH1D("W_pion","W p #pi^{+} only",bins,  w_min, w_max);
TH1D *Q2_pion = new TH1D("Q2_pion","Q^{2} p #pi^{+} only",bins, q2_min, q2_max);

TH2D *WvsQ2_single_pi = new TH2D("WvsQ2_single_pi","W vs Q^{2} #pi^{+}", bins, w_min, w_max, bins, q2_min, q2_max);
TH1D *W_single_pi = new TH1D("W_single_pi","W #pi^{+}",bins,  w_min, w_max);
TH1D *Q2_single_pi = new TH1D("Q2_single_pi","Q^{2} #pi^{+}",bins, q2_min, q2_max);

void Fill_proton_WQ2(double W, double Q2){
	WvsQ2_proton->Fill(W,Q2);
	W_proton->Fill(W);
	Q2_proton->Fill(Q2);
}

void Fill_single_pi_WQ2(double W, double Q2){
	WvsQ2_single_pi->Fill(W,Q2);
	W_single_pi->Fill(W);
	Q2_single_pi->Fill(Q2);
}

void WvsQ2_Fill(double E_prime, double W, double Q2, double xb){
	E_prime_hist->Fill(E_prime);
	WvsQ2_hist->Fill(W,Q2);
	W_hist->Fill(W);
	Q2_hist->Fill(Q2);
	Q2_vs_xb->Fill(xb,Q2);
}

void Fill_e_proton_only_WQ2(double W, double Q2){
	WvsQ2_e_proton_only_found->Fill(W,Q2);
	W_e_proton_only_found->Fill(W);
	Q2_e_proton_only_found->Fill(Q2);
}

void Fill_pion_WQ2(double W, double Q2){
	WvsQ2_pion->Fill(W,Q2);
	W_pion->Fill(W);
	Q2_pion->Fill(Q2);
}

void WvsQ2_Write(){
	WvsQ2_hist->SetXTitle("W (GeV)");
	WvsQ2_hist->SetYTitle("Q^{2} (GeV^{2})");
	WvsQ2_hist->Write();

	W_hist->SetXTitle("W (GeV)");
	W_hist->Write();

	Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
	Q2_hist->Write();

	E_prime_hist->SetXTitle("Energy (GeV)");
	E_prime_hist->Write();

	Q2_vs_xb->SetXTitle("x_{b}");
	Q2_vs_xb->SetYTitle("Q^{2}");
	Q2_vs_xb->Write();

	WvsQ2_proton->SetXTitle("W (GeV)");
	WvsQ2_proton->SetYTitle("Q^{2} (GeV^{2})");
	WvsQ2_proton->Write();

	W_proton->SetXTitle("W (GeV)");
	W_proton->Write();

	Q2_proton->SetXTitle("Q^{2} (GeV^{2})");
	Q2_proton->Write();

	//WvsQ2_e_proton_only_found->SetXTitle("W (GeV)");
	//WvsQ2_e_proton_only_found->SetYTitle("Q^{2} (GeV^{2})");
	//WvsQ2_e_proton_only_found->Write();

	//W_e_proton_only_found->SetXTitle("W (GeV)");
	//W_e_proton_only_found->Write();

	//Q2_e_proton_only_found->SetXTitle("Q^{2} (GeV^{2})");
	//Q2_e_proton_only_found->Write();

	WvsQ2_pion->SetXTitle("W (GeV)");
	WvsQ2_pion->SetYTitle("Q^{2} (GeV^{2})");
	WvsQ2_pion->Write();

	W_pion->SetXTitle("W (GeV)");
	W_pion->Write();

	Q2_pion->SetXTitle("Q^{2} (GeV^{2})");
	Q2_pion->Write();


	WvsQ2_single_pi->SetXTitle("W (GeV)");
	WvsQ2_single_pi->SetYTitle("Q^{2} (GeV^{2})");
	WvsQ2_single_pi->Write();

	W_single_pi->SetXTitle("W (GeV)");
	W_single_pi->Write();

	Q2_single_pi->SetXTitle("Q^{2} (GeV^{2})");
	Q2_single_pi->Write();
}
#endif
