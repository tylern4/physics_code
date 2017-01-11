/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef HISTOGRAM_H_GUARD
#define HISTOGRAM_H_GUARD
#include "fid_hists.hpp"
#include "missing_mass_hists.hpp"
#include "delta_t_hist.hpp" 
//#include "WvsQ2_hists.hpp"
//#include "momentum_hists.hpp"
#include "cc_hist.hpp"
#include "ec_hist.hpp"

class Histogram {
	public:
		Histogram ();
		~Histogram();
		// W and Q^2
		void Fill_proton_WQ2(double W, double Q2);
		void Fill_single_pi_WQ2(double W, double Q2);
		void Fill_single_proton_WQ2(double W, double Q2);
		void WvsQ2_Fill(double E_prime, double W, double Q2, double xb);
		void Fill_pion_WQ2(double W, double Q2);
		void WvsQ2_Write();

		// P and E
		void MomVsBeta_Fill_pos(double P, double Beta);
		void MomVsBeta_Fill_neg(double P, double Beta);
		void Fill_proton_ID_P(double p, double beta);
		void Fill_Pi_ID_P(double p,double beta);
		void Fill_proton_Pi_ID_P(double p,double beta);
		void MomVsBeta_Fill(double Energy, double P, double Beta);
		void MomVsBeta_Write();

	private:
		// W and Q^2
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
		
		TH2D *WvsQ2_proton = new TH2D("WvsQ2_proton","W vs Q^{2} P", bins, w_min, w_max, bins, q2_min, q2_max);
		TH1D *W_proton = new TH1D("W_proton","W P",bins,  w_min, w_max);
		TH1D *Q2_proton = new TH1D("Q2_proton","Q^{2} P",bins, q2_min, q2_max);
		
		TH2D *WvsQ2_pion = new TH2D("WvsQ2_pion","W vs Q^{2} #pi^{+} only", bins, w_min, w_max, bins, q2_min, q2_max);
		TH1D *W_pion = new TH1D("W_pion","W #pi^{+} only",bins,  w_min, w_max);
		TH1D *Q2_pion = new TH1D("Q2_pion","Q^{2} #pi^{+} only",bins, q2_min, q2_max);
		
		TH2D *WvsQ2_single_pi = new TH2D("WvsQ2_single_pi","W vs Q^{2} #pi^{+}", bins, w_min, w_max, bins, q2_min, q2_max);
		TH1D *W_single_pi = new TH1D("W_single_pi","W #pi^{+}",bins,  w_min, w_max);
		TH1D *Q2_single_pi = new TH1D("Q2_single_pi","Q^{2} #pi^{+}",bins, q2_min, q2_max);
		
		TH2D *WvsQ2_single_proton = new TH2D("WvsQ2_single_proton","W vs Q^{2} P", bins, w_min, w_max, bins, q2_min, q2_max);
		TH1D *W_single_proton = new TH1D("W_single_proton","W P",bins,  w_min, w_max);
		TH1D *Q2_single_proton = new TH1D("Q2_single_proton","Q^{2} P",bins, q2_min, q2_max);
		// W and Q^2

		// P and E
		int bins_pvb = 500;
		float p_min = 0;
		float p_max = 2.5;
		float b_min = 0.1;
		float b_max = 1.2;
		
		TH2D *MomVsBeta_hist = new TH2D("MomVsBeta","Momentum versus #beta", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
		TH2D *MomVsBeta_hist_pos = new TH2D("MomVsBeta_pos","Momentum versus #beta Positive", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
		TH2D *MomVsBeta_hist_neg = new TH2D("MomVsBeta_neg","Momentum versus #beta Negative", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
		TH1D *Mom = new TH1D("Momentum","Momentum",bins,0,2.5);
		TH1D *Energy_hist = new TH1D("Energy_hist","Energy_hist",bins,0.0,2.5);
		
		TH2D *MomVsBeta_proton_ID = new TH2D("MomVsBeta_proton_ID","Momentum versus #beta p", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
		TH2D *MomVsBeta_Pi_ID = new TH2D("MomVsBeta_Pi_ID","Momentum versus #beta #pi^{+}", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
		
		TH2D *MomVsBeta_proton_Pi_ID = new TH2D("MomVsBeta_proton_Pi_ID","Momentum versus #beta P #pi^{+}", bins_pvb, p_min, p_max, bins_pvb, b_min, b_max);
		// P and E




};

Histogram::Histogram () {
	makeHists_delta_t();
	makeHists_CC();
	makeHists_fid();
}

Histogram::~Histogram () {
}

// W and Q^2
void Histogram::Fill_proton_WQ2(double W, double Q2){
	WvsQ2_proton->Fill(W,Q2);
	W_proton->Fill(W);
	Q2_proton->Fill(Q2);
}

void Histogram::Fill_single_pi_WQ2(double W, double Q2){
	WvsQ2_single_pi->Fill(W,Q2);
	W_single_pi->Fill(W);
	Q2_single_pi->Fill(Q2);
}

void Histogram::Fill_single_proton_WQ2(double W, double Q2){
	WvsQ2_single_proton->Fill(W,Q2);
	W_single_proton->Fill(W);
	Q2_single_proton->Fill(Q2);
}

void Histogram::WvsQ2_Fill(double E_prime, double W, double Q2, double xb){
	E_prime_hist->Fill(E_prime);
	WvsQ2_hist->Fill(W,Q2);
	W_hist->Fill(W);
	Q2_hist->Fill(Q2);
	Q2_vs_xb->Fill(xb,Q2);
}

void Histogram::Fill_pion_WQ2(double W, double Q2){
	WvsQ2_pion->Fill(W,Q2);
	W_pion->Fill(W);
	Q2_pion->Fill(Q2);
}

void Histogram::WvsQ2_Write(){
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

	WvsQ2_single_proton->SetXTitle("W (GeV)");
	WvsQ2_single_proton->SetYTitle("Q^{2} (GeV^{2})");
	WvsQ2_single_proton->Write();

	W_single_proton->SetXTitle("W (GeV)");
	W_single_proton->Write();

	Q2_single_proton->SetXTitle("Q^{2} (GeV^{2})");
	Q2_single_proton->Write();
}
// W and Q^2

//P and E
void Histogram::MomVsBeta_Fill_pos(double P, double Beta){
	MomVsBeta_hist_pos->Fill(P,Beta);
}

void Histogram::MomVsBeta_Fill_neg(double P, double Beta){
	MomVsBeta_hist_neg->Fill(P,Beta);
}

void Histogram::Fill_proton_ID_P(double p, double beta){
	MomVsBeta_proton_ID->Fill(p,beta);
}

void Histogram::Fill_Pi_ID_P(double p,double beta){
	MomVsBeta_Pi_ID->Fill(p,beta);
}

void Histogram::Fill_proton_Pi_ID_P(double p,double beta){
	MomVsBeta_proton_Pi_ID->Fill(p,beta);
}
void Histogram::MomVsBeta_Fill(double Energy, double P, double Beta){
	Energy_hist->Fill(Energy);
	MomVsBeta_hist->Fill(P,Beta);
	Mom->Fill(P);
}

void Histogram::MomVsBeta_Write(){
	MomVsBeta_hist->SetXTitle("Momentum (GeV)");
	MomVsBeta_hist->SetYTitle("#beta");
	MomVsBeta_hist_pos->SetXTitle("Momentum (GeV)");
	MomVsBeta_hist_pos->SetYTitle("#beta");
	MomVsBeta_hist_neg->SetXTitle("Momentum (GeV)");
	MomVsBeta_hist_neg->SetYTitle("#beta");
	Mom->SetXTitle("Momentum (GeV)");

	MomVsBeta_proton_ID->SetXTitle("Momentum (GeV)");
	MomVsBeta_proton_ID->SetYTitle("#beta");
	MomVsBeta_proton_ID->Write();

	MomVsBeta_Pi_ID->SetXTitle("Momentum (GeV)");
	MomVsBeta_Pi_ID->SetYTitle("#beta");
	MomVsBeta_Pi_ID->Write();

	MomVsBeta_proton_Pi_ID->SetXTitle("Momentum (GeV)");
	MomVsBeta_proton_Pi_ID->SetYTitle("#beta");
	MomVsBeta_proton_Pi_ID->Write();

	Energy_hist->Write();
	MomVsBeta_hist->Write();
	MomVsBeta_hist_pos->Write();
	MomVsBeta_hist_neg->Write();
	Mom->Write();
}







#endif
