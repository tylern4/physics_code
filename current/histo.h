/************************************************************************/
/*									
/*									
/*  Created by Nick Tyler					
/*	University Of South Carolina			
/************************************************************************/

#ifndef HISTO_H_GUARD
#define HISTO_H_GUARD
//#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

//
//Histogram declarations, fills, and write
//
//
TH1D *PxHist = new TH1D("PxHist", "PxHist", 100, -4, 4);
TH1D *PyHist = new TH1D("PyHist", "PyHist", 100, -4, 4);
TH1D *PzHist = new TH1D("PzHist", "PzHist", 100, -4, 4);

TH1D *XHist = new TH1D("XHist", "XHist", 100, -10, 10);
TH1D *YHist = new TH1D("YHist", "YHist", 100, -10, 10);
TH1D *ZHist = new TH1D("ZHist", "ZHist", 100, -10, 10);

TH1I *PartID = new TH1I("PartID", "PartID",423,-211,211);

TH2D *WvsQ2_hist = new TH2D("WvsQ2_hist","W vs Q^{2}", 1000, 0.0, 3.25, 1000, 0, 3.25);
TH2D *Q2vsW_hist = new TH2D("Q2vsW_hist","Q^{2} vs W", 1000, 0.0, 3.25, 1000, 0, 3.25);
TH1D *W_hist = new TH1D("W","W",100, 0.0, 3.25);
TH1D *Q2_hist = new TH1D("Q2","Q^{2}",100, 0.0, 3.35);
TH1D *E_prime_hist = new TH1D("E_prime","Scattered Electron Energy",100,0.0,6.0);

TH2D *MomVsBeta_hist = new TH2D("MomVsBeta","Momentum Vs #beta", 100, 0, 5.0, 100, 0.0, 1.5);
TH1D *Mom = new TH1D("Momentum","Momentum",100,0,5.0);
TH1D *Energy_hist = new TH1D("Energy_hist","Energy_hist",500,0.0,6.0);

double Px, Py, Pz, P;
double x,y,z;
int ID;
double W, Q2, E_prime; 
double Beta, Energy;

void WvsQ2_Fill(){
	E_prime_hist->Fill(E_prime);
	WvsQ2_hist->Fill(W,Q2);
	Q2vsW_hist->Fill(Q2,W);
	W_hist->Fill(W);
	Q2_hist->Fill(Q2);
}
void WvsQ2_Write(){
	WvsQ2_hist->SetXTitle("W (GeV/c)");
	WvsQ2_hist->SetYTitle("Q^{2} (GeV/c^{2})");

	Q2vsW_hist->SetXTitle("Q^{2} (GeV/c^{2})");
	Q2vsW_hist->SetYTitle("W (GeV/c)");

	W_hist->SetXTitle("W (GeV/c)");
	Q2_hist->SetXTitle("Q^{2} (GeV/c^{2})");
	E_prime_hist->SetXTitle("Energy (GeV)");

	E_prime_hist->Write();
	WvsQ2_hist->Write();
	Q2vsW_hist->Write();
	W_hist->Write();
	Q2_hist->Write();
}

void MomVsBeta_Fill(){
	Energy_hist->Fill(Energy);
	MomVsBeta_hist->Fill(P,Beta);
	Mom->Fill(P);
}
void MomVsBeta_Write(){
	MomVsBeta_hist->SetXTitle("Momentum (GeV/c)");
	MomVsBeta_hist->SetYTitle("#beta");
	Mom->SetXTitle("Momentum (GeV/c)");

	Energy_hist->Write();
	//MomVsBeta_hist->Write();
	Mom->Write();
}

void WriteHists(){
	PxHist->Write();
	PyHist->Write();
	PzHist->Write();

	XHist->Write();
	YHist->Write();
	ZHist->Write();

	PartID->Write();
}

void FillHist(){
	PxHist->Fill(Px);
	PyHist->Fill(Py);
	PzHist->Fill(Pz);

	XHist->Fill(x);
	YHist->Fill(y);
	ZHist->Fill(z);

	PartID->Fill(ID);
}

#endif
