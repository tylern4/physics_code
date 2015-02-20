#ifndef HISTO_H_GUARD
#define HISTO_H_GUARD
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

//
//Histogram declarations and write
//
//
TH1D *PxHist = new TH1D("PxHist", "PxHist", 100, -4, 4);
TH1D *PyHist = new TH1D("PyHist", "PyHist", 100, -4, 4);
TH1D *PzHist = new TH1D("PzHist", "PzHist", 100, -4, 4);

TH1D *XHist = new TH1D("XHist", "XHist", 100, -10, 10);
TH1D *YHist = new TH1D("YHist", "YHist", 100, -10, 10);
TH1D *ZHist = new TH1D("ZHist", "ZHist", 100, -10, 10);

TH1I *PartID = new TH1I("PartID", "PartID",423,-211,211);

TH2D *WvsQ2_hist = new TH2D("WvsQ2_hist","WvsQ2_hist", 1000, 0.0, 3.25, 1000, 0, 3.25);
TH2D *Q2vsW_hist = new TH2D("Q2vsW_hist","Q2vsW_hist", 1000, 0.0, 3.25, 1000, 0, 3.25);
TH1D *W_hist = new TH1D("W","W",100, 0.0, 3.25);
TH1D *Q2_hist = new TH1D("Q2","Q2",100, 0.0, 3.35);
TH1D *E_prime_hist = new TH1D("E_prime","E_prime",100,0.0,5.0);

TH2D *MomVsBeta_hist = new TH2D("MomVsBeta","MomVsBeta", 100,0,5.0,100,0.5,1.5);
TH1D *Mom = new TH1D("Momentum","Momentum",100,0,5.0);

double Px, Py, Pz, P;
double x,y,z;
int ID;
double W, Q2, E_prime; 
double Beta;

void WvsQ2_Fill(){
	E_prime_hist->Fill(E_prime);
	WvsQ2_hist->Fill(W,Q2);
	Q2vsW_hist->Fill(Q2,W);
	W_hist->Fill(W);
	Q2_hist->Fill(Q2);
}
void WvsQ2_Write(){
	E_prime_hist->Write();
	WvsQ2_hist->Write();
	Q2vsW_hist->Write();
	W_hist->Write();
	Q2_hist->Write();
}

void MomVsBeta_Fill(){
	MomVsBeta_hist->Fill(P,Beta);
	Mom->Fill(P);
}
void MomVsBeta_Write(){
	MomVsBeta_hist->Write();
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
