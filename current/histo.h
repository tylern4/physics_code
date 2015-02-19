#ifndef HISTO_H_GUARD
#define HISTO_H_GUARD
////#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

//
//Histogram declarations and write
//
//
TH1F *PxHist = new TH1F("PxHist", "PxHist", 100, -4, 4);
TH1F *PyHist = new TH1F("PyHist", "PyHist", 100, -4, 4);
TH1F *PzHist = new TH1F("PzHist", "PzHist", 100, -4, 4);

TH1F *XHist = new TH1F("XHist", "XHist", 100, -10, 10);
TH1F *YHist = new TH1F("YHist", "YHist", 100, -10, 10);
TH1F *ZHist = new TH1F("ZHist", "ZHist", 100, -10, 10);

TH1I *PartID = new TH1I("PartID", "PartID",423,-211,211);

TH2F *WvsQ2_hist = new TH2F("WvsQ2_hist","WvsQ2_hist", 1000, 0.0, 3.25, 1000, 0, 3.25);
TH2F *Q2vsW_hist = new TH2F("Q2vsW_hist","Q2vsW_hist", 1000, 0.0, 3.25, 1000, 0, 3.25);
TH1F *W_hist = new TH1F("W","W",100, 0.0, 3.25);
TH1F *Q2_hist = new TH1F("Q2","Q2",100, 0.0, 3.35);
TH1F *E_prime_hist = new TH1F("E_prime","E_prime",100,0.0,5.0);

TH2F *MomVsBeta_hist = new TH2F("MomVsBeta","MomVsBeta", 100,0,5.0,100,0.0,2.0);
TH1F *Mom = new TH1F("Momentum","Momentum",100,0,5.0);

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