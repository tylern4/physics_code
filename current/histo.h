////////////////////#include <omp.h>
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

<<<<<<< HEAD
TH2F *WvsQ2_hist = new TH2F("WvsQ2_hist","WvsQ2_hist", 1000, 0.0, 3.0, 1000, 2.5, 3.2);
=======
TH2F *WvsQ2_hist = new TH2F("WvsQ2_hist","WvsQ2_hist", 1000, 0.0, 3.0, 1000, 2.5, 4.0);
TH1F *W_hist = new TH1F("W","W",100, 0.0, 3.0);
TH1F *Q2_hist = new TH1F("Q2","Q2",100, 0.0, 3.0);
>>>>>>> FETCH_HEAD

Double_t Px, Py, Pz;
Double_t x,y,z;
Int_t ID;
Double_t W, Q2;



void WvsQ2_Fill(){
	WvsQ2_hist->Fill(Q2,W);
	W_hist->Fill(W);
	Q2_hist->Fill(Q2);
}

void WvsQ2_Write(){
	WvsQ2_hist->Write();
	W_hist->Write();
	Q2_hist->Write();
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