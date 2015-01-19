#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

//
//Histogram declarations and write
//
//
TH1F *PxHist = new TH1F("PxHist", "PxHist", 100, 0, 4);
TH1F *PyHist = new TH1F("PyHist", "PyHist", 100, 0, 4);
TH1F *PzHist = new TH1F("PzHist", "PzHist", 100, 0, 4);

TH1F *XHist = new TH1F("XHist", "XHist", 100, -10, 10);
TH1F *YHist = new TH1F("YHist", "YHist", 100, -10, 10);
TH1F *ZHist = new TH1F("ZHist", "ZHist", 100, -10, 10);

TH1F *EHist = new TH1F("EHist", "EHist", 1000, 0, 25000);

TH2D *EHistX = new TH2D("EHistX", "EHistX", 1000, 0, 25000, 100, 0, 4);
TH2D *EHistY = new TH2D("EHistY", "EHistY", 1000, 0, 25000, 100, 0, 4);
TH2D *EHistZ = new TH2D("EHistZ", "EHistZ", 1000, 0, 25000, 100, 0, 5);

TH1F *Q2Hist = new TH1F("Q2Hist","Q2Hist",100,0,50000);


Double_t E = 0;
Double_t Px, Py, Pz;
Double_t x,y,z,Q2;


void WriteHists(){

	PxHist->Write();
	PyHist->Write();
	PzHist->Write();
	XHist->Write();
	YHist->Write();
	ZHist->Write();
	EHist->Write();
	EHistX->Write();
	EHistY->Write();
	EHistZ->Write();

	Q2Hist->Write();

}

void FillHist(){
	PxHist->Fill(Px);
	PyHist->Fill(Py);
	PzHist->Fill(Pz);
	EHist->Fill(E);

	XHist->Fill(x);
	YHist->Fill(y);
	ZHist->Fill(z);


	EHistX->Fill(E,x);
	EHistY->Fill(E,y);
	EHistZ->Fill(E,z);

	Q2Hist->Fill(Q2);

}