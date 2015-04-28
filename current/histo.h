/************************************************************************/
/*									
/*									
/*  Created by Nick Tyler					
/*	University Of South Carolina			
/************************************************************************/

#ifndef HISTO_H_GUARD
#define HISTO_H_GUARD
#include <omp.h>
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

TH1I *PartID = new TH1I("PartID", "PartID",10,0,10);

TH2D *WvsQ2_hist = new TH2D("WvsQ2_hist","W vs Q^{2}", 500, 0.0, 3.25, 500, 0, 3.25);
TH2D *Q2vsW_hist = new TH2D("Q2vsW_hist","Q^{2} vs W", 500, 0.0, 3.25, 500, 0, 3.25);
TH1D *W_hist = new TH1D("W","W",100, 0.0, 3.25);
TH1D *Q2_hist = new TH1D("Q2","Q^{2}",100, 0.0, 3.35);
TH1D *E_prime_hist = new TH1D("E_prime","Scattered Electron Energy",100,0.0,6.0);
TH2D *Q2_vs_xb = new TH2D("Q2_vs_xb","Q^{2} vs x_{b}",500,0.1,0.6,500,1.0,3.5);

TH2D *MomVsBeta_hist = new TH2D("MomVsBeta","Momentum Vs #beta", 100, 0, 5.0, 100, 0.0, 1.5);
TH1D *Mom = new TH1D("Momentum","Momentum",100,0,5.0);
TH1D *Energy_hist = new TH1D("Energy_hist","Energy_hist",500,0.0,6.0);

double Px, Py, Pz, P;
//double x,y,z;
int ID;
double W, Q2, E_prime, xb; 
double Beta, Energy;

void WvsQ2_Fill(){
	E_prime_hist->Fill(E_prime);
	WvsQ2_hist->Fill(W,Q2);
	Q2vsW_hist->Fill(Q2,W);
	W_hist->Fill(W);
	Q2_hist->Fill(Q2);
	Q2_vs_xb->Fill(xb,Q2);
}
void WvsQ2_Write(){
	WvsQ2_hist->SetXTitle("W (GeV)");
	WvsQ2_hist->SetYTitle("Q^{2} (GeV^{2})");

	Q2vsW_hist->SetXTitle("Q^{2} (GeV^{2})");
	Q2vsW_hist->SetYTitle("W (GeV)");

	W_hist->SetXTitle("W (GeV)");
	Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
	E_prime_hist->SetXTitle("Energy (GeV)");

	Q2_vs_xb->SetXTitle("x_{b}");
	Q2_vs_xb->SetYTitle("Q^{2}");

	E_prime_hist->Write();
	//WvsQ2_hist->Write();
	Q2vsW_hist->Write();
	W_hist->Write();
	Q2_hist->Write();
	Q2_vs_xb->Write();

	PartID->Write();
}

void MomVsBeta_Fill(){
	Energy_hist->Fill(Energy);
	MomVsBeta_hist->Fill(P,Beta);
	Mom->Fill(P);
}
void MomVsBeta_Write(){
	MomVsBeta_hist->SetXTitle("Momentum (GeV)");
	MomVsBeta_hist->SetYTitle("#beta");
	Mom->SetXTitle("Momentum (GeV)");

	Energy_hist->Write();
	MomVsBeta_hist->Write();
	Mom->Write();
}

void WriteHists(){
	//PxHist->Write();
	//PyHist->Write();
	//PzHist->Write();

	//XHist->Write();
	//YHist->Write();
	//ZHist->Write();

	PartID->GetXaxis()->SetBinLabel(2,"proton");
	PartID->GetXaxis()->SetBinLabel(3,"nuetron");
	PartID->GetXaxis()->SetBinLabel(4,"#pi^{+}");
	PartID->GetXaxis()->SetBinLabel(5,"#pi^{-}");
	PartID->GetXaxis()->SetBinLabel(6,"#pi^{0}");
	PartID->GetXaxis()->SetBinLabel(7,"#gamma");
	PartID->GetXaxis()->SetBinLabel(1,"Electron");
	PartID->SetFillColor(kBlue);
	PartID->Write();
}

void FillHist(){
	//PxHist->Fill(Px);
	//PyHist->Fill(Py);
	//PzHist->Fill(Pz);

	//XHist->Fill(x);
	//YHist->Fill(y);
	//ZHist->Fill(z);
	switch (ID) {
		case 2212:
			PartID->Fill(1);
			break;
		case 2112:
			PartID->Fill(2);
			break;
		case 211:
			PartID->Fill(3);
			break;
		case -211:
			PartID->Fill(4);
			break;
		case 111:
			PartID->Fill(5);
			break;
		case 22:
			PartID->Fill(6);
			break;
		case 11:
			PartID->Fill(0);
			break;
		default:
			PartID->Fill(10);
			break;
	}

	
}

//Overload of above so that you can input which ID to fill
void FillHist(int Particle_ID){
	//PxHist->Fill(Px);
	//PyHist->Fill(Py);
	//PzHist->Fill(Pz);

	//XHist->Fill(x);
	//YHist->Fill(y);
	//ZHist->Fill(z);
	switch (Particle_ID) {
		case 2212:
			PartID->Fill(1);
			break;
		case 2112:
			PartID->Fill(2);
			break;
		case 211:
			PartID->Fill(3);
			break;
		case -211:
			PartID->Fill(4);
			break;
		case 111:
			PartID->Fill(5);
			break;
		case 22:
			PartID->Fill(6);
			break;
		case 11:
			PartID->Fill(0);
			break;
		default:
			PartID->Fill(10);
			break;
	}

	
}

#endif
