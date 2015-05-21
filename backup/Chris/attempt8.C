//Chris McLauchlin
//This is an attempt to read and plot Root Files!
//Using Nick's root files. 

#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;
int main()	
{

	TFile *myFile;
	TTree *myTree;
	TCanvas *c1 = new TCanvas("c1","c1",0,0,500,500);

	myFile = new TFile("nt23495_02.root","READ");
	myTree = (TTree*)myFile->Get("h10");
	
	
	Int_t nentries= (Int_t) myTree->GetEntries();
	Float_t p[1000];
	myTree->SetBranchAddress("p",&p);


	TH1F *plot = new TH1F("plot","momentum?",1000,0,5);

	for(Int_t i=0;i<nentries;i++)
	{
		myTree->GetEntry(i);
		for (int ii = 0; ii < 3; ii++)
		{
			plot->Fill(p[ii]);	
		}
		
	}
	c1->cd();
	plot->Draw();
	c1->Print("momentum.ps");

	return 0;
}
