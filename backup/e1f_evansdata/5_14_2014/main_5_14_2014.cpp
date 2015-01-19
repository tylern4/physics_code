/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler, Canisius College								*/
/*						University Of South Carolina					*/
/************************************************************************/

#define PI 3.14159265;

#include "TTree.h"
#include "TROOT.h"
#include "TH2.h"
#include <TLorentzVector.h>
#include <TFile.h>
#include "TStyle.h"
#include "TCanvas.h"
#include <fstream>
#include "TF1.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TSystem.h"
#include "main.h"
#include "TMath.h"

using namespace std;


//  dataHandeler
//
//	hopefully this works
//

void dataHandeler(char *fin="all.lis", char *RootFile="outFile.root", Int_t MaxEvents=0, Int_t dEvents=10000){
	gROOT->Reset();
	Int_t ii;
	Int_t nEvents;	
	Int_t TotEvents = 0;
	Double_t E = 0;
	Double_t Px, Py, Pz;
//
	TH1F *PxHist = new TH1F("PxHist", "PxHist", 100, 0, 4);
	TH1F *PyHist = new TH1F("PyHist", "PyHist", 100, 0, 4);
	TH1F *PzHist = new TH1F("PzHist", "PzHist", 100, 0, 4);
	TH1F *EHist = new TH1F("EHist", "EHist", 1000, 0, 25000);
	TH1F *EHist2 = new TH1F("EHist2", "EHist2", 1000, 0, 25000);

	TLorentzVector *_e0, *_p0, *_e1;//, *_p1;

	_e0 = new TLorentzVector();
	_p0 = new TLorentzVector();
	_e1 = new TLorentzVector();
	_e0->SetPxPyPzE(0,0,E1F_E0,E1F_E0);
	_p0->SetPxPyPzE(0,0,0,MASS_P);

	TFile *myFile;
	TFile *rootOutFile;
	TTree *myTree;
	Int_t ncols=0;
	Int_t nfiles = 0;
	char rootFile[500];

	rootOutFile = new TFile(RootFile,"RECREATE");
	

	cout << "Analyzing file " << fin << endl;


	FILE *in1 = fopen(fin,"r");
	if (in1 == NULL) perror ("Error opening file");

	while (1){

		ncols = fscanf(in1,"%s",rootFile); 

		if (ncols<0) break;
		myFile = new TFile(rootFile, "READ");
		myTree = (TTree *)myFile->Get("h10clone/h10");

		getBranches(myTree);

		nEvents = (Int_t)myTree->GetEntries();


		ii = 0; 

		while(ii<nEvents){

			myTree->GetEntry(ii);
			for (int j = 0; j < gpart; j++)
			{
				if(id[j] == 11){

					E =((Square(p[j])/(2*MASS_E))+ Square(MASS_E));
					Px = cx[j]*p[j];
					Py = cy[j]*p[j];
					Pz = cz[j]*p[j];
//
				
					PxHist->Fill(Px);
					PyHist->Fill(Py);
					PzHist->Fill(Pz);
					EHist->Fill(E);

					_e1->SetPxPyPzE(Px,Py,Pz,E);

					EHist2->Fill(_e1->E());

				}

			}




			ii++; 		  // increment event counter
			TotEvents++; // increment total event counter 
		}

		myTree->Delete(); 						// delete Tree object
		myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
		nfiles++; 								// increment file counter
		
	}
	rootOutFile->cd();
	PxHist->Write();
	PyHist->Write();
	PzHist->Write();
	EHist->Write();
	EHist2->Write();
	rootOutFile->Write();
	rootOutFile->Close();
	fclose(in1); 							// close file with input file list
	cout<<TotEvents<<" events in "<<nfiles<< " files."<<endl; // print out stats

}









int main(int argc, char **argv){

	float time1 = clock();
	gSystem->Load("libTree");

	char  infilename[128];
	char  outfilename[128];


	sprintf(infilename,"%s",argv[1]);
	sprintf(outfilename,"%s",argv[2]);



	dataHandeler(infilename,outfilename);

	 

 	//time to complete function
 	float time2 = clock();
  	float minutes = 0;
	float seconds = 0;
  	minutes = (time2 - time1)/1000000;
  	minutes = (minutes)/60;
  	seconds = fmod(minutes,1);
	minutes = minutes-seconds;
	seconds = seconds*60;


  	if (minutes==0) cout<<endl<<"Completed in "<<seconds<<" seconds."<<endl<<endl;
  	else cout<<endl<<"Completed in "<<minutes<<" minutes and "<<seconds<<" seconds."<<endl<<endl;


	return 0;
}


