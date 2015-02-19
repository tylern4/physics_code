#ifndef DATAHANDELER_H_GUARD
#define DATAHANDELER_H_GUARD
#include "TTree.h"
#include <math.h>
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
#include "TMath.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


//  dataHandeler
//
//	hopefully this works
//
void dataHandeler(char *fin="all.lis", char *RootFile_output="outFile.root", Int_t MaxEvents=0, Int_t dEvents=10000){
	gROOT->Reset();
	Int_t current_event;
	Int_t num_of_events;
	Int_t total_events = 0;

	TLorentzVector *_e0, *_p0, *_e1;//, *_p1;

	_e0 = new TLorentzVector();
	_p0 = new TLorentzVector();
	_e1 = new TLorentzVector();
	_e0->SetPxPyPzE(0,0,E1D_E0,E1D_E0);
	_p0->SetPxPyPzE(0,0,0,MASS_P);

	TFile *myFile;
	TFile *RootOutputFile;
	TTree *myTree;
	Int_t number_cols=0;
	Int_t number_files = 0;
	char rootFile[500];

	RootOutputFile = new TFile(RootFile_output,"RECREATE");

	cout << "Analyzing file " << fin << endl;


	FILE *input_file = fopen(fin,"r");
	if (input_file == NULL) perror ("Error opening file");

	while (1){

		number_cols = fscanf(input_file,"%s",rootFile); 

		if (number_cols<0) break;
		myFile = new TFile(rootFile, "READ");

		myTree = (TTree *)myFile->Get("h10");


		getBranches(myTree);

		num_of_events = (Int_t)myTree->GetEntries();


		current_event = 0; 

		while(current_event<num_of_events){

			myTree->GetEntry(current_event);

			////////////if (current_event%10000 == 0)	cout<<current_event<<"/"<<num_of_events<<endl;

			#pragma omp parallel for
			for(int event_number = 0; event_number < gpart; event_number++)
			{

				Px = cx[event_number]*p[event_number];
				Py = cy[event_number]*p[event_number];
				Pz = cz[event_number]*p[event_number];

				x = vx[event_number];
				y = vy[event_number];
				z = vz[event_number];

				ID = id[event_number];


				FillHist();

			}

			current_event++; 		  	// increment event counter
			total_events++; 					// increment total event counter 
		}

		myTree->Delete(); 						// delete Tree object
		myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
		number_files++; 						// increment file counter

	}
	RootOutputFile->cd();
	WriteHists();
	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
	cout<<total_events<<" events in "<<number_files<< " files."<<endl; // print out stats
}

//Function to go through data files and calculate W and Q2
//Fill in W vs Q2 hist and save to output root file
//
void WvsQ2(char *fin, char *RootFile_output){
	gROOT->Reset();
	Int_t current_event;
	Int_t num_of_events;
	Int_t total_events = 0;

	TFile *myFile;
	TFile *RootOutputFile;
	TTree *myTree;
	Int_t number_cols=0;
	Int_t number_files = 0;
	char rootFile[500];

	RootOutputFile = new TFile(RootFile_output,"RECREATE");

	cout << "Analyzing file " << fin << endl;


	FILE *input_file = fopen(fin,"r");
	if (input_file == NULL) perror ("Error opening file");

	while (1){

		number_cols = fscanf(input_file,"%s",rootFile); 

		if (number_cols<0) break;
		myFile = new TFile(rootFile, "READ");

		myTree = (TTree *)myFile->Get("h10");


		getBranches(myTree);

		num_of_events = (Int_t)myTree->GetEntries();


		current_event = 0; 

		while(current_event<num_of_events){

			myTree->GetEntry(current_event);

			///////////if (current_event%10000 == 0)	cout<<current_event<<"/"<<num_of_events<<endl;

			#pragma omp parallel for
			for(int event_number = 0; event_number < gpart; event_number++){
				P = P_calc(p[event_number],cx[event_number],cy[event_number],cz[event_number]);
				Beta = b[event_number];
				//If event is a scattered electron find it's energy
				//	and from it's energy find Q2 and W and then fill
				//	histograms for those variables
				if (id[event_number] == ELECTRON){
					E_prime = E_calc(p[event_number],cx[event_number],cy[event_number],cz[event_number]);
					Q2 = Q2_calc(cz[event_number],E_prime);
					W = W_calc(E_prime);

					WvsQ2_Fill();
				}
				MomVsBeta_Fill();
			}

			current_event++; 		  	// increment event counter
			total_events++; 					// increment total event counter 
		}

		myTree->Delete(); 						// delete Tree object
		myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
		number_files++; 						// increment file counter

	}
	RootOutputFile->cd();
	WvsQ2_Write();
	MomVsBeta_Write();
	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
	cout<<total_events<<" events in "<<number_files<< " files."<<endl; // print out stats
}
#endif
