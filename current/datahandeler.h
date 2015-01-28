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
#include "TMath.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


//  dataHandeler
//
//	hopefully this works
//

void dataHandeler(char *fin="all.lis", char *RootFile="outFile.root", Int_t MaxEvents=0, Int_t dEvents=10000){
	gROOT->Reset();
	Int_t current_event_number;
	Int_t num_of_events;
	Int_t total_events = 0;

	TLorentzVector *_e0, *_p0, *_e1;//, *_p1;

	_e0 = new TLorentzVector();
	_p0 = new TLorentzVector();
	_e1 = new TLorentzVector();
	_e0->SetPxPyPzE(0,0,E1F_E0,E1F_E0);
	_p0->SetPxPyPzE(0,0,0,MASS_P);

	TFile *myFile;
	TFile *rootOutFile;
	TTree *myTree;
	Int_t number_cols=0;
	Int_t number_files = 0;
	char rootFile[500];

	rootOutFile = new TFile(RootFile,"RECREATE");

	cout << "Analyzing file " << fin << endl;


	FILE *in1 = fopen(fin,"r");
	if (in1 == NULL) perror ("Error opening file");

	while (1){

		number_cols = fscanf(in1,"%s",rootFile); 

		if (number_cols<0) break;
		myFile = new TFile(rootFile, "READ");

		myTree = (TTree *)myFile->Get("h10");


		getBranches(myTree);

		num_of_events = (Int_t)myTree->GetEntries();


		current_event_number = 0; 

		while(current_event_number<num_of_events){

			myTree->GetEntry(current_event_number);

			////////////if (current_event_number%10000 == 0)	cout<<current_event_number<<"/"<<num_of_events<<endl;

			#pragma omp parallel for
			for(int j = 0; j < gpart; j++)
			{

				Px = cx[j]*p[j];
				Py = cy[j]*p[j];
				Pz = cz[j]*p[j];

				x = vx[j];
				y = vy[j];
				z = vz[j];

				ID = id[j];


				FillHist();

			}

			current_event_number++; 		  	// increment event counter
			total_events++; 					// increment total event counter 
		}

		myTree->Delete(); 						// delete Tree object
		myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
		number_files++; 						// increment file counter

	}
	rootOutFile->cd();
	WriteHists();
	rootOutFile->Write();
	rootOutFile->Close();
	fclose(in1); 														// close file with input file list
	cout<<total_events<<" events in "<<number_files<< " files."<<endl; // print out stats
}

void count_after_cut(char *fin="all.lis", char *RootFile="outFile.root"){
	gROOT->Reset();
	//Int_t current_event_number;
	//Int_t num_of_events;
	//Int_t total_events = 0;

	//TLorentzVector *_e0, *_p0, *_e1;//, *_p1;

	//_e0 = new TLorentzVector();
	//_p0 = new TLorentzVector();
	//_e1 = new TLorentzVector();
	//_e0->SetPxPyPzE(0,0,E1F_E0,E1F_E0);
	//_p0->SetPxPyPzE(0,0,0,MASS_P);

	//TFile *myFile;
	//TFile *rootOutFile;
	//TTree *myTree;
	//Int_t number_cols=0;
	//Int_t number_files = 0;
	//char rootFile[500];

	//rootOutFile = new TFile(RootFile,"RECREATE");

	//cout << "Analyzing file " << fin << endl;


	//FILE *in1 = fopen(fin,"r");
	//if (in1 == NULL) perror ("Error opening file");

	//while (1){

		//number_cols = fscanf(in1,"%s",rootFile); 

		//if (number_cols<0) break;
		//myFile = new TFile(rootFile, "READ");

		//myTree = (TTree *)myFile->Get("h10");


		//getBranches(myTree);

		//num_of_events = (Int_t)myTree->GetEntries();


		//current_event_number = 0; 

		//while(current_event_number<num_of_events){

			//myTree->GetEntry(current_event_number);

			////////////if (current_event_number%10000 == 0)	cout<<current_event_number<<"/"<<num_of_events<<endl;

			//#pragma omp parallel for
			//for(int j = 0; j < gpart; j++){

				//Px = cx[j]*p[j];
				//Py = cy[j]*p[j];
				//Pz = cz[j]*p[j];

				//x = vx[j];
				//y = vy[j];
				//z = vz[j];

				//ID = id[j];


				//FillHist();

			//}

			//current_event_number++; 		  	// increment event counter
			//total_events++; 					// increment total event counter 
		//}

		//myTree->Delete(); 						// delete Tree object
		//myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
		//number_files++; 						// increment file counter

	//}
	//rootOutFile->cd();
	//WriteHists();
	//rootOutFile->Write();
	//rootOutFile->Close();
	//fclose(in1); 														// close file with input file list
	//cout<<total_events<<" events in "<<number_files<< " files."<<endl; // print out stats
}