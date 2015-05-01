/************************************************************************/
/*									
/*									
/*  Created by Nick Tyler					
/*	University Of South Carolina			
/************************************************************************/

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
#include <vector>
#include "TF1.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "physics.hpp"

using namespace std;

//Function to go through data files and calculate W and Q2
//Fill in W vs Q2 hist and save to output root file
//
void dataHandeler(char *fin, char *RootFile_output){
	gROOT->Reset();
	int current_event;
	int num_of_events;
	int total_events = 0;

	int num_elec = 0, num_pip =0;
	int files_in_lis = 2486;


	TFile *myFile;
	TFile *RootOutputFile;
	TTree *myTree;
	int number_cols = 0;
	int number_files = 0;
	char rootFile[500];
	
	RootOutputFile = new TFile(RootFile_output,"RECREATE");

	cout << blue <<"Analyzing file " << green << fin << def << bgdef << endl;

	FILE *input_file = fopen(fin,"r");
	if (input_file == NULL) perror ("Error opening file");

	while (1){

		number_cols = fscanf(input_file,"%s",rootFile);

		if (number_cols<0) break;

		loadbar(number_files,files_in_lis);

		myFile = new TFile(rootFile, "READ");

		myTree = (TTree *)myFile->Get("h10");

		getBranches(myTree);

		num_of_events = (Int_t)myTree->GetEntries();

		current_event = 0;

		while(current_event<num_of_events){

			std::vector<int> Part_ID;
			myTree->GetEntry(current_event);

			#pragma omp parallel for
			for(int event_number = 0; event_number < gpart; event_number++){
				Part_ID.push_back(id[event_number]);
			}
			
			current_event++; 		  	// increment event counter
			total_events++; 			// increment total event counter
		}
		myTree->Delete(); 						// delete Tree object
		myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
		number_files++; 						// increment file counter

	}
	RootOutputFile->cd();
	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
	cout<<blue<<total_events<<" events in "<<number_files<< " files."<<def<<endl; // print out stats
	//cout << "Number of electrons " << num_elec <<  " Number PIP " << num_pip << " ratio:" << (double)num_pip/(double)num_elec << endl;
}
#endif
