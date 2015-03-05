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
#include "TF1.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
//#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "physics.hpp"

using namespace std;

//Function to go through data files and calculate W and Q2
//Fill in W vs Q2 hist and save to output root file
//
void WvsQ2(char *fin, char *RootFile_output){
	gROOT->Reset();
	int current_event;
	int num_of_events;
	int total_events = 0;
	int files_in_lis = 2466;

	TFile *myFile;
	TFile *RootOutputFile;
	TTree *myTree;
	int number_cols = 0;
	int number_files = 0;
	char rootFile[500];

	RootOutputFile = new TFile(RootFile_output,"RECREATE");

	cout << "Analyzing file " << fin << endl;

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

			myTree->GetEntry(current_event);

			// Changed id to id[0] because scattered elctron should be first particle (i.e. id[0])
			if (id[0] == ELECTRON && gpart > 1 && stat[0] > 0 && q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0){
				E_prime = E_calc(p[0],cx[0],cy[0],cz[0]);
				Q2 = Q2_calc(cz[0],E_prime);
				W = W_calc(E_prime);
				WvsQ2_Fill();
			}


			#pragma omp parallel for
			for(int event_number = 0; event_number < gpart; event_number++){
				ID = id[event_number];
				P = P_calc(p[event_number],cx[event_number],cy[event_number],cz[event_number]);
				Beta = b[event_number];
				Energy = E_calc(p[event_number],cx[event_number],cy[event_number],cz[event_number],id[event_number]);
				//	If event is a scattered electron find it's energy
				//	and from it's energy find Q2 and W and then fill
				//	histograms for those variables

				/*//if (id[event_number] == ELECTRON){
				// Changed id to id[0] because scattered elctron should be first particle (i.e. id[0])
				if (id[0] == ELECTRON && gpart > 1 && stat[0] > 0 && q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0){
					E_prime = E_calc(p[event_number],cx[event_number],cy[event_number],cz[event_number]);

					Q2 = Q2_calc(cz[event_number],E_prime);
					W = W_calc(E_prime);
					if(W != W ) {
						cout << "here" << endl;
					} else {
						WvsQ2_Fill();
						//MomVsBeta_Fill();
					}

				}*/
				FillHist();
				MomVsBeta_Fill();
			}
			current_event++; 		  	// increment event counter
			total_events++; 			// increment total event counter
		}
		myTree->Delete(); 						// delete Tree object
		myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
		number_files++; 						// increment file counter

	}
	RootOutputFile->cd();
	WvsQ2_Write();
	MomVsBeta_Write();
	WriteHists();
	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
	cout<<total_events<<" events in "<<number_files<< " files."<<endl; // print out stats
}
#endif
