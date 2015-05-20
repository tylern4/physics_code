/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef WVSQ2_H_GUARD
#define WVSQ2_H_GUARD
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
#include "main.h"
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

	int num_elec = 0, num_pip =0;
	int files_in_lis = 2486;
	//int files_in_lis = 153;

	TFile *myFile;
	TFile *RootOutputFile;
	TTree *myTree;
	int number_cols = 0;
	int number_files = 0;
	char rootFile[500];

	TVector3 e_mu_prime_3;
	TLorentzVector e_mu_prime;
	TLorentzVector e_mu(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0);

	TVector3 Particle3(0.0,0.0,0.0);
	TLorentzVector Particle4(0.0,0.0,0.0,0.0);

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

			myTree->GetEntry(current_event);

			// Check to see whether the first particle is an Electron
			// Changed id to id[0] because scattered elctron should be first particle (i.e. id[0])
			if (id[0] == ELECTRON && gpart > 1 && stat[0] > 0 && q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0 /****/ && b[0] < 1 ){
				//Setup scattered electron 4 vector
				e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cy[0]);	
				e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

				//Get energy of scattered elctron from 4 vector and calculate Q2 and W
				E_prime = e_mu_prime.E();
				Q2 = Q2_calc(e_mu,e_mu_prime);
				W = W_calc(e_mu,e_mu_prime);

				//This is my testing to see if momentum calculations are equal
				P = e_mu_prime.P();
				P1 = p[0];
				PminusP_Fill();

				xb = xb_calc(Q2,E_prime);

				WvsQ2_Fill();
				
/*
check beta vs P for a few different cases
	only E
	only P
	only pi+

Check how P is calculated
	fill Tlorentz and get P that way
*/
				#pragma omp parallel for
				for(int event_number = 1; event_number < gpart; event_number++){
					Particle3.SetXYZ(p[event_number]*cx[event_number], p[event_number]*cy[event_number], p[event_number]*cz[event_number]);
					Particle4.SetVectM(Particle3,Get_Mass(id[event_number]));
					Energy = Particle4.E();
					Beta = b[event_number];
					MomVsBeta_Fill();

					if(id[event_number] == PIP && q[event_number] == 1 /***/ && b[event_number] < 1) {
						Fill_e_pi_found(W_calc(e_mu,e_mu_prime),Q2_calc(e_mu,e_mu_prime),Particle4.P(),b[event_number]);

						for (int event_number_1 = 1; event_number_1 < gpart; event_number_1++){
							if(id[event_number_1] == PROTON && q[event_number_1] == 1 /****/ && b[event_number_1] < 1) {
								Fill_e_proton_pi_found(W_calc(e_mu,e_mu_prime),Q2_calc(e_mu,e_mu_prime),Particle4.P(),b[event_number]);
							}
						}
					} else if (id[event_number] == PROTON && q[event_number] == 1 /***/ && b[event_number] < 1 ){
						Fill_e_proton_found(W_calc(e_mu,e_mu_prime),Q2_calc(e_mu,e_mu_prime),Particle4.P(),b[event_number]);
					} 

				} 
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
	Write_found_hists();
	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
	cout<<blue<<endl<<total_events<<" events in "<<number_files<< " files."<<def<<endl; // print out stats
}
#endif
