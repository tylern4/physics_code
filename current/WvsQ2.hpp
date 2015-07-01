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

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;

	//TVector3 e_mu_prime_3;
	//TLorentzVector e_mu_prime;
	TLorentzVector e_mu(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0);

	//TVector3 Particle3(0.0,0.0,0.0);
	//TLorentzVector Particle4(0.0,0.0,0.0,0.0);

	RootOutputFile = new TFile(RootFile_output,"RECREATE"); 

	cout << blue <<"Analyzing file " << green << fin << def << bgdef << endl;

	FILE *input_file = fopen(fin,"r");
	if (input_file == NULL) perror ("Error opening file");

	TChain chain("h10");
	while (1){
		number_cols = fscanf(input_file,"%s",rootFile);
		if (number_cols<0) break;
		chain.Add(rootFile);
	}

	getBranches(&chain);
	num_of_events = (int)chain.GetEntries();

	//#pragma omp parallel for
	for (int current_event = 0; current_event <= num_of_events; current_event++) {
		TVector3 e_mu_prime_3;
		TLorentzVector e_mu_prime;
		TLorentzVector total(0.0,0.0,0.0,0.0);

		loadbar(current_event,num_of_events);
		chain.GetEntry(current_event);

		// Check to see whether the first particle is an Electron
		// Changed id to id[0] because scattered elctron should be first particle (i.e. id[0])
		if ((id[0] == ELECTRON || id[0] == 0) && gpart > 1 && stat[0] > 0 && (int)q[0] == -1 && cc[0] > 0 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0 ){
			//Setup scattered electron 4 vector
			e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

			//Get energy of scattered elctron from 4 vector and calculate Q2 and W
			WvsQ2_Fill(e_mu_prime.E(),W_calc(e_mu, e_mu_prime),Q2_calc(e_mu, e_mu_prime),xb_calc(Q2_calc(e_mu,e_mu_prime), e_mu_prime.E() ) );
			total += e_mu_prime;

			//#pragma omp parallel for
			for(int event_number = 0; event_number <= gpart; event_number++){
				//Get particles 3 and 4 vector for current event.
				TVector3 Particle3(0.0,0.0,0.0);
				TLorentzVector Particle4(0.0,0.0,0.0,0.0);
				Particle3.SetXYZ(p[event_number]*cx[event_number], p[event_number]*cy[event_number], p[event_number]*cz[event_number]);
				Particle4.SetVectM(Particle3,Get_Mass(id[event_number]));

				MomVsBeta_Fill(Particle4.E(),Particle4.P(),b[event_number]);
				total += Particle4;

				//If Pi+
				if(id[event_number] == PIP && (int)q[event_number] == 1 /*&& sc[event_number] > 0 && dc[event_number] > 0*/) {
					Fill_e_pi_found(W_calc(e_mu,e_mu_prime),Q2_calc(e_mu,e_mu_prime),Particle4.P(),b[event_number]);
					//If Pi+ and Proton
					//#pragma omp parallel for
					for (int event_number_1 = 0; event_number_1 <= gpart; event_number_1++){
						if(id[event_number_1] == PROTON && (int)q[event_number_1] == 1 /*&& sc[event_number_1] > 0 && dc[event_number_1] > 0*/) {
							Fill_e_proton_pi_found(W_calc(e_mu,e_mu_prime),Q2_calc(e_mu,e_mu_prime),Particle4.P(),b[event_number]);
							//total += Particle4;
						}
					}
				//If Proton	
				} else if (id[event_number] == PROTON && (int)q[event_number] == 1 /*&& sc[event_number] > 0 && dc[event_number] > 0*/){
					Fill_e_proton_found(W_calc(e_mu,e_mu_prime),Q2_calc(e_mu,e_mu_prime),Particle4.P(),b[event_number]);
				} 
			} 
			//Fill_Missing_Mass((total-e_mu).M());
		}
	}
	chain.Reset();						// delete Tree object

	RootOutputFile->cd();
	WvsQ2_Write();
	MomVsBeta_Write();
	Write_found_hists();
	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
	cout << endl << blue << "Completed " << num_of_events << " in " << def;
}
#endif
