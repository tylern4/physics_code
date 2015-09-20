/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
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
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "physics.hpp"
#include "THnSparse.h"
#include "TRandom.h"
#include "TH3.h"
#include "delta_t.hpp"
#include "delta_t_cut.hpp"
#include <thread> 

// Mashing together W vs Q2 and Delta T cuts into one file
// Saving the old files in a new folder to refer back to.
//

void dataHandeler(char *fin, char *RootFile_output){

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;
	bool cuts;

	TVector3 e_mu_prime_3;
	TLorentzVector e_mu_prime;
	TLorentzVector e_mu(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0);

	TVector3 Particle3(0.0,0.0,0.0);
	TLorentzVector Particle4(0.0,0.0,0.0,0.0);

	RootOutputFile = new TFile(RootFile_output,"RECREATE"); 

	TChain chain("h10");
	cout << blue <<"Analyzing file " << green << fin << def << bgdef << endl;

	FILE *input_file = fopen(fin,"r");
	if (input_file == NULL) perror ("Error opening file");

	while (1){
		number_cols = fscanf(input_file,"%s",rootFile);

		if (number_cols<0) break;
		chain.Add(rootFile);
	}

	getBranches(&chain);

	num_of_events = (int)chain.GetEntries();

	for (int current_event = 0; current_event < num_of_events; current_event++) {
		loadbar(current_event,num_of_events);
		chain.GetEntry(current_event);

		if (id[0] == ELECTRON && gpart > 0 && stat[0] > 0 && (int)q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0){
			//Setup scattered electron 4 vector
			e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
			cuts = ( b[0] != 0 );

			if(cuts){
				delta_t_cut();
				WvsQ2(e_mu,e_mu_prime);
				TLorentzVector gamma_mu = (e_mu - e_mu_prime);
				missing_mass(gamma_mu);
					
				/*std::thread thread1(WvsQ2,e_mu,e_mu_prime);
				std::thread thread2(delta_t_cut);
				thread1.detach();
				thread2.detach(); */
			}
		}
	}

	//end stuff
	chain.Reset();						// delete Tree object

	RootOutputFile->cd();
	//write stuff
	//Can probble make this into a single write function in main.h or write functions in histo.h
	//WvsQ2
	TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
	WvsQ2_folder->cd();
	WvsQ2_Write();

	TDirectory *MomVsBeta_folder = RootOutputFile->mkdir("Momentum vs beta");
	MomVsBeta_folder->cd();
	MomVsBeta_Write();

	//Delta_t Write
	TDirectory *delta_t_folder = RootOutputFile->mkdir("Delta t");
	delta_t_folder->cd();
	delta_t_Write();

	//Extra Write
	TDirectory *extras_folder = RootOutputFile->mkdir("Missing Mass");
	extras_folder->cd();
	Write_Missing_Mass();


	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
}
#endif
