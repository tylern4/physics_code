/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef DELTA_T_CUT_H_GUARD
#define DELTA_T_CUT_H_GUARD
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
#include "delta_t.hpp"
#include "delta_t_hist.hpp"
#include "TProof.h"

using namespace std;

//Function to go through data files and calculate W and Q2
//Fill in W vs Q2 hist and save to output root file
//
void delta_t_cut(char *fin, char *RootFile_output){

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;
	double electron_vertex, delta_t_P, delta_t_PIP;

	//TVector3 e_mu_prime_3;
	//TLorentzVector e_mu_prime;
	TLorentzVector e_mu(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0);

	//TVector3 Particle3(0.0,0.0,0.0);
	//TLorentzVector Particle4(0.0,0.0,0.0,0.0);

	RootOutputFile = new TFile(RootFile_output,"UPDATE"); 
	//RootOutputFile = new TFile(RootFile_output,"RECREATE");

	cout << blue <<"Analyzing file " << green << fin << def << bgdef << endl;

	FILE *input_file = fopen(fin,"r");
	if (input_file == NULL) perror ("Error opening file");

	TChain chain("h10");
	//TProof *plite = TProof::Open("");

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

			electron_vertex = vertex_time(sc_t[sc[0]-1], sc_r[sc[0]-1], 1.0);

			for(int event_number = 0; event_number < gpart; event_number++){
				//Get particles 3 and 4 vector for current event.
				TVector3 Particle3(0.0,0.0,0.0);
				TLorentzVector Particle4(0.0,0.0,0.0,0.0);
				Particle3.SetXYZ(p[event_number]*cx[event_number], p[event_number]*cy[event_number], p[event_number]*cz[event_number]);
				Particle4.SetVectM(Particle3,Get_Mass(id[event_number]));
				delta_t_P = delta_t(electron_vertex, MASS_P, p[event_number], sc_t[sc[event_number]-1], sc_r[sc[event_number]-1]);
				delta_t_PIP = delta_t(electron_vertex, MASS_PIP, p[event_number], sc_t[sc[event_number]-1], sc_r[sc[event_number]-1]);
				if (Particle4.P() != 0 && (int)q[event_number] == 1)
				{
					delta_t_Fill(Particle4.P(),delta_t_P,3);
					delta_t_Fill(Particle4.P(), delta_t_PIP, 4);

					//If Pi+
					if(id[event_number] == PROTON && (int)q[event_number] == 1) {
						delta_t_Fill(Particle4.P(), delta_t_P, 1);
					//If Proton	
					} else if (id[event_number] == PIP && (int)q[event_number] == 1){
						delta_t_Fill(Particle4.P(), delta_t_PIP, 2);
					} 
				}
			} 
		}
	}
	chain.Reset();
	TDirectory *delta_t_folder = RootOutputFile->mkdir("Delta_t");
	delta_t_folder->cd();
	delta_t_Write();
	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
	cout << endl << blue << "Completed " << num_of_events << " in " << def;
}
#endif
