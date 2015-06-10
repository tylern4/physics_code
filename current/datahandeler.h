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
//#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "physics.hpp"
#include "THnSparse.h"
#include "TRandom.h"
#include "TH3.h"

using namespace std;

//Function to go through data files and calculate W and Q2
//Fill in W vs Q2 hist and save to output root file
//
void dataHandeler(char *fin, char *RootFile_output){

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;

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
//start stuff
	const Int_t ndims = 5;
   	Int_t bins[ndims] = {500, 500, 500, 500, 500};
   	Double_t xmin[ndims] = {0., 0., 0., 0, -M_PI};
   	Double_t xmax[ndims] = {3.25, 10., 4.4, M_PI/2.0, M_PI};
   	Double_t x[ndims];

	const Int_t ndims_N = 5;
   	Int_t bins_N[ndims_N] = {500, 500, 500, 500, 500};
   	Double_t xmin_N[ndims_N] = {0, -M_PI, 0, 0, 0};
   	Double_t xmax_N[ndims_N] = {M_PI/2.0, M_PI, 1.2, 4.0, 2};
   	Double_t x_N[ndims_N];


	THnSparse* Electron_NSparse = new THnSparseD("Electron_NSparse", "Histogram", ndims, bins, xmin, xmax);
	THnSparse* Particle_NSparse = new THnSparseD("Particle_NSparse", "Histogram", ndims_N, bins_N, xmin_N, xmax_N);
	//THnSparse* testNsparse = new THnSparseF("testNsparse", "Sparse Histogram", 5, 100000, 0, 100000);

	/*
	SOOOOOOOO i think what happens is that you make an nsparce with arrays as everything instead of the usual number except ndims which is the size of arrays
	basically an of size ndims where each column is a th1d
	fill it by filling a single row in an array
	array_to_fill[colum1,colum2,colum3,....,ndims];
	Electron_NSparse->Fill(array_to_fill);
	get values back in same way but with give bin content, puts the output into another array
	Electron_NSparse->GetBin(bin_number,output_array);



	make an arrary with array[ID,W,Q2,xb,E]
	take for loop from WvsQ2 program and fill the array

	do the writes and try a project(W,Q2) ie TH2D* WvsQ2 = Electron_NSparse->Projection(2,3);
	*/
	
	Electron_NSparse->GetAxis(0)->SetTitle(" W ");
	Electron_NSparse->GetAxis(1)->SetTitle(" Q2 ");
	Electron_NSparse->GetAxis(2)->SetTitle(" P ");
	Electron_NSparse->GetAxis(3)->SetTitle(" #theta ");
	Electron_NSparse->GetAxis(4)->SetTitle(" #phi ");

	Particle_NSparse->GetAxis(0)->SetTitle(" #theta ");
	Particle_NSparse->GetAxis(1)->SetTitle(" #phi ");
	Particle_NSparse->GetAxis(2)->SetTitle(" #beta ");
	Particle_NSparse->GetAxis(3)->SetTitle(" P ");
	Particle_NSparse->GetAxis(4)->SetTitle(" Perp ");

	for (int current_event = 0; current_event <= num_of_events; current_event++) {
		loadbar(current_event,num_of_events);
		chain.GetEntry(current_event);

		// Check to see whether the first particle is an Electron
		// Changed id to id[0] because scattered elctron should be first particle (i.e. id[0])
		if (id[0] == ELECTRON && gpart > 1 && stat[0] > 0 && (int)q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0){
			//Setup scattered electron 4 vector
			e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
			//e_mu_prime.SetVectM(e_mu_prime_3, m[0]);

			//Get energy of scattered elctron from 4 vector and calculate Q2 and W
			x[0] = W_calc(e_mu, e_mu_prime);
			x[1] = Q2_calc(e_mu, e_mu_prime);
			x[2] = e_mu_prime.P();
			x[3] = e_mu_prime.Theta();
			x[4] = e_mu_prime.Phi();

			//#pragma omp parallel for
			for(int event_number = 1; event_number < gpart; event_number++){
				//Get particles 3 and 4 vector for current event.
				Particle3.SetXYZ(p[event_number]*cx[event_number], p[event_number]*cy[event_number], p[event_number]*cz[event_number]);
				Particle4.SetVectM(Particle3,Get_Mass(id[event_number]));
				//Particle4.SetVectM(Particle3,m[event_number]);

				x_N[0] = Particle4.Theta();
				x_N[1] = Particle4.Phi(); 
				x_N[2] = b[event_number];
				x_N[3] = Particle4.P();
				x_N[4] = Particle4.Perp();

				Particle_NSparse->Fill(x_N);
			}

			Electron_NSparse->Fill(x);
		}
		//Electron_NSparse->Fill(x);
	}

	//TH2D* h2proj_1 = Electron_NSparse->Projection(1,0);
	for (int j = 0; j < ndims; ++j){
		for (int jj = 0; jj < ndims; ++jj){
			if(j != jj) TH2D* h2proj = Electron_NSparse->Projection(j,jj);
		}
	}
	for (int k = 0; k < ndims_N; ++k){
		for (int kk = 0; kk < ndims_N; ++kk){
			if(k != kk) TH2D* h2proj = Particle_NSparse->Projection(k,kk);
		}
	}

//end stuff
	chain.Reset();						// delete Tree object

	RootOutputFile->cd();
//write stuff
	Electron_NSparse->Write();
	Particle_NSparse->Write();
	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
}
#endif
