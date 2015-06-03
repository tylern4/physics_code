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
	const Int_t ndims = 2;
   	Int_t bins[ndims] = {500, 500};
   	Double_t xmin[ndims] = {0., 0.};
   	Double_t xmax[ndims] = {3.25, 10.};
	THnSparse* hs = new THnSparseD("hs", "Sparse Histogram", ndims, bins, xmin, xmax);
	//THnSparse* testNsparse = new THnSparseF("testNsparse", "Sparse Histogram", 5, 100000, 0, 100000);

	/*
	SOOOOOOOO i think what happens is that you make an nsparce with arrays as everything instead of the usual number except ndims which is the size of arrays
	basically an of size ndims where each column is a th1d
	fill it by filling a single row in an array
	array_to_fill[colum1,colum2,colum3,....,ndims];
	hs->Fill(array_to_fill);
	get values back in same way but with give bin content, puts the output into another array
	hs->GetBin(bin_number,output_array);



	make an arrary with array[ID,W,Q2,xb,E]
	take for loop from WvsQ2 program and fill the array

	do the writes and try a project(W,Q2) ie TH2D* WvsQ2 = hs->Projection(2,3);
	*/

	Double_t x[ndims];
	/*long stupid = 10000000;
	for (long i = 0; i < stupid; ++i) {
		for (Int_t d = 0; d < ndims; ++d) {
        	switch (d) {
         	case 0: x[d] = gRandom->Gaus()*2 + 3.; break;
         	case 1:
         	case 2: 
         	case 3: x[d] = (x[d-1]*x[d-1] - 1.5)/1.5 + (0.5*gRandom->Rndm()); break;
         	default: x[d] = sin(gRandom->Gaus()*i/1000.) + 1.;
        	}
        }
        loadbar(i, stupid);
        hs->Fill(x);
    }*/
int stupid;
cout << num_of_events << endl;
	for (int current_event = 0; current_event <= num_of_events; current_event++) {
		loadbar(current_event,num_of_events);
		chain.GetEntry(current_event);

		// Check to see whether the first particle is an Electron
		// Changed id to id[0] because scattered elctron should be first particle (i.e. id[0])
		if (id[0] == ELECTRON && gpart > 1 && stat[0] > 0 && (int)q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0 /*** && b[0] <= 1 /**/){
			//Setup scattered electron 4 vector
			e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
			//e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
			e_mu_prime.SetVectM(e_mu_prime_3, m[0]);

			//Get energy of scattered elctron from 4 vector and calculate Q2 and W
			x[0] = W_calc(e_mu, e_mu_prime);
			x[1] = Q2_calc(e_mu, e_mu_prime);
			stupid++;
		}
		hs->Fill(x);
	}
   TH2D* h2proj = hs->Projection(1,0);
cout << stupid << endl;
   //h2proj->Write();
    


//end stuff
	chain.Reset();						// delete Tree object

	RootOutputFile->cd();
//write stuff
	hs->Write();
	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
}
#endif
