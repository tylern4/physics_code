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
////#include <omp.h>
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
	_e0->SetPxPyPzE(0,0,E1F_E0,E1F_E0);
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


// Count_after_cut *probably should be renamed?
//
//	Look at my channel and find the files with the best events.
//  Output of this should go to a function from Ye's goldenrun.C program to make
//	the appropriate graphs.
//

void count_after_cut(char *fin="all.lis", char *RootFile_output="outFile.root"){
	gROOT->Reset();
	Int_t current_event = 0, num_of_events = 0, total_events = 0, number_cols = 0;

	TFile *myFile;
	TFile *RootOutputFile;
	TTree *myTree;
	Int_t number_files = 0;
	char rootFile[500];

	//Define variables //change these for my channel
	float totalQ=0, qcurr, qprev, deltaq, q_temp;
	TLorentzVector ni_4vec(0, 0, 0, 0);
	TVector3 ni_3vec(0, 0, 0);
	TLorentzVector n_4vec(0,0,0,0);
	TVector3 n_3vec(0,0,0);
	TVector3 q_3vec(0,0,0);
	ni_4vec.SetVectM(ni_3vec, 0.9396);
	TLorentzVector ei_4vec(0, 0, 0, 0);
	TVector3 ei_3vec(0, 0, 2.039);
	ei_4vec.SetVectM(ei_3vec, 0.000511);
	TVector3 ef_3vec(0, 0, 0);
	TLorentzVector ef_4vec(0, 0, 0, 0);
	TVector3 pionf_3vec(0, 0, 0);
	TLorentzVector pionf_4vec(0, 0, 0, 0);
	TVector3 pf_3vec(0, 0, 0);
	TLorentzVector pf_4vec(0, 0, 0, 0);
	TVector3 protonf_3vec(0, 0, 0);
	TLorentzVector protonf_4vec(0, 0, 0, 0);
	TLorentzVector w_4vec(0,0,0,0);
	TLorentzVector w_n_q_4vec(0,0,0,0);
	TLorentzVector ZERO_4vec(0,0,0,0);
	TLorentzVector pionf_CM_4vec(0,0,0,0);
	//TLorentzVector cmq_4vec(0,0,0,0);
	//TLorentzVector cmPip_4vec;
	TVector3 b_3vec(0,0,0);
	TLorentzVector m_4vec(0,0,0,0);
	//End of variables

	RootOutputFile = new TFile(RootFile_output,"RECREATE"); //Open rootfile for output if it's not there create it.

	FILE *input_file = fopen(fin,"r"); //open the *.lis input file.

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
			//if (current_event%10000 == 0)	cout<<current_event<<"/"<<num_of_events<<endl;
			//Ye's line 131
			myTree->GetEntry(current_event); //Should load all the variables properly
			//Total Farday cup charge
			
			qcurr = q_l;

			if (q_l > 0.)
			{
				if (qcurr > qprev)
				{
					deltaq = qcurr - qprev;
					totalQ += deltaq;
					cout<<"qcurr="<<qcurr<<"qprev="<<qprev<<"deltaq="<<deltaq<<endl;
				}
				qprev = qcurr;
			}

/*			#pragma omp parallel for
			for(int j = 0; j < gpart; j++){

				/*
				Here is where I need to make some additions from Ye's code.
				Copy over the calculaions from count_after_cut.C and add plots from plot_golden.C
				



			} */

			current_event++; 		  	// increment event counter
			total_events++; 			// increment total event counter 
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