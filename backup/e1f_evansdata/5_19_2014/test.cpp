/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler, Canisius College								*/
/*						University Of South Carolina					*/
/************************************************************************/

#define PI 3.14159265;

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
#include "test.h"
#include "TMath.h"
#include "histo.h"
//
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
	Int_t ii;
	Int_t nEvents;	
	Int_t TotEvents = 0;

	TLorentzVector *_e0, *_p0, *_e1;//, *_p1;

	_e0 = new TLorentzVector();
	_p0 = new TLorentzVector();
	_e1 = new TLorentzVector();
	_e0->SetPxPyPzE(0,0,E1F_E0,E1F_E0);
	_p0->SetPxPyPzE(0,0,0,MASS_P);

	TFile *myFile;
	TFile *rootOutFile;
	TTree *myTree;
	Int_t ncols=0;
	Int_t nfiles = 0;
	char rootFile[500];

	rootOutFile = new TFile(RootFile,"RECREATE");
	



	cout << "Analyzing file " << fin << endl;


	FILE *in1 = fopen(fin,"r");
	if (in1 == NULL) perror ("Error opening file");

	while (1){

		ncols = fscanf(in1,"%s",rootFile); 

		if (ncols<0) break;
		myFile = new TFile(rootFile, "READ");
		myTree = (TTree *)myFile->Get("h10clone/h10");

		getBranches(myTree);

		nEvents = (Int_t)myTree->GetEntries();


		ii = 0; 

		while(ii<nEvents){

			myTree->GetEntry(ii);

			#pragma omp parallel for
			for(int j = 0; j < gpart; j++)
			{
				if(id[j] == 11){

					E =((Square(p[j])/(2*MASS_E))+ Square(MASS_E));
					Px = cx[j]*p[j];
					Py = cy[j]*p[j];
					Pz = cz[j]*p[j];

					x = vx[j];
					y = vy[j];
					z = vz[j];

					FillHist();

					_e1->SetPxPyPzE(Px,Py,Pz,E);


				}

			}




			ii++; 		  // increment event counter
			TotEvents++; // increment total event counter 
		}

		myTree->Delete(); 						// delete Tree object
		myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
		nfiles++; 								// increment file counter
		
	}
	rootOutFile->cd();
	WriteHists();
	rootOutFile->Write();
	rootOutFile->Close();
	fclose(in1); 							// close file with input file list
	cout<<TotEvents<<" events in "<<nfiles<< " files."<<endl; // print out stats

}



int main(int argc, char **argv){

	float time1 = clock();
	gSystem->Load("libTree");

	char  infilename[128];
	char  outfilename[128];


	sprintf(infilename,"%s",argv[1]);
	sprintf(outfilename,"%s",argv[2]);



	dataHandeler(infilename,outfilename);



	Timer(time1);
	return 0;
}