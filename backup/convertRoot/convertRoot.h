#include "branches.h"

using namespace std;

void enoughArguments(int argc){
	if (argc < 3)
	{		
		cout<<"To use:"<<endl;
		cout<<"./convertRoot input.lis output.root"<<endl<<endl;
		exit;
	}	
}

void convert(char *fin, char *RootFile){
	Int_t ii;
	Int_t nEvents;	
	Int_t TotEvents = 0;
	
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
///
		if (ncols<0) break;
		myFile = new TFile(rootFile, "READ");
		myTree = (TTree *)myFile->Get("h10");
///
		getBranches(myTree); //change to be assocaited with branches from my root file

		nEvents = (Int_t)myTree->GetEntries();


		ii = 0; 

		while(ii<nEvents){

			myTree->GetEntry(ii);

			#pragma omp parallel for
			for(int j = 0; j < gpart; j++)
			{
				cout<<
			}

			ii++; 		  // increment event counter
			TotEvents++; // increment total event counter 
		}

		myTree->Delete(); 						// delete Tree object
		myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
		nfiles++; 								// increment file counter
		
	}
	rootOutFile->cd();
	rootOutFile->Write();
	rootOutFile->Close();
	cout<<TotEvents<<" events in "<<nfiles<< " files."<<endl; // print out stats

}