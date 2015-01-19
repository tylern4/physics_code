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
#include "TMath.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "convertRoot.h"

using namespace std;
#ifndef __CINT__
int main(int argc, char **argv){
	char  infilename[128];
	char  outfilename[128];

	gSystem->Load("libTree");
	gROOT->Reset();

	enoughArguments(argc);
	sprintf(infilename,"%s",argv[1]);
	sprintf(outfilename,"%s",argv[2]);

	convert(infilename,outfilename);



	return 0;
}
#endif
