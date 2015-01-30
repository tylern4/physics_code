/************************************************************************/
/*									*/
/*									*/
/*  Created by Nick Tyler						*/
/*	University Of South Carolina					*/
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
#include "main.h"
#include "TMath.h"
#include "histo.h"
#include "datahandeler.h"
//#include "golden_run.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>


using namespace std;

int main(int argc, char **argv){
	
	//float time1 = clock();
	gSystem->Load("libTree");

	char  infilename[128];
	char  outfilename[128];


	sprintf(infilename,"%s",argv[1]);
	sprintf(outfilename,"%s",argv[2]);

	dataHandeler(infilename,outfilename);
	count_after_cut(infilename,outfilename);
	//golden_run(infilename,outfilename);

	//Timer(time1);
	return 0;
}
