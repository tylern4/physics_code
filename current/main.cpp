/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

//#define PI 3.14159265;

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
//#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <cstring>
#include "time.h"
//mine
#include "WvsQ2_hists.hpp"
#include "extra_hists.hpp"
#include "momentum_hists.hpp"
#include "delta_t_hist.hpp"
#include "physics.hpp"
#include "WvsQ2.hpp"
#include "datahandeler.h"
#include "TStopwatch.h"
#include "delta_t_cut.hpp"

using namespace std;

int main(int argc, char **argv){
	TStopwatch *Watch = new TStopwatch;
	Watch->Start();

	//bad work around until I fix using strings in datahandeler/WvsQ2
	char  infilename_bad[128];
	char  outfilename_bad[128];

	string infilename = (char*)argv[1];
	string outfilename = outputFileName(argc,argv);

	//bad work around until I fix using strings in datahandeler/WvsQ2
	sprintf(infilename_bad,"%s",infilename.c_str());
	sprintf(outfilename_bad,"%s",outfilename.c_str());
	dataHandeler(infilename_bad,outfilename_bad);
	//WvsQ2(infilename_bad,outfilename_bad);
	//delta_t_cut(infilename_bad,outfilename_bad);

	Watch->Stop();
	cout << red << Watch->RealTime() << "sec" << def << endl;

	return 0;
}
