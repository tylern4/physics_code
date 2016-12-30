/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

//Only My Includes. All others in main.h
#include "main.h"
#include "TStopwatch.h"
#include "classes.hpp"
#include "delta_t_hist.hpp" 
#include "WvsQ2_hists.hpp"
#include "momentum_hists.hpp"
#include "missing_mass.hpp"
#include "missing_mass_hists.hpp"
#include "physics.hpp"
#include "delta_t_cut.hpp"
#include "ec_cut.hpp"
#include "cc_cut.hpp"
#include "fid_hists.hpp"
#include "datahandeler.hpp"

using namespace std;

int main(int argc, char **argv){
	TStopwatch *Watch = new TStopwatch;
	Watch->Start();
	gStyle->SetOptFit(1111);
	makeHists_delta_t();
	makeHists_CC();
	makeHists_fid();

	if (argc == 3) {
		char* infilename = argv[1];
		char*  outfilename = argv[2];
		dataHandeler(infilename,outfilename,true);
	}

	if (argc == 4) {
		char* infilename = argv[1];
		char*  outfilename = argv[2];
		dataHandeler(infilename,outfilename,false);
	}

	Watch->Stop();
	cout << red << Watch->RealTime() << "sec" << def << endl;

	return 0;
}