/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

//#define PI 3.14159265;
//Only My Includes. All others in main.h
#include "main.h"
#include "WvsQ2_hists.hpp"
#include "missing_mass_hists.hpp"
#include "momentum_hists.hpp"
//#include "delta_t_hist.hpp"
#include "physics.hpp"
#include "WvsQ2.hpp"
#include "TStopwatch.h"
#include "delta_t_cut.hpp"
#include "missing_mass.hpp"
#include "datahandeler.h"

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

	Watch->Stop();
	cout << red << Watch->RealTime() << "sec" << def << endl;

	return 0;
}