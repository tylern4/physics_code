/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

//#define PI 3.14159265;
//Only My Includes. All others in main.h
#include "main.h"
#include "classes.hpp"
#include "TStopwatch.h"
#include "physics.hpp"
#include "delta_t_cut.hpp"
#include "missing_mass.hpp"
#include "skim.hpp"

using namespace std;

int main(int argc, char **argv){
	TStopwatch *Watch = new TStopwatch;
	Watch->Start();
	gStyle->SetOptFit(1111);
	//bad work around until I fix using strings in datahandeler/WvsQ2
	char  infilename_bad[128];
	char  outfilename_bad[128];

	string infilename = (char*)argv[1];
	string outfilename = outputFileName(argc,argv);

	//bad work around until I fix using strings in datahandeler/WvsQ2
	sprintf(infilename_bad,"%s",infilename.c_str());
	sprintf(outfilename_bad,"%s",outfilename.c_str());
	skim(infilename_bad,outfilename_bad);
	//dataHandeler(infilename_bad,outfilename_bad);

	Watch->Stop();
	//cout << red << Watch->RealTime() << "sec" << def << endl;

	return 0;
}