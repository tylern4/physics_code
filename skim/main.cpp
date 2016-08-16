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
#include "delta_t.hpp"
#include "missing_mass.hpp"
#include "skim.hpp"

using namespace std;

int main(int argc, char **argv){
	TStopwatch *Watch = new TStopwatch;
	Watch->Start();
	gStyle->SetOptFit(1111);

	/*if (argc == 3) {
		char* infilename = argv[1];
		char*  outfilename = argv[2];
		skim(infilename,outfilename);
	}*/

	if (argc == 3) {
		char* infilename = argv[1];
		char*  outfilename = argv[2];
		double mean = 0.947655;
		double sigma = 0.00983651;
		skim(infilename,outfilename,mean,sigma);
	}

	Watch->Stop();
	//cout << red << Watch->RealTime() << "sec" << def << endl;

	return 0;
}