/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

//#define PI 3.14159265;
//Only My Includes. All others in main.h
#include "main.h"
#include "TStopwatch.h"
#include "classes.hpp"
#include "missing_mass.hpp"
#include "missing_mass_hists.hpp"
#include "physics.hpp"

#include "json.hpp"
using json = nlohmann::json;

#include "datahandeler.h"

using namespace std;

int main(int argc, char **argv){
	TStopwatch *Watch = new TStopwatch;
	Watch->Start();
	gStyle->SetOptFit(1111);

	if (argc == 3) {
		char* infilename = argv[1];
		char*  outfilename = argv[2];
		dataHandeler(infilename,outfilename);
	}

	Watch->Stop();
	cout << red << Watch->RealTime() << "sec" << def << endl;

	return 0;
}