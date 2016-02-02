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
	char  infilename_bad[128];

	string infilename = (char*)argv[1];
	sprintf(infilename_bad,"%s",infilename.c_str());

	make_db(infilename_bad);

	Watch->Stop();
	cout << red << Watch->RealTime() << "sec" << def << endl;

	return 0;
}