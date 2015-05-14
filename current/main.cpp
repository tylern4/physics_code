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
#include "histo.h"
//#include "datahandeler.h" //dataHandeler() 
//#include "count_after_cut.h" //count_after_cut()
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <cstring>
//#include "stupid.h" //PrintEverything()
#include "time.h"
#include "physics.hpp"
#include "WvsQ2.hpp" //WvsQ2()

using namespace std;

int main(int argc, char **argv){
	
	gSystem->Load("libTree");
	char  infilename[128];
	char  outfilename[128];
	string minutes, hour;
	//char* minutes, hour;

	
	sprintf(infilename,"%s",argv[1]);
	
	//Either name outputfile or use time to name output file
	if(argc == 3 ) {
		sprintf(outfilename,"%s",argv[2]);
	} else {
		time_t currentTime;
		time(&currentTime); 
  		struct tm *localTime = localtime(&currentTime);  // Convert the current time to the local time;
  		if(localTime->tm_min < 10) {
  			minutes = "0" + to_string(localTime->tm_min);
  		}
  		else {
  			minutes = to_string(localTime->tm_min);
  		}

  		if(localTime->tm_hour > 10) {
  			hour = to_string(localTime->tm_hour - 12);
  		}
  		else {
  			hour = to_string(localTime->tm_hour);
  		}
  		string time = "outputFiles/release_" + to_string(localTime->tm_mon+1) + "-" + to_string(localTime->tm_mday) + "_"
	 	+ hour + "-" + minutes + ".root";
		sprintf(outfilename,"%s",time.c_str());
	}

	//dataHandeler(infilename,outfilename);
	WvsQ2(infilename,outfilename);

	return 0;
}
