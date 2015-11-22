/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD
#include "TMath.h"
#include "TTree.h"
#include "TROOT.h"
#include <stdlib.h>
#include <stdio.h>
#include <TLorentzVector.h>
#include "color.hpp"
#include <string.h>
#include <string>
#include <cstring>
#include "time.h"
//#include <omp.h>
#include <string.h>
#include <string>
#include <cstring>
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

using namespace std;

static const int MAX_PARTS = 10000;

//static const float PI = TMath::Pi();
//static const float D2R = PI/180.0;
static const float E1D_E0 = 4.802; //GeV
//static const float E1D_E0 = 2.03939; //GeV ///This is for Ye's data. Actually E1E_E0

static const float SOL = 29.9792458;
//misc. constants
static const float FSC = 0.00729735253;
static const float NA = 6.02214129E23; //Avigadro's number
static const float QE = 1.60217646E-19;  //Charge or electron

//particle codes, usually PDG codes, but always those used in BOS
static const int PROTON = 2212;
static const int NEUTRON = 2112;
static const int PIP = 211;
static const int PIM = -211;
static const int PI0 = 111;
static const int KP = 321;
static const int KM = -321;
static const int PHOTON = 22;
static const int ELECTRON = 11;

//PDG particle masses in GeV/c2
static const float MASS_P = 0.93827203;
static const float MASS_N = 0.93956556;
static const float MASS_E = 0.000511;
static const float MASS_PIP = 0.13957018;
static const float MASS_PIM = 0.13957018;
static const float MASS_PI0 = 0.1349766;
static const float MASS_KP = 0.493677;
static const float MASS_KM = 0.493677;
static const float MASS_G = 0.0;
static const float MASS_OMEGA = 0.78265;

//static const double NaN = (0.0/0.0);
static const double NaN = std::nan("1");

UChar_t npart;
//UChar_t evstat;
UInt_t  evntid;
Char_t  evtype;
Char_t  evntclas;
Char_t  evthel;
Int_t   evntclas2;
Float_t q_l;
Float_t t_l;
Float_t tr_time;
Float_t rf_time1;
Float_t rf_time2;

Int_t   gpart;
Short_t id[MAX_PARTS];   //[gpart]
Char_t  stat[MAX_PARTS];   //[gpart]
UChar_t dc[MAX_PARTS];   //[gpart]
UChar_t cc[MAX_PARTS];   //[gpart]
UChar_t sc[MAX_PARTS];   //[gpart]
UChar_t ec[MAX_PARTS];   //[gpart]
UChar_t lec[MAX_PARTS];   //[gpart]
Float_t p[MAX_PARTS];   //[gpart]
//Float_t m[MAX_PARTS];   //[gpart]
Char_t  q[MAX_PARTS];   //[gpart]
Float_t b[MAX_PARTS];   //[gpart]
Float_t cx[MAX_PARTS];   //[gpart]
Float_t cy[MAX_PARTS];   //[gpart]
Float_t cz[MAX_PARTS];   //[gpart]
Float_t vx[MAX_PARTS];   //[gpart]
Float_t vy[MAX_PARTS];   //[gpart]
Float_t vz[MAX_PARTS];   //[gpart]

//Each of the folowing has multiple parts 
 
Int_t   dc_part;		
UChar_t dc_sect[MAX_PARTS];   //[dc_part]
UChar_t dc_trk[MAX_PARTS];   //[dc_part]
Char_t  dc_stat[MAX_PARTS];   //[dc_part]
Float_t dc_xsc[MAX_PARTS];   //[dc_part]
Float_t dc_ysc[MAX_PARTS];   //[dc_part]
Float_t dc_zsc[MAX_PARTS];   //[dc_part]
Float_t dc_cxsc[MAX_PARTS];   //[dc_part]
Float_t dc_cysc[MAX_PARTS];   //[dc_part]
Float_t dc_czsc[MAX_PARTS];   //[dc_part]
Float_t dc_xec[MAX_PARTS];   //[dc_part]
Float_t dc_yec[MAX_PARTS];   //[dc_part]
Float_t dc_zec[MAX_PARTS];   //[dc_part]
Float_t dc_thcc[MAX_PARTS];   //[dc_part]
Float_t dc_c2[MAX_PARTS];   //[dc_part]

Int_t   ec_part;
UShort_t ec_stat[MAX_PARTS];   //[ec_part]
UChar_t ec_sect[MAX_PARTS];   //[ec_part]
Int_t   ec_whol[MAX_PARTS];   //[ec_part]
//Int_t   ec_inst[MAX_PARTS];   //[ec_part]
//Int_t   ec_oust[MAX_PARTS];   //[ec_part]
Float_t etot[MAX_PARTS];   //[ec_part]
Float_t ec_ei[MAX_PARTS];   //[ec_part]
Float_t ec_eo[MAX_PARTS];   //[ec_part]
Float_t ec_t[MAX_PARTS];   //[ec_part]
Float_t ec_r[MAX_PARTS];   //[ec_part]
Float_t ech_x[MAX_PARTS];   //[ec_part]
Float_t ech_y[MAX_PARTS];   //[ec_part]
Float_t ech_z[MAX_PARTS];   //[ec_part]
//Float_t ec_m2[MAX_PARTS];   //[ec_part]
//Float_t ec_m3[MAX_PARTS];   //[ec_part]
//Float_t ec_m4[MAX_PARTS];   //[ec_part]
Float_t ec_c2[MAX_PARTS];   //[ec_part]

Int_t   sc_part;
UChar_t sc_sect[MAX_PARTS];   //[sc_part]
UChar_t sc_hit[MAX_PARTS];   //[sc_part]
UChar_t sc_pd[MAX_PARTS];   //[sc_part]
UChar_t sc_stat[MAX_PARTS];   //[sc_part]
Float_t edep[MAX_PARTS];   //[sc_part]
Float_t sc_t[MAX_PARTS];   //[sc_part]
Float_t sc_r[MAX_PARTS];   //[sc_part]
Float_t sc_c2[MAX_PARTS];   //[sc_part]

Int_t   cc_part;
UChar_t cc_sect[MAX_PARTS];   //[cc_part]
UChar_t cc_hit[MAX_PARTS];   //[cc_part]
Int_t   cc_segm[MAX_PARTS];   //[cc_part]
UShort_t	nphe[MAX_PARTS];   //[cc_part]
Float_t cc_t[MAX_PARTS];   //[cc_part]
Float_t cc_r[MAX_PARTS];   //[cc_part]
Float_t cc_c2[MAX_PARTS];   //[cc_part]

/*Int_t   lac_part;
Int_t   lec_sect[MAX_PARTS];   //[lac_part]
Int_t   lec_hit[MAX_PARTS];   //[lac_part]
Int_t   lec_stat[MAX_PARTS];   //[lac_part]
Float_t lec_etot[MAX_PARTS];   //[lac_part]
Float_t lec_t[MAX_PARTS];   //[lac_part]
Float_t lec_r[MAX_PARTS];   //[lac_part]
Float_t lec_x[MAX_PARTS];   //[lac_part]
Float_t lec_y[MAX_PARTS];   //[lac_part]
Float_t lec_z[MAX_PARTS];   //[lac_part]
Float_t lec_c2[MAX_PARTS];   //[lac_part]
*/

//Color outputs
Color::Modifier red(Color::FG_RED);
Color::Modifier blue(Color::FG_BLUE);
Color::Modifier def(Color::FG_DEFAULT);
Color::Modifier green(Color::FG_GREEN);
Color::Modifier bgred(Color::BG_RED);
Color::Modifier bgblue(Color::BG_BLUE);
Color::Modifier bgdef(Color::BG_DEFAULT);
Color::Modifier bggreen(Color::BG_GREEN);


void loadbar(long x, long n){

	int w = 50;
    if ( (x != n) && (x % (n/100+1) != 0) ) return;
 
    float ratio  =  x/(float)n;
    int   c      =  ratio * w;
 
    cout << blue << " [";
    for (int x=0; x<c; x++) cout <<green << "=" << def;
    cout << green << ">" << def;
    for (int x=c; x<w; x++) cout << " ";
    cout << blue << (int)(ratio*100) << "%]\r" << def << flush;
}

double Square(double a){
	return a*a;
}

void getBranches(TTree* myTree){

	myTree->SetBranchAddress("npart", &npart);			//number of final particles
	//myTree->SetBranchAddress("evstat", &evstat);
	myTree->SetBranchAddress("evntid", &evntid);		//event number
	myTree->SetBranchAddress("evntclas", &evntclas);
	myTree->SetBranchAddress("q_l", &q_l);
	myTree->SetBranchAddress("t_l", &t_l);
	myTree->SetBranchAddress("tr_time", &tr_time);
	myTree->SetBranchAddress("gpart", &gpart);			//number of particles in a single event (geometric particles)
	myTree->SetBranchAddress("id", &id);				//particle ID of i'th element id[i]  
	myTree->SetBranchAddress("stat", &stat);
	myTree->SetBranchAddress("dc", &dc);
	myTree->SetBranchAddress("cc", &cc);
	myTree->SetBranchAddress("sc", &sc);
	myTree->SetBranchAddress("ec", &ec);
	myTree->SetBranchAddress("lec", &lec);
	myTree->SetBranchAddress("p", &p);					//momentum of i'th particle p[i] (GeV/C)
	//myTree->SetBranchAddress("m", &m);					//mass of i'th particle m[i] (GeV/C)
	myTree->SetBranchAddress("q", &q);					//charge of i'th particle q[i] (charge in e's 1,0,-1)
	myTree->SetBranchAddress("b", &b);					//Velocity of i'th particle b[i] (in terms of c) ie. Beta
	myTree->SetBranchAddress("cx", &cx);				//X direction cosine at origin
	myTree->SetBranchAddress("cy", &cy);				//Y direction cosine at origin			
	myTree->SetBranchAddress("cz", &cz);				//Z direction cosine at origin
	myTree->SetBranchAddress("vx", &vx);				//X coordinate of vertex (cm)
	myTree->SetBranchAddress("vy", &vy);				//y coordinate of vertex (cm)
	myTree->SetBranchAddress("vz", &vz);				//z coordinate of vertex (cm)
	myTree->SetBranchAddress("dc_part", &dc_part);
	myTree->SetBranchAddress("dc_sect", &dc_sect);
	myTree->SetBranchAddress("dc_trk", &dc_trk);
	myTree->SetBranchAddress("dc_stat", &dc_stat);
	myTree->SetBranchAddress("dc_xsc", &dc_xsc);
	myTree->SetBranchAddress("dc_ysc", &dc_ysc);
	myTree->SetBranchAddress("dc_zsc", &dc_zsc);
	myTree->SetBranchAddress("dc_cxsc", &dc_cxsc);
	myTree->SetBranchAddress("dc_cysc", &dc_cysc);
	myTree->SetBranchAddress("dc_czsc", &dc_czsc);
	myTree->SetBranchAddress("dc_c2", &dc_c2);
	myTree->SetBranchAddress("ec_part", &ec_part);
	myTree->SetBranchAddress("ec_stat", &ec_stat);
	myTree->SetBranchAddress("ec_sect", &ec_sect);
	myTree->SetBranchAddress("ec_whol", &ec_whol);
	//myTree->SetBranchAddress("ec_inst", &ec_inst);
	//myTree->SetBranchAddress("ec_oust", &ec_oust);
	myTree->SetBranchAddress("etot", &etot);
	myTree->SetBranchAddress("ec_ei", &ec_ei);
	myTree->SetBranchAddress("ec_eo", &ec_eo);
	myTree->SetBranchAddress("ec_t", &ec_t);
	myTree->SetBranchAddress("ec_r", &ec_r);
	myTree->SetBranchAddress("ech_x", &ech_x);
	myTree->SetBranchAddress("ech_y", &ech_y);
	myTree->SetBranchAddress("ech_z", &ech_z);
	//myTree->SetBranchAddress("ec_m2", &ec_m2);
	//myTree->SetBranchAddress("ec_m3", &ec_m3);
	//myTree->SetBranchAddress("ec_m4", &ec_m4);
	myTree->SetBranchAddress("ec_c2", &ec_c2);
	myTree->SetBranchAddress("sc_part", &sc_part);
	myTree->SetBranchAddress("sc_sect", &sc_sect);
	myTree->SetBranchAddress("sc_hit", &sc_hit);
	myTree->SetBranchAddress("sc_pd", &sc_pd);
	myTree->SetBranchAddress("sc_stat", &sc_stat);
	myTree->SetBranchAddress("edep", &edep);
	myTree->SetBranchAddress("sc_t", &sc_t);
	myTree->SetBranchAddress("sc_r", &sc_r);
	myTree->SetBranchAddress("sc_c2", &sc_c2);
	myTree->SetBranchAddress("cc_part", &cc_part);
	myTree->SetBranchAddress("cc_sect", &cc_sect);
	myTree->SetBranchAddress("cc_hit", &cc_hit);
	myTree->SetBranchAddress("cc_segm", &cc_segm);
	myTree->SetBranchAddress("nphe", &nphe);
	myTree->SetBranchAddress("cc_t", &cc_t);
	myTree->SetBranchAddress("cc_r", &cc_r);
	myTree->SetBranchAddress("cc_c2", &cc_c2);
	/*myTree->SetBranchAddress("lac_part", &lac_part);
	myTree->SetBranchAddress("lec_sect", &lec_sect);
	myTree->SetBranchAddress("lec_hit", &lec_hit);
	myTree->SetBranchAddress("lec_stat", &lec_stat);
	myTree->SetBranchAddress("lec_etot", &lec_etot);
	myTree->SetBranchAddress("lec_t", &lec_t);
	myTree->SetBranchAddress("lec_r", &lec_r);
	myTree->SetBranchAddress("lec_x", &lec_x);
	myTree->SetBranchAddress("lec_y", &lec_y);
	myTree->SetBranchAddress("lec_z", &lec_z);
	myTree->SetBranchAddress("lec_c2", &lec_c2); */

	myTree->SetBranchStatus("*",1);
} 

string outputFileName(int argc, char **argv){
	string minutes, hour;
	if(argc == 3 ) {
		return (char*)argv[2];
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
  		string time = "outputFiles/release_" + to_string(localTime->tm_mon+1) 
  						+ "-" + to_string(localTime->tm_mday) + "_" + hour + "-" + minutes + ".root";
		return time;
	}
}

#endif
