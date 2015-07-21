/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef WVSQ2_H_GUARD
#define WVSQ2_H_GUARD
#include "TTree.h"
#include <math.h>
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
#include "TMath.h"
//#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "physics.hpp"

using namespace std;

//Function to go through data files and calculate W and Q2
//Fill in W vs Q2 hist and save to output root file
//
void WvsQ2(TLorentzVector e_mu, TLorentzVector e_mu_prime){
	
	WvsQ2_Fill(e_mu_prime.E(),W_calc(e_mu, e_mu_prime),Q2_calc(e_mu, e_mu_prime),xb_calc(Q2_calc(e_mu,e_mu_prime), e_mu_prime.E() ) );

	for(int event_number = 0; event_number < gpart; event_number++){
		//Get particles 3 and 4 vector for current event.
		TVector3 Particle3(0.0,0.0,0.0);
		TLorentzVector Particle4(0.0,0.0,0.0,0.0);
		Particle3.SetXYZ(p[event_number]*cx[event_number], p[event_number]*cy[event_number], p[event_number]*cz[event_number]);
		Particle4.SetVectM(Particle3,Get_Mass(id[event_number]));
		if (Particle4.P() != 0 ){
			MomVsBeta_Fill(Particle4.E(),Particle4.P(),b[event_number]);
			if (q[event_number] == 1){
				MomVsBeta_Fill_pos(Particle4.P(),b[event_number]);
			} else if(q[event_number] == -1) {
				MomVsBeta_Fill_neg(Particle4.P(),b[event_number]);
			}

			if(id[event_number] == PIP && (int)q[event_number] == 1 /*&& sc[event_number] > 0 && dc[event_number] > 0*/) {
				Fill_e_pi_found(W_calc(e_mu,e_mu_prime),Q2_calc(e_mu,e_mu_prime),Particle4.P(),b[event_number]);
			} else if (id[event_number] == PROTON && (int)q[event_number] == 1 /*&& sc[event_number] > 0 && dc[event_number] > 0*/){
				Fill_e_proton_found(W_calc(e_mu,e_mu_prime),Q2_calc(e_mu,e_mu_prime),Particle4.P(),b[event_number]);
			} 
		} 
	}
}
#endif
