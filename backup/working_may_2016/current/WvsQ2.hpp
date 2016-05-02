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

//Function to go through data files and calculate W and Q2
//Fill in W vs Q2 hist and save to output root file
//
void WvsQ2(TLorentzVector e_mu, TLorentzVector e_mu_prime){
	TVector3 Particle3(0.0,0.0,0.0);
	TLorentzVector Particle4(0.0,0.0,0.0,0.0);
	int numOfPis = 0, numOfProtons = 0;
	double W, Q2;
	W = W_calc(e_mu, e_mu_prime);
	Q2 = Q2_calc(e_mu, e_mu_prime);

	WvsQ2_Fill(e_mu_prime.E(),W,Q2,xb_calc(Q2, e_mu_prime.E()));

	for(int event_number = 0; event_number < gpart; event_number++){
		//Get particles 3 and 4 vector for current event.
		Particle3.SetXYZ(p[event_number]*cx[event_number], p[event_number]*cy[event_number], p[event_number]*cz[event_number]);
		Particle4.SetVectM(Particle3,Get_Mass(id[event_number]));

		if (Particle4.P() != 0 ){
			MomVsBeta_Fill(Particle4.E(),Particle4.P(),b[event_number]);
			if (q[event_number] == 1){
				MomVsBeta_Fill_pos(Particle4.P(),b[event_number]);
			} else if(q[event_number] == -1) {
				MomVsBeta_Fill_neg(Particle4.P(),b[event_number]);
			}

			if(id[event_number] == PIP && (int)q[event_number] == 1 ) {
				Fill_e_pi_found_WQ2(W,Q2);
				Fill_e_pi_found_P(Particle4.P(),b[event_number]);
				numOfPis++;
			} else if (id[event_number] == PROTON && (int)q[event_number] == 1 ){
				Fill_e_proton_found_WQ2(W,Q2);
				Fill_e_proton_found_P(Particle4.P(),b[event_number]);
				numOfProtons++;
			} 
		}
	}
	if (numOfProtons == 1 && gpart == 2) Fill_e_proton_only_WQ2(W,Q2);
	if (numOfPis == 1 && numOfProtons == 1 && gpart == 3) Fill_e_proton_pi_only_WQ2(W,Q2);
}
#endif
