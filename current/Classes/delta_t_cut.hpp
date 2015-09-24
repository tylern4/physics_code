/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef DELTA_T_CUT_H_GUARD
#define DELTA_T_CUT_H_GUARD
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
#include "delta_t.hpp"
#include "delta_t_hist.hpp"
#include "TProof.h"
#include "TProfile.h"
#include "delta_t.cpp"
//Function to go through data files and calculate W and Q2
//Fill in W vs Q2 hist and save to output root file
//
void delta_t_cut(){
	D_T cuts;
	TVector3 Particle3(0.0,0.0,0.0);
	TLorentzVector Particle4(0.0,0.0,0.0,0.0);
	double delta_t_P, delta_t_PIP, delta_t_ELECTRON;
	double electron_vertex = cuts.vertex_time(sc_t[sc[0]-1], sc_r[sc[0]-1], 1.0);

	for(int event_number = 0; event_number < gpart; event_number++){
		//Get particles 3 and 4 vector for current event.
		Particle3.SetXYZ(p[event_number]*cx[event_number], p[event_number]*cy[event_number], p[event_number]*cz[event_number]);
		Particle4.SetVectM(Particle3,Get_Mass(id[event_number]));

		delta_t_P = cuts.delta_t(electron_vertex, MASS_P, p[event_number], sc_t[sc[event_number]-1], sc_r[sc[event_number]-1]);
		delta_t_PIP = cuts.delta_t(electron_vertex, MASS_PIP, p[event_number], sc_t[sc[event_number]-1], sc_r[sc[event_number]-1]);
		delta_t_ELECTRON = cuts.delta_t(electron_vertex, MASS_E, p[event_number], sc_t[sc[event_number]-1], sc_r[sc[event_number]-1]);

		if (event_number == 0 && id[0] == ELECTRON && gpart > 0 && stat[0] > 0 && (int)q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0) {
			delta_t_Fill(Particle4.P(), delta_t_PIP, 8);
		}

		if (Particle4.P() != 0 && (int)q[event_number] == 1) {
			delta_t_Fill(Particle4.P(),delta_t_P,3);

			delta_t_Fill(Particle4.P(), delta_t_PIP,4);

			delta_t_Fill(Particle4.P(), delta_t_ELECTRON, 5); 

			//If Pi+
			if(id[event_number] == PROTON && (int)q[event_number] == 1) {
				delta_t_Fill(Particle4.P(), delta_t_P, 1);
				delta_t_Fill(Particle4.P(), delta_t_ELECTRON, 7);
			//If Proton	
			} else if (id[event_number] == PIP && (int)q[event_number] == 1){
				delta_t_Fill(Particle4.P(), delta_t_PIP, 2);
				delta_t_Fill(Particle4.P(), delta_t_ELECTRON, 6);
			} 
		}
	} 

}
#endif
