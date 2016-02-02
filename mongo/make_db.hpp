/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef MAKE_DB_H_GUARD
#define MAKE_DB_H_GUARD
#include "main.h"

void make_db(char *fin){

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	bool cuts, electron_cuts;

	TVector3 e_mu_prime_3;
	TLorentzVector e_mu_prime;
	TLorentzVector e_mu(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0);

	TVector3 Particle3(0.0,0.0,0.0);
	TLorentzVector Particle4(0.0,0.0,0.0,0.0);

	TChain chain("h10");
	cout << blue <<"Analyzing file " << green << fin << def << bgdef << endl;
	chain.AddFile(fin);

	getBranches(&chain);

	num_of_events = (int)chain.GetEntries();

	for (int current_event = 0; current_event < num_of_events; current_event++) {
		chain.GetEntry(current_event);

		if (id[0] == ELECTRON && gpart > 0 && stat[0] > 0 && (int)q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0){
			//Setup scattered electron 4 vector
			e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

			W = W_calc(e_mu,e_mu_prime);
			Q2 = Q2_calc(e_mu,e_mu_prime);

			add_db(W, Q2);
		}
	}

	chain.Reset();						// delete Tree object


}
#endif
