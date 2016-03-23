/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef SKIM_H_GUARD
#define SKIM_H_GUARD
#include "main.h"
// Mashing together W vs Q2 and Delta T cuts into one file
// Saving the old files in a new folder to refer back to.
//
void skim(char *fin, char *RootFile_output){

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;
	bool electron_cuts, strict_cuts;

	Delta_T dt;
	D_time dt_cuts;
	MissingMass MM;
	double MissMass, W, Q2;
	double dt_proton[MAX_PARTS], dt_pip[MAX_PARTS];
	int num_of_pis;

	TVector3 e_mu_prime_3;
	TLorentzVector e_mu_prime;
	TLorentzVector e_mu(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0);

	TVector3 Particle3(0.0,0.0,0.0);
	TLorentzVector Particle4(0.0,0.0,0.0,0.0);

	RootOutputFile = new TFile(RootFile_output,"RECREATE"); 

	TChain chain("h10");
	cout << blue <<"Analyzing file " << green << fin << def << bgdef << endl;
	chain.AddFile(fin);

	getBranches(&chain);

	num_of_events = (int)chain.GetEntries();

	TTree *skim = chain.CloneTree(0);
	TBranch *MissMass_branch = skim->Branch("MissingMass",&MissMass);
	TBranch *W_branch = skim->Branch("W",&W);
	TBranch *Q2_branch = skim->Branch("Q2",&Q2);
	TBranch *Delta_t_proton_branch = skim->Branch("dt_proton",dt_proton);
	TBranch *Delta_t_pip_branch = skim->Branch("dt_pip",dt_pip);

	for (int current_event = 0; current_event < num_of_events; current_event++) {
		chain.GetEntry(current_event);
		electron_cuts = true;
		//electron cuts
		electron_cuts &= (id[0] == ELECTRON); //First particle is electron
		electron_cuts &= (gpart > 0); //Number of good particles is greater than 0
		electron_cuts &= (stat[0] > 0); //First Particle hit stat
		electron_cuts &= ((int)q[0] == -1); //First particle is negative Q
		electron_cuts &= (sc[0] > 0); //First Particle hit sc
		electron_cuts &= (dc[0] > 0); // ``` ``` ``` dc
		electron_cuts &= (ec[0] > 0); // ``` ``` ``` ec
		electron_cuts &= (dc_stat[dc[0]-1] > 0); //??

		//Extra cuts
		strict_cuts = true;
		//strict_cuts &= (gpart == 2);

		if (electron_cuts && strict_cuts){
			//Setup scattered electron 4 vector
			e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
			dt_cuts.SetVertexTimes(sc_t[sc[0]-1],sc_r[sc[0]-1]);

			//strict_cuts = (b[0] != 0);
			//strict_cuts &= (id[1] == PIP);

			if(strict_cuts){
				W = W_calc(e_mu,e_mu_prime);
				Q2 = Q2_calc(e_mu,e_mu_prime);
				for(int part_num = 1; part_num < gpart; part_num++){
					
					dt_cuts.SetTandP(sc_t[sc[part_num]-1],sc_r[sc[part_num]-1],p[part_num]);
					dt_cuts = dt_cuts.delta_t_calc();
					Particle3.SetXYZ(p[part_num]*cx[part_num],p[part_num]*cy[part_num],p[part_num]*cz[part_num]);
					Particle4.SetVectM(Particle3, Get_Mass(id[part_num]));
					
					if (Particle4.P() != 0 && (int)q[part_num] == 1) {
						dt_proton[part_num] = dt_cuts.proton_time;
						dt_pip[part_num] = dt_cuts.pip_time;
					} else {
						dt_proton[part_num] = NaN;
						dt_pip[part_num] = NaN;
					}
					
					num_of_pis = 0;
					if(id[part_num] == PIP){
						num_of_pis++;
						TLorentzVector gamma_mu = (e_mu - e_mu_prime);
						MM.MissingMassPxPyPz(p[part_num]*cx[part_num],p[part_num]*cy[part_num],p[part_num]*cz[part_num]);
						MM = MM.missing_mass(gamma_mu);
					}
				}
				//If the number of pions is 1 MissMass = calcualted mass : else NaNs
				MissMass = (num_of_pis == 1) ? MM.mass : NaN;

				skim->Fill(); //Fill the banks after the skim
			}
		}
	}

	//
	//end stuff
	chain.Reset();						// delete Tree object

	RootOutputFile->cd();

	RootOutputFile->Write();
	RootOutputFile->Close();

}
#endif
