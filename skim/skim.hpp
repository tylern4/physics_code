/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef SKIM_H_GUARD
#define SKIM_H_GUARD
#include "main.h"

void skim(char *fin, char *RootFile_output){

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;
	bool electron_cuts;

	Double_t W, Q2;
	Int_t MyID[MAX_PARTS];

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
	TBranch *W_branch = skim->Branch("W",&W);
	TBranch *Q2_branch = skim->Branch("Q2",&Q2);
	TBranch *ID_branch = skim->Branch("MyID",&MyID);


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

		if (electron_cuts){
			W = W_calc(e_mu,e_mu_prime);
			Q2 = Q2_calc(e_mu,e_mu_prime);
			MyID[0] = ELECTRON;
			skim->Fill(); //Fill the banks after the skim
		}
	}

	//
	//end stuff
	chain.Reset();						// delete Tree object

	RootOutputFile->cd();
	RootOutputFile->Write();
	RootOutputFile->Close();

}

void skim(char* fin, char* RootFile_output, double mean, double sigma){

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;
	bool electron_cuts, MM_cut;
	//From missing_mass.hpp :: missing_mass_calc()
	MissingMass MissingMassNeutron;

	Float_t W, Q2, MM;
	Int_t MyID[MAX_PARTS];
	Double_t dt_proton[MAX_PARTS], dt_pip[MAX_PARTS];
	Int_t num_of_pis;

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
	TBranch *W_branch = skim->Branch("W",&W);
	TBranch *Q2_branch = skim->Branch("Q2",&Q2);
	TBranch *MM_branch = skim->Branch("MM",&MM);
	TBranch *ID_branch = skim->Branch("MyID",MyID);
	TBranch *DeltaT_P_branch = skim->Branch("DeltaT_P",dt_proton);
	TBranch *DeltaT_Pip_branch = skim->Branch("DeltaT_Pip",dt_pip);
	TBranch *NumPI_branch = skim->Branch("NumPI",&num_of_pis);


	for (int current_event = 0; current_event < num_of_events; current_event++) {
		chain.GetEntry(current_event);
		electron_cuts = true;
		//electron cuts
		electron_cuts &= (id[0] == ELECTRON); //First particle is electron
		electron_cuts &= (gpart > 0); //Number of good particles is greater than 0
		electron_cuts &= (stat[0] > 0); //First Particle hit stat
		electron_cuts &= ((int)q[0] == -1); //First particle is negative Q
		electron_cuts &= (sc[0] > 0); //First Particle hit sc
		electron_cuts &= (dc[0] > 0); // ``` ``` ``` d
		electron_cuts &= (ec[0] > 0); // ``` ``` ``` ec
		electron_cuts &= (dc_stat[dc[0]-1] > 0); //??

		e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
		e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

		delta_t_array(dt_proton,MASS_P);
		delta_t_array(dt_pip,MASS_PIP);

		for(int part_num = 1; part_num < gpart; part_num++){
			num_of_pis = 0;
			if(id[part_num] == PIP){
				num_of_pis++;
				TLorentzVector gamma_mu = (e_mu - e_mu_prime);
				MissingMassNeutron.MissingMassPxPyPz(p[part_num]*cx[part_num],p[part_num]*cy[part_num],p[part_num]*cz[part_num]);
				MissingMassNeutron = MissingMassNeutron.missing_mass(gamma_mu);
			}
		}
		MM = (MissingMassNeutron.mass >=0 ) ? MissingMassNeutron.mass : NaN;
		MM = (num_of_pis == 1) ? MissingMassNeutron.mass : NaN;

		MM_cut = true;
		MM_cut &= (MM == MM); //removes NaN

		MM_cut &= (MM <= mean + 10 * sigma);
		MM_cut &= (MM >= mean - 10 * sigma);

		if (electron_cuts){ //&& MM_cut
			W = W_calc(e_mu,e_mu_prime);
			Q2 = Q2_calc(e_mu,e_mu_prime);
			MyID[0] = ELECTRON;
			skim->Fill(); //Fill the banks after the skim
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
