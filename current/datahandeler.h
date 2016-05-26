/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef DATAHANDELER_H_GUARD
#define DATAHANDELER_H_GUARD
#include "main.h"
// Mashing together W vs Q2 and Delta T cuts into one file
// Saving the old files in a new folder to refer back to.
//
void dataHandeler(char *fin, char *RootFile_output, bool first_run){

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;
	bool cuts, electron_cuts;
	ofstream cut_outputs;
	cut_outputs.open("outputFiles/cut_outputs.csv");
	cut_outputs << "Cut,Mean,Sigma" << endl;

	//From missing_mass.hpp :: missing_mass_calc()
	MissingMass MissingMassNeutron;

	int num_of_pis,num_of_proton;

	TVector3 e_mu_prime_3;
	TLorentzVector e_mu_prime;
	TLorentzVector e_mu(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0);

	TVector3 Particle3(0.0,0.0,0.0);
	TLorentzVector Particle4(0.0,0.0,0.0,0.0);
	//double W,Q2;
	//End declrare variables

	//Open outputfile
	RootOutputFile = new TFile(RootFile_output,"RECREATE"); 

	//Load chain from branch h10
	TChain chain("h10");
	cout << blue <<"Analyzing file " << green << fin << def << bgdef << endl;

	//Try to open lis file throw an error if it won't open
	FILE *input_file = fopen(fin,"r");
	if (input_file == NULL) perror ("Error opening file");

	//Go through lis file adding each files h10 branch to the chain
	while (1){
		number_cols = fscanf(input_file,"%s",rootFile);
		if (number_cols<0) break;
		chain.Add(rootFile);
	}
	//get branches from the chain
	getBranches(&chain);
	if(!first_run) getWQ2branch(&chain);

	num_of_events = (int)chain.GetEntries();

	for (int current_event = 0; current_event < num_of_events; current_event++) {
		//update loadbar and get current event
		loadbar(current_event,num_of_events);
		chain.GetEntry(current_event);


		//reset electron cut bool
		electron_cuts = true;

		//electron cuts
		electron_cuts &= (id[0] == ELECTRON); //First particle is electron
		electron_cuts &= (gpart > 0); //Number of good particles is greater than 0
		electron_cuts &= (stat[0] > 0); //First Particle hit stat
		electron_cuts &= ((int)q[0] == -1); //First particle is negative Q
		electron_cuts &= (sc[0] > 0); //First Particle hit sc
		electron_cuts &= (dc[0] > 0); // ``` ``` ``` dc
		electron_cuts &= (ec[0] > 0); // ``` ``` ``` ec
		electron_cuts &= (dc_stat[dc[0]-1] > 0);

		if(electron_cuts){
			//Setup scattered electron 4 vector
			e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
			//Set the vertex time (time of electron hit) 
			delta_t_cut();

			if(first_run){	
				W = W_calc(e_mu, e_mu_prime);
				Q2 = Q2_calc(e_mu, e_mu_prime);
			}
			WvsQ2_Fill(e_mu_prime.E(),W,Q2,xb_calc(Q2, e_mu_prime.E()));
			num_of_proton = num_of_pis = 0;

			for(int part_num = 1; part_num < gpart; part_num++){
				if (p[part_num] == 0) continue;
				if (id[part_num] == PIP) Fill_pion_WQ2(W,Q2);
				if (id[part_num] == PROTON) Fill_proton_WQ2(W,Q2);

				if (id[part_num] == PIP) Fill_Pi_ID_P(p[part_num],b[part_num]);
				if (id[part_num] == PROTON) Fill_proton_ID_P(p[part_num],b[part_num]);

				if (id[part_num] == PIP || id[part_num] == PROTON) Fill_proton_Pi_ID_P(p[part_num],b[part_num]);

				Fill_Mass(m[part_num]);
				Particle3.SetXYZ(p[part_num]*cx[part_num],p[part_num]*cy[part_num],p[part_num]*cz[part_num]);
				/////Particle4.SetVectM(Particle3, Get_Mass(id[part_num]));
				Particle4.SetVectM(Particle3, m[part_num]);
				MomVsBeta_Fill(Particle4.E(),p[part_num],b[part_num]);
				if (q[part_num] == 1){
					MomVsBeta_Fill_pos(p[part_num],b[part_num]);
				} else if(q[part_num] == -1) {
					MomVsBeta_Fill_neg(p[part_num],b[part_num]);
				}
				if(id[part_num] == PROTON) num_of_proton++;
				if(id[part_num] == PIP){
					num_of_pis++;
					TLorentzVector gamma_mu = (e_mu - e_mu_prime);
					MissingMassNeutron.MissingMassPxPyPz(p[part_num]*cx[part_num],p[part_num]*cy[part_num],p[part_num]*cz[part_num]);
					MissingMassNeutron = MissingMassNeutron.missing_mass(gamma_mu);
				}
			}
	
			if(num_of_pis == 1) {
				Fill_Missing_Mass(MissingMassNeutron.mass);
				Fill_Missing_Mass_square(Square(MissingMassNeutron.mass));
				Fill_single_pi_WQ2(W,Q2);
			}
			if(num_of_proton == 1) Fill_single_proton_WQ2(W,Q2);
		}
	}

	// Start of cuts
	Cuts MissingMassNeutron_cut;
	double fit_range_min = 0.88;
	double fit_range_max = 1.0;
	MissingMassNeutron_cut.FitGaus(Missing_Mass,fit_range_min,fit_range_max);

	cut_outputs << "MM_N";
	cut_outputs << "," << MissingMassNeutron_cut.mean;
	cut_outputs << "," << MissingMassNeutron_cut.sigma << endl;

	Cuts MissingMassSquare_cut;
	fit_range_min = 0.5;
	fit_range_max = 1.1;
	MissingMassSquare_cut.FitGaus(Missing_Mass_square,fit_range_min,fit_range_max);

	cut_outputs << "MM_N_2";
	cut_outputs << "," << MissingMassSquare_cut.mean;
	cut_outputs << "," << MissingMassSquare_cut.sigma << endl;



	//
	//end stuff
	chain.Reset();						// delete Tree object

	RootOutputFile->cd();
	TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
	WvsQ2_folder->cd();
	WvsQ2_Write();

	TDirectory *MomVsBeta_folder = RootOutputFile->mkdir("Momentum vs beta");
	MomVsBeta_folder->cd();
	MomVsBeta_Write();

	//Missing Mass Write
	TDirectory *MissMass = RootOutputFile->mkdir("Missing_Mass");
	MissMass->cd();
	Write_Missing_Mass();

	//Missing Mass Write
	TDirectory *DeltaT = RootOutputFile->mkdir("Delta_T");
	DeltaT->cd();
	delta_t_Write();


	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
	cut_outputs.close();

}
#endif
