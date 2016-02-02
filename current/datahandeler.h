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
void dataHandeler(char *fin, char *RootFile_output){

	std::vector<double> W_vec, Q2_vec, MM_vec;

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;
	bool cuts, electron_cuts;

	//From delta_t.hpp :: vertex_time() delta_t()
	Delta_T dt;
	//From delta_t_cut.hpp :: delta_t_calc()
	D_time dt_cuts;
	//From missing_mass.hpp :: missing_mass_calc()
	MissingMass MM;
	int num_of_pis;

	TVector3 e_mu_prime_3;
	TLorentzVector e_mu_prime;
	TLorentzVector e_mu(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0);

	TVector3 Particle3(0.0,0.0,0.0);
	TLorentzVector Particle4(0.0,0.0,0.0,0.0);
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

	num_of_events = (int)chain.GetEntries();
///rewite from here down:
	///two loops:
	/// 1) calc and fill original hists
	/// 2) use calc from first loop and cuts to fill cut hists
		///Look at making std::vec for W, Q2, mm_cuts, delta_t_cuts
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

		if (electron_cuts){
			//Setup scattered electron 4 vector
			e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

			cuts = ( b[0] != 0 ); //beta doesn't equal zero (speed doesn't equal c)

			//Set the vertex time (time of electron hit) 
			dt_cuts.SetVertexTimes(sc_t[sc[0]-1],sc_r[sc[0]-1]);

			if(cuts){
				W_vec.push_back(W_calc(e_mu, e_mu_prime));
				Q2_vec.push_back(Q2_calc(e_mu, e_mu_prime));

				//Part of WvsQ2.hpp 
				WvsQ2(e_mu,e_mu_prime);
				for(int part_num = 1; part_num < gpart; part_num++){

					dt_cuts.SetTandP(sc_t[sc[part_num]-1],sc_r[sc[part_num]-1],p[part_num]);
					dt_cuts = dt_cuts.delta_t_calc();
					Particle3.SetXYZ(p[part_num]*cx[part_num],p[part_num]*cy[part_num],p[part_num]*cz[part_num]);
					Particle4.SetVectM(Particle3, Get_Mass(id[part_num]));

					if (Particle4.P() != 0 && (int)q[part_num] == 1) {
						delta_t_Fill(Particle4.P(), dt_cuts.proton_time, 3);
						delta_t_Fill(Particle4.P(), dt_cuts.pip_time, 4);
						delta_t_Fill(Particle4.P(), dt_cuts.electron_time, 5); 
			
						//If Pi+
						if(id[part_num] == PROTON && (int)q[part_num] == 1) {
							delta_t_Fill(Particle4.P(), dt_cuts.proton_time, 1);
							delta_t_Fill(Particle4.P(), dt_cuts.electron_time, 7);
						//If Proton	
						} else if (id[part_num] == PIP && (int)q[part_num] == 1){
							delta_t_Fill(Particle4.P(), dt_cuts.pip_time, 2);
							delta_t_Fill(Particle4.P(), dt_cuts.electron_time, 6);
						} 
					}

					num_of_pis = 0;
					if(id[part_num] == PIP){
						num_of_pis++;
						TLorentzVector gamma_mu = (e_mu - e_mu_prime);
						MM.MissingMassPxPyPz(p[part_num]*cx[part_num],p[part_num]*cy[part_num],p[part_num]*cz[part_num]);
						MM = MM.missing_mass(gamma_mu);
					}
				}
				if(num_of_pis == 1) {
					Fill_Missing_Mass(MM.mass);
					MM_vec.push_back(MM.mass);
				}
			}
		}
	}

	// Start of cuts
	Cuts mm_cut;
	double *par;
	//par[0] = MASS_N;
	//mm_cut.CutFit(Missing_Mass,0.9,1.0, par);

/*
	while (1){
		number_cols = fscanf(input_file,"%s",rootFile);

		if (number_cols<0) break;
		chain.Add(rootFile);
	}

	getBranches(&chain);

	num_of_events = (int)chain.GetEntries();

	for (int current_event = 0; current_event < num_of_events; current_event++) {
		loadbar(current_event,num_of_events);
		chain.GetEntry(current_event);

		if (id[0] == ELECTRON && gpart > 0 && stat[0] > 0 && (int)q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0){
			//Setup scattered electron 4 vector
			e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
			cuts = ( b[0] != 0 );

			if(cuts){
				delta_t_cut();
				WvsQ2(e_mu,e_mu_prime);
				TLorentzVector gamma_mu = (e_mu - e_mu_prime);
				missing_mass(gamma_mu);
				text_output << current_event <<"\t"<< W_calc(e_mu, e_mu_prime) <<"\t"<< Q2_calc(e_mu, e_mu_prime) << endl;
			}
		}
	}
*/

	//
	//end stuff
	chain.Reset();						// delete Tree object

	RootOutputFile->cd();
	//write stuff
	//Can probble make this into a single write function in main.h or write functions in histo.h
	//WvsQ2
	TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W_vs_Q2");
	WvsQ2_folder->cd();
	WvsQ2_Write();

	TDirectory *MomVsBeta_folder = RootOutputFile->mkdir("Momentum_vs_beta");
	MomVsBeta_folder->cd();
	MomVsBeta_Write();

	//Delta_t Write
	TDirectory *delta_t_folder = RootOutputFile->mkdir("Delta_t");
	delta_t_folder->cd();
	delta_t_Write();

	//Extra Write
	TDirectory *extras_folder = RootOutputFile->mkdir("Missing_Mass");
	extras_folder->cd();
	Write_Missing_Mass();


	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list

}
#endif
