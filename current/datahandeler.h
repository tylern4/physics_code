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

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;
	bool cuts, electron_cuts;
	ofstream cut_outputs;
	cut_outputs.open ("outputFiles/cut_outputs.csv");
	cut_outputs << "Cut,Mean,Sigma" << endl;

	//From delta_t.hpp :: vertex_time() delta_t()
	Delta_T delta_t;
	//From delta_t_cut.hpp :: delta_t_calc()
	D_time delta_time;
	double delta_t_P, delta_t_PIP, delta_t_ELECTRON;
	//From missing_mass.hpp :: missing_mass_calc()
	MissingMass MissingMassNeutron;

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
			delta_time.SetVertexTimes(sc_t[sc[0]-1],sc_r[sc[0]-1]);
			for(int part_num = 1; part_num < gpart; part_num++){
				Fill_Mass(m[part_num]);
				delta_time.SetTandP(sc_t[sc[part_num]-1],sc_r[sc[part_num]-1],p[part_num]);
				delta_time = delta_time.delta_t_calc();
				Particle3.SetXYZ(p[part_num]*cx[part_num],p[part_num]*cy[part_num],p[part_num]*cz[part_num]);
				Particle4.SetVectM(Particle3, Get_Mass(id[part_num]));
	
				if (Particle4.P() != 0) {  // && (int)q[part_num] == 1
					delta_t_Fill(Particle4.P(), delta_time.proton_time, 3);
					delta_t_Fill(Particle4.P(), delta_time.pip_time, 4);
					delta_t_Fill(Particle4.P(), delta_time.electron_time, 5);
	
					//If Pi+
					if(id[part_num] == PROTON) { //&& (int)q[part_num] == 1
						delta_t_Fill(Particle4.P(), delta_time.proton_time, 1);
						delta_t_Fill(Particle4.P(), delta_time.electron_time, 7);
					//If Proton	
					} else if (id[part_num] == PIP){ //&& (int)q[part_num] == 1
						delta_t_Fill(Particle4.P(), delta_time.pip_time, 2);
						delta_t_Fill(Particle4.P(), delta_time.electron_time, 6);
					} 
				}
			}

			for(int part_num = 1; part_num < gpart; part_num++){
				num_of_pis = 0;
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
			}
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
