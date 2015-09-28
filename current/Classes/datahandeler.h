/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef DATAHANDELER_H_GUARD
#define DATAHANDELER_H_GUARD
// Mashing together W vs Q2 and Delta T cuts into one file
// Saving the old files in a new folder to refer back to.
//
void dataHandeler(char *fin, char *RootFile_output){

	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;
	bool cuts;

	Delta_T dt;
	D_time dt_cuts;
	MissingMass MM;

	TVector3 e_mu_prime_3;
	TLorentzVector e_mu_prime;
	TLorentzVector e_mu(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0);

	TVector3 Particle3(0.0,0.0,0.0);
	TLorentzVector Particle4(0.0,0.0,0.0,0.0);

	RootOutputFile = new TFile(RootFile_output,"RECREATE"); 

	TChain chain("h10");
	cout << blue <<"Analyzing file " << green << fin << def << bgdef << endl;

	FILE *input_file = fopen(fin,"r");
	ofstream text_output;
	text_output.open("outputFiles/output.txt");

	if (input_file == NULL) perror ("Error opening file");

	while (1){
		number_cols = fscanf(input_file,"%s",rootFile);

		if (number_cols<0) break;
		chain.Add(rootFile);
	}

	getBranches(&chain);

	num_of_events = (int)chain.GetEntries();
///rewite from here down:
	///two loops:
	/// 1) calc and fill original hists
	/// 2) use calc from first loop and cuts to fill cut hists
		///Look at making std::vec for W, Q2, mm_cuts, delta_t_cuts
	for (int current_event = 0; current_event < num_of_events; current_event++) {
		loadbar(current_event,num_of_events);
		chain.GetEntry(current_event);

		if (id[0] == ELECTRON && gpart > 0 && stat[0] > 0 && (int)q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0){
			//Setup scattered electron 4 vector
			e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
			cuts = ( b[0] != 0 );

			if(cuts){
				//dt_cuts = delta_t_calc();
				dt_cuts.SetVertexTimes(sc_t[sc[0]-1],sc_r[sc[0]-1]);
				WvsQ2(e_mu,e_mu_prime);
				for(int part_num = 1; part_num < gpart; part_num++){
					dt_cuts.SetTandP(sc_t[sc[part_num]-1],sc_r[sc[part_num]-1],p[part_num]);
					dt_cuts = dt_cuts.delta_t_calc();


					//if (event_number == 0 && id[0] == ELECTRON && gpart > 0 && stat[0] > 0 && (int)q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0) {
					//	delta_t_Fill(Particle4.P(), delta_t_PIP, 8);
					//}
					//Fill delta T hists here now
					/*

						delta_t_Fill(Particle4.P(), delta_t_PIP, 8);
				
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

					*/
					if(id[part_num] == PIP){
						TLorentzVector gamma_mu = (e_mu - e_mu_prime);
						MissingMass MM;
						MM.MissingMassPxPyPz(p[part_num]*cx[part_num],p[part_num]*cy[part_num],p[0]*cz[part_num]);
						MM = MM.missing_mass(gamma_mu);
						Fill_Missing_Mass(MM.mass);
					}
				}
				//text_output << current_event <<"\t"<< W_calc(e_mu, e_mu_prime) <<"\t"<< Q2_calc(e_mu, e_mu_prime) << endl;
			}
		}
	}

	// Start of cuts
	Cuts mm_cut;
	double *par;
	//par[0] = MASS_N;
	mm_cut.CutFit(Missing_Mass,0.9,0.98, par);
	cout << "Mean: " << mm_cut.mean << "\tSigma: " << mm_cut.sigma << endl; 
	cout << blue << mm_cut.mean + 3 * mm_cut.sigma <<"\t" << mm_cut.mean - 3 * mm_cut.sigma << def << endl;
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
	text_output.close();
}
#endif
