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
{
	TFile *RootOutputFile;
	int number_cols = 0;
	char rootFile[500];
	int num_of_events, total_events;
	bool cuts, electron_cuts;

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

{	//electron cuts
		electron_cuts &= (id[0] == ELECTRON); //First particle is electron
		electron_cuts &= (gpart > 0); //Number of good particles is greater than 0
		electron_cuts &= (stat[0] > 0); //First Particle hit stat
		electron_cuts &= ((int)q[0] == -1); //First particle is negative Q
		electron_cuts &= (sc[0] > 0); //First Particle hit sc
		electron_cuts &= (dc[0] > 0); // ``` ``` ``` dc
		electron_cuts &= (ec[0] > 0); // ``` ``` ``` ec
		electron_cuts &= (dc_stat[dc[0]-1] > 0);
}

		if(electron_cuts){
			//Setup scattered electron 4 vector
			e_mu_prime_3.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

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

	//
	//end stuff
	chain.Reset();						// delete Tree object

	RootOutputFile->cd();

	//Missing Mass Write
	TDirectory *MissMass = RootOutputFile->mkdir("Missing_Mass");
	MissMass->cd();
	Write_Missing_Mass();

	

	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list

}
#endif
