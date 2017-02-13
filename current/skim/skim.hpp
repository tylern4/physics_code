/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler
 */
/*	University Of South Carolina
 */
/************************************************************************/

#ifndef SKIM_H_GUARD
#define SKIM_H_GUARD
#include "main.h"

void skim(char *fin, char *RootFile_output) {

  TFile *RootOutputFile;
  int number_cols = 0;
  char rootFile[500];
  int num_of_events, total_events, num_of_pis;
  bool electron_cuts, MM_cut, has_neutron;

  MissingMass *MM_neutron = new MissingMass();
  MM_neutron->Set_target_mass(MASS_P);
  MM_neutron->Set_target_PxPyPz(0);
  Delta_T *delta_t = new Delta_T();

  Float_t W, Q2, MM;
  std::vector<bool> is_proton, is_pip, is_electron, is_pim;
  std::vector<double> dt_proton, dt_pip;

  TVector3 e_mu_prime_3;
  TLorentzVector e_mu_prime;
  TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)), E1D_E0);

  TVector3 Particle3(0.0, 0.0, 0.0);
  TLorentzVector Particle4(0.0, 0.0, 0.0, 0.0);

  RootOutputFile = new TFile(RootFile_output, "RECREATE");

  TChain chain("h10");
  cout << blue << "Analyzing file " << green << fin << def << bgdef << endl;
  chain.AddFile(fin);

  getBranches(&chain);

  num_of_events = (int)chain.GetEntries();

  TTree *skim = chain.CloneTree(0);
  TBranch *W_branch = skim->Branch("W", &W);
  TBranch *Q2_branch = skim->Branch("Q2", &Q2);
  TBranch *MM_branch = skim->Branch("MM", &MM);

  TBranch *is_Electron = skim->Branch("is_electron", &is_electron);
  TBranch *is_Proton = skim->Branch("is_proton", &is_proton);
  TBranch *is_Pip = skim->Branch("is_pip", &is_pip);
  TBranch *is_Pim = skim->Branch("is_pim", &is_pim);

  TBranch *DeltaT_P_branch =
      skim->Branch("DeltaT_P", "vector<double>", &dt_proton);
  TBranch *DeltaT_Pip_branch =
      skim->Branch("DeltaT_Pip", "vector<double>", &dt_pip);
  TBranch *NumPI_branch = skim->Branch("NumPI", &num_of_pis);
  TBranch *Neutron_branch = skim->Branch("has_neutron", &has_neutron);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain.GetEntry(current_event);
    is_proton = std::vector<bool>(gpart, false);
    is_electron = std::vector<bool>(gpart, false);
    is_pip = std::vector<bool>(gpart, false);
    is_pim = std::vector<bool>(gpart, false);

    electron_cuts = true;
    // electron cuts
    electron_cuts &= (id[0] == ELECTRON); // First particle is electron
    electron_cuts &= (gpart > 0); // Number of good particles is greater than 0
    electron_cuts &= (stat[0] > 0);            // First Particle hit stat
    electron_cuts &= ((int)q[0] == -1);        // First particle is negative Q
    electron_cuts &= (sc[0] > 0);              // First Particle hit sc
    electron_cuts &= (dc[0] > 0);              // ``` ``` ``` d
    electron_cuts &= (ec[0] > 0);              // ``` ``` ``` ec
    electron_cuts &= (dc_stat[dc[0] - 1] > 0); //??

    e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

    dt_proton = delta_t->delta_t_array(MASS_P, gpart);
    dt_pip = delta_t->delta_t_array(MASS_PIP, gpart);

    for (int part_num = 1; part_num < gpart; part_num++) {
      num_of_pis = 0;
      if (dt_pip.at(part_num) >= Pip_Neg_fit(p[part_num]) &&
          dt_pip.at(part_num) <= Pip_Pos_fit(p[part_num])) {
        is_pip.at(part_num) = true;
        num_of_pis++;
        TLorentzVector gamma_mu = (e_mu - e_mu_prime);
        MM_neutron->Set_PxPyPz(p[part_num] * cx[part_num],
                               p[part_num] * cy[part_num],
                               p[part_num] * cz[part_num]);
        MM = MM_neutron->missing_mass(gamma_mu);
      }
      if (dt_proton.at(part_num) >= Proton_Neg_fit(p[part_num]) &&
          dt_proton.at(part_num) <= Proton_Pos_fit(p[part_num]) &&
          q[part_num] == 1) {

        is_proton.at(part_num) = true;
      }
      if (dt_pip.at(part_num) >= Pip_Neg_fit(p[part_num]) &&
          dt_pip.at(part_num) <= Pip_Pos_fit(p[part_num])) {
        is_pim.at(part_num) = true;
      }
    }

    has_neutron = between_mm(MM);

    if (electron_cuts && has_neutron) {
      W = W_calc(e_mu, e_mu_prime);
      Q2 = Q2_calc(e_mu, e_mu_prime);
      is_electron.at(0) = true;
      skim->Fill(); // Fill the banks after the skim
    }
  }
  //
  // end stuff
  chain.Reset(); // delete Tree object

  RootOutputFile->cd();
  RootOutputFile->Write();
  RootOutputFile->Close();
}

#endif
