/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#ifndef SKIM_H_GUARD
#define SKIM_H_GUARD
#include "main.h"

bool is_proton_dt(double dt, double p) {
  bool neg = (dt > (p * 0.4172 - 1.509));
  bool pos = (dt < (p * -0.3473 + 1.037));
  return (pos && neg);
}

bool is_pip_dt(double dt, double p) {
  bool pos = (dt < (p * -0.0118 + 0.9657));
  bool neg = (dt > (p * -0.1204 - 0.8795));
  return (pos && neg);
}

void skim(char *fin, char *RootFile_output) {
  TFile *RootOutputFile;
  int number_cols = 0;
  char rootFile[500];
  int num_of_events, total_events, num_of_pis;
  bool electron_cuts, MM_cut, has_neutron;

  MissingMass *MM_neutron = new MissingMass();
  MM_neutron->Set_target_mass(MASS_P);
  MM_neutron->Set_target_PxPyPz(0);

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
  cout << BLUE << "Analyzing file " << GREEN << fin << DEF << endl;
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

  TBranch *DeltaT_P_branch = skim->Branch("DeltaT_P", "vector<double>", &dt_proton);
  TBranch *DeltaT_Pip_branch = skim->Branch("DeltaT_Pip", "vector<double>", &dt_pip);
  TBranch *NumPI_branch = skim->Branch("NumPI", &num_of_pis);
  TBranch *Neutron_branch = skim->Branch("has_neutron", &has_neutron);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain.GetEntry(current_event);
    /*********/
    int n_pip = 0;
    int n_other = 0;
    for (int x = 0; x < gpart; x++)
      if (id[x] == PIP)
        n_pip++;
      else
        n_other++;

    if (n_pip != 1 || n_other > 1) continue;
    /*********/
    is_proton = std::vector<bool>(gpart, false);
    is_electron = std::vector<bool>(gpart, false);
    is_pip = std::vector<bool>(gpart, false);
    is_pim = std::vector<bool>(gpart, false);

    electron_cuts = true;
    // electron cuts
    electron_cuts &= (id[0] == ELECTRON);       // First particle is electron
    electron_cuts &= (gpart > 0);               // Number of good particles is greater than 0
    electron_cuts &= (stat[0] > 0);             // First Particle hit stat
    electron_cuts &= ((int)q[0] == -1);         // First particle is negative Q
    electron_cuts &= (sc[0] > 0);               // First Particle hit sc
    electron_cuts &= (dc[0] > 0);               // ``` ``` ``` d
    electron_cuts &= (ec[0] > 0);               // ``` ``` ``` ec
    electron_cuts &= (dc_stat[dc[0] - 1] > 0);  //??
    if (electron_cuts) {
      electron_cuts &= (etot[ec[0] - 1] / p[0]) < 0.4;
      electron_cuts &= (etot[ec[0] - 1] / p[0]) > 0.2;
    }

    e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
    Delta_T *delta_t = new Delta_T();
    dt_proton = delta_t->delta_t_array(MASS_P, gpart);
    dt_pip = delta_t->delta_t_array(MASS_PIP, gpart);

    for (int part_num = 1; part_num < gpart; part_num++) {
      num_of_pis = 0;
      // Hard code of values to use for cut
      if (q[part_num] == 1) {
        is_proton.at(part_num) = is_proton_dt(dt_proton.at(part_num), p[part_num]);
        if (is_pip_dt(dt_pip.at(part_num), p[part_num])) {
          is_pip.at(part_num) = true;
          num_of_pis++;
          TLorentzVector gamma_mu = (e_mu - e_mu_prime);
          MM_neutron->Set_PxPyPz(p[part_num] * cx[part_num], p[part_num] * cy[part_num], p[part_num] * cz[part_num]);
          MM_neutron->missing_mass(gamma_mu);
          MM = MM_neutron->Get_MM();
        }
      }
    }

    has_neutron = false;  // between_mm(MM);

    // Simple cut for missing mass
    if (num_of_pis == 1) has_neutron = (MM > 0.5 && MM < 1.5);

    if (electron_cuts && has_neutron) {
      W = physics::W_calc(e_mu, e_mu_prime);
      Q2 = physics::Q2_calc(e_mu, e_mu_prime);
      is_electron.at(0) = true;
      skim->Fill();  // Fill the banks after the skim
    }
  }
  //
  // end stuff
  chain.Reset();  // delete Tree object

  RootOutputFile->cd();
  RootOutputFile->Write();
  RootOutputFile->Close();
}

#endif
