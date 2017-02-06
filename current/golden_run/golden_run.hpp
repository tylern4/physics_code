/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DATAHANDELER_H_GUARD
#define DATAHANDELER_H_GUARD
#include "main.h"
// Mashing together W vs Q2 and Delta T cuts into one file
// Saving the old files in a new folder to refer back to.
//
void golden_run(char *fin, char *fout) {
  ifstream input(fin);
  ofstream golden_run(fout);
  golden_run << "run_num,file_num,num_of_events,total_q" << endl;
  string file_num, run_num, line;
  int n_evnt, num_of_events = 0;
  double total_q = 0.0, curr_q = 0.0, prev_q = 0.0, delta_q = 0.0;
  bool electron_cuts;
  TVector3 e_mu_prime_3;
  TLorentzVector e_mu_prime;
  TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)), E1D_E0);
  TLorentzVector ZERO(0.0, 0.0, 0.0, 0.0);

  while (input.is_open() && input.good()) {
    getline(input, line);
    if (!line.empty()) {
      run_num = line.substr(line.find("/r2") + 2, 5);
      file_num = line.substr(line.find("_") + 1, 2);

      TChain chain("h10");
      chain.AddFile(line.c_str());
      getBranches(&chain);
      n_evnt = (int)chain.GetEntries();
      num_of_events = 0;
      total_q = 0.0;

      for (int current_event = 0; current_event < n_evnt; current_event++) {
        loadbar(current_event + 1, n_evnt);
        chain.GetEntry(current_event);
        electron_cuts = true;
        // electron cuts
        electron_cuts &= ((int)id[0] == ELECTRON); // First particle is electron
        electron_cuts &=
            ((int)gpart > 0); // Number of good particles is greater than 0
        electron_cuts &= ((int)stat[0] > 0); // First Particle hit stat
        electron_cuts &= ((int)q[0] == -1);  // First particle is negative Q
        electron_cuts &= ((int)sc[0] > 0);   // First Particle hit sc
        electron_cuts &= ((int)dc[0] > 0);   // ``` ``` ``` dc
        electron_cuts &= ((int)dc_stat[dc[0] - 1] > 0);

        if (electron_cuts) {
          e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
          e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
          curr_q = q_l;
          if (curr_q > 0.0 && curr_q > prev_q) {
            delta_q = curr_q - prev_q;
            total_q += delta_q;
          }
          prev_q = curr_q;
          for (int part_num = 1; part_num < gpart; part_num++) {
            if (e_mu_prime != ZERO)
              num_of_events++;
          }
        }
      }

      if (num_of_events != 0 && total_q != 0)
        golden_run << run_num << "," << file_num << "," << num_of_events << ","
                   << total_q << endl;
      chain.Clear();
    }
  }
  input.close();
  golden_run.close();
}
#endif
