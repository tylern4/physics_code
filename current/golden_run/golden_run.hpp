/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef GOLDEN_RUN_H
#define GOLDEN_RUN_H
#include "main.h"

// Finds the files with golden runs

void golden_run(char *fin, char *fout) {
  ifstream input(fin);
  ofstream golden_run(fout);
  golden_run << "run_num,file_num,num_of_events,total_q" << endl;
  string file_num, run_num, line;
  int n_evnt, num_of_events = 0;
  double total_q = 0.0, curr_q = 0.0, prev_q = 0.0, delta_q = 0.0;
  bool electron_cuts;
  TVector3 e_mu_prime_3, pi_3vec;
  TLorentzVector e_mu_prime, pi_4vec;
  TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)), E1D_E0);
  TLorentzVector ZERO(0.0, 0.0, 0.0, 0.0);
  int line_num = 0;
  const char *progress = "-\\|/";

  while (input.is_open() && input.good()) {
    getline(input, line);
    if (!line.empty()) {
      run_num = line.substr(line.find("/h10") + 6, 5);
      // cout << run_num << endl;
      file_num = line.substr(line.find("/h10") + 12, 2);
      // cout << file_num << endl;
      TChain chain("h10");
      chain.AddFile(line.c_str());
      getBranches(&chain);
      n_evnt = (int)chain.GetEntries();
      num_of_events = 0;
      total_q = 0.0;

      for (int current_event = 0; current_event < n_evnt; current_event++) {
        cout << BLUE << "\t[ " << progress[(current_event / 100) % 4] << " ]\r" << DEF << flush;
        chain.GetEntry(current_event);
        electron_cuts = true;
        // electron cuts
        electron_cuts &= ((int)id[0] == ELECTRON);  // First particle is electron
        electron_cuts &= ((int)gpart > 0);          // Number of good particles is greater than 0
        electron_cuts &= ((int)stat[0] > 0);        // First Particle hit stat
        electron_cuts &= ((int)q[0] == -1);         // First particle is negative Q
        electron_cuts &= ((int)sc[0] > 0);          // First Particle hit sc
        electron_cuts &= ((int)dc[0] > 0);          // ``` ``` ``` dc
        electron_cuts &= ((int)dc_stat[dc[0] - 1] > 0);

        if (electron_cuts) {
          e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
          e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
          curr_q = q_l;
          if (curr_q > 0.0) {
            if (curr_q > prev_q) {
              delta_q = curr_q - prev_q;
              total_q += delta_q;
            }
            prev_q = curr_q;
          }
          for (int part_num = 1; part_num < gpart; part_num++) {
            if (id[part_num] == PIP) {
              pi_3vec.SetXYZ(p[part_num] * cx[part_num], p[part_num] * cy[part_num],
                             p[part_num] * cz[part_num]);
              pi_4vec.SetVectM(pi_3vec, MASS_PIP);
            }
            if (e_mu_prime != ZERO && pi_4vec != ZERO) num_of_events++;
          }
        }
      }

      if (num_of_events != 0 && total_q != 0)
        golden_run << run_num << "," << file_num << "," << num_of_events << "," << total_q << endl;
      chain.Clear();
    }
  }
  input.close();
  golden_run.close();
}
#endif
