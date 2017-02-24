/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include <iostream>
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "src/physics.hpp"
#include "src/constants.hpp"
#include "src/branches.hpp"
#include <TLorentzVector.h>

class H10 {
private:
  double Square(double a) { return a * a; }

public:
  std::vector<double> W_vec;
  std::vector<double> Q2_vec;
  std::vector<float> p_vec;
  std::vector<float> b_vec;
  // W and Q^2
  int bins = 500;
  double w_min = 0;
  double w_max = 3.25;
  double q2_min = 0;
  double q2_max = 10;

  TH2D *WvsQ2_hist = new TH2D("WvsQ2_hist", "W vs Q^{2}", bins, w_min, w_max,
                              bins, q2_min, q2_max);

  H10() {}

  ~H10() {}
  int num_of_events;
  bool electron_cuts;

  void pass_chain(TTree &chain) {
    getBranches(&chain);
    int num_of_events = (int)chain.GetEntries();
    std::cout << num_of_events << std::endl;
  }

  void loop(TTree &chain) {
    TVector3 e_mu_prime_3;
    TLorentzVector e_mu_prime;
    TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)),
                        E1D_E0);

    getBranches(&chain);
    int num_of_events = (int)chain.GetEntries();
    for (int current_event = 0; current_event < num_of_events;
         current_event++) {
      chain.GetEntry(current_event);
      // reset electron cut bool
      electron_cuts = true;
      // electron cuts
      electron_cuts &= (ec[0] > 0); // ``` ``` ``` ec
      // if (electron_cuts) hists->EC_fill(etot[ec[0]-1],p[0]);
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
        W = W_calc(e_mu, e_mu_prime);
        W_vec.push_back(W);
        Q2 = Q2_calc(e_mu, e_mu_prime);
        Q2_vec.push_back(Q2);
        WvsQ2_hist->Fill(W, Q2);
      }
      for (int part_num = 1; part_num < gpart; part_num++) {
        if (p[part_num] == 0)
          continue;
        p_vec.push_back((float)p[part_num]);
        b_vec.push_back((float)b[part_num]);
      }
    }
  }
};
