/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include <iostream>
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "../src/constants.hpp"
#include "../src/physics.hpp"
#include "../src/classes.hpp"
#include <TLorentzVector.h>

class h10_data {
public:
  std::vector<double> W_vec;
  std::vector<double> Q2_vec;
  std::vector<double> MM_vec;
  std::vector<std::vector<std::string>> particle_list;

  h10_data() {}

  ~h10_data() {}
};

class H10 {
private:
  double Square(double a) { return a * a; }

public:
  //  std::vector<double> W_vec;
  //  std::vector<double> Q2_vec;
  //  std::vector<double> MM_vec;
  //  std::vector<std::vector<std::string>> particle_list;

  H10() {}

  ~H10() {}
  int num_of_events;
  bool electron_cuts;

  void pass_chain(TTree &chain) {
    getBranches(&chain);
    int num_of_events = (int)chain.GetEntries();
    std::cout << num_of_events << std::endl;
  }

  h10_data dataHandeler(TTree &chain) {
    h10_data out;
    // From missing_mass.hpp :: missing_mass_calc()
    MissingMass *MM_neutron = new MissingMass();
    MM_neutron->Set_target_mass(MASS_P);
    MM_neutron->Set_target_PxPyPz(0);

    TCanvas *c1 = new TCanvas("c1", "c1", 100, 100);
    auto number_cols = 0;
    char rootFile[500];
    int num_of_events, total_events;
    bool cuts, electron_cuts;
    double e_E;

    int num_of_pis, num_of_proton;

    TVector3 e_mu_prime_3;
    TLorentzVector e_mu_prime;
    TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)),
                        E1D_E0);

    TVector3 Particle3(0.0, 0.0, 0.0);
    TLorentzVector Particle4(0.0, 0.0, 0.0, 0.0);
    double theta, phi;
    int sector;

    getBranches(&chain);

    num_of_events = (int)chain.GetEntries();

    for (int current_event = 0; current_event < num_of_events;
         current_event++) {
      std::vector<std::string> particle_vector;
      chain.GetEntry(current_event);

      // reset electron cut bool
      electron_cuts = true;
      // electron cuts
      electron_cuts &= (ec[0] > 0);              // ``` ``` ``` ec
      electron_cuts &= ((int)id[0] == ELECTRON); // First particle is electron
      electron_cuts &= (p[0] > MIN_P_CUT);       // Minimum Momentum cut
      electron_cuts &=
          ((int)gpart > 0); // Number of good particles is greater than 0
      // electron_cuts &= ((int)stat[0] > 0); // First Particle hit stat
      // std::cout << "stat " << (int)stat[0] << std::endl;
      electron_cuts &= ((int)q[0] == -1); // First particle is negative Q
      electron_cuts &= ((int)sc[0] > 0);  // First Particle hit sc
      electron_cuts &= ((int)dc[0] > 0);  // ``` ``` ``` dc
      electron_cuts &= ((int)dc_stat[dc[0] - 1] > 0);

      // Setup scattered electron 4 vector
      e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

      out.W_vec.push_back(W_calc(e_mu, e_mu_prime));
      out.Q2_vec.push_back(Q2_calc(e_mu, e_mu_prime));
      for (int part_num = 0; part_num < gpart; part_num++) {
        particle_vector.push_back(std::to_string(id[part_num]));
      }
      out.particle_list.push_back(particle_vector);
      /*
      num_of_proton = num_of_pis = 0;
        for (int part_num = 1; part_num < gpart; part_num++) {
          if (p[part_num] == 0)
            continue;

          Particle3.SetXYZ(p[part_num] * cx[part_num], p[part_num] *
        cy[part_num],
                           p[part_num] * cz[part_num]);
          Particle4.SetVectM(Particle3, Get_Mass(id[part_num]));

          if (id[part_num] == PIP) {
            num_of_pis++;
            TLorentzVector gamma_mu = (e_mu - e_mu_prime);

            MM_neutron->Set_PxPyPz(p[part_num] * cx[part_num],
                                   p[part_num] * cy[part_num],
                                   p[part_num] * cz[part_num]);
            out.MM_vec.push_back(MM_neutron->missing_mass(gamma_mu));
          } else {
            out.MM_vec.push_back(-999.0);
          }
        }
        */
    }

    chain.Reset();
    return out;
  }
};
