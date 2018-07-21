/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "skim.hpp"

Skim::Skim(std::string input) {
  fin = input;
  chain = new TChain("h10");
  chain->AddFile(fin.c_str());

  fout = fin.substr(0, fin.size() - 5) + "_skim.root";
  RootOutputFile = new TFile(fout.c_str(), "RECREATE");

  e_mu.SetPxPyPzE(0.0, 0.0, sqrt((E1D_E0 * E1D_E0) - (MASS_E * MASS_E)), E1D_E0);
}
Skim::~Skim() {}

void Skim::Process() {
  int num_of_events;
  bool electron_cuts;
  int num_proton, num_pi;
  std::cout << BLUE << "Analyzing file " << GREEN << fin << DEF << std::endl;
  getBranches(chain);
  num_of_events = (int)chain->GetEntries();
  TTree *skim = chain->CloneTree(0);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    num_proton = 0;
    num_pi = 0;

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
      electron_cuts &= sf_cut((etot[ec[0] - 1] / p[0]), p[0]);
      electron_cuts &= ((int)nphe[cc[0] - 1] > 30);
    }

    e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

    Delta_T *dt = new Delta_T(sc_t[sc[0] - 1], sc_r[sc[0] - 1]);
    std::vector<double> dt_proton = dt->delta_t_array(MASS_P, gpart);
    std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, gpart);
    for (auto proton_dt : dt_proton)
      if (abs(proton_dt) < 0.5) num_proton++;
    for (auto pi_dt : dt_pi)
      if (abs(pi_dt) < 0.5) num_pi++;

    if (electron_cuts && (num_pi == 1 || num_proton == 1)) {
      skim->Fill();  // Fill the banks after the skim
    }
    delete dt;
  }
  chain->Reset();  // delete Tree object
  delete chain;

  RootOutputFile->cd();
  RootOutputFile->Write();
  RootOutputFile->Close();
  delete RootOutputFile;
}

double Skim::sf_top_fit(double P) {
  double par[3] = {0.363901, -0.00992778, 5.84749e-06};
  double x[1] = {P};
  return func::ec_fit_func(x, par);
}
double Skim::sf_bot_fit(double P) {
  double par[3] = {0.103964, 0.0524214, -3.64355e-05};
  double x[1] = {P};
  return func::ec_fit_func(x, par);
}
bool Skim::sf_cut(double sf, double P) { return ((sf > sf_bot_fit(P)) && (sf < sf_top_fit(P))); }

double Skim::dt_P_bot_fit(double P) {
  double par[2] = {-1.509, 0.4172};
  double x[1] = {P};
  return func::dt_fit(x, par);
}
double Skim::dt_P_top_fit(double P) {
  double par[2] = {1.307, -0.3473};
  double x[1] = {P};
  return func::dt_fit(x, par);
}
bool Skim::dt_P_cut(double dt, double P) { return ((dt > dt_P_bot_fit(P)) && (dt < dt_P_top_fit(P))); }
double Skim::dt_Pip_bot_fit(double P) {
  double par[2] = {-0.9285, -0.04094};
  double x[1] = {P};
  return func::dt_fit(x, par);
}
double Skim::dt_Pip_top_fit(double P) {
  double par[2] = {0.9845, -0.05473};
  double x[1] = {P};
  return func::dt_fit(x, par);
}
bool Skim::dt_Pip_cut(double dt, double P) { return ((dt > dt_Pip_bot_fit(P)) && (dt < dt_Pip_top_fit(P))); }
