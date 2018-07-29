/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "skim.hpp"

Skim::Skim(std::vector<std::string> input, std::string output) {
  fin = input;
  fout = output;
  chain = new TChain("h10");
  for (auto f : fin) chain->AddFile(f.c_str());

  RootOutputFile = new TFile(fout.c_str(), "RECREATE");

  e_mu.SetPxPyPzE(0.0, 0.0, sqrt((E1D_E0 * E1D_E0) - (MASS_E * MASS_E)), E1D_E0);
  MM_neutron = new MissingMass(MASS_P, 0.0);
}
Skim::~Skim() {}

void Skim::Basic() {
  int num_of_events;
  bool electron_cuts, mm_cut;
  int num_proton, num_pi;
  std::cout << BLUE << "Skim file " << GREEN << fout << DEF << std::endl;
  getBranches(chain);
  num_of_events = (int)chain->GetEntries();
  TTree *skim = chain->CloneTree(0);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    Cuts *check = new Cuts();

    // Not Used
    check->Set_electron_id((int)id[0]);  // First particle is electron
    // electron cuts
    check->Set_charge((int)q[0]);
    check->Set_ec_cut(ec[0] > 0);  // ``` ``` ``` ec
    check->Set_gpart((int)gpart);  // Number of good particles is greater than 0
    check->Set_cc_cut((int)cc[0] > 0);
    check->Set_stat_cut((int)stat[0] > 0);  // First Particle hit stat
    check->Set_sc_cut((int)sc[0] > 0);
    check->Set_dc_cut((int)dc[0] > 0);
    check->Set_dc_stat_cut((int)dc_stat[dc[0] - 1] > 0);

    if (check->isElecctron()) {
      skim->Fill();  // Fill the banks after the skim
    }
    delete check;
  }
  chain->Reset();  // delete Tree object
  delete chain;

  RootOutputFile->cd();
  RootOutputFile->Write();
  RootOutputFile->Close();
  delete RootOutputFile;
}

void Skim::Strict() {
  int num_of_events;
  bool electron_cuts, mm_cut;
  int num_proton, num_pip;
  int num_PPIP = 0;
  std::cout << BLUE << "Skim file " << GREEN << fout << DEF << std::endl;
  getBranches(chain);
  num_of_events = (int)chain->GetEntries();
  TTree *skim = chain->CloneTree(0);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (gpart > 4) continue;
    Cuts *check = new Cuts();

    num_proton = 0;
    num_pip = 0;
    // num_PPIP = 0;
    // Not Used
    check->Set_electron_id((int)id[0]);  // First particle is electron
    // electron cuts
    check->Set_charge((int)q[0]);
    check->Set_ec_cut(ec[0] > 0);  // ``` ``` ``` ec
    check->Set_gpart((int)gpart);  // Number of good particles is greater than 0
    check->Set_cc_cut((int)cc[0] > 0);
    check->Set_stat_cut((int)stat[0] > 0);  // First Particle hit stat
    check->Set_sc_cut((int)sc[0] > 0);
    check->Set_dc_cut((int)dc[0] > 0);
    check->Set_dc_stat_cut((int)dc_stat[dc[0] - 1] > 0);

    // isStrictElecctron
    check->Set_p((double)p[0]);
    check->Set_Sf((double)etot[ec[0] - 1] / (double)p[0]);
    check->Set_num_phe((int)nphe[cc[0] - 1]);
    check->Set_BeamPosition((double)dc_vx[dc[0] - 1], (double)dc_vy[dc[0] - 1], (double)dc_vz[dc[0] - 1]);

    if (!check->isElecctron()) continue;
    e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
    TLorentzVector gamma_mu = (e_mu - e_mu_prime);

    Delta_T *dt = new Delta_T(sc_t[sc[0] - 1], sc_r[sc[0] - 1]);
    std::vector<double> dt_proton = dt->delta_t_array(MASS_P, gpart);
    std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, gpart);

    mm_cut = true;
    for (int part_num = 0; part_num < gpart; part_num++) {
      if (q[part_num] == NEGATIVE) continue;
      particle_3.SetXYZ(p[part_num] * cx[part_num], p[part_num] * cy[part_num], p[part_num] * cz[part_num]);
      if (check->dt_P_cut(dt_proton[part_num], p[part_num]) && check->dt_Pip_cut(dt_pi[part_num], p[part_num]))
        num_PPIP++;

      if (check->dt_P_cut(dt_proton[part_num], p[part_num])) {
        num_proton++;
        particle.SetVectM(particle_3, MASS_P);
      }
      if (check->dt_Pip_cut(dt_pi[part_num], p[part_num])) {
        num_pip++;
        particle.SetVectM(particle_3, MASS_PIP);
        MM_neutron->Set_4Vec(particle);
        MM_neutron->missing_mass(gamma_mu);
      }
    }
    /* TODO:
    Here's the problem:
      What if I have two pions?
      What if I have more than 3 particles?
    */
    mm_cut &= (MM_neutron->Get_MM() < 1.5);
    mm_cut &= (MM_neutron->Get_MM() > 0.5);

    if (check->isElecctron() && num_pip >= 1) {
      skim->Fill();  // Fill the banks after the skim}
    }
    // delete dt;
    delete check;
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
