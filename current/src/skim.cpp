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
  data = std::make_shared<Branches>(chain);
  RootOutputFile = new TFile(fout.c_str(), "RECREATE");
  if (getenv("BEAM_E") != NULL) {
    BEAM_ENERGY = atof(getenv("BEAM_E"));
    std::cout << RED << "Beam energy set to: " << BEAM_ENERGY << DEF << std::endl;
  }

  e_mu.SetPxPyPzE(0.0, 0.0, sqrt((BEAM_ENERGY * BEAM_ENERGY) - (MASS_E * MASS_E)), BEAM_ENERGY);
  MM_neutron = new MissingMass(MASS_P, 0.0);
}
Skim::~Skim() {}

void Skim::Basic() {
  int num_of_events;
  bool electron_cuts, mm_cut;
  int bad, num_pi;
  std::cout << BLUE << "Skim file " << GREEN << fout << DEF << std::endl;

  num_of_events = (int)chain->GetEntries();
  TTree *skim = chain->CloneTree(0);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    num_pi = 0;
    bad = 0;
    chain->GetEntry(current_event);
    Cuts *check = new Cuts();

    // Not Used
    check->Set_electron_id(data->id(0));  // First particle is electron
    // electron cuts
    check->Set_charge(data->q(0));
    check->Set_ec_cut(data->ec(0) > 0);  // ``` ``` ``` ec
    check->Set_gpart(data->gpart());     // Number of good particles is greater than 0
    check->Set_cc_cut(data->cc(0) > 0);
    check->Set_stat_cut(data->stat(0) > 0);  // First Particle hit stat
    check->Set_sc_cut(data->sc(0) > 0);
    check->Set_dc_cut(data->dc(0) > 0);
    check->Set_dc_stat_cut(data->dc_stat(data->dc(0) - 1) > 0);

    Delta_T *dt = new Delta_T(data->sc_t(data->sc(0) - 1), data->sc_r(data->sc(0) - 1));
    std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, data.get());
    delete dt;

    for (int part_num = 1; part_num < data->gpart(); part_num++)
      if (data->q(part_num) == POSITIVE && check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num)))
        num_pi++;
      else
        bad++;

    if (check->isElecctron() && num_pi == 1 && bad == 0) skim->Fill();  // Fill the banks after the skim

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
  std::cout << BLUE << "Strict Skim file " << GREEN << fout << DEF << std::endl;

  num_of_events = (int)chain->GetEntries();
  TTree *skim = chain->CloneTree(0);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (data->gpart() > 5) continue;
    auto check = std::make_shared<Cuts>();
    double theta = 0.0;
    double phi = 0.0;
    double sector = 0.0;

    check->Set_electron_id(data->id(0));  // First particle is electron
    // electron cuts
    check->Set_charge(data->q(0));
    check->Set_ec_cut(data->ec(0) > 0);  // ``` ``` ``` ec
    check->Set_gpart(data->gpart());     // Number of good particles is greater than 0
    check->Set_cc_cut(data->cc(0) > 0);
    check->Set_stat_cut(data->stat(0) > 0);  // First Particle hit stat
    check->Set_sc_cut(data->sc(0) > 0);
    check->Set_dc_cut(data->dc(0) > 0);
    check->Set_dc_stat_cut(data->dc_stat(data->dc(0) - 1) > 0);

    theta = physics::theta_calc(data->cz(0));
    phi = physics::phi_calc(data->cx(0), data->cy(0));
    sector = physics::get_sector(phi);
    check->Set_elec_fid(theta, phi, sector);

    // isStrictElecctron
    check->Set_p(data->p(0));
    check->Set_Sf(data->etot(data->ec(0) - 1) / data->p(0));
    check->Set_num_phe(data->nphe(data->cc(0) - 1));
    check->Set_BeamPosition(data->dc_vx(data->dc(0) - 1), data->dc_vy(data->dc(0) - 1), data->dc_vz(data->dc(0) - 1));

    e_mu_prime_3.SetXYZ(data->p(0) * data->cx(0), data->p(0) * data->cy(0), data->p(0) * data->cz(0));
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
    double W = physics::W_calc(e_mu, e_mu_prime);
    double Q2 = physics::Q2_calc(e_mu, e_mu_prime);

    bool cuts = true;

    cuts &= check->Beam_cut();
    cuts &= check->isElecctron();
    cuts &= (Q2 >= 1.0);
    cuts &= (W >= 0.8);
    cuts &= (W <= 2.0);

    if (cuts) {
      skim->Fill();  // Fill the banks after the skim
    }
  }
  chain->Reset();  // delete Tree object
  delete chain;
  RootOutputFile->cd();
  RootOutputFile->Write();
  RootOutputFile->Close();
  delete RootOutputFile;
}

void Skim::Final() {
  int num_of_events;
  bool electron_cuts, mm_cut;
  bool cuts;
  int total_events;
  int num_of_pips, num_of_pims, num_of_proton, PID;
  double theta;
  double phi;
  int sector;
  MissingMass *MM_neutron = new MissingMass(MASS_P, 0.0);
  std::cout << BLUE << "Skim file " << GREEN << fout << DEF << std::endl;

  num_of_events = (int)chain->GetEntries();
  TTree *skim = chain->CloneTree(0);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (data->gpart() > 4) continue;
    Cuts *check = new Cuts();

    // electron cuts
    check->Set_charge(data->q(0));
    check->Set_ec_cut(data->ec(0) > 0);   // ``` ``` ``` ec
    check->Set_electron_id(data->id(0));  // First particle is electron
    check->Set_gpart(data->gpart());      // Number of good particles is greater than 0
    check->Set_cc_cut(data->cc(0) > 0);
    check->Set_stat_cut(data->stat(0) > 0);  // First Particle hit stat
    check->Set_sc_cut(data->sc(0) > 0);
    check->Set_dc_cut(data->dc(0) > 0);
    check->Set_dc_stat_cut(data->dc_stat(data->dc(0) - 1) > 0);
    check->Set_p(data->p(0));
    check->Set_Sf(data->etot(data->ec(0) - 1) / data->p(0));
    check->Set_num_phe(data->nphe(data->cc(0) - 1));
    // Beam position cut
    check->Set_BeamPosition(data->dc_vx(data->dc(0) - 1), data->dc_vy(data->dc(0) - 1), data->dc_vz(data->dc(0) - 1));

    theta = physics::theta_calc(data->cz(0));
    phi = physics::phi_calc(data->cx(0), data->cy(0));
    sector = physics::get_sector(phi);
    check->Set_elec_fid(theta, phi, sector);

    if (!check->isStrictElecctron()) {
      delete check;
      continue;
    }
    // Setup scattered electron 4 vector
    TLorentzVector e_mu_prime = physics::fourVec(data->p(0), data->cx(0), data->cy(0), data->cz(0), MASS_E);

    Delta_T *dt = new Delta_T(data->sc_t(data->sc(0) - 1), data->sc_r(data->sc(0) - 1));
    std::vector<double> dt_proton = dt->delta_t_array(MASS_P, data.get());
    std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, data.get());
    delete dt;

    // W = physics::W_calc(*e_mu, e_mu_prime);
    // Q2 = physics::Q2_calc(*e_mu, e_mu_prime);

    TLorentzVector gamma_mu = (e_mu - e_mu_prime);

    num_of_proton = num_of_pips = num_of_pims = 0;
    for (int part_num = 1; part_num < data->gpart(); part_num++) {
      PID = -99;
      if (data->q(part_num) == POSITIVE) {
        if (check->dt_P_cut(dt_proton.at(part_num), data->p(part_num)))
          PID = PROTON;
        else if (check->dt_Pip_cut(dt_proton.at(part_num), data->p(part_num)))
          PID = PIP;
      } else if (data->q(part_num) == NEGATIVE) {
        if (check->dt_Pip_cut(dt_proton.at(part_num), data->p(part_num))) PID = PIM;
      }

      TLorentzVector particle =
          physics::fourVec(data->p(part_num), data->cx(part_num), data->cy(part_num), data->cz(part_num), PID);

      if (data->q(part_num) == POSITIVE) {
        if (check->dt_P_cut(dt_proton.at(part_num), data->p(part_num))) {
          num_of_proton++;
        } else if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
          num_of_pips++;

          MM_neutron->Set_4Vec(particle);
          MM_neutron->missing_mass(gamma_mu);
        }
      } else if (data->q(part_num) == NEGATIVE && check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
        num_of_pims++;
      }
    }

    bool mm_cut = true;
    mm_cut &= (MM_neutron->Get_MM() < 1.1);
    mm_cut &= (MM_neutron->Get_MM() > 0.8);

    if (num_of_pips == 1 && mm_cut && num_of_proton == 0 && data->gpart() < 3) {
      skim->Fill();
    }

    delete check;
  }
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
