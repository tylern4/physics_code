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
  std::cout << BLUE << "Basic Skim " << GREEN << fout << DEF << std::endl;

  int num_of_events = (int)chain->GetEntries();
  TTree *skim = chain->CloneTree(0);
  Branches *data = new Branches(chain);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    auto check = std::make_shared<Cuts>(data);
    if (check->isElecctron()) skim->Fill();
  }

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
    cuts &= check->isStrictElecctron();
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
  std::cout << BLUE << "Event selection of file " << GREEN << fout << DEF << std::endl;

  int num_of_events = (int)chain->GetEntries();
  TTree *skim = chain->CloneTree(0);
  Branches *data = new Branches(chain);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (data->ec_ei(0) < 0.01 || data->ec_eo(0) < 0.01) continue;
    auto check = std::make_shared<Cuts>(data);
    if (!check->isElecctron()) continue;

    auto event = std::make_shared<Reaction>();
    event->SetElec(data->p(0), data->cx(0), data->cy(0), data->cz(0));

    if (check->isElecctron()) {
      auto dt = std::make_shared<Delta_T>(data->sc_t(0), data->sc_r(0));
      std::vector<double> dt_proton = dt->delta_t_array(MASS_P, data);
      std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, data);

      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        if (data->q(part_num) == POSITIVE) {
          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
            event->SetPip(data->p(part_num), data->cx(part_num), data->cy(part_num), data->cz(part_num));
          } else if (check->dt_P_cut(dt_proton.at(part_num), data->p(part_num))) {
            event->SetProton(data->p(part_num), data->cx(part_num), data->cy(part_num), data->cz(part_num));
          }
        } else if (data->q(part_num) == NEGATIVE) {
          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
            event->SetPim(data->p(part_num), data->cx(part_num), data->cy(part_num), data->cz(part_num));
          }
        }
      }

      bool mm_cut = true;
      mm_cut &= (event->MM() < 1.00091);
      mm_cut &= (event->MM() > 0.911698);
      if (event->SinglePip() && mm_cut) skim->Fill();
    }
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
