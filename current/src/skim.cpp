/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "skim.hpp"

Skim::Skim(std::vector<std::string> input, std::string output) {
  fin = std::move(input);
  fout = std::move(output);
  chain = new TChain("h10");

  for (auto f : fin) chain->AddFile(f.c_str());

  RootOutputFile = new TFile(fout.c_str(), "RECREATE");
  RootOutputFile->SetCompressionSettings(404);
  if (getenv("BEAM_E") != NULL) {
    BEAM_ENERGY = strtod(getenv("BEAM_E"),(char **)NULL);
    std::cout << RED << "Beam energy set to: " << BEAM_ENERGY << DEF << std::endl;
  }
  e_mu = std::make_shared<TLorentzVector>();
  e_mu->SetPxPyPzE(0.0, 0.0, sqrt((BEAM_ENERGY * BEAM_ENERGY) - (MASS_E * MASS_E)), BEAM_ENERGY);
  MM_neutron = std::make_shared<MissingMass>(MASS_P, 0.0);
}
Skim::~Skim() {}

float Skim::Basic() {
  std::cout << BLUE << "Basic Skim " << GREEN << fout << DEF << std::endl;

  int num_of_events = (int)chain->GetEntries();
  TTree *skim = chain->CloneTree(0);
  auto data = std::make_shared<Branches>(chain);
  int total = 0;
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    auto check = std::make_unique<Cuts>(data);
    if (check->isElecctron()) {
      total++;
      skim->Fill();
    }
  }

  delete chain;
  RootOutputFile->cd();
  RootOutputFile->Write();
  RootOutputFile->Close();
  delete RootOutputFile;

  return (total/(float)num_of_events);
}

void Skim::Strict() {
  int num_of_events;
  std::cout << BLUE << "Strict Skim file " << GREEN << fout << DEF << std::endl;

  num_of_events = (int)chain->GetEntries();
  auto skim = chain->CloneTree(0);
  auto data = std::make_shared<Branches>(chain);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (data->gpart() >= 4) continue;
    auto check = std::make_unique<Cuts>(data);
    auto event = std::make_unique<Reaction>(data);

    bool cuts = true;
    cuts &= check->isElecctron();
    cuts &= check->Beam_cut();
    cuts &= (event->Q2() >= 1.0);
    cuts &= (event->W() >= 0.8);
    cuts &= (event->W() <= 2.0);
    if (cuts) skim->Fill();  // Fill the banks after the skim
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
  auto data = std::make_shared<Branches>(chain);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (data->ec_ei(0) < 0.01 || data->ec_eo(0) < 0.01) continue;
    auto check = std::make_unique<Cuts>(data);
    if (!check->isElecctron()) continue;

    auto event = std::make_unique<Reaction>(data);

    auto dt = std::make_unique<Delta_T>(data->sc_t(0), data->sc_r(0));
    std::vector<double> dt_proton = dt->delta_t_array(MASS_P, data);
    std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, data);

    for (int part_num = 1; part_num < data->gpart(); part_num++) {
      if (data->q(part_num) == POSITIVE) {
        if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
          event->SetPip(part_num);
        } else if (check->dt_P_cut(dt_proton.at(part_num), data->p(part_num))) {
          event->SetProton(part_num);
        }
      } else if (data->q(part_num) == NEGATIVE) {
        if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
          event->SetPim(part_num);
        }
      }
    }

    bool mm_cut = true;
    mm_cut &= (event->MM() < 1.5);
    mm_cut &= (event->MM() > 0.7);
    if ((event->SinglePip() || event->NeutronPip()) && mm_cut) skim->Fill();
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
