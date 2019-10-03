/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "skim.hpp"

Skim::Skim(const std::shared_ptr<TChain> &chain) : _chain(chain) {}

Skim::~Skim() {}

std::shared_ptr<TTree> Skim::Basic() {
  int num_of_events = (int)_chain->GetEntries();
  std::shared_ptr<TTree> skim(_chain->CloneTree(0));
  auto data = std::make_shared<Branches>(_chain);
  int total = 0;
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    _chain->GetEntry(current_event);
    auto check = std::make_unique<Cuts>(data);
    if (check->isElecctron()) {
      total++;
      skim->Fill();
    }
  }
  return skim;
}

/*
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

  RootOutputFile->cd();
  RootOutputFile->Write();
  RootOutputFile->Close();
}
*/

std::shared_ptr<TTree> Skim::Final() {
  int num_of_events = (int)_chain->GetEntries();
  std::shared_ptr<TTree> skim(_chain->CloneTree(0));
  auto data = std::make_shared<Branches>(_chain);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    _chain->GetEntry(current_event);
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

    if (event->SinglePip() || event->NeutronPip() || event->PPi0()) skim->Fill();
  }

  return skim;
}
