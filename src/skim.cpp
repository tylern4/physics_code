/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "skim.hpp"

Skim::Skim(const std::shared_ptr<TChain> &chain) : _chain(chain) {}

Skim::~Skim() {}

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

    auto dt = std::make_unique<Delta_T>(data);

    for (int part_num = 1; part_num < data->gpart(); part_num++) {
      if (data->q(part_num) == POSITIVE) {
        if (check->dt_Pip_cut(part_num)) {
          event->SetPip(part_num);
        } else if (check->dt_P_cut(part_num)) {
          event->SetProton(part_num);
        }
      } else if (data->q(part_num) == NEGATIVE) {
        if (check->dt_Pip_cut(part_num)) {
          event->SetPim(part_num);
        }
      }
    }

    if (event->SinglePip() || event->NeutronPip() || event->PPi0()) skim->Fill();
  }

  return skim;
}
