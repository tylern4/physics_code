/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#ifndef SKIM_H_GUARD
#define SKIM_H_GUARD
#include <iostream>
#include "TChain.h"
#include "branches.hpp"
#include "color.hpp"
#include "cuts.hpp"
#include "delta_t.hpp"
#include "reaction.hpp"

class Skim {
 private:
  std::shared_ptr<TChain> _chain;

 public:
  Skim(const std::shared_ptr<TChain>& _chain);
  ~Skim();

  template <class CutType>
  std::shared_ptr<TTree> Basic() {
    int num_of_events = (int)_chain->GetEntries();
    std::shared_ptr<TTree> skim(_chain->CloneTree(0));
    auto data = std::make_shared<Branches>(_chain);
    for (int current_event = 0; current_event < num_of_events; current_event++) {
      _chain->GetEntry(current_event);
      auto check = std::make_unique<CutType>(data);
      if (check->check_banks()) {
        skim->Fill();
      }
    }
    return skim;
  }

  template <class CutType>
  std::shared_ptr<TTree> Elastic() {
    int num_of_events = (int)_chain->GetEntries();
    std::shared_ptr<TTree> skim(_chain->CloneTree(0));
    auto data = std::make_shared<Branches>(_chain);
    float beam_energy;

    if (std::is_same<CutType, e1f_Cuts>::value) {
      beam_energy = E1F_E0;
    } else if (std::is_same<CutType, e16_Cuts>::value) {
      beam_energy = E16_E0;
    } else {
      beam_energy = E1D_E0;
    }

    for (int current_event = 0; current_event < num_of_events; current_event++) {
      _chain->GetEntry(current_event);
      auto check = std::make_unique<CutType>(data);
      if (!check->isElecctron()) continue;

      auto event = std::make_shared<Reaction>(data, beam_energy);
      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        if (check->Pip(part_num)) {
          event->SetPip(part_num);
        } else if (check->Prot(part_num)) {
          event->SetProton(part_num);
        } else if (check->Pim(part_num)) {
          event->SetPim(part_num);
        } else
          event->SetOther(part_num);
      }
      if (event->elastic()) skim->Fill();
    }
    return skim;
  }

  std::shared_ptr<TTree> Final();
};
#endif
