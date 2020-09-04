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
      if (check->Sanity()) {
        skim->Fill();
      }
    }
    return skim;
  }
  // void Strict();
  std::shared_ptr<TTree> Final();
};
#endif
