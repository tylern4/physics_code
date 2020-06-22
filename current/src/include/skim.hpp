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

  std::shared_ptr<TTree> Basic();
  // void Strict();
  std::shared_ptr<TTree> Final();
};
#endif
