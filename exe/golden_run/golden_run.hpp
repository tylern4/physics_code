/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef GOLDEN_RUN_H
#define GOLDEN_RUN_H
#include "branches.hpp"
#include "cuts.hpp"
#include "main.h"

// Finds the files with golden runs
#define COUNTER 10000
std::string golden_run(const std::vector<std::string>& fins) {
  std::string golden_run;
  std::string file_num, run_num, line;
  int n_evnt = 0;
  double total_q = 0.0, curr_q = 0.0, prev_q = 0.0, delta_q = 0.0;
  bool electron_cuts;
  LorentzVector e_mu(0.0, 0.0, E1D_E0, MASS_E);
  LorentzVector ZERO(0.0, 0.0, 0.0, 0.0);
  int line_num = 0;
  // const char* progress = "-\\|/";

  for (auto& fin : fins) {
    line_num++;
    run_num = fin.substr(fin.find("/h10") + 6, 5);
    file_num = fin.substr(fin.find("/h10") + 12, 2);

    // run_num = fin.substr(fin.find("run") + 4, 5);
    // file_num = fin.substr(fin.find("pass1") + 7, 2);

    // run_num = fin.substr(fin.find("skim_r") + 6, 5);
    // file_num = fin.substr(fin.find("skim_r") + 12, 2);

    // std::cout << run_num << "\n";
    // std::cout << file_num << "\n";

    auto chain = std::make_shared<TChain>("h10");
    chain->Add(fin.c_str());
    auto data = std::make_shared<Branches>(chain);

    size_t num_of_events = chain->GetEntries();
    total_q = 0.0;

    for (int current_event = 0; current_event < num_of_events; current_event++) {
      chain->GetEntry(current_event);
      auto check = std::make_shared<e1d_Cuts>(data);

      if (check->check_banks()) {
        n_evnt++;
        curr_q = data->q_l();
        if (curr_q > 0.0) {
          if (curr_q > prev_q) {
            delta_q = curr_q - prev_q;
            total_q += delta_q;
          }
          prev_q = curr_q;
        }
      }
    }
    if (n_evnt != 0 && total_q != 0)
      golden_run += run_num + "," + file_num + "," + num_of_events + "," + total_q + "\n";
  }
  return golden_run;
}
#endif
