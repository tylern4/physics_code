/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef GOLDEN_RUN_H
#define GOLDEN_RUN_H
#include "branches.hpp"
#include "cuts.hpp"
#include "main.h"

std::string golden_run(const std::string& fin) {
  std::string golden_run;
  std::string file_num, run_num;
  size_t n_evnt = 0;

  // run_num = fin.substr(fin.find("/h10") + 6, 5);
  // file_num = fin.substr(fin.find("/h10") + 12, 2);

  // run_num = fin.substr(fin.find("run") + 4, 5);
  // file_num = fin.substr(fin.find("pass1") + 7, 2);
  auto file_name = fin.substr(fin.find_last_of("/\\") + 1);

  run_num = file_name.substr(file_name.find("_r2") + 2, 5);
  file_num = file_name.substr(file_name.find("_r2") + 8, 2);
  // std::cout << file_name << std::endl;
  // std::cout << run_num << "\t";
  // std::cout << file_num << "\n";

  auto chain = std::make_shared<TChain>("h10");
  chain->Add(fin.c_str());
  auto data = std::make_shared<Branches>(chain);

  size_t num_of_events = chain->GetEntries();
  float total_q = 0.0;
  float qcurr;
  float qprev;
  float deltaq;
  float q_temp;

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    int npip = 0;
    chain->GetEntry(current_event);
    auto check = std::make_unique<e1d_Cuts>(data);

    if (data->gpart() < 0) continue;
    q_temp = data->q_l();
    qcurr = q_temp;
    // cout<<"q_l="<<q_l<<endl;
    if (q_temp <= 0) continue;

    if (q_temp > 0.) {
      // cout<<"q_l"<<q_temp<<"qcurr"<<qcurr<<endl;
      if ((qcurr > qprev)) {
        deltaq = qcurr - qprev;
        total_q += deltaq;
        // cout << setw(10) << "qcurr= " << qcurr << setw(10) << " qprev= " << qprev << setw(10) << " deltaq= " <<
        // deltaq
        //      << endl;
        // cout<<"q_l"<<q_temp<<"qcurr"<<qcurr<<"qprev"<<qprev<<"deltaq"<<deltaq<<endl;
      }
      qprev = qcurr;
    }

    if (check->isElectron()) {
      auto event = std::make_shared<Reaction>(data);
      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        if (check->Pip(part_num)) {
          npip++;
          event->SetPip(part_num);
        } else if (check->Prot(part_num)) {
          event->SetProton(part_num);
        } else if (check->Pim(part_num)) {
          event->SetPim(part_num);
        } else
          event->SetOther(part_num);
      }
      if (npip >= 1) n_evnt++;
    }
  }
  golden_run += run_num + "," + file_num + "," + n_evnt + "," + num_of_events + "," + total_q + "," +
                num_of_events / total_q + "\n";
  std::cout << golden_run;
  return golden_run;
}
#endif
