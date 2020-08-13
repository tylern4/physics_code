/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include <algorithm>
#include <future>
#include <iostream>
#include <thread>
#include <vector>
#include "TStopwatch.h"
#include "branches.hpp"
#include "constants.hpp"
#include "datahandeler.hpp"
#include "physics.hpp"

int main(int argc, char** argv) {
  auto start = std::chrono::high_resolution_clock::now();
  std::vector<std::string> infilenames;
  std::string outfilename;
  if (argc >= 2) {
    outfilename = argv[1];
    for (size_t i = 2; i < argc; i++) infilenames.push_back(argv[i]);
  } else {
    return 1;
  }
  auto output_root = std::make_shared<TFile>("mm_for_dt.root", "RECREATE");
  auto mm_hist = std::make_shared<TH1D>("mm_for_dt", "mm_for_dt", 500, 0, 3);
  output_root->cd();
  std::cout.imbue(std::locale(""));

  std::ofstream csv_output;
  csv_output.open(outfilename);
  csv_output << "sec,theta,phi,charge,vertex,p,sc_t,sc_r\n";

  auto chain = std::make_shared<TChain>("h10");
  for (auto& f : infilenames) chain->Add(f.c_str());
  size_t events = chain->GetEntries();
  auto data = std::make_shared<Branches>(chain);
  int progress = (events / 100);

  std::cout << " [" << std::flush;
  for (int current_event = 0; current_event < events; current_event++) {
    if (current_event % progress == 0) std::cout << "%" << std::flush;
    chain->GetEntry(current_event);
    auto check = std::make_unique<Cuts>(data);
    if (!check->isStrictElecctron()) continue;
    if (!check->Fid_cut()) continue;

    auto event = std::make_shared<Reaction>(data);

    for (int part_num = 1; part_num < data->gpart(); part_num++) {
      if (check->Protish(part_num)) {
        event->SetProton(part_num);
      }
    }
    if (event->elastic()) {
      mm_hist->Fill(event->MM2());
      auto dt = std::make_shared<Delta_T>(data);
      csv_output << *dt;
    }
  }
  std::cout << "]\n";

  output_root->Write();

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout.imbue(std::locale(""));
  std::cout << RED << events << " events\t" << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;
  return 0;
}
