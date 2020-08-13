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
#include "glob_files.hpp"
#include "main.h"
#include "physics.hpp"

int main(int argc, char** argv) {
  ROOT::EnableThreadSafety();
  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
  std::string outfilename;
  if (argc >= 2) {
    outfilename = argv[1];
    for (int i = 2; i < argc; i++) infilenames[i % NUM_THREADS].push_back(argv[i]);
  } else {
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();

  auto hist = std::make_shared<mcHistogram>(outfilename);
  auto mom_corr = std::make_shared<MomCorr>();
  std::cout.imbue(std::locale(""));
  size_t events = 0;

  std::future<size_t> t[NUM_THREADS];

  for (size_t i = 0; i < NUM_THREADS; i++) {
    t[i] = std::async(run_e1f_file, infilenames.at(i), hist, mom_corr, i);
  }
  for (size_t i = 0; i < NUM_THREADS; i++) {
    events += t[i].get();
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout.imbue(std::locale(""));
  hist->Write();
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;
  return 0;
}
