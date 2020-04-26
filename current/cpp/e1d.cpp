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
#include "main.h"
#include "physics.hpp"

size_t run_file(const std::vector<std::string>& in, const std::shared_ptr<Histogram>& hists,
                const std::shared_ptr<MomCorr>& mom_corr, int thread_id) {
  auto dh = std::make_unique<DataHandeler>(in, hists, mom_corr);
  dh->setLoadBar(false);
  if (thread_id == 0) dh->setLoadBar(true);
  size_t tot = 0;
  tot += dh->Run();
  return tot;
}

int main(int argc, char** argv) {
  auto start = std::chrono::high_resolution_clock::now();
  ROOT::EnableThreadSafety();
  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
  std::string outfilename;

  if (argc >= 2) {
    outfilename = argv[1];
    for (size_t i = 2; i < argc; i++) infilenames[i % NUM_THREADS].push_back(argv[i]);
  } else {
    return 1;
  }

  // Create the shared objects for hists and momentum corrections
  auto hist = std::make_shared<Histogram>();
  auto mom_corr = std::make_shared<MomCorr>();

  std::cout.imbue(std::locale(""));
  size_t events = 0;

  // Make an array of futures which return an integer
  std::future<size_t> threads[NUM_THREADS];

  // For the number of threads run the list of file
  for (size_t i = 0; i < NUM_THREADS; i++) {
    threads[i] = std::async(run_file, infilenames.at(i), hist, mom_corr, i);
  }

  // Let the threads run and get the results when they're done
  for (size_t i = 0; i < NUM_THREADS; i++) {
    events += threads[i].get();
  }

  std::cout.imbue(std::locale(""));
  ROOT::EnableImplicitMT(2);
  // Write the histograms to file
  hist->Write(outfilename);

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << events << " events\t" << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;

  std::cout << "Amdahl " << NUM_THREADS << "," << elapsed_full.count() << std::endl;

  return 0;
}
