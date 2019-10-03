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

size_t run_file(std::vector<std::string> in, std::shared_ptr<Histogram> hists, int thread_id) {
  auto dh = std::make_unique<DataHandeler>(in, hists);
  dh->setLoadBar(false);
  if (thread_id == 0) dh->setLoadBar(true);
  size_t tot = 0;
  tot += dh->Run();
  return tot;
}

int main(int argc, char **argv) {
  ROOT::EnableThreadSafety();
  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
  std::string outfilename;
  if (argc >= 2) {
    outfilename = argv[1];
    for (size_t i = 2; i < argc; i++) infilenames[i % NUM_THREADS].push_back(argv[i]);
  } else {
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();

  auto hist = std::make_shared<Histogram>();
  std::cout.imbue(std::locale(""));
  size_t events = 0;

  std::future<size_t> threads[NUM_THREADS];

  for (size_t i = 0; i < NUM_THREADS; i++) {
    threads[i] = std::async(run_file, infilenames.at(i), hist, i);
  }

  for (size_t i = 0; i < NUM_THREADS; i++) {
    events += threads[i].get();
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout.imbue(std::locale(""));
  ROOT::EnableImplicitMT(NUM_THREADS);
  hist->Write(outfilename);
  std::cout << RED << events << " events\t" << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;
  return 0;
}
