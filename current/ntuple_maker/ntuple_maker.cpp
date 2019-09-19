/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include <algorithm>
#include <future>
#include <iostream>
#include <thread>
#include <vector>
#include "ntuple.hpp"

size_t run_file(std::vector<std::string> in, std::string out_path, int thread_id) {
  std::string outf = out_path + "/ntuple_" + std::to_string(thread_id) + ".root";
  auto chain = std::make_shared<TChain>("h10");
  for (auto& f : in) chain->Add(f.c_str());
  auto dh = std::make_unique<Ntuple>(outf);
  size_t tot = dh->Run(chain);
  return tot;
}

int main(int argc, char* argv[]) {
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
  std::cout.imbue(std::locale(""));
  size_t events = 0;

  std::future<size_t> t[NUM_THREADS];
  for (size_t i = 0; i < NUM_THREADS; i++) {
    t[i] = std::async(run_file, infilenames.at(i), outfilename, i);
  }
  for (size_t i = 0; i < NUM_THREADS; i++) {
    events += t[i].get();
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;

  return 0;
}
