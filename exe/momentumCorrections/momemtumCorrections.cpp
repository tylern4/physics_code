/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "momemtumCorrections.hpp"
#include <future>
#include <thread>

int main(int argc, char **argv) {
  ROOT::EnableThreadSafety();
  auto start = std::chrono::high_resolution_clock::now();

  if (argc < 3) {
    exit(1);
  }
  std::string outfilename = argv[1];
  auto syncFile = std::make_shared<SyncFile>(outfilename);
  syncFile->write("e_p,e_theta,e_phi,p_p,p_theta,p_phi,W_uncorr,Q2_uncorr,sector,type\n");

  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
  for (int i = 2; i < argc; i++) infilenames[i % NUM_THREADS].push_back(argv[i]);

  std::future<size_t> mom_corrections[NUM_THREADS];
  auto mom_corr = nullptr;  // std::make_shared<MomCorr>();
  int i = 0;
  for (auto &fils : infilenames) mom_corrections[i++] = std::async(mom_correction_csv, fils, syncFile, mom_corr);
  size_t total = 0;
  for (auto &num : mom_corrections) {
    total += num.get();
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout.imbue(std::locale(""));
  std::cout << RED << NUM_THREADS << " " << total / elapsed_full.count() << " Hz\t";
  std::cout << DEF << std::endl;
  return 0;
}
