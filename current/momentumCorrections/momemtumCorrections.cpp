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
  if (argc < 3) {
    exit(1);
  }
  std::string outfilename = argv[1];
  std::ofstream myfile;
  myfile.open(outfilename);
  std::vector<std::string> out;
  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
  for (int i = 2; i < argc; i++) infilenames[i % NUM_THREADS].push_back(argv[i]);

  std::future<std::string> mom_corrections[NUM_THREADS];
  int i = 0;
  for (auto &fils : infilenames) mom_corrections[i++] = std::async(mom_correction_csv, fils);

  myfile << "e_px,e_py,e_pz,p_px,p_py,p_pz,W_uncorr,Q2_uncorr,sector" << std::endl;
  for (auto &o : mom_corrections) myfile << o.get();

  return 0;
}
