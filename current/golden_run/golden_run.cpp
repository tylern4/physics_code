/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "golden_run.hpp"
#include <future>
#include <thread>
#include "branches.hpp"
#include "constants.hpp"

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

  std::future<std::string> golden_runs[NUM_THREADS];
  int i = 0;
  for (auto &fils : infilenames) golden_runs[i++] = std::async(golden_run, fils);

  myfile << "run_num,file_num,num_of_events,total_q" << std::endl;
  for (auto &o : golden_runs) myfile << o.get();

  return 0;
}
