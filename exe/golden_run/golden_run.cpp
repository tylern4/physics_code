/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "golden_run.hpp"
#include <future>
#include <iostream>
#include <thread>
#include "branches.hpp"
#include "constants.hpp"

int main(int argc, char **argv) {
  ROOT::EnableThreadSafety();

  if (argc < 3) {
    exit(1);
  }

  std::string outfilename = argv[1];
  std::fstream myfile;
  myfile.open(outfilename, std::ios::out | std::ios::binary);
  if (!myfile.is_open()) exit(1);

  myfile << "run_num,file_num,num_of_events,total_q" << std::endl;
  for (int i = 2; i < argc; i++) myfile << golden_run(argv[i]);

  return 0;
}
