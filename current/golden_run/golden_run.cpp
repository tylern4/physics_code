/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "golden_run.hpp"
#include "TStopwatch.h"
#include "branches.hpp"
#include "constants.hpp"
#include "main.h"

int main(int argc, char **argv) {
  TStopwatch *Watch = new TStopwatch;
  Watch->Start();
  gStyle->SetOptFit(1111);

  if (argc >= 3) {
    std::string outfilename = argv[1];
    std::vector<std::string> infilenames;
    for (int i = 2; i < argc; i++) infilenames.push_back(argv[i]);
    golden_run(infilenames, outfilename);
  }

  Watch->Stop();
  std::cout << RED << Watch->RealTime() << "sec" << DEF << std::endl;

  return 0;
}
