/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler
 */
/*	University Of South Carolina
 */
/************************************************************************/

// Only My Includes. All others in main.h
#include "../src/classes.hpp"
#include "../src/glob_files.hpp"
#include "../src/physics.hpp"
#include "TStopwatch.h"
#include "skim.hpp"

using namespace std;

int main(int argc, char **argv) {
  gROOT->SetBatch(true);
  TStopwatch *Watch = new TStopwatch;
  Watch->Start();
  gStyle->SetOptFit(1111);

  std::vector<std::string> files;
  if (argc >= 2) {
    files = glob(argv[1]);
  } else if (argc < 2) {
    std::cerr << RED << "Error: \n";
    std::cerr << BOLDRED << "\tNeed input file and output file\n";
    std::cerr << RESET << "Usage:\n\t";
    std::cerr << BOLDWHITE << argv[0] << " '/path/to/files/to/skim_*.root'\n\n";
    std::cerr << RESET << std::endl;
    return 1;
  }

  int num_files = files.size();
#pragma omp parallel for private(files)
  for (int i = 0; i < num_files; i++) {
    skim(files[i].c_str());
  }

  return 0;
}
