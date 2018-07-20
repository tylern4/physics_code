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
#include "../src/skim.hpp"
#include "TStopwatch.h"

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
  int i = 0;
  for (i = 0; i < num_files; i++) {
    Skim *s = new Skim(files.at(i).c_str());
    s->Process();
  }

  return 0;
}
