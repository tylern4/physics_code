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

  std::vector<std::string> infile;
  std::string outfile;
  if (argc == 2) {
    infile = glob(argv[1]);
    outfile = infile[0].substr(0, infile[0].size() - 8) + "_skim.root";
  } else if (argc == 3) {
    infile = glob(argv[1]);
    outfile = argv[2];
  } else {
    std::cerr << RED << "Error: \n";
    std::cerr << BOLDRED << "\tNeed input file and output file\n";
    std::cerr << RESET << "Usage:\n\t";
    std::cerr << BOLDWHITE << argv[0] << " '/path/to/files/to/skim_*.root'\n\n";
    std::cerr << RESET << std::endl;
    return 1;
  }

  Skim *s = new Skim(infile, outfile);
  s->Process();

  return 0;
}
