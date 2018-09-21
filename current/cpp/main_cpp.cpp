/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "../src/branches.hpp"
#include "../src/classes.hpp"
#include "../src/constants.hpp"
#include "../src/datahandeler.hpp"
#include "../src/glob_files.hpp"
#include "../src/physics.hpp"
#include "TStopwatch.h"
#include "main.h"

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
    std::cerr << BOLDWHITE << argv[0] << " infile.root outfile.root\n\n";
    std::cerr << RESET << std::endl;
    return 1;
  }

  std::string outfilename;
  if (argc == 2) {
    outfilename = "out.root";
  } else if (argc == 3) {
    outfilename = argv[2];
  }
  DataHandeler *dh[files.size()];
  Histogram *hist[files.size()];

  int i = 0;
#pragma omp parallel for private(i, hist, dh)
  for (i = 0; i < files.size(); i++) {
    std::cout << i << "\t" << files.at(i) << std::endl;
    // hist[i] = new Histogram();
    // dh[i] = new DataHandeler(files.at(i), hist[i]);
  }

  Watch->Stop();
  cout << RED << Watch->RealTime() << "sec" << DEF << endl;

  return 0;
}
