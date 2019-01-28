/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include <memory>
#include "TStopwatch.h"
#include "branches.hpp"
#include "constants.hpp"
#include "datahandeler.hpp"
#include "glob_files.hpp"
#include "main.h"
#include "physics.hpp"

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

  DataHandeler *dh = new DataHandeler();
  Histogram *hist = new Histogram();

  if (files.size() > 1) {
    for (int i = 0; i < files.size(); i++) {
      loadbar(i, files.size() - 1);
      dh->Run(files.at(i), hist);
    }
    hist->Write(outfilename);
  } else {
    dh->Run(files.at(0), hist);
    hist->Write(outfilename, true);
  }
  Watch->Stop();
  cout << RED << Watch->RealTime() << "sec" << DEF << endl;

  return 0;
}
