/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

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

  auto dh = std::make_unique<DataHandeler>();
  auto hist = std::make_shared<Histogram>();
  auto start = std::chrono::high_resolution_clock::now();
  size_t events = 0;
  if (files.size() > 1) {
    for (int i = 0; i < files.size(); i++) {
      loadbar(i, files.size() - 1);
      events += dh->Run(files.at(i), hist);
    }
    hist->Write(outfilename, true);
    std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
    std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
    std::cout.imbue(std::locale(""));
    std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;
  } else {
    dh->Run(files.at(0), hist);
    hist->Write(outfilename, true);
  }

  return 0;
}
