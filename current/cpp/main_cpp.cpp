/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

//#include <memory>
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
    try {
      files = glob(argv[1]);
    } catch (...) {
      files = glob("/Users/tylern/Data/e1d/new_cook/skim*.root");
    }
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

  auto hist = std::make_shared<Histogram>();
  auto dh = std::make_unique<DataHandeler>(files, hist);
  dh->setLoadBar(true);
  auto start = std::chrono::high_resolution_clock::now();
  size_t events = 0;
  events += dh->Run();
  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout.imbue(std::locale(""));
  std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;
  hist->Write(outfilename, true);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;
  return 0;
}
