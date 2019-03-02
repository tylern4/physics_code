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

  auto dh = std::make_shared<DataHandeler>();
  auto hist = std::make_shared<Histogram>();
  auto Watch = std::make_shared<TStopwatch>();
  Watch->Start();
  size_t events = 0;
  if (files.size() > 1) {
    for (int i = 0; i < files.size(); i++) {
      loadbar(i, files.size() - 1);
      events += dh->Run(files.at(i), hist);
    }
    Watch->Stop();
    cout << BOLDGREEN << "\n\n" << events / Watch->RealTime() << "Hz" << DEF << endl;
    hist->Write(outfilename, true);
    cout << RED << Watch->RealTime() << "sec" << DEF << endl;
    cout << BOLDYELLOW << "\n\n" << events / Watch->RealTime() << "Hz" << DEF << endl;

    return 0;
  }
  /*
  else {
    dh->Run(files.at(0), hist);
    hist->Write(outfilename, true);
    return 0;
  }
  */

  return 1;
}
