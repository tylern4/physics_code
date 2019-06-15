/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"
#include "yeilds.hpp"

using namespace std;

int main(int argc, char **argv) {
  std::string outfilename = "test.csv";
  std::vector<std::string> infilename;
  if (argc >= 2) {
    for (int i = 1; i < argc; i++) infilename.push_back(argv[i]);
  } else {
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();
  std::cout.imbue(std::locale(""));
  int events = 0;
  int j = 0;
  Yeilds *dh = new Yeilds(outfilename);

  for (int i = 0; i < infilename.size(); i++) {
    events += dh->Run(infilename.at(i));

    std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
    std::cout << RED << elapsed_full.count();
    std::cout << BOLDYELLOW << "\t\t" << events / elapsed_full.count() << " Hz\r\r" << DEF << std::flush;
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;

  return 0;
}
