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
  std::vector<std::string> infilename;
  std::string outfilename;
  if (argc >= 2) {
    outfilename = argv[1];
    for (int i = 2; i < argc; i++) infilename.push_back(argv[i]);
  } else {
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();
  std::cout.imbue(std::locale(""));
  size_t events = 0;
  int j = 0;
  Yeilds *dh = new Yeilds(outfilename, true);

  for (auto &in : infilename) {
    events += dh->RunNtuple(in);
    std::chrono::duration<double> elapsed_temp = (std::chrono::high_resolution_clock::now() - start);
    std::cout << BOLDYELLOW << "\t" << events / elapsed_temp.count() << " Hz\r\r" << DEF << std::flush;
  }

  delete dh;
  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;

  return 0;
}
