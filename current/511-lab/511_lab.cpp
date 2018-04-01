/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "main.h"
#include "../src/classes.hpp"
#include "../src/constants.hpp"
#include "../src/branches.hpp"
#include "TStopwatch.h"
#include "../src/physics.hpp"
#include "511_lab.hpp"

int main(int argc, char **argv) {
  TStopwatch *Watch = new TStopwatch;
  Watch->Start();

  if (argc == 2) {
    char *infilename = argv[1];
    make_electron_csv(infilename);
  } else if (argc == 3) {
    char *infilename = argv[1];
    char *outfilename = argv[2];
    analyze(infilename, outfilename);
  } else {
    std::cerr << RED << "Wrong" << DEF << std::endl;
  }

  Watch->Stop();
  std::cout << RED << Watch->RealTime() << "sec" << DEF << std::endl;

  return 0;
}
