/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "511_lab.hpp"
#include "TStopwatch.h"
#include "branches.hpp"
#include "classes.hpp"
#include "constants.hpp"
#include "main.h"

int main(int argc, char **argv) {
  TStopwatch *Watch = new TStopwatch;
  Watch->Start();

  if (argc == 2) {
    std::string infilename = argv[1];
    make_electron_csv(infilename);
    make_mm_csv(infilename);
  } else if (argc == 1) {
    analyze_wq2("511_lab_E_data.root", "out_wq2.root");
    analyze_MM("511_lab_E_PIP_data.root", "out_mm.root");
  } else {
    std::cerr << RED << "Wrong" << DEF << std::endl;
  }

  Watch->Stop();
  std::cout << RED << Watch->RealTime() << "sec" << DEF << std::endl;

  return 0;
}
