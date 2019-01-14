/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "TStopwatch.h"
#include "branches.hpp"
#include "constants.hpp"
#include "csv_maker.hpp"
#include "main.h"
#include "physics.hpp"

using namespace std;

int main(int argc, char **argv) {
  TStopwatch *Watch = new TStopwatch;
  Watch->Start();
  gStyle->SetOptFit(1111);

  if (argc == 3) {
    char *infilename = argv[1];
    char *outfilename = argv[2];
    make_electron_csv(infilename, outfilename);
  }

  Watch->Stop();
  cout << RED << Watch->RealTime() << "sec" << DEF << endl;

  return 0;
}
