/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "classes.hpp"
#include "main.h"
#include "constants.h"
#include "TStopwatch.h"
#include "physics.hpp"
#include "delta_t_cut.hpp"
#include "datahandeler.hpp"

using namespace std;

int main(int argc, char **argv) {
  TStopwatch *Watch = new TStopwatch;
  Watch->Start();
  gStyle->SetOptFit(1111);

  if (argc == 3) {
    char *infilename = argv[1];
    char *outfilename = argv[2];
    dataHandeler(infilename, outfilename, true);
  }

  if (argc == 4) {
    char *infilename = argv[1];
    char *outfilename = argv[2];
    dataHandeler(infilename, outfilename, false);
  }

  Watch->Stop();
  cout << red << Watch->RealTime() << "sec" << def << endl;

  return 0;
}
