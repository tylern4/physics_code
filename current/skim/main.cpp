/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler
 */
/*	University Of South Carolina
 */
/************************************************************************/

// Only My Includes. All others in main.h
//#include "main.h"
#include "../src/classes.hpp"
#include "TStopwatch.h"
#include "../src/physics.hpp"
#include "../src/fit_functions.hpp"
#include "../src/missing_mass_gaussians.hpp"
#include "skim.hpp"

using namespace std;

int main(int argc, char **argv) {
  TStopwatch *Watch = new TStopwatch;
  Watch->Start();
  gStyle->SetOptFit(1111);

  if (argc == 3) {
    char *infilename = argv[1];
    char *outfilename = argv[2];
    skim(infilename, outfilename);
  }

  Watch->Stop();
  // cout << red << Watch->RealTime() << "sec" << def << endl;

  return 0;
}
