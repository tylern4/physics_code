/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler
 */
/*	University Of South Carolina
 */
/************************************************************************/

// Only My Includes. All others in main.h
#include "../src/classes.hpp"
//#include "../src/fit_functions.hpp"
//#include "../src/missing_mass_gaussians.hpp"
#include "../src/physics.hpp"
#include "TStopwatch.h"
#include "skim.hpp"

using namespace std;

int main(int argc, char **argv) {
  if (argc == 3) {
    char *infilename = argv[1];
    char *outfilename = argv[2];
    skim(infilename, outfilename);
  }

  return 0;
}
