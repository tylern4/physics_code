/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef MAIN_H
#define MAIN_H
#include <TFile.h>
#include <TFileCollection.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include "../src/color.hpp"
#include "../src/constants.hpp"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "time.h"

#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)
#define Square(a) (a * a)

using namespace std;

// ofstream cut_outputs;

void loadbar(long x, long n) {
  int w = 50;
  if ((x != n) && (x % (n / 100 + 1) != 0)) return;

  double ratio = x / (double)n;
  int c = ratio * w;

  cout << BLUE << " [";
  for (int x = 0; x < c; x++) cout << GREEN << "=" << DEF;
  cout << GREEN << ">" << DEF;
  for (int x = c; x < w; x++) cout << " ";
  cout << BLUE << (int)(ratio * 100) << "%]\r" << DEF << flush;
}

// double Square(double a) { return a * a; }

#endif
