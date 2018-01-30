/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef MAIN_H
#define MAIN_H
#include <stdlib.h>
#include <stdio.h>
#include <TFileCollection.h>
#include <TFile.h>
#include "TMath.h"
#include <TGraph.h>
#include "THnSparse.h"
#include "TTree.h"
#include "TROOT.h"
#include <TLorentzVector.h>
#include <string.h>
#include <string>
#include <cstring>
#include "time.h"
#include <string>
#include <cstring>
#include "TTree.h"
#include "TROOT.h"
#include "TH2.h"
#include <TFile.h>
#include "TStyle.h"
#include "TCanvas.h"
#include <fstream>
#include "TF1.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TSystem.h"
#include <vector>
#include <fstream>
#include "../src/color.cpp"
#include "../src/constants.hpp"

using namespace std;

ofstream cut_outputs;

// Color outputs
Color::Modifier red(Color::FG_RED);
Color::Modifier blue(Color::FG_BLUE);
Color::Modifier def(Color::FG_DEFAULT);
Color::Modifier green(Color::FG_GREEN);
Color::Modifier bgred(Color::BG_RED);
Color::Modifier bgblue(Color::BG_BLUE);
Color::Modifier bgdef(Color::BG_DEFAULT);
Color::Modifier bggreen(Color::BG_GREEN);

void loadbar(long x, long n) {

  int w = 50;
  if ((x != n) && (x % (n / 100 + 1) != 0))
    return;

  double ratio = x / (double)n;
  int c = ratio * w;

  cout << blue << " [";
  for (int x = 0; x < c; x++)
    cout << green << "=" << def;
  cout << green << ">" << def;
  for (int x = c; x < w; x++)
    cout << " ";
  cout << blue << (int)(ratio * 100) << "%]\r" << def << flush;
}

double Square(double a) { return a * a; }

#endif
