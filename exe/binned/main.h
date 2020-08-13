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
#include "TChain.h"
#include "TF1.h"
#include "TH2.h"
#include "THn.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "time.h"

#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)
#define Square(a) (a * a)

#endif
