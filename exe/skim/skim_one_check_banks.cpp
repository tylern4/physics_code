/**************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

// Only My Includes. All others in main.h
#include <future>
#include <thread>
#include "TFile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "glob_files.hpp"
#include "physics.hpp"
#include "skim.hpp"

std::shared_ptr<TTree> run_file(const std::string& in) {
  std::shared_ptr<TChain> h10 = std::make_shared<TChain>("h10");
  h10->Add(in.c_str());
  auto s = std::make_unique<Skim>(h10);
  return s->Basic<Cuts>();
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << argv[0] << " infile.root";
    exit(1);
  }
  auto start = std::chrono::high_resolution_clock::now();

  std::string in_file = argv[1];

  auto skimmed_tree = run_file(in_file);

  auto filename = "skim.root";
  if (argc == 3) {
    filename = argv[2];
  }

  auto outFile = std::make_unique<TFile>(filename, "RECREATE");
  outFile->cd();
  skimmed_tree->Write();
  outFile->Write();
  outFile->Close();

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout.imbue(std::locale(""));
  std::cout << RED << elapsed_full.count() << " sec\t";
  std::cout << DEF << std::endl;

  return 0;
}
