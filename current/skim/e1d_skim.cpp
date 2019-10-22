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

std::shared_ptr<TTree> run_file(const std::vector<std::string>& in) {
  std::shared_ptr<TChain> h10 = std::make_shared<TChain>("h10");
  for (auto& f : in) h10->Add(f.c_str());
  auto s = std::make_unique<Skim>(h10);
  return s->Final();
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << argv[0] << " infiles*.root";
    return 1;
  }

  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);

  if (argc > 2) {
    for (int i = 2; i < argc; i++) infilenames[i % NUM_THREADS].push_back(argv[i]);
  } else {
    return 1;
  }

  ROOT::EnableThreadSafety();
  auto start = std::chrono::high_resolution_clock::now();

  std::future<std::shared_ptr<TTree>> threads[NUM_THREADS];
  std::vector<std::shared_ptr<TTree>> skimmed_trees;

  for (size_t i = 0; i < NUM_THREADS; i++) threads[i] = std::async(run_file, infilenames.at(i));
  for (size_t i = 0; i < NUM_THREADS; i++) skimmed_trees.push_back(threads[i].get());

  int i = 0;
  for (auto&& f : skimmed_trees) {
    auto outFile = std::make_unique<TFile>(Form("%s_%d.root", "h10_skim", ++i), "RECREATE");
    outFile->cd();
    f->Write();
    outFile->Write();
    outFile->Close();
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout.imbue(std::locale(""));
  std::cout << RED << elapsed_full.count() << " sec\t";
  std::cout << DEF << std::endl;

  return 0;
}
