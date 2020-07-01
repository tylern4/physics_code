/**************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

// Only My Includes. All others in main.h
#include <algorithm>
#include <future>
#include <iostream>
#include <thread>
#include "TFile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "glob_files.hpp"
#include "physics.hpp"
#include "skim.hpp"

struct fileNames {
  std::string folder;
  int run_num;
  friend std::ostream& operator<<(std::ostream& os, const fileNames& fn) {
    os << fn.folder << " " << fn.run_num;
    return os;
  }
};

std::vector<fileNames> getFileNames(const std::vector<std::string>& inputNames) {
  std::vector<int> run_nums;
  int _temp_run = 0;
  std::string folder;
  for (auto&& in : inputNames) {
    int folder_end = in.find_last_of("/") + 1;
    folder = in.substr(0, folder_end - 1);
    auto name = in.substr(folder_end, in.length());
    int run_start = name.find_first_of("_") + 1;
    int run_num = stol(name.substr(run_start, run_start + 7));
    run_nums.push_back(run_num);
  }

  auto ip = std::unique(run_nums.begin(), run_nums.end());
  run_nums.resize(std::distance(run_nums.begin(), ip));
  std::sort(run_nums.begin(), run_nums.end());

  std::vector<fileNames> x;
  for (auto&& r : run_nums) {
    fileNames z;
    z.run_num = r;
    z.folder = folder;
    x.push_back(z);
  }

  return x;
}

struct tree_run {
  std::shared_ptr<TTree> tree;
  int run_num;
  std::string folder;
};

tree_run run_file(const fileNames& in) {
  std::shared_ptr<TChain> h10 = std::make_shared<TChain>("h10");
  h10->Add(Form("%s/run_%d_pass1.a*.root", in.folder.c_str(), in.run_num));
  auto s = std::make_unique<Skim>(h10);
  tree_run out;
  out.tree = s->Basic<e1f_Cuts>();
  out.run_num = in.run_num;
  out.folder = in.folder;
  return out;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << argv[0] << " infiles*.root";
    exit(1);
  }

  std::vector<std::string> inputs;
  if (argc > 2) {
    for (int i = 2; i < argc; i++) inputs.push_back(argv[i]);
  }
  auto file_names = getFileNames(inputs);

  auto number_of_threads = std::min(size_t(NUM_THREADS), file_names.size());

  ROOT::EnableThreadSafety();
  auto start = std::chrono::high_resolution_clock::now();

  std::future<tree_run> threads[number_of_threads];
  std::vector<tree_run> skimmed_trees;

  int total_num = 0;
  for (auto&& fn : file_names) {
    std::cout << "Start: " << fn.run_num << "\tThread: " << (total_num % number_of_threads) << std::endl;
    threads[total_num++ % number_of_threads] = std::async(run_file, fn);

    if (total_num % number_of_threads == 0) {
      for (size_t i = 0; i < number_of_threads; i++) {
        auto f = threads[i].get();
        std::cout << Form("%s/skim/%s_%d.root", f.folder.c_str(), "e1f_skim", f.run_num) << std::endl;
        auto outFile =
            std::make_unique<TFile>(Form("%s/skim/%s_%d.root", f.folder.c_str(), "e1f_skim", f.run_num), "RECREATE");
        outFile->cd();
        f.tree->Write();
        outFile->Write();
        outFile->Close();
      }
    }
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout.imbue(std::locale(""));
  std::cout << RED << elapsed_full.count() << " sec\t";
  std::cout << DEF << std::endl;

  return 0;
}
