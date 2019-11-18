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

struct fileNames {
  std::string folder;
  int run_num;
  std::vector<std::string> file_num;
};

std::vector<fileNames> getFileNames(const std::vector<std::string>& inputNames) {
  std::vector<int> run_nums;
  int _temp_run = 0;

  for (auto&& in : inputNames) {
    int folder_end = in.find_last_of("/") + 1;
    int run_num = stol(in.substr(folder_end + 5, 5));
    run_nums.push_back(run_num);
  }

  auto ip = std::unique(run_nums.begin(), run_nums.end());
  run_nums.resize(std::distance(run_nums.begin(), ip));
  std::sort(run_nums.begin(), run_nums.end());

  std::vector<fileNames> x;
  for (auto&& r : run_nums) {
    fileNames z;
    for (auto&& in : inputNames) {
      int folder_end = in.find_last_of("/") + 1;
      int run_end = in.find_last_of("_") + 1;
      int run_num = stoi(in.substr(folder_end + 5, 5));
      if (r != run_num) continue;
      z.run_num = stol(in.substr(folder_end + 5, 5));
      z.folder = in.substr(0, folder_end);
      z.file_num.push_back(in.substr(run_end, 2));
    }
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
  for (auto& f : in.file_num) h10->Add(Form("%s/h10_r%d_%s.root", in.folder.c_str(), in.run_num, f.c_str()));
  auto s = std::make_unique<Skim>(h10);
  tree_run out;
  out.tree = s->Basic();
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

  ROOT::EnableThreadSafety();
  auto start = std::chrono::high_resolution_clock::now();

  std::future<tree_run> threads[NUM_THREADS];
  std::vector<tree_run> skimmed_trees;

  int total_num = 0;
  for (auto&& fn : file_names) {
    std::cout << "Start: " << fn.run_num << "\tThread: " << (total_num % NUM_THREADS) << std::endl;
    threads[total_num % NUM_THREADS] = std::async(run_file, fn);
    total_num++;
  }

  for (size_t i = 0; i < total_num; i++) {
    skimmed_trees.push_back(threads[i].get());
  }

  int i = 0;
  for (auto&& f : skimmed_trees) {
    TThread::Lock();
    std::cout << Form("%sskim/%s_%d.root", f.folder.c_str(), "h10_skim", f.run_num) << std::endl;
    auto outFile =
        std::make_unique<TFile>(Form("%sskim/%s_%d.root", f.folder.c_str(), "h10_skim", f.run_num), "RECREATE");
    outFile->cd();
    f.tree->Write();
    outFile->Write();
    outFile->Close();
    TThread::UnLock();
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout.imbue(std::locale(""));
  std::cout << RED << elapsed_full.count() << " sec\t";
  std::cout << DEF << std::endl;

  return 0;
}
