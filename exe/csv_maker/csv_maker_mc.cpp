/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include <algorithm>
#include "branches.hpp"
#include "clipp.h"
#include "constants.hpp"
#include "glob_files.hpp"
#include "physics.hpp"
#include "yeilds.hpp"

using namespace std;

int main(int argc, char **argv) {
  ROOT::EnableThreadSafety();
  std::cout.imbue(std::locale(""));

  std::string e1d_string = "";
  std::string e1f_string = "";
  std::string outputfile = "test_mc.csv";
  size_t num_cut = 0;
  bool print_help = false;
  std::vector<std::vector<std::string>> infilenames_e1d(NUM_THREADS);
  std::vector<std::vector<std::string>> infilenames_e1f(NUM_THREADS);

  auto cli = (clipp::option("-h", "--help").set(print_help) % "print help",
              clipp::option("-e1d") & clipp::value("e1d input file path", e1d_string),
              clipp::option("-e1f") & clipp::value("e1f input file path", e1f_string),
              clipp::option("-cut") & clipp::value("Cut the number of files", num_cut),
              clipp::option("-o") & clipp::value("output filename", outputfile));

  auto good_args = clipp::parse(argc, argv, cli);
  if (!good_args) {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  } else if ((e1d_string == "") && (e1f_string == "")) {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  }
  auto start = std::chrono::high_resolution_clock::now();

  std::vector<std::string> e1d_files = glob(e1d_string);
  std::vector<std::string> e1f_files = glob(e1f_string);
  if (num_cut == 0) num_cut = e1d_files.size();
  if (num_cut == 0) num_cut = e1f_files.size();
  if (num_cut == 0) exit(4);

  int i = 0;
  for (auto &&f : e1d_files) {
    if (i <= num_cut) infilenames_e1d[i++ % NUM_THREADS].push_back(f);
  }

  i = 0;
  for (auto &&f : e1f_files) {
    if (i <= num_cut) infilenames_e1f[i++ % NUM_THREADS].push_back(f);
  }

  PRINT_TIMEING(start, "Read in files: ");
  std::string name = outputfile + "_e1d.csv";  //+ "_" + to_string(num)
  auto csv_file = std::make_shared<SyncFile>(name);

  auto e1dworker = [&csv_file](auto &&fls, auto &&num) mutable {
    auto dh = std::make_shared<mcYeilds>(csv_file);
    size_t total = 0;

    auto data_chain = std::make_shared<TChain>("h10");
    for (auto &&f : fls) {
      data_chain->Add(f.c_str());
    }

    size_t num_events = (size_t)data_chain->GetEntries();

    total += dh->Run<e1d_Cuts>(data_chain, "mc_rec", num);

    auto mc_chain = std::make_shared<TChain>("h10");
    for (auto &&f : fls) {
      mc_chain->Add(f.c_str());
    }

    total += dh->RunMC<e1d_Cuts>(mc_chain, num);

    return total;
  };

  // auto e1fworker = [&outputfile](auto &&fls, auto &&num) mutable {
  //   std::string name = outputfile + "_" + to_string(num) + "_e1f.csv";
  //   auto csv_file = std::make_shared<SyncFile>(name);
  //   auto dh = std::make_shared<mcYeilds>(csv_file);
  //   size_t total = 0;
  //   auto chain = std::make_shared<TChain>("h10");
  //   for (auto &&f : fls) {
  //     chain->Add(f.c_str());
  //   }
  //   total += dh->RunMC<e1f_Cuts>(chain, num);
  //   total += dh->Run<e1f_Cuts>(chain, "mc_rec", num);

  //   return total;
  // };

  PRINT_TIMEING(start, "Make lambdas: ");

  size_t events = 0;
  if (e1d_files.size() > 0) {
    std::future<size_t> threads_e1d[NUM_THREADS];
    for (size_t i = 0; i < NUM_THREADS; i++) {
      PRINT_TIMEING(start, "Make thread " << i << ": ");
      threads_e1d[i] = std::async(e1dworker, infilenames_e1d.at(i), i);
    }

    for (size_t i = 0; i < NUM_THREADS; i++) {
      PRINT_TIMEING(start, "Get thread " << i << ": ");
      events += threads_e1d[i].get();
    }
  }
  // Write out to csv file
  csv_file->writeToFile();

  // if (e1f_files.size() > 0) {
  //   std::future<size_t> threads_e1f[NUM_THREADS];
  //   for (size_t i = 0; i < NUM_THREADS; i++) {
  //     threads_e1f[i] = std::async(e1fworker, infilenames_e1f.at(i), i);
  //   }

  //   for (size_t i = 0; i < NUM_THREADS; i++) {
  //     events += threads_e1f[i].get();
  //   }
  // }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;

  return 0;
}
