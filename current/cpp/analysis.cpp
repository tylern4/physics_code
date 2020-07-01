/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include <algorithm>
#include <future>
#include <iostream>
#include <thread>
#include <vector>
#include "TStopwatch.h"
#include "branches.hpp"
#include "clipp.h"
#include "constants.hpp"
#include "datahandeler.hpp"
#include "glob_files.hpp"
#include "main.h"
#include "physics.hpp"

int main(int argc, char** argv) {
  ROOT::EnableThreadSafety();

  std::string e1d_string = "";
  std::string e1f_string = "";
  std::string output_file = "";
  std::vector<std::vector<std::string>> infilenames_e1d(NUM_THREADS);
  std::vector<std::vector<std::string>> infilenames_e1f(NUM_THREADS);

  auto cli = (clipp::option("-e1d") & clipp::value("e1d input file path", e1d_string),
              clipp::option("-e1f") & clipp::value("e1f input file path", e1f_string),
              clipp::option("-o") & clipp::value("output filename", output_file));

  if (!clipp::parse(argc, argv, cli)) {
    std::cout << clipp::make_man_page(cli, argv[0]);
  }

  auto e1d_files = glob(e1d_string);
  auto e1f_files = glob(e1f_string);

  int i = 0;
  for (auto&& f : e1d_files) {
    infilenames_e1d[i++ % NUM_THREADS].push_back(f);
  }
  i = 0;
  for (auto&& f : e1f_files) {
    infilenames_e1f[i++ % NUM_THREADS].push_back(f);
  }

  // Create the shared objects for hists and momentum corrections
  auto hist = std::make_shared<Histogram>();
  auto mom_corr = std::make_shared<MomCorr>();

  std::cout.imbue(std::locale(""));
  size_t events = 0;

  // Make an array of futures which return an integer
  std::future<size_t> threads_e1d[NUM_THREADS];

  auto start = std::chrono::high_resolution_clock::now();
  // For the number of threads run the list of file
  for (size_t i = 0; i < NUM_THREADS; i++) {
    threads_e1d[i] = std::async(run_e1d_file, infilenames_e1d.at(i), hist, mom_corr, i);
  }
  // Let the threads run and get the results when they're done
  for (size_t i = 0; i < NUM_THREADS; i++) {
    events += threads_e1d[i].get();
  }

  // Make an array of futures which return an integer
  std::future<size_t> threads_e1f[NUM_THREADS];

  // For the number of threads run the list of file
  for (size_t i = 0; i < NUM_THREADS; i++) {
    threads_e1f[i] = std::async(run_e1f_file, infilenames_e1f.at(i), hist, mom_corr, i);
  }
  // Let the threads run and get the results when they're done
  for (size_t i = 0; i < NUM_THREADS; i++) {
    events += threads_e1f[i].get();
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);

  std::cout.imbue(std::locale(""));
  ROOT::EnableImplicitMT(2);
  // Write the histograms to file
  hist->Write(output_file);

  std::cout << RED << events << " events\t" << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;
  return 0;
}
