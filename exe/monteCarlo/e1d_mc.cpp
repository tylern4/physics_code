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
  std::string file_string = "";
  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
  std::string outfilename = "e1d_mc.root";
  bool print_help = false;

  auto cli = (clipp::option("-h", "--help").set(print_help) % "print help",
              clipp::option("-i") & clipp::value("e1d input file path", file_string),
              clipp::option("-o") & clipp::value("output filename", outfilename));

  auto good_args = clipp::parse(argc, argv, cli);
  if (!good_args) {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  } else if (file_string == "") {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  }

  std::vector<std::string> e1d_files = glob(file_string);

  int i = 0;
  for (auto&& f : e1d_files) {
    infilenames[i++ % NUM_THREADS].push_back(f);
  }

  auto start = std::chrono::high_resolution_clock::now();

  auto hist = std::make_shared<mcHistogram>(outfilename);
  auto mom_corr = nullptr;  // std::make_shared<MomCorr>();
  std::cout.imbue(std::locale(""));
  size_t events = 0;

  std::future<size_t> t[NUM_THREADS];

  for (size_t i = 0; i < NUM_THREADS; i++) {
    t[i] = std::async(run_e1d_file, infilenames.at(i), hist, mom_corr, i);
  }
  for (size_t i = 0; i < NUM_THREADS; i++) {
    events += t[i].get();
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout.imbue(std::locale(""));
  hist->Write();
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;
  return 0;
}
