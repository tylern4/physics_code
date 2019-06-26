/**************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

// Only My Includes. All others in main.h
#include <future>
#include <thread>
#include "TStopwatch.h"
#include "clipp.h"
#include "glob_files.hpp"
#include "physics.hpp"
#include "skim.hpp"

float run_file(const std::vector<std::string>& in, std::string outfile, int thread_id) {
  if (in.size() < 1) return 0.0;
  float tot = 0;
  outfile = outfile.substr(0, outfile.find(".root"));
  outfile += "_" + std::to_string(thread_id) + ".root";
  auto s = std::make_unique<Skim>(in, outfile);
  tot = s->Basic();
  return tot;
}

int main(int argc, char** argv) {
  short n_threads = 16;
  if (argc > 2) {
    n_threads = atoi(argv[1]);
  } else {
    std::cerr << argv[0] << " n_threads outfile.root infiles*.root";
    return 1;
  }
  std::vector<std::vector<std::string>> infilenames(n_threads);
  std::string outfilename;

  if (argc >= 3) {
    outfilename = argv[2];
    for (int i = 3; i < argc; i++) infilenames[i % n_threads].push_back(argv[i]);
  } else {
    return 1;
  }

  ROOT::EnableThreadSafety();
  auto start = std::chrono::high_resolution_clock::now();
  // auto s = std::make_unique<Skim>(infile, outfile);
  float per = 0;
  std::future<float> threads[n_threads];
  for (size_t i = 0; i < n_threads; i++) {
    threads[i] = std::async(run_file, infilenames.at(i), outfilename, i);
  }
  for (size_t i = 0; i < n_threads; i++) {
    per += threads[i].get();
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout.imbue(std::locale(""));
  std::cout << RED << elapsed_full.count() << " sec\t";
  std::cout << BOLDYELLOW << "(" << (per / n_threads) * 100 << ") % converted" << DEF << std::endl;

  return 0;
}
