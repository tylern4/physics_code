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
  std::string e1d_string = "";
  std::string e1f_string = "";
  std::string outputfile = "test_mc.csv";
  bool print_help = false;
  std::vector<std::vector<std::string>> infilenames_e1d(NUM_THREADS);
  std::vector<std::vector<std::string>> infilenames_e1f(NUM_THREADS);

  auto cli = (clipp::option("-h", "--help").set(print_help) % "print help",
              clipp::option("-e1d") & clipp::value("e1d input file path", e1d_string),
              clipp::option("-e1f") & clipp::value("e1f input file path", e1f_string),
              clipp::option("-o") & clipp::value("output filename", outputfile));

  if (!clipp::parse(argc, argv, cli)) {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  } else if ((!(e1d_string == "") && !(e1f_string == ""))) {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  }

  std::vector<std::string> e1d_files = glob(e1d_string);

  auto e1f_files = glob(e1f_string);
  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
  int i = 0;
  for (auto &&f : e1d_files) {
    infilenames[i++ % NUM_THREADS].push_back(f);
  }

  auto start = std::chrono::high_resolution_clock::now();
  std::cout.imbue(std::locale(""));

  size_t events = 0;
  auto csv_file = std::make_shared<SyncFile>(outputfile);
  auto dh = std::make_unique<mcYeilds>(csv_file);

  auto e1dworker = [&dh](auto &&fls) mutable {
    size_t total = 0;
    for (auto &&f : fls) {
      total += dh->RunMC<e1d_Cuts>(f);
      total += dh->Run<e1d_Cuts>(f, "mc_rec");
      std::cout << "  " << total << "\r\r" << std::flush;
    }
    return total;
  };

  std::future<size_t> threads[NUM_THREADS];
  for (size_t i = 0; i < NUM_THREADS; i++) {
    threads[i] = std::async(e1dworker, infilenames.at(i));
  }

  for (size_t i = 0; i < NUM_THREADS; i++) {
    events += threads[i].get();
  }

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;

  return 0;
}
