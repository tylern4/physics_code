/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include <algorithm>
#include "branches.hpp"
#include "clipp.h"
#include "constants.hpp"
#include "csv_data.hpp"
#include "glob_files.hpp"
#include "physics.hpp"
#include "yeilds.hpp"

using namespace std;

int main(int argc, char **argv) {
  ROOT::EnableThreadSafety();
  std::cout.imbue(std::locale(""));

  std::string e1d_string = "";
  std::string outputfile = "test_mc.csv";
  bool print_help = false;
  std::vector<std::vector<std::string>> infilenames_e1d(NUM_THREADS);

  auto cli = (clipp::option("-h", "--help").set(print_help) % "print help",
              clipp::option("-e1d") & clipp::value("e1d input file path", e1d_string),
              clipp::option("-o") & clipp::value("output filename", outputfile));

  auto good_args = clipp::parse(argc, argv, cli);
  if (!good_args) {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  } else if (e1d_string == "") {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  }
  auto start = std::chrono::high_resolution_clock::now();

  std::vector<std::string> e1d_files = glob(e1d_string);

  int i = 0;
  for (auto &&f : e1d_files) {
    infilenames_e1d[i++ % NUM_THREADS].push_back(f);
  }

  // PRINT_TIMEING(start, "Read in files: ");
  std::string name = outputfile;  //+ "_" + to_string(num)
  auto csv_file = std::make_shared<SyncFile>(name);
  csv_file->write(fid_data_pip::header());
  // Worker function
  auto e1dworker = [&csv_file](auto &&fls, auto &&num) mutable {
    auto dh = std::make_shared<Yeilds>(csv_file);
    size_t total = 0;
    auto data_chain = std::make_shared<TChain>("h10");
    for (auto &&f : fls) {
      data_chain->Add(f.c_str());
    }
    size_t num_events = (size_t)data_chain->GetEntries();
    total += dh->Run_fidFitsPip<e1d_Cuts>(data_chain, num);
    // total += dh->Run_fidFitsPip<Cuts>(data_chain, num);

    return total;
  };

  size_t events = 0;
  if (e1d_files.size() > 0) {
    std::future<size_t> threads_e1d[NUM_THREADS];
    for (size_t i = 0; i < NUM_THREADS; i++) {
      threads_e1d[i] = std::async(e1dworker, infilenames_e1d.at(i), i);
    }

    for (size_t i = 0; i < NUM_THREADS; i++) {
      events += threads_e1d[i].get();
    }
  }
  // Write out to csv file
  csv_file->writeToFile();

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;

  return 0;
}
