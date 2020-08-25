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
  std::string outputfile = "test.csv";
  bool print_help = false;
  std::vector<std::vector<std::string>> infilenames_e1d(NUM_THREADS);
  std::vector<std::vector<std::string>> infilenames_e1f(NUM_THREADS);

  auto cli = (clipp::option("-h", "--help").set(print_help) % "print help",
              clipp::option("-e1d") & clipp::value("e1d input file path", e1d_string),
              clipp::option("-e1f") & clipp::value("e1f input file path", e1f_string),
              clipp::option("-o") & clipp::value("output filename", outputfile));
  auto good_args = clipp::parse(argc, argv, cli);
  if (!good_args) {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  } else if ((e1d_string == "") && (e1f_string == "")) {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  }

  std::vector<std::string> e1d_files = glob(e1d_string);
  std::vector<std::string> e1f_files = glob(e1f_string);

  int i = 0;
  for (auto &&f : e1d_files) {
    infilenames_e1d[i++ % NUM_THREADS].push_back(f);
  }

  i = 0;
  for (auto &&f : e1f_files) {
    infilenames_e1f[i++ % NUM_THREADS].push_back(f);
  }

  auto start = std::chrono::high_resolution_clock::now();
  std::cout.imbue(std::locale(""));

  size_t events = 0;

  auto e1dworker = [&outputfile](auto &&fls, auto &&num) mutable {
    std::string name = outputfile + "_" + to_string(num) + "_e1d.csv";
    auto csv_file = std::make_shared<SyncFile>(name);
    auto dh = std::make_unique<Yeilds>(csv_file);
    size_t total = 0;
    for (auto &&f : fls) {
      total += dh->Run<e1d_Cuts>(f, "rec");
      std::cout << "  " << total << "\r\r" << std::flush;
    }
    return total;
  };

  auto e1fworker = [&outputfile](auto &&fls, auto &&num) mutable {
    std::string name = outputfile + "_" + to_string(num) + "_e1f.csv";
    auto csv_file = std::make_shared<SyncFile>(name);
    auto dh = std::make_unique<Yeilds>(csv_file);
    size_t total = 0;
    for (auto &&f : fls) {
      total += dh->Run<e1f_Cuts>(f, "rec");
      std::cout << "  " << total << "\r\r" << std::flush;
    }
    return total;
  };

  if (e1d_files.size() > 0) {
    std::future<size_t> threads[NUM_THREADS];
    for (size_t i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::async(e1dworker, infilenames_e1d.at(i), i);
    }

    for (size_t i = 0; i < NUM_THREADS; i++) {
      events += threads[i].get();
    }
  }

  if (e1f_files.size() > 0) {
    std::future<size_t> threads_e1f[NUM_THREADS];
    for (size_t i = 0; i < NUM_THREADS; i++) {
      threads_e1f[i] = std::async(e1fworker, infilenames_e1f.at(i), i);
    }

    for (size_t i = 0; i < NUM_THREADS; i++) {
      events += threads_e1f[i].get();
    }
  }
  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;

  return 0;
}
