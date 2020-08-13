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

  if (!clipp::parse(argc, argv, cli)) {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  } else if (e1d_string == "" || e1f_string == "") {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(2);
  }

  auto e1d_files = glob(e1d_string);
  auto e1f_files = glob(e1f_string);

  auto start = std::chrono::high_resolution_clock::now();
  std::cout.imbue(std::locale(""));
  int events = 0;
  int j = 0;
  Yeilds *dh = new Yeilds(outputfile);
  dh->WriteHeader();

  auto e1dworker = [start, events, dh](auto &&f) mutable {
    events += dh->Run<e1d_Cuts>(f, "rec");
    std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
    std::cout << BOLDYELLOW << " " << events / elapsed_full.count() << " Hz\r\r" << DEF << std::flush;
  };

  auto e1fworker = [start, events, dh](auto &&f) mutable {
    events += dh->Run<e1f_Cuts>(f, "rec");
    std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
    std::cout << BOLDYELLOW << " " << events / elapsed_full.count() << " Hz\r\r" << DEF << std::flush;
  };

#ifdef DOCKER
  std::for_each(std::execution::par, e1d_files.begin(), e1d_files.end(), e1dworker);
  std::for_each(std::execution::par, e1f_files.begin(), e1f_files.end(), e1fworker);
#else
  std::for_each(e1d_files.begin(), e1d_files.end(), e1dworker);
  std::for_each(e1f_files.begin(), e1f_files.end(), e1fworker);
#endif

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;

  return 0;
}
