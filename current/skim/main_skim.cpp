/**************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

// Only My Includes. All others in main.h
#include "TStopwatch.h"
#include "clipp.h"
#include "glob_files.hpp"
#include "physics.hpp"
#include "skim.hpp"

using namespace std;

int main(int argc, char **argv) {
  bool _mc = false;
  bool _basic = false;
  bool _strict = false;
  bool _final = false;
  bool print_help = false;
  std::vector<std::string> infile;
  std::string outfile = "skim.root";
  auto cli =
      (clipp::option("-h", "--help").set(print_help) % "print help", clipp::option("-mc", "--MC").set(_mc) % "mc Skim",
       clipp::option("-b", "--basic").set(_basic) % "basic skim",
       clipp::option("-s", "--strict").set(_strict) % "strict skim",
       clipp::option("-f", "--final").set(_final) % "final skim",
       clipp::value("skim.root", outfile), clipp::values("inputFile.root", infile));

  clipp::parse(argc, argv, cli);
  if (print_help) {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(-1);
  }

  auto start = std::chrono::high_resolution_clock::now();
  auto s = std::make_unique<Skim>(infile, outfile);
  float per = 0;
  if (_basic)
    per = s->Basic();
  else if (_strict)
    s->Strict();
  else if (_final)
    s->Final();
  else
    per = s->Basic();


  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout.imbue(std::locale(""));
  std::cout << RED << elapsed_full.count() << " sec\t";
  std::cout << BOLDYELLOW << "(" << per*100 << ") % converted" << DEF << std::endl;


  return 0;
}
