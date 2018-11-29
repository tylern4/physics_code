/**************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

// Only My Includes. All others in main.h
#include "../src/classes.hpp"
#include "../src/clipp.h"
#include "../src/glob_files.hpp"
#include "../src/physics.hpp"
#include "../src/skim.hpp"
#include "TStopwatch.h"

using namespace std;

int main(int argc, char **argv) {
  bool _mc = false;
  bool _basic = false;
  bool _strict = false;
  bool _final = false;
  bool print_help = false;
  std::vector<std::string> infile;
  std::string outfile;
  auto cli =
      (clipp::option("-h", "--help").set(print_help) % "print help", clipp::option("-mc", "--MC").set(_mc) % "mc Skim",
       clipp::option("-b", "--basic").set(_basic) % "basic skim",
       clipp::option("-s", "--strict").set(_strict) % "strict skim",
       clipp::option("-f", "--final").set(_mc) % "final skim", clipp::value("skim.root", outfile),
       clipp::values("inputFile.root", infile));

  clipp::parse(argc, argv, cli);
  if (print_help) {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(-1);
  }

  Skim *s = new Skim(infile, outfile);

  if (_basic)
    s->Basic();
  else if (_strict)
    s->Strict();
  else if (_final)
    s->Final();
  else
    s->Basic();

  return 0;
}
