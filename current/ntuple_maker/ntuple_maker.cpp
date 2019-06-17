/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include <algorithm>
#include <iostream>
#include <thread>
#include <vector>
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"
#include "yeilds.hpp"

using namespace std;
size_t run_file(std::string in, int thread_id) {
  std::string outf = in.substr(in.find("h10_"));
  outf = "out/ntuple_" + std::to_string(thread_id) + "_" + outf;
  auto chain = std::make_unique<TChain>("h10");
  chain->Add(in.c_str());
  auto dh = std::make_unique<Yeilds>(outf.c_str(), true);
  return dh->RunNtuple(std::move(chain));
}

void call_from_thread(std::string in) { std::cout << in << std::endl; }

int main(int argc, char** argv) {
  ROOT::EnableThreadSafety();
  std::vector<std::string> infilename;
  std::string outfilename;
  if (argc >= 2) {
    outfilename = argv[1];
    for (int i = 2; i < argc; i++) infilename.push_back(argv[i]);
  } else {
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();
  std::cout.imbue(std::locale(""));
  size_t events = 0;
  int j = 0;

for (size_t i = 0; i < count; i++) {

}
  std::thread t[NUM_THREADS];
  for (size_t i = 0; i < NUM_THREADS; i++) {
    t[i] = std::thread(run_file, infilename.at(i), i);
  }
  for (size_t i = 0; i < NUM_THREADS; i++) {
    t[i].join();
  }

  // events += run_file(in);

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << "\n\n" << events / elapsed_full.count() << " Hz" << DEF << std::endl;

  return 0;
}
