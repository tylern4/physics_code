/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/
#include "datahandeler.hpp"
#include <algorithm>

std::vector<std::vector<std::string>> DataHandler::DistributeFiles(const std::vector<std::string>& files,
                                                                   size_t nthreads) {
  std::vector<std::vector<std::string>> out(nthreads);
  if (nthreads == 0) return out;
  if (nthreads == 1) {
    out[0] = files;
    return out;
  }

  // Longest-processing-time-first: biggest file to the least loaded thread.
  std::vector<size_t> order(files.size());
  for (size_t i = 0; i < order.size(); i++) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&files](size_t a, size_t b) {
    return FileSize(files[a]) > FileSize(files[b]);
  });

  std::vector<unsigned long long> load(nthreads, 0);
  for (size_t idx : order) {
    size_t least = 0;
    for (size_t t = 1; t < nthreads; t++)
      if (load[t] < load[least]) least = t;
    out[least].push_back(files[idx]);
    load[least] += FileSize(files[idx]);
  }
  return out;
}

DataHandler::DataHandler() = default;

DataHandler::DataHandler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists)
    : _input_files(fin), _hists(hists) {
  _chain = std::make_shared<TChain>("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain);
}

DataHandler::DataHandler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists,
                         const std::shared_ptr<MomCorr>& mom_corr)
    : _input_files(fin), _hists(hists), _mom_corr(mom_corr) {
  _chain = std::make_shared<TChain>("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain);
}

DataHandler::DataHandler(const std::string& fin, const std::shared_ptr<Histogram>& hists,
                         const std::shared_ptr<MomCorr>& mom_corr)
    : _hists(hists), _mom_corr(mom_corr) {
  _chain = std::make_shared<TChain>("h10");
  _chain->Add(fin.c_str());
  _data = std::make_shared<Branches>(_chain);
}

DataHandler::~DataHandler() = default;

void DataHandler::loadbar(long x, long n) {
  int w = 50;
  if (x >= n) return;

  float ratio = x / (float)n;
  int c = ratio * w;

  std::cerr << BLUE << " [";
  for (int x = 0; x < c; x++) std::cerr << GREEN << "=" << DEF;
  std::cerr << GREEN << ">" << DEF;
  for (int x = c; x < w; x++) std::cerr << " ";
  std::cerr << BLUE << (int)(ratio * 100) << "%]\r" << DEF << std::flush;
}

void DataHandler::setLoadBar(bool load) { _loadbar = load; }

mcHandler::mcHandler(const std::vector<std::string>& fin, const std::shared_ptr<mcHistogram>& hists) {
  _input_files = fin;
  _hists = hists;
  _mc_hists = hists;
  _chain = std::make_shared<TChain>("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain, true);
}

mcHandler::mcHandler(const std::vector<std::string>& fin, const std::shared_ptr<mcHistogram>& hists,
                     const std::shared_ptr<MomCorr>& mom_corr) {
  _input_files = fin;
  _hists = hists;
  _mc_hists = hists;
  _mom_corr = mom_corr;
  _chain = std::make_shared<TChain>("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain, true);
}
