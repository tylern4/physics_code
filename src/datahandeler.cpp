/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/
#include "datahandeler.hpp"

DataHandeler::DataHandeler() = default;

DataHandeler::DataHandeler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists)
    : _input_files(fin), _hists(hists) {
  _chain = std::make_shared<TChain>("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain);
}

DataHandeler::DataHandeler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists,
                           const std::shared_ptr<MomCorr>& mom_corr)
    : _input_files(fin), _hists(hists), _mom_corr(mom_corr) {
  _chain = std::make_shared<TChain>("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain);
}

DataHandeler::~DataHandeler() = default;

void DataHandeler::loadbar(long x, long n) {
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

void DataHandeler::setLoadBar(bool load) { _loadbar = load; }

mcHandeler::mcHandeler(const std::vector<std::string>& fin, const std::shared_ptr<mcHistogram>& hists) {
  _input_files = fin;
  _hists = hists;
  _mc_hists = hists;
  _chain = std::make_shared<TChain>("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain, true);
}

mcHandeler::mcHandeler(const std::vector<std::string>& fin, const std::shared_ptr<mcHistogram>& hists,
                       const std::shared_ptr<MomCorr>& mom_corr) {
  _input_files = fin;
  _hists = hists;
  _mc_hists = hists;
  _mom_corr = mom_corr;
  _chain = std::make_shared<TChain>("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain, true);
}
