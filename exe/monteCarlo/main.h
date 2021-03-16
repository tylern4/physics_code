/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef MAIN_H
#define MAIN_H
#include <iostream>
#include "color.hpp"
#include "datahandeler.hpp"
#include "histogram.hpp"

#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)
#define Square(a) (a * a)

void loadbar(long x, long n) {
  int w = 50;
  if ((x != n) && (x % (n / 100 + 1) != 0)) return;

  double ratio = x / (double)n;
  int c = ratio * w;

  std::cout << BLUE << " [";
  for (int x = 0; x < c; x++) std::cout << GREEN << "=" << DEF;
  std::cout << GREEN << ">" << DEF;
  for (int x = c; x < w; x++) std::cout << " ";
  std::cout << BLUE << (int)(ratio * 100) << "%]\r" << DEF << std::flush;
}

size_t run_e1d_file(const std::vector<std::string>& in, const std::shared_ptr<mcHistogram>& hists,
                    const std::shared_ptr<MomCorr>& mom_corr, int thread_id) {
  auto dh = std::make_unique<mcHandler>(in, hists, mom_corr);
  dh->setLoadBar(false);
  if (thread_id == 0) dh->setLoadBar(true);
  size_t tot = 0;
  tot += dh->Run<e1d_mcCuts>();
  return tot;
}

size_t run_e1f_file(const std::vector<std::string>& in, const std::shared_ptr<mcHistogram>& hists,
                    const std::shared_ptr<MomCorr>& mom_corr, int thread_id) {
  auto dh = std::make_unique<mcHandler>(in, hists, mom_corr);
  dh->setLoadBar(false);
  if (thread_id == 0) dh->setLoadBar(true);
  size_t tot = 0;
  tot += dh->Run<e1f_Cuts>();
  return tot;
}

#endif
