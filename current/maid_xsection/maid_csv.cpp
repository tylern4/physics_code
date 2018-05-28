#include "maid.h"
#include <exception>
#include <iostream>

int main() {
  float W_bins[] = {1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
  float Q2_bins[] = {0.00001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
  double dsigma = 0;

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 10; j++) {
      dsigma = maid_dsigma( // beam_energy
          5.5,
          // W
          W_bins[i],
          // Q^2
          Q2_bins[j],
          // cos(theta*)
          0.5,
          // phi* (degrees)
          1,
          // helicity
          0,
          // model_opt
          5,
          // channel_opt
          3,
          // resonance_opt
          0);
      std::cout << dsigma << ",";
    }
    std::cout << std::endl;
  }
}
