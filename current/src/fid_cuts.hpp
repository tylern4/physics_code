#ifndef CUT_FID_H
#define CUT_FID_H

#include "physics.hpp"

extern "C" {
bool fidu_e_sub_(float *pel, float *thetael, float *phiel, bool *in_fid_reg);
bool hadronfid_(bool *result, float *theta, float *phi, int *s);
}

bool fid_e(float p, float cosz, float cosx, float cosy) {
  float phiel_s =
      physics::phi_calc(cosx, cosy) < -30 ? physics::phi_calc(cosx, cosy) + 360 : physics::phi_calc(cosx, cosy);
  float thetael = physics::theta_calc(cosz);
  bool in_fid_reg = true;
  /*
  float pshift = 0.14;
  float c1 = 12.0;
  float c2 = 18.5;
  float c3 = 0.25;
  float c4 = 25.0;
  float factor = 0.416667;


  float thetacut = c1 + c2 / ((p + pshift));
  float expon = c3 * pow(p, factor);

  float del_phie = c4 * pow(sin((thetael - thetacut) * D2R), expon);
  if (abs(phiel_s) <= del_phie && thetael >= thetacut)
    in_fid_reg = true;
  else
    in_fid_reg = false;

  std::cout << fidu_e_sub_(&p, &thetael, &phiel_s, &in_fid_reg) << "\t" << in_fid_reg << std::endl;
*/
  return fidu_e_sub_(&p, &thetael, &phiel_s, &in_fid_reg);
}

bool hadronfid(bool result, float theta, float phi, int s) { return hadronfid_(&result, &theta, &phi, &s); }

#endif
