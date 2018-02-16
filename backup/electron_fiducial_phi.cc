/*
  These two functions provide the lower and upper bounds for a
  fiducial cut on the electron's lab frame phi angle given the theta
  angle and momentum.
*/

double fiducial_phi_lo(double theta_e, double theta_e_min, double k, double m) {
  return -37.14 * pow(sin((theta_e - theta_e_min) * 0.01745),
                      k + m / theta_e + 1500. / (theta_e * theta_e));
}

double fiducial_phi_hi(double theta_e, double theta_e_min, double k, double m) {
  return 37.14 * pow(sin((theta_e - theta_e_min) * 0.01745),
                     k + m / theta_e + 1500. / (theta_e * theta_e));
}

/*
  If you define the electron's lab frame 4-momentum as the "emom"
  variable, then you can do the following in your event loop:

  NOTE: Assumes you have defined the usual h10 variables using their
  C++ names, such as cc_segm
*/
{
  double e_E = e_mu_prime.Energy();
  double costheta_e = e_mu_prime.CosTheta();
  double theta_e = e_mu_prime.Theta();
  double phi_e = e_mu_prime.Phi();
  double sector_width = 3.141592653589793 / 3.0;
  int esector_index;
  for (esector_index = -3; (esector_index < 4); ++esector_index) {
    if ((-0.5 * sector_width + esector_index * sector_width) <= phi_e &&
        phi_e <= (0.5 * sector_width + esector_index * sector_width)) {
      break;
    }
  }

  int esector;
  if ((esector_index < 0)) {
    esector = (esector_index + 7);
  } else {
    esector = (esector_index + 1);
  }
  // phi angle relative to the sector:
  double phi_e_rel = phi_e - sector_width * esector_index;
  double efid_theta_min = 9.5 + 17.0 / (e_mu_prime.P() + 0.17);
  double efid_k = 0.705 + 1.1 * e_mu_prime.P();
  double efid_m = -63.5 + (-30.0 * e_mu_prime.P());

  /*
    and then the functions can be used to cut like this, assuming you
    use radians for your theta and phi angles:
  */

  double efid_phi_hi =
      deg_to_rad *
      fiducial_phi_hi(theta_e / D2R, efid_theta_min, efid_k, efid_m);
  double efid_phi_lo =
      deg_to_rad *
      fiducial_phi_lo(theta_e / D2R, efid_theta_min, efid_k, efid_m);

  bool efid_passes_cut = efid_phi_lo <= phi_e_rel && phi_e_rel <= efid_phi_hi;
}
