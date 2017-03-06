/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include <iostream>
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "../src/constants.hpp"
#include "../src/physics.hpp"
#include "../src/classes.hpp"
#include <TLorentzVector.h>

class H10 {
private:
  double Square(double a) { return a * a; }

public:
  std::vector<double> W_vec;
  std::vector<double> Q2_vec;
  std::vector<float> p_vec;
  std::vector<float> b_vec;
  // W and Q^2
  int bins = 500;
  double w_min = 0;
  double w_max = 3.25;
  double q2_min = 0;
  double q2_max = 10;

  // TH2D *WvsQ2_hist = new TH2D("WvsQ2_hist", "W vs Q^{2}", bins, w_min, w_max,
  //                            bins, q2_min, q2_max);

  H10() {}

  ~H10() {}
  int num_of_events;
  bool electron_cuts;

  void pass_chain(TTree &chain) {
    getBranches(&chain);
    int num_of_events = (int)chain.GetEntries();
    std::cout << num_of_events << std::endl;
  }

  Histogram dataHandeler(TTree &chain) {
    bool first_run = true;
    // My Classes
    Histogram hists;
    // From missing_mass.hpp :: missing_mass_calc()
    MissingMass *MM_neutron = new MissingMass();
    MM_neutron->Set_target_mass(MASS_P);
    MM_neutron->Set_target_PxPyPz(0);

    TCanvas *c1 = new TCanvas("c1", "c1", 100, 100);
    auto number_cols = 0;
    char rootFile[500];
    int num_of_events, total_events;
    bool cuts, electron_cuts;
    double e_E;

    std::vector<bool> pip_vec(MAX_PARTS, false);
    std::vector<bool> pim_vec(MAX_PARTS, false);
    std::vector<bool> proton_vec(MAX_PARTS, false);
    std::vector<bool> elec_vec(MAX_PARTS, false);

    int num_of_pis, num_of_proton;

    TVector3 e_mu_prime_3;
    TLorentzVector e_mu_prime;
    TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)),
                        E1D_E0);

    TVector3 Particle3(0.0, 0.0, 0.0);
    TLorentzVector Particle4(0.0, 0.0, 0.0, 0.0);
    double theta, phi;
    int sector;

    getBranches(&chain);
    if (!first_run)
      getMorebranchs(&chain);

    num_of_events = (int)chain.GetEntries();

    //#pragma omp parallel for
    for (int current_event = 0; current_event < num_of_events;
         current_event++) {

      chain.GetEntry(current_event);

      // reset electron cut bool
      electron_cuts = true;
      // electron cuts
      electron_cuts &= (ec[0] > 0); // ``` ``` ``` ec
      if (electron_cuts)
        hists.EC_fill(etot[ec[0] - 1], p[0]);
      electron_cuts &= ((int)id[0] == ELECTRON); // First particle is electron
      electron_cuts &=
          ((int)gpart > 0); // Number of good particles is greater than 0
      // electron_cuts &= ((int)stat[0] > 0); // First Particle hit stat
      // std::cout << "stat " << (int)stat[0] << std::endl;
      electron_cuts &= ((int)q[0] == -1); // First particle is negative Q
      electron_cuts &= ((int)sc[0] > 0);  // First Particle hit sc
      electron_cuts &= ((int)dc[0] > 0);  // ``` ``` ``` dc
      electron_cuts &= ((int)dc_stat[dc[0] - 1] > 0);

      if (electron_cuts && cc[0] > 0) {
        int cc_sector = cc_sect[cc[0] - 1];
        int cc_segment = (cc_segm[0] % 1000) / 10;
        int cc_pmt = cc_segm[0] / 1000 - 1;
        int cc_nphe = nphe[cc[0] - 1];
        hists.CC_fill(cc_sector, cc_segment, cc_pmt, cc_nphe);
      }

      if (electron_cuts) {
        if (first_run) {
          is_electron = &elec_vec;
          is_electron->at(0) = true;
          for (int part_num = 1; part_num < gpart; part_num++) {
            is_pip = &pip_vec;
            is_pim = &pim_vec;
            is_proton = &proton_vec;
            is_pip->at(part_num) = (id[part_num] == PIP);
            is_proton->at(part_num) = (id[part_num] == PROTON);
            is_pim->at(part_num) = (id[part_num] == PIM);
          }
        }

        // Setup scattered electron 4 vector
        e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
        e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
        // Set the vertex time (time of electron hit)
        Delta_T *delta_t = new Delta_T();
        delta_t->delta_t_cut(&hists, first_run);

        theta = theta_calc(cz[0]);
        // phi = center_phi_calc(cx[0],cy[0]);
        phi = phi_calc(cx[0], cy[0]);
        sector = get_sector(phi);

        // Fill_fid(theta,phi,get_sector(phi_calc(cx[0],cy[0])));

        hists.Fill_fid(theta, phi, sector);

        if (first_run) {
          W = W_calc(e_mu, e_mu_prime);
          Q2 = Q2_calc(e_mu, e_mu_prime);
          e_E = e_mu_prime.E();
        }

        hists.WvsQ2_Fill(e_E, W, Q2, xb_calc(Q2, e_E));
        num_of_proton = num_of_pis = 0;

        //#pragma omp parallel for
        for (int part_num = 1; part_num < gpart; part_num++) {
          if (p[part_num] == 0)
            continue;
          // if(is_proton->at(part_num) == is_pip->at(part_num)) continue;

          hists.Fill_Mass(m[part_num]);
          Particle3.SetXYZ(p[part_num] * cx[part_num],
                           p[part_num] * cy[part_num],
                           p[part_num] * cz[part_num]);
          Particle4.SetVectM(Particle3, Get_Mass(id[part_num]));

          hists.MomVsBeta_Fill(Particle4.E(), p[part_num], b[part_num]);
          if (q[part_num] == 1) {
            hists.MomVsBeta_Fill_pos(p[part_num], b[part_num]);
            if (is_proton->at(part_num) && (id[part_num] == PROTON)) {
              num_of_proton++;
              hists.Fill_proton_WQ2(W, Q2);
              hists.Fill_proton_ID_P(p[part_num], b[part_num]);
            }
            if (is_pip->at(part_num) && (id[part_num] == PIP)) {
              num_of_pis++;
              hists.Fill_pion_WQ2(W, Q2);
              hists.Fill_Pi_ID_P(p[part_num], b[part_num]);
              TLorentzVector gamma_mu = (e_mu - e_mu_prime);
              if (first_run) {
                MM_neutron->Set_PxPyPz(p[part_num] * cx[part_num],
                                       p[part_num] * cy[part_num],
                                       p[part_num] * cz[part_num]);
                MM = MM_neutron->missing_mass(gamma_mu);
              }
              hists.Fill_Missing_Mass(MM);
              hists.Fill_Missing_Mass_square(Square(MM));
            }

            if ((is_pip->at(part_num) && (id[part_num] == PIP)) ||
                (is_proton->at(part_num) && (id[part_num] == PROTON))) {
              hists.Fill_proton_Pi_ID_P(p[part_num], b[part_num]);
            }
          } else if (q[part_num] == -1) {
            hists.MomVsBeta_Fill_neg(p[part_num], b[part_num]);
          }
        }

        if (num_of_pis == 1)
          hists.Fill_single_pi_WQ2(W, Q2);
        if (num_of_proton == 1)
          hists.Fill_single_proton_WQ2(W, Q2);
      }
    }

    chain.Reset();
    return hists;
  }

  void deltat_slicefit(TH2D delta_t_mass_P, TH2D delta_t_mass_PIP) {
    Header *fit_functions = new Header("../src/fit_functions.hpp", "FF");
    TCanvas *c1 = new TCanvas;
    TF1 *peak = new TF1("peak", "gaus", -1, 1);
    //[0]*exp(-[1]*x) +
    std::string func = "[0]*exp(-[1]*x) + [2]*x*x + [3]*x + [4]";
    delta_t_mass_P.FitSlicesY(peak, 0, -1, 10, "QRG5");
    TH1D *delta_t_mass_P_0 = (TH1D *)gDirectory->Get("delta_t_mass_P_0");
    TH1D *delta_t_mass_P_1 = (TH1D *)gDirectory->Get("delta_t_mass_P_1");
    TH1D *delta_t_mass_P_2 = (TH1D *)gDirectory->Get("delta_t_mass_P_2");
    double x[500];
    double y_plus[500];
    double y_minus[500];
    int num = 0;
    for (int i = 0; i < 500; i++) {
      if (delta_t_mass_P_1->GetBinContent(i) != 0) {
        // Get momentum from bin center
        x[num] = (double)delta_t_mass_P_1->GetBinCenter(i);
        // mean + 3sigma
        y_plus[num] = (double)delta_t_mass_P_1->GetBinContent(i) +
                      N_SIGMA * (double)delta_t_mass_P_2->GetBinContent(i);
        // mean - 3simga
        y_minus[num] = (double)delta_t_mass_P_1->GetBinContent(i) -
                       N_SIGMA * (double)delta_t_mass_P_2->GetBinContent(i);
        num++;
      }
    }

    TGraph *P = new TGraph(num, x, y_plus);
    TGraph *M = new TGraph(num, x, y_minus);
    TF1 *Proton_Pos_fit = new TF1("Proton_Pos_fit", func.c_str());
    TF1 *Proton_Neg_fit = new TF1("Proton_Neg_fit", func.c_str());
    P->Fit(Proton_Pos_fit, "QRG5", "", 0.2, 2);
    P->Write();
    M->Fit(Proton_Neg_fit, "QRG5", "", 0.2, 2);
    M->Write();
    Proton_Pos_fit->Write();
    Proton_Neg_fit->Write();
    P->Draw("Same");
    M->Draw("Same");
    Proton_Pos_fit->Draw("Same");
    Proton_Neg_fit->Draw("Same");

    fit_functions->NewFunction();
    fit_functions->Set_RetrunType("double");
    fit_functions->Set_FuncName("Proton_Pos_fit");
    fit_functions->Set_FuncInputs("double x");
    fit_functions->Set_Function(Proton_Pos_fit->GetExpFormula("P"));
    fit_functions->WriteFunction();

    fit_functions->NewFunction();
    fit_functions->Set_RetrunType("double");
    fit_functions->Set_FuncName("Proton_Neg_fit");
    fit_functions->Set_FuncInputs("double x");
    fit_functions->Set_Function(Proton_Neg_fit->GetExpFormula("P"));
    fit_functions->WriteFunction();

    delta_t_mass_PIP.FitSlicesY(peak, 0, -1, 10, "QRG5");
    TH1D *delta_t_mass_PIP_0 = (TH1D *)gDirectory->Get("delta_t_mass_PIP_0");
    TH1D *delta_t_mass_PIP_1 = (TH1D *)gDirectory->Get("delta_t_mass_PIP_1");
    TH1D *delta_t_mass_PIP_2 = (TH1D *)gDirectory->Get("delta_t_mass_PIP_2");
    double x_pip[500];
    double y_plus_pip[500];
    double y_minus_pip[500];
    num = 0;
    for (int i = 0; i < 500; i++) {
      if (delta_t_mass_PIP_1->GetBinContent(i) != 0) {
        // Get momentum from bin center
        x_pip[num] = (double)delta_t_mass_PIP_1->GetBinCenter(i);
        // mean + 3sigma
        y_plus_pip[num] =
            (double)delta_t_mass_PIP_1->GetBinContent(i) +
            N_SIGMA * (double)delta_t_mass_PIP_2->GetBinContent(i);
        // mean - 3simga
        y_minus_pip[num] =
            (double)delta_t_mass_PIP_1->GetBinContent(i) -
            N_SIGMA * (double)delta_t_mass_PIP_2->GetBinContent(i);
        num++;
      }
    }

    TGraph *P_pip = new TGraph(num, x_pip, y_plus_pip);
    TGraph *M_pip = new TGraph(num, x_pip, y_minus_pip);
    TF1 *Pip_Pos_fit = new TF1("Pip_Pos_fit", func.c_str());
    TF1 *Pip_Neg_fit = new TF1("Pip_Neg_fit", func.c_str());
    P_pip->Fit(Pip_Pos_fit, "QRG5", "", 0.1, 1.75);
    P_pip->Write();
    M_pip->Fit(Pip_Neg_fit, "QRG5", "", 0.1, 1.75);
    M_pip->Write();
    Pip_Pos_fit->Write();
    Pip_Neg_fit->Write();
    P_pip->Draw("Same");
    M_pip->Draw("Same");
    Pip_Pos_fit->Draw("Same");
    Pip_Neg_fit->Draw("Same");

    fit_functions->NewFunction();
    fit_functions->Set_RetrunType("double");
    fit_functions->Set_FuncName("Pip_Pos_fit");
    fit_functions->Set_FuncInputs("double x");
    fit_functions->Set_Function(Pip_Pos_fit->GetExpFormula("P"));
    fit_functions->WriteFunction();

    fit_functions->NewFunction();
    fit_functions->Set_RetrunType("double");
    fit_functions->Set_FuncName("Pip_Neg_fit");
    fit_functions->Set_FuncInputs("double x");
    fit_functions->Set_Function(Pip_Neg_fit->GetExpFormula("P"));
    fit_functions->WriteFunction();

    fit_functions->NewFunction();
    fit_functions->Set_RetrunType("bool");
    fit_functions->Set_FuncName("Between_Pip_fit");
    fit_functions->Set_FuncInputs("double dt, double p");
    fit_functions->AddLine("bool between = true");
    fit_functions->AddLine("between &= (dt >= Pip_Neg_fit(p))");
    fit_functions->AddLine("between &= (dt <= Pip_Pos_fit(p))");
    fit_functions->Set_Function("between");
    fit_functions->WriteFunction();

    fit_functions->NewFunction();
    fit_functions->Set_RetrunType("bool");
    fit_functions->Set_FuncName("Between_Proton_fit");
    fit_functions->Set_FuncInputs("double dt, double p");
    fit_functions->AddLine("bool between = true");
    fit_functions->AddLine("between &= (dt >= Proton_Neg_fit(p))");
    fit_functions->AddLine("between &= (dt <= Proton_Pos_fit(p))");
    fit_functions->Set_Function("between");
    fit_functions->WriteFunction();
    delete fit_functions;
  }
};
