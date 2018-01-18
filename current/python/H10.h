/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include <iostream>
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
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
      electron_cuts &= (p[0] > MIN_P_CUT);       // Minimum Momentum cut
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
        hists.CC_fill(cc_sector, cc_segment, cc_pmt, cc_nphe, 0);
        hists.Fill_Beam_Position((double)dc_vx[dc[0] - 1],
                                 (double)dc_vy[dc[0] - 1],
                                 (double)dc_vz[dc[0] - 1]);
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

        hists.Fill_electron_fid(theta, phi, sector);

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

  void missing_mass_fit(TH1D Missing_Mass) {
    // Start of cuts
    TCanvas *c2 = new TCanvas;
    Fits MM_neutron_cut;
    double fit_range_min = 0.88;
    double fit_range_max = 1.0;
    MM_neutron_cut.FitGaus(&Missing_Mass, fit_range_min, fit_range_max);

    Header *MM_header = new Header("../src/missing_mass_gaussians.hpp", "MM");
    MM_header->WriteGaussian("mm", 1, MM_neutron_cut.mean,
                             MM_neutron_cut.sigma);

    delete MM_header;
    std::cout << "Here" << std::endl;
    // return Missing_Mass;
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

  void Skim(char *fin, char *RootFile_output) {
    TFile *RootOutputFile;
    int number_cols = 0;
    char rootFile[500];
    int num_of_events, total_events, num_of_pis;
    bool electron_cuts, MM_cut, has_neutron;

    MissingMass *MM_neutron = new MissingMass();
    MM_neutron->Set_target_mass(MASS_P);
    MM_neutron->Set_target_PxPyPz(0);
    Delta_T *delta_t = new Delta_T();

    Float_t W, Q2, MM;
    std::vector<bool> is_proton, is_pip, is_electron, is_pim;
    std::vector<double> dt_proton, dt_pip;

    TVector3 e_mu_prime_3;
    TLorentzVector e_mu_prime;
    TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)),
                        E1D_E0);

    TVector3 Particle3(0.0, 0.0, 0.0);
    TLorentzVector Particle4(0.0, 0.0, 0.0, 0.0);

    RootOutputFile = new TFile(RootFile_output, "RECREATE");

    TChain chain("h10");
    cout << "Analyzing file " << fin << endl;
    chain.AddFile(fin);

    getBranches(&chain);

    num_of_events = (int)chain.GetEntries();

    TTree *skim = chain.CloneTree(0);
    TBranch *W_branch = skim->Branch("W", &W);
    TBranch *Q2_branch = skim->Branch("Q2", &Q2);
    TBranch *MM_branch = skim->Branch("MM", &MM);

    TBranch *is_Electron = skim->Branch("is_electron", &is_electron);
    TBranch *is_Proton = skim->Branch("is_proton", &is_proton);
    TBranch *is_Pip = skim->Branch("is_pip", &is_pip);
    TBranch *is_Pim = skim->Branch("is_pim", &is_pim);

    TBranch *DeltaT_P_branch =
        skim->Branch("DeltaT_P", "vector<double>", &dt_proton);
    TBranch *DeltaT_Pip_branch =
        skim->Branch("DeltaT_Pip", "vector<double>", &dt_pip);
    TBranch *NumPI_branch = skim->Branch("NumPI", &num_of_pis);
    TBranch *Neutron_branch = skim->Branch("has_neutron", &has_neutron);

    for (int current_event = 0; current_event < num_of_events;
         current_event++) {
      chain.GetEntry(current_event);
      is_proton = std::vector<bool>(gpart, false);
      is_electron = std::vector<bool>(gpart, false);
      is_pip = std::vector<bool>(gpart, false);
      is_pim = std::vector<bool>(gpart, false);

      electron_cuts = true;
      // electron cuts
      electron_cuts &= (id[0] == ELECTRON); // First particle is electron
      electron_cuts &=
          (gpart > 0); // Number of good particles is greater than 0
      electron_cuts &= (stat[0] > 0);            // First Particle hit stat
      electron_cuts &= ((int)q[0] == -1);        // First particle is negative Q
      electron_cuts &= (sc[0] > 0);              // First Particle hit sc
      electron_cuts &= (dc[0] > 0);              // ``` ``` ``` d
      electron_cuts &= (ec[0] > 0);              // ``` ``` ``` ec
      electron_cuts &= (dc_stat[dc[0] - 1] > 0); //??

      e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

      dt_proton = delta_t->delta_t_array(MASS_P, gpart);
      dt_pip = delta_t->delta_t_array(MASS_PIP, gpart);

      for (int part_num = 1; part_num < gpart; part_num++) {
        num_of_pis = 0;
        if (dt_pip.at(part_num) >= -20 && dt_pip.at(part_num) <= 20) {
          is_pip.at(part_num) = true;
          num_of_pis++;
          TLorentzVector gamma_mu = (e_mu - e_mu_prime);
          MM_neutron->Set_PxPyPz(p[part_num] * cx[part_num],
                                 p[part_num] * cy[part_num],
                                 p[part_num] * cz[part_num]);
          MM = MM_neutron->missing_mass(gamma_mu);
        }
        if (dt_proton.at(part_num) >= -20 && dt_proton.at(part_num) <= 20 &&
            q[part_num] == 1) {
          is_proton.at(part_num) = true;
        }
        if (dt_pip.at(part_num) >= -20 && dt_pip.at(part_num) <= 20) {
          is_pim.at(part_num) = true;
        }
      }

      has_neutron = true; // between_mm(MM);

      if (electron_cuts && has_neutron) {
        W = W_calc(e_mu, e_mu_prime);
        Q2 = Q2_calc(e_mu, e_mu_prime);
        is_electron.at(0) = true;
        skim->Fill(); // Fill the banks after the skim
      }
    }
    //
    // end stuff
    chain.Reset(); // delete Tree object

    RootOutputFile->cd();
    RootOutputFile->Write();
    RootOutputFile->Close();
  }
};
