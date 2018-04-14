// A basic analysis program looking at
// a electron beam on a proton target
//
// Follow the TODO portions to get the analysis working properly
//

#define Square(x) ((x) * (x))

static const double BEAM = 4.81726;       // Beam energy in GeV
static const double MASS_P = 0.93827203;  // Mass in GeV
static const double MASS_E = 0.000511;    // Mass in GeV

TLorentzVector e_mu(0.0, 0.0, TMath::Sqrt((BEAM * BEAM) - (MASS_E * MASS_E)), BEAM);
TLorentzVector p_mu(0.0, 0.0, 0.0, MASS_P);

double Q2_calc(TLorentzVector q_mu) { return -q_mu.Mag2(); }

double W_calc(TLorentzVector q_mu) { return (p_mu + q_mu).Mag(); }

void analyze() {
  const char* fin = "511_lab_E_data_small.root";
  const char* fout = "WvsQ2.root";
  TH1D* W_hist = new TH1D("W_hist", "W_hist", 500, 0.0, 3.0);
  TH1D* Q2_hist = new TH1D("Q2_hist", "Q2_hist", 500, 0.0, 4.0);
  TH2D* W_vs_Q2_hist = new TH2D("W_vs_Q2_hist", "W_vs_Q2_hist", 500, 0.0, 3.0, 500, 0.0, 4.0);

  // Load chain from branch lab
  TFile* OutputFile = new TFile(fout, "RECREATE");
  TChain chain("lab");
  chain.Add(fin);
  double e_p, e_cx, e_cy, e_cz;
  chain.SetBranchAddress("p", &e_p);
  chain.SetBranchAddress("cx", &e_cx);
  chain.SetBranchAddress("cy", &e_cy);
  chain.SetBranchAddress("cz", &e_cz);

  // Create 4 vectors for the scattered electron
  TVector3 e_mu_prime_3;
  TLorentzVector e_mu_prime;
  TLorentzVector q_mu;
  double W, Q2;

  int num_of_events = (int)chain.GetEntries();
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain.GetEntry(current_event);
    e_mu_prime_3.SetXYZ(e_p * e_cx, e_p * e_cy, e_p * e_cz);
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
    q_mu = (e_mu - e_mu_prime);
    Q2 = Q2_calc(q_mu);
    W = W_calc(q_mu);
    W_hist->Fill(W);
    Q2_hist->Fill(Q2);
    W_vs_Q2_hist->Fill(W, Q2);
  }
  //
  // end stuff
  chain.Reset();
  OutputFile->cd();
  W_hist->Write();
  Q2_hist->Write();
  W_vs_Q2_hist->Write();
  OutputFile->Close();
}
