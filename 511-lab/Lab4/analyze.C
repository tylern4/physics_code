// A basic analysis program looking at
// a electron beam on a proton target
//
// Follow the TODO portions to get the analysis working properly
//

#define Square(x) ((x) * (x))

static const double BEAM = 4.81726;         // Beam energy in GeV
static const double MASS_P = 0.93827203;    // Mass in GeV
static const double MASS_E = 0.000511;      // Mass in GeV
static const double MASS_PIP = 0.13957018;  // Mass in GeV

TCanvas c1("c1", "Plots", 1280, 720);

double Breit(double *x, double *par) { return par[2] * TMath::BreitWigner(x[0], par[0], par[1]); }

double missing_mass(TLorentzVector gamma_mu, TLorentzVector pip_mu) {
  TVector3 target_3;
  TLorentzVector target;
  // Set target vector
  target_3.SetXYZ(0.0, 0.0, 0.0);
  target.SetVectM(target_3, MASS_P);
  TLorentzVector reaction;
  reaction = (gamma_mu + target - pip_mu);

  return reaction.M();
}

void analyze() {
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  double e_p, e_cx, e_cy, e_cz;
  double pip_p, pip_cx, pip_cy, pip_cz;
  TH1D *MM = new TH1D("mm", "mm", 500, 0, 3);

  TF1 *bw = new TF1("bw", Breit, 0, 2, 3);

  // Load chain from branch h10
  TChain lab("lab");
  lab.Add("511_lab_E_PIP_data.root");

  lab.SetBranchAddress("e_p", &e_p);
  lab.SetBranchAddress("e_cx", &e_cx);
  lab.SetBranchAddress("e_cy", &e_cy);
  lab.SetBranchAddress("e_cz", &e_cz);
  lab.SetBranchAddress("pip_p", &pip_p);
  lab.SetBranchAddress("pip_cx", &pip_cx);
  lab.SetBranchAddress("pip_cy", &pip_cy);
  lab.SetBranchAddress("pip_cz", &pip_cz);

  int num_of_events = (int)lab.GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    lab.GetEntry(current_event);
    // Setup scattered electron 4 vector
    TVector3 e_mu_prime_3;
    TLorentzVector e_mu_prime;
    TVector3 pip_mu_prime_3;
    TLorentzVector pip_mu_prime;
    TLorentzVector e_mu(0.0, 0.0, sqrt(Square(BEAM) - Square(MASS_E)), BEAM);

    e_mu_prime_3.SetXYZ(e_p * e_cx, e_p * e_cy, e_p * e_cz);
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

    pip_mu_prime_3.SetXYZ(pip_p * pip_cx, pip_p * pip_cy, pip_p * pip_cz);
    pip_mu_prime.SetVectM(pip_mu_prime_3, MASS_PIP);

    TLorentzVector gamma_mu;
    gamma_mu = e_mu - e_mu_prime;
    double mm = missing_mass(gamma_mu, pip_mu_prime);
    MM->Fill(mm);
  }
  //
  // end stuff
  lab.Reset();
  c1.cd();
  bw->SetParName(0, "Mean");
  bw->SetParName(1, "Width");
  bw->SetParName(2, "Const");
  bw->SetParameter(0, 1);
  bw->SetParameter(1, 1);
  bw->SetParameter(2, 1000);
  bw->SetParLimits(2, 100.0, 10000.0);
  MM->Fit("bw", "M", "", 0.5, 1.1);
  MM->Draw();
}
