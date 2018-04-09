#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)

double poissonf(double *x, double *par) { return par[0] * TMath::Poisson(x[0], par[1]); }

const char *fin[8] = {"0.5ms_1000s", "0.5ms_100s", "0.5ms_10s", "50ms_10000s",
                      "5ms_1000s",   "5ms_100s",   "5ms_10s",   "500ms_50000s.n0"};

TCanvas c1("c1", "canvas", 1280, 720);

const int bins[8] = {8, 8, 8, 65, 17, 16, 13, 125};
const int min_val[8] = {0, 0, 0, 25, 0, 0, 0, 425};
const int max_val[8] = {7, 5, 5, 90, 17, 16, 13, 675};
const float dt[8] = {0.5, 0.5, 0.5, 50, 5, 5, 5, 500};
const float T[6] = {1000, 100, 10, 1000, 100, 10};

float process(int file) {
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TChain *chain = new TChain("Sr90");
  char file_name[50];
  sprintf(file_name, "%s.root", fin[file]);
  chain->Add(file_name);

  int _x1;
  chain->SetBranchAddress("x1", &_x1);
  int num_of_events = (int)chain->GetEntries();
  TH1I *hist = new TH1I(fin[file], fin[file], bins[file], min_val[file], max_val[file]);
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    hist->Fill(_x1);
  }

  chain->Reset();
  c1.cd((int)file + 1);
  hist->Draw();
  if (file == 3 || file == 7) {
    hist->Fit("gaus", "Q+", "", min_val[file], max_val[file]);
    return -99;
  } else {
    TF1 *pois = new TF1("pois", poissonf, 0, 10, 2);
    pois->SetParName(0, "Const");
    pois->SetParName(1, "#mu");
    pois->SetParameter(0, 1);
    pois->SetParameter(1, 1);
    hist->Fit("pois", "Q+");
    return pois->GetParameter("#mu");
  }
}

void run() {
  std::vector<float> r_vec;
  float temp;
  gStyle->SetOptStat(111);
  gStyle->SetOptFit(111);
  c1.Divide(4, 2);
  for (int i = 0; i < 8; i++) {
    temp = process(i);
    if (temp != -99) r_vec.push_back(temp / dt[i]);
  }
  for (auto const &value : r_vec) std::cout << value << std::endl;
}
