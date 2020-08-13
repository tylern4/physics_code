#include "TTree.h"

void convert_data_to_root() {
  // Name the tree to fill
  TTree *tree = new TTree("data", "data");
  // Tell how many columns there are and what you want them to be named and the
  // input file
  tree->ReadFile("../Lab2/data.dat", "x:y:yerr", ',');
  // Save to a root file
  // tree->SaveAs("data.root");

  // Setup variables
  float x = 0.0;
  float y = 0.0;
  float yerr = 0.0;
  std::vector<float> xVec;
  std::vector<float> yVec;
  std::vector<float> yerrVec;

  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("yerr", &yerr);

  int n = (int)tree->GetEntries();
  // Loop over the events in the file and print out the variables
  for (int i = 0; i < n; i++) {
    tree->GetEntry(i);
    xVec.push_back(x);
    yVec.push_back(y);
    yerrVec.push_back(yerr);
    std::cout << x << ",";
    std::cout << y << ",";
    std::cout << yerr << std::endl;
  }

  TGraphErrors *gr =
      new TGraphErrors(xVec.size(), &xVec[0], &yVec[0], 0, &yerrVec[0]);

  gr->Fit("pol2");
  gr->Draw("A*");
}
