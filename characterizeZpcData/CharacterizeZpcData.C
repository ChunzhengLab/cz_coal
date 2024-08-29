#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TBits.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include <iostream>
#include <fstream>

struct zpcInfo
{
  int nevent;
  int nparton;
  int ID[100000];
  float Px[100000];
  float Py[100000];
  float Pz[100000];
  float Mass[100000];
  float X[100000];
  float Y[100000];
  float Z[100000];
  float Time[100000];
};

int color[9] = {kRed, kBlue, kGreen, kYellow, kMagenta, kCyan, kOrange, kBlack, kViolet};

void CharacterizeZpcData()
{
  std::cout << "CharacterizeZpcData" << std::endl;
  zpcInfo zpc;
  TFile *file = new TFile("zpc-1.root");
  if (!file) {
    std::cerr << "Could not open file" << std::endl;
    return;
  }
  TTree *tree = (TTree*)file->Get("AMPT");
  if (!tree) {
    std::cerr << "Could not open tree" << std::endl;
    return;
  }

  tree->SetBranchAddress("Event", &zpc);
  tree->SetBranchAddress("ID", zpc.ID);
  tree->SetBranchAddress("Px", zpc.Px);
  tree->SetBranchAddress("Py", zpc.Py);
  tree->SetBranchAddress("Pz", zpc.Pz);
  tree->SetBranchAddress("Mass", zpc.Mass);
  tree->SetBranchAddress("X", zpc.X);
  tree->SetBranchAddress("Y", zpc.Y);
  tree->SetBranchAddress("Z", zpc.Z);
  tree->SetBranchAddress("Time", zpc.Time);
  std::cout << "Number of entries: " << tree->GetEntries() << std::endl;
  
  float maxX = 0, maxY = 0, maxZ = 0;
  float minX = 0, minY = 0, minZ = 0;
  float maxPx = 0, maxPy = 0, maxPz = 0;
  float minPx = 0, minPy = 0, minPz = 0;
  float maxMass = 0, minMass = 0;
  float maxTime = 0, minTime = 0;

  TH1D* h_z_PzHiger2000 = new TH1D("h_z_PzHiger2000", "z distribution of partons with Pz > 2000", 200, -10000, 10000);
  TH1D* h_Pz[10];
  for (int iBit = 0; iBit < 9; iBit++) {
    // 0.0001 = pow(10, -4), 0.001 = pow(10, -3), 0.01 = pow(10, -2), 0.1 = pow(10, -1), 1 = pow(10, 0), 10 = pow(10, 1), 100 = pow(10, 2), 1000 = pow(10, 3), 10000 = pow(10, 4)
    h_Pz[iBit] = new TH1D(Form("h_Pz_%d", iBit), Form("Pz distribution of partons with Z > %f", pow(10, iBit - 4)), 200, 0, 2000);
  }

  int nEvents = tree->GetEntries();
  for (int iEvent = 0; iEvent < nEvents; iEvent++) {
    tree->GetEntry(iEvent);
    int nPartons = zpc.nparton;
    for (int iParton = 0; iParton < nPartons; iParton++) {
      if (abs(zpc.X[iParton]) > maxX) maxX = abs(zpc.X[iParton]);
      if (abs(zpc.Y[iParton]) > maxY) maxY = abs(zpc.Y[iParton]);
      if (abs(zpc.Z[iParton]) > maxZ) maxZ = abs(zpc.Z[iParton]);
      if (abs(zpc.Px[iParton]) > maxPx) maxPx = abs(zpc.Px[iParton]);
      if (abs(zpc.Py[iParton]) > maxPy) maxPy = abs(zpc.Py[iParton]);
      if (abs(zpc.Pz[iParton]) > maxPz) maxPz = abs(zpc.Pz[iParton]);
      if (abs(zpc.Mass[iParton]) > maxMass) maxMass = abs(zpc.Mass[iParton]);
      if (abs(zpc.Time[iParton]) > maxTime) maxTime = abs(zpc.Time[iParton]);

      // if (abs(zpc.Pz[iParton]) > 200) {
      //   std::cout << "====================="<< std::endl;
      //   std::cout << "We got a crazy parton!!!" << std::endl;
      //   std::cout << "Event: " << zpc.nevent << std::endl;
      //   std::cout << "Parton: " << iParton << std::endl;
      //   std::cout << "ID: " << zpc.ID[iParton] << std::endl;
      //   std::cout << "Pz: " << zpc.Pz[iParton] << std::endl;
      //   std::cout << "X: " << zpc.X[iParton] << std::endl;
      //   std::cout << "Y: " << zpc.Y[iParton] << std::endl;
      //   std::cout << "Z: " << zpc.Z[iParton] << std::endl;
      //   std::cout << "Time: " << zpc.Time[iParton] << std::endl;
      //   std::cout << "====================="<< std::endl;
      //   h_z_PzHiger2000->Fill(zpc.Z[iParton]);
      // }

      TBits nBits(9);
      nBits.SetBitNumber(1, abs(zpc.Z[iParton]) > 0.001 && abs(zpc.Z[iParton]) < 0.01);
      nBits.SetBitNumber(2, abs(zpc.Z[iParton]) > 0.01 && abs(zpc.Z[iParton]) < 0.1);
      nBits.SetBitNumber(3, abs(zpc.Z[iParton]) > 0.1 && abs(zpc.Z[iParton]) < 1);
      nBits.SetBitNumber(4, abs(zpc.Z[iParton]) > 1 && abs(zpc.Z[iParton]) < 10);
      nBits.SetBitNumber(5, abs(zpc.Z[iParton]) > 10 && abs(zpc.Z[iParton]) < 100);
      nBits.SetBitNumber(6, abs(zpc.Z[iParton]) > 100 && abs(zpc.Z[iParton]) < 1000);
      nBits.SetBitNumber(7, abs(zpc.Z[iParton]) > 1000 && abs(zpc.Z[iParton]) < 10000);
      nBits.SetBitNumber(8, abs(zpc.Z[iParton]) > 10000 && abs(zpc.Z[iParton]) < 100000);

      for (int iBit = 1; iBit < 9; iBit++) {
        if (nBits.TestBitNumber(iBit)) {
          h_Pz[iBit]->Fill(abs(zpc.Pz[iParton]));
        }
      }
    }
  }

  for (int iBit = 1; iBit < 9; iBit++) {
    h_Pz[iBit]->Scale(1.0 / h_Pz[iBit]->Integral());
    h_Pz[iBit]->SetLineColor(color[iBit]);
  }

  TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  for (int iBit = 1; iBit < 8; iBit++) {
    legend->AddEntry(h_Pz[iBit], Form("Z > %f and Z < %f", pow(10, iBit - 4), pow(10, iBit - 3)), "l");
  }

  //写到一个txt文件中
  std::ofstream outfile;
  outfile.open("CharacterizeZpcData.txt");
  outfile << "X: [" << minX << ", " << maxX << "]" << std::endl;
  outfile << "Y: [" << minY << ", " << maxY << "]" << std::endl;
  outfile << "Z: [" << minZ << ", " << maxZ << "]" << std::endl;
  outfile << "Px: [" << minPx << ", " << maxPx << "]" << std::endl;
  outfile << "Py: [" << minPy << ", " << maxPy << "]" << std::endl;
  outfile << "Pz: [" << minPz << ", " << maxPz << "]" << std::endl;
  outfile << "Mass: [" << minMass << ", " << maxMass << "]" << std::endl;
  outfile << "Time: [" << minTime << ", " << maxTime << "]" << std::endl;
  outfile.close();

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->SetLogy();
  c1->SetLogx();
  for (int iBit = 1; iBit < 8; iBit++) {
    h_Pz[iBit]->GetYaxis()->SetTitle("Normalized counts");
    h_Pz[iBit]->GetXaxis()->SetTitle("Pz (GeV/c)");
    h_Pz[iBit]->Draw("same L");
  }
  legend->Draw();

}