#include <iostream>
#include "TGraph2D.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;

struct AmptEvent {
  int nevent;
  int nparton;
  int   ID[30000];   //[nparton]
  float Px[30000];   //[nparton]
  float Py[30000];   //[nparton]
  float Pz[30000];   //[nparton]
  float X[30000];   //[nparton]
  float Y[30000];   //[nparton]
  float Z[30000];   //[nparton]
};


int BinPDG(int pdg) {
  // 第一个bin是u
  if (pdg == 2) return 0;
  // 第二个bin是ubar
  if (pdg == -2) return 1;
  // 第三个bin是d
  if (pdg == 1) return 2;
  // 第四个bin是dbar
  if (pdg == -1) return 3;
  // 第五个bin是s
  if (pdg == 3) return 4;
  // 第六个bin是sbar
  if (pdg == -3) return 5;
  // 第七个bin是c
  if (pdg == 4) return 6;
  // 第八个bin是cbar
  if (pdg == -4) return 7;
  // 第九个bin是b
  if (pdg == 5) return 8;
  // 第十个bin是bbar
  if (pdg == -5) return 9;
  // 第十一个bin是其他的正粒子
  if (pdg > 0) return 10;
  // 第十二个bin是其他的反粒子
  if (pdg < 0) return 11;
  return -1;
}

void GetAMPTPartonDist() {
  TFile* file = new TFile("../zpc-1.root", "READ");
  if (!file) {
    std::cerr << "Cannot open file" << std::endl;
    return;
  }
  TTree* tree = (TTree*)file->Get("AMPT");

  AmptEvent amptEvent;
  tree->SetBranchAddress("Event", &amptEvent);
  tree->SetBranchAddress("ID", amptEvent.ID);
  tree->SetBranchAddress("Px", amptEvent.Px);
  tree->SetBranchAddress("Py", amptEvent.Py);
  tree->SetBranchAddress("Pz", amptEvent.Pz);
  tree->SetBranchAddress("X", amptEvent.X);
  tree->SetBranchAddress("Y", amptEvent.Y);
  tree->SetBranchAddress("Z", amptEvent.Z);

  TH1I* h_nparton = new TH1I("h_nparton", "Number of partons per event; Number of partons; Number of events", 50, 5000, 20000);
  TH3D* h3_xyz = new TH3D("h3_xyz", "x-y-z distribution; x(fm); y(fm); z(fm)", 300, -30, 30, 300, -30, 30, 300, -30, 30);
  TH3D* h3_pxpypz = new TH3D("h3_pxpypz", "px-py-pz distribution; px(GeV/c); py(GeV/c); pz(GeV/c)", 300, -30, 30, 300, -30, 30, 300, -30, 30);
  TH1I* h_pdg = new TH1I("h_pdg", "PDG distribution", 12, 0, 12);
  h_pdg->GetXaxis()->SetBinLabel(1, "u");
  h_pdg->GetXaxis()->SetBinLabel(2, "ubar");
  h_pdg->GetXaxis()->SetBinLabel(3, "d");
  h_pdg->GetXaxis()->SetBinLabel(4, "dbar");
  h_pdg->GetXaxis()->SetBinLabel(5, "s");
  h_pdg->GetXaxis()->SetBinLabel(6, "sbar");
  h_pdg->GetXaxis()->SetBinLabel(7, "c");
  h_pdg->GetXaxis()->SetBinLabel(8, "cbar");
  h_pdg->GetXaxis()->SetBinLabel(9, "b");
  h_pdg->GetXaxis()->SetBinLabel(10, "bbar");
  h_pdg->GetXaxis()->SetBinLabel(11, "other");
  h_pdg->GetXaxis()->SetBinLabel(12, "otherbar");

  int nEvents = tree->GetEntries();
  for (int iEvent = 0; iEvent < nEvents; iEvent++) {
    tree->GetEntry(iEvent);
    int serial = amptEvent.nevent;
    cout<<"Event: "<<serial<<endl;
    int nPartons = amptEvent.nparton;
    h_nparton->Fill(nPartons);
    for (int i = 0; i < amptEvent.nparton; i++) {
      h3_xyz->Fill(amptEvent.X[i], amptEvent.Y[i], amptEvent.Z[i]);
      h3_pxpypz->Fill(amptEvent.Px[i], amptEvent.Py[i], amptEvent.Pz[i]);
      h_pdg->Fill(BinPDG(amptEvent.ID[i]));
    }
  }

  TH2D* h2_xy = (TH2D*)h3_xyz->Project3D("xy");
  h2_xy->SetName("h2_xy");
  h2_xy->SetTitle("x-y distribution; x(fm); y(fm)");

  TH1D* h_z = (TH1D*)h3_xyz->Project3D("z");
  h_z->SetName("h_z");
  h_z->SetTitle("z distribution; z(fm)");

  TH2D* h2_pxpy = (TH2D*)h3_pxpypz->Project3D("xy");
  h2_pxpy->SetName("h2_pxpy");
  h2_pxpy->SetTitle("px-py distribution; px(GeV/c); py(GeV/c)");

  TH1D* h_pz = (TH1D*)h3_pxpypz->Project3D("z");
  h_pz->SetName("h_pz");
  h_pz->SetTitle("pz distribution; pz(GeV/c)");

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  c1->Divide(2,2);
  c1->cd(1);
  h2_xy->Draw("colz");
  c1->cd(2);
  h_z->Draw();
  c1->cd(3);
  h2_pxpy->Draw("colz");
  c1->cd(4);
  h_pz->Draw();

  TCanvas* c2 = new TCanvas("c2", "c2", 800, 800);
  h_pdg->Draw();

  TCanvas* c3 = new TCanvas("c3", "c3", 800, 800);
  h_nparton->Draw();

  TFile* outFile = new TFile("AMPTPartonDist.root", "RECREATE");
  outFile->cd();
  h_nparton->Write();
  // h3_xyz->Write();
  h2_xy->Write();
  h_z->Write();
  // h3_pxpypz->Write();
  h2_pxpy->Write();
  h_pz->Write();
  h_pdg->Write();
  outFile->Close();
}