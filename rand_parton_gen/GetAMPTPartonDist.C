#include <iostream>
#include "TGraph2D.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TCanvas.h"

using namespace std;

float phi(float x , float y) {
  float phi = atan2(y , x);
  if (phi < 0) phi += 2 * M_PI;
  return phi;
}

float dphi(float phi1 , float phi2) {
  float dphi = phi1 - phi2;
  while (dphi < -0.5 * M_PI) dphi += 2 * M_PI;
  while (dphi > 1.5 * M_PI) dphi -= 2 * M_PI;
  return dphi;
}

float eta(float x , float y , float z) {
  float r = sqrt(x * x + y * y);
  return 0.5 * log((r + z) / (r - z));
}

float deta(float eta1 , float eta2) {
  return eta1 - eta2;
}

int pairType(int pdg1 , int pdg2) {
  //uu
  if (pdg1 == 2 && pdg2 == 2) return 0;
  //ud
  if (pdg1 == 2 && pdg2 == 1) return 1;
  //us
  if (pdg1 == 2 && pdg2 == 3) return 2;
  //uubar
  if (pdg1 == 2 && pdg2 == -2) return 3;
  //udbar
  if (pdg1 == 2 && pdg2 == -1) return 4;
  //usbar
  if (pdg1 == 2 && pdg2 == -3) return 5;
  //dd
  if (pdg1 == 1 && pdg2 == 1) return 6;
  //ds
  if (pdg1 == 1 && pdg2 == 3) return 7;
  //ddbar
  if (pdg1 == 1 && pdg2 == -1) return 8;
  //dsbar
  if (pdg1 == 1 && pdg2 == -3) return 9;
  //ss
  if (pdg1 == 3 && pdg2 == 3) return 10;
  //sbarbar
  if (pdg1 == -3 && pdg2 == -3) return 11;
  return -1;
}


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
  TFile* file = new TFile("../test/zpc-1.root" , "READ");
  if (!file) {
    std::cerr << "Cannot open file" << std::endl;
    return;
  }
  TTree* tree = (TTree*)file->Get("AMPT");

  AmptEvent amptEvent;
  tree->SetBranchAddress("Event" , &amptEvent);
  tree->SetBranchAddress("ID" , amptEvent.ID);
  tree->SetBranchAddress("Px" , amptEvent.Px);
  tree->SetBranchAddress("Py" , amptEvent.Py);
  tree->SetBranchAddress("Pz" , amptEvent.Pz);
  tree->SetBranchAddress("X" , amptEvent.X);
  tree->SetBranchAddress("Y" , amptEvent.Y);
  tree->SetBranchAddress("Z" , amptEvent.Z);

  TH1I* h_nparton = new TH1I("h_nparton" , "Number of partons per event; Number of partons; Number of events" , 50 , 5000 , 20000);
  TH3D* h3_xyz = new TH3D("h3_xyz" , "x-y-z distribution; x(fm); y(fm); z(fm)" , 300 , -30 , 30 , 300 , -30 , 30 , 300 , -30 , 30);
  TH3D* h3_pxpypz = new TH3D("h3_pxpypz" , "px-py-pz distribution; px(GeV/c); py(GeV/c); pz(GeV/c)" , 300 , -30 , 30 , 300 , -30 , 30 , 300 , -30 , 30);
  TH1I* h_pdg = new TH1I("h_pdg" , "PDG distribution" , 12 , 0 , 12);
  h_pdg->GetXaxis()->SetBinLabel(1 , "u");
  h_pdg->GetXaxis()->SetBinLabel(2 , "ubar");
  h_pdg->GetXaxis()->SetBinLabel(3 , "d");
  h_pdg->GetXaxis()->SetBinLabel(4 , "dbar");
  h_pdg->GetXaxis()->SetBinLabel(5 , "s");
  h_pdg->GetXaxis()->SetBinLabel(6 , "sbar");
  h_pdg->GetXaxis()->SetBinLabel(7 , "c");
  h_pdg->GetXaxis()->SetBinLabel(8 , "cbar");
  h_pdg->GetXaxis()->SetBinLabel(9 , "b");
  h_pdg->GetXaxis()->SetBinLabel(10 , "bbar");
  h_pdg->GetXaxis()->SetBinLabel(11 , "other");
  h_pdg->GetXaxis()->SetBinLabel(12 , "otherbar");

  TH1D* h_dphip_uu = new TH1D("h2_dphip_uu" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphip_ud = new TH1D("h2_dphip_ud" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphip_us = new TH1D("h2_dphip_us" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphip_uubar = new TH1D("h2_dphip_uubar" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphip_udbar = new TH1D("h2_dphip_udbar" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphip_usbar = new TH1D("h2_dphip_usbar" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphip_dd = new TH1D("h2_dphip_dd" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphip_ds = new TH1D("h2_dphip_ds" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphip_ddbar = new TH1D("h2_dphip_ddbar" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphip_dsbar = new TH1D("h2_dphip_dsbar" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphip_ss = new TH1D("h2_dphip_ss" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphip_sbarbar = new TH1D("h2_dphip_sbarbar" , "dphip distribution; dphip" , 30 , -0.5 * M_PI , 1.5 * M_PI);

  TH1D* h_dphis_uu = new TH1D("h2_dphis_uu" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphis_ud = new TH1D("h2_dphis_ud" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphis_us = new TH1D("h2_dphis_us" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphis_uubar = new TH1D("h2_dphis_uubar" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphis_udbar = new TH1D("h2_dphis_udbar" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphis_usbar = new TH1D("h2_dphis_usbar" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphis_dd = new TH1D("h2_dphis_dd" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphis_ds = new TH1D("h2_dphis_ds" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphis_ddbar = new TH1D("h2_dphis_ddbar" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphis_dsbar = new TH1D("h2_dphis_dsbar" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphis_ss = new TH1D("h2_dphis_ss" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);
  TH1D* h_dphis_sbarbar = new TH1D("h2_dphis_sbarbar" , "dphis distribution; dphis" , 30 , -0.5 * M_PI , 1.5 * M_PI);

  TProfile* p_delta_p = new TProfile("p_delta_p" , ";;delta" , 12 , 0 , 12);
  p_delta_p->GetXaxis()->SetBinLabel(1 , "uu");
  p_delta_p->GetXaxis()->SetBinLabel(2 , "ud");
  p_delta_p->GetXaxis()->SetBinLabel(3 , "us");
  p_delta_p->GetXaxis()->SetBinLabel(4 , "dd");
  p_delta_p->GetXaxis()->SetBinLabel(5 , "ds");
  p_delta_p->GetXaxis()->SetBinLabel(6 , "ss");
  p_delta_p->GetXaxis()->SetBinLabel(7 , "uubar");
  p_delta_p->GetXaxis()->SetBinLabel(8 , "udbar");
  p_delta_p->GetXaxis()->SetBinLabel(9 , "usbar");
  p_delta_p->GetXaxis()->SetBinLabel(10 , "ddbar");
  p_delta_p->GetXaxis()->SetBinLabel(11 , "dsbar");
  p_delta_p->GetXaxis()->SetBinLabel(12 , "sbarbar");

  TProfile* p_delta_s = new TProfile("p_delta_s" , ";;delta" , 12 , 0 , 12);
  p_delta_s->GetXaxis()->SetBinLabel(1 , "uu");
  p_delta_s->GetXaxis()->SetBinLabel(2 , "ud");
  p_delta_s->GetXaxis()->SetBinLabel(3 , "us");
  p_delta_s->GetXaxis()->SetBinLabel(4 , "dd");
  p_delta_s->GetXaxis()->SetBinLabel(5 , "ds");
  p_delta_s->GetXaxis()->SetBinLabel(6 , "ss");
  p_delta_s->GetXaxis()->SetBinLabel(7 , "uubar");
  p_delta_s->GetXaxis()->SetBinLabel(8 , "udbar");
  p_delta_s->GetXaxis()->SetBinLabel(9 , "usbar");
  p_delta_s->GetXaxis()->SetBinLabel(10 , "ddbar");
  p_delta_s->GetXaxis()->SetBinLabel(11 , "dsbar");
  p_delta_s->GetXaxis()->SetBinLabel(12 , "sbarbar");


  int nEvents = tree->GetEntries();
  for (int iEvent = 0; iEvent < nEvents; iEvent++) {
    tree->GetEntry(iEvent);
    int serial = amptEvent.nevent;
    // if (serial % 100 == 0) {
    //   std::cout << "Processing event " << serial << std::endl;
    // }
    std::cout << "Processing event " << serial << std::endl;
    int nPartons = amptEvent.nparton;
    h_nparton->Fill(nPartons);

    for (int i = 0; i < amptEvent.nparton; i++) {
      h3_xyz->Fill(amptEvent.X[i] , amptEvent.Y[i] , amptEvent.Z[i]);
      h3_pxpypz->Fill(amptEvent.Px[i] , amptEvent.Py[i] , amptEvent.Pz[i]);
      h_pdg->Fill(BinPDG(amptEvent.ID[i]));

      for (int j = i + 1; j < amptEvent.nparton; j++) {
        float eta1 = eta(amptEvent.Px[i] , amptEvent.Py[i] , amptEvent.Pz[i]);
        if (eta1 > 0.8 || eta1 < -0.8) continue;
        float eta2 = eta(amptEvent.Px[j] , amptEvent.Py[j] , amptEvent.Pz[j]);
        if (eta2 > 0.8 || eta2 < -0.8) continue;
        float deta_ = deta(eta1 , eta2);
        if (deta_ > 1.6 || deta_ < -1.6) continue;

        float dphis = dphi(phi(amptEvent.X[i] , amptEvent.Y[i]) , phi(amptEvent.X[j] , amptEvent.Y[j]));
        float cosdphis = cos(dphis);

        float dphip = dphi(phi(amptEvent.Px[i] , amptEvent.Py[i]) , phi(amptEvent.Px[j] , amptEvent.Py[j]));
        float cosdphip = cos(dphip);

        switch (pairType(amptEvent.ID[i] , amptEvent.ID[j])) {
        case 0:
          h_dphip_uu->Fill(dphip);
          h_dphis_uu->Fill(dphis);
          p_delta_p->Fill(0. , cosdphip);
          p_delta_s->Fill(0. , cosdphis);
          break;
        case 1:
          h_dphip_ud->Fill(dphip);
          h_dphis_ud->Fill(dphis);
          p_delta_p->Fill(1. , cosdphip);
          p_delta_s->Fill(1. , cosdphis);
          break;
        case 2:
          h_dphip_us->Fill(dphip);
          h_dphis_us->Fill(dphis);
          p_delta_p->Fill(2. , cosdphip);
          p_delta_s->Fill(2. , cosdphis);
          break;
        case 3:
          h_dphip_uubar->Fill(dphip);
          h_dphis_uubar->Fill(dphis);
          p_delta_p->Fill(6. , cosdphip);
          p_delta_s->Fill(6. , cosdphis);
          break;
        case 4:
          h_dphip_udbar->Fill(dphip);
          h_dphis_udbar->Fill(dphis);
          p_delta_p->Fill(7. , cosdphip);
          p_delta_s->Fill(7. , cosdphis);
          break;
        case 5:
          h_dphip_usbar->Fill(dphip);
          h_dphis_usbar->Fill(dphis);
          p_delta_p->Fill(8. , cosdphip);
          p_delta_s->Fill(8. , cosdphis);
      
          break;
        case 6:
          h_dphip_dd->Fill(dphip);
          h_dphis_dd->Fill(dphis);
          p_delta_p->Fill(3. , cosdphip);
          p_delta_s->Fill(3. , cosdphis);
          break;
        case 7:
          h_dphip_ds->Fill(dphip);
          h_dphis_ds->Fill(dphis);
          p_delta_p->Fill(4. , cosdphip);
          p_delta_s->Fill(4. , cosdphis);
          break;
        case 8:
          h_dphip_ddbar->Fill(dphip);
          h_dphis_ddbar->Fill(dphis);
          p_delta_p->Fill(9. , cosdphip);
          p_delta_s->Fill(9. , cosdphis);
          break;
        case 9:
          h_dphip_dsbar->Fill(dphip);
          h_dphis_dsbar->Fill(dphis);
          p_delta_p->Fill(10. , cosdphip);
          p_delta_s->Fill(10. , cosdphis);
          break;
        case 10:
          h_dphip_ss->Fill(dphip);
          h_dphis_ss->Fill(dphis);
          p_delta_p->Fill(5. , cosdphip);
          p_delta_s->Fill(5. , cosdphis);
          break;
        case 11:
          h_dphip_sbarbar->Fill(dphip);
          h_dphis_sbarbar->Fill(dphis);
          p_delta_p->Fill(11. , cosdphip);
          p_delta_s->Fill(11. , cosdphis);
          break;
        default:
          break;
        }

      }
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

  TCanvas* c1 = new TCanvas("c1" , "c1" , 800 , 800);
  c1->Divide(2 , 2);
  c1->cd(1);
  h2_xy->Draw("colz");
  c1->cd(2);
  h_z->Draw();
  c1->cd(3);
  h2_pxpy->Draw("colz");
  c1->cd(4);
  h_pz->Draw();

  TCanvas* c2 = new TCanvas("c2" , "c2" , 800 , 800);
  h_pdg->Draw();

  TCanvas* c3 = new TCanvas("c3" , "c3" , 800 , 800);
  h_nparton->Draw();

  TFile* outFile = new TFile("AMPTPartonDist.root" , "RECREATE");
  outFile->cd();
  h_nparton->Write();
  // h3_xyz->Write();
  h2_xy->Write();
  h_z->Write();
  // h3_pxpypz->Write();
  h2_pxpy->Write();
  h_pz->Write();
  h_pdg->Write();

  h_dphip_uu->Scale((float)h_dphip_uu->GetNbinsX() / h_dphip_uu->Integral());
  h_dphip_ud->Scale((float)h_dphip_ud->GetNbinsX() / h_dphip_ud->Integral());
  h_dphip_us->Scale((float)h_dphip_us->GetNbinsX() / h_dphip_us->Integral());
  h_dphip_uubar->Scale((float)h_dphip_uubar->GetNbinsX() / h_dphip_uubar->Integral());
  h_dphip_udbar->Scale((float)h_dphip_udbar->GetNbinsX() / h_dphip_udbar->Integral());
  h_dphip_usbar->Scale((float)h_dphip_usbar->GetNbinsX() / h_dphip_usbar->Integral());
  h_dphip_dd->Scale((float)h_dphip_dd->GetNbinsX() / h_dphip_dd->Integral());
  h_dphip_ds->Scale((float)h_dphip_ds->GetNbinsX() / h_dphip_ds->Integral());
  h_dphip_ddbar->Scale((float)h_dphip_ddbar->GetNbinsX() / h_dphip_ddbar->Integral());
  h_dphip_dsbar->Scale((float)h_dphip_dsbar->GetNbinsX() / h_dphip_dsbar->Integral());
  h_dphip_ss->Scale((float)h_dphip_ss->GetNbinsX() / h_dphip_ss->Integral());
  h_dphip_sbarbar->Scale((float)h_dphip_sbarbar->GetNbinsX() / h_dphip_sbarbar->Integral());

  h_dphip_uu->Write();
  h_dphip_ud->Write();
  h_dphip_us->Write();
  h_dphip_uubar->Write();
  h_dphip_udbar->Write();
  h_dphip_usbar->Write();
  h_dphip_dd->Write();
  h_dphip_ds->Write();
  h_dphip_ddbar->Write();
  h_dphip_dsbar->Write();
  h_dphip_ss->Write();
  h_dphip_sbarbar->Write();
  p_delta_p->Write();

  h_dphis_uu->Scale((float)h_dphis_uu->GetNbinsX() / h_dphis_uu->Integral());
  h_dphis_ud->Scale((float)h_dphis_ud->GetNbinsX() / h_dphis_ud->Integral());
  h_dphis_us->Scale((float)h_dphis_us->GetNbinsX() / h_dphis_us->Integral());
  h_dphis_uubar->Scale((float)h_dphis_uubar->GetNbinsX() / h_dphis_uubar->Integral());
  h_dphis_udbar->Scale((float)h_dphis_udbar->GetNbinsX() / h_dphis_udbar->Integral());
  h_dphis_usbar->Scale((float)h_dphis_usbar->GetNbinsX() / h_dphis_usbar->Integral());
  h_dphis_dd->Scale((float)h_dphis_dd->GetNbinsX() / h_dphis_dd->Integral());
  h_dphis_ds->Scale((float)h_dphis_ds->GetNbinsX() / h_dphis_ds->Integral());
  h_dphis_ddbar->Scale((float)h_dphis_ddbar->GetNbinsX() / h_dphis_ddbar->Integral());
  h_dphis_dsbar->Scale((float)h_dphis_dsbar->GetNbinsX() / h_dphis_dsbar->Integral());
  h_dphis_ss->Scale((float)h_dphis_ss->GetNbinsX() / h_dphis_ss->Integral());
  h_dphis_sbarbar->Scale((float)h_dphis_sbarbar->GetNbinsX() / h_dphis_sbarbar->Integral());
  

  h_dphis_uu->Write();
  h_dphis_ud->Write();
  h_dphis_us->Write();
  h_dphis_uubar->Write();
  h_dphis_udbar->Write();
  h_dphis_usbar->Write();
  h_dphis_dd->Write();
  h_dphis_ds->Write();
  h_dphis_ddbar->Write();
  h_dphis_dsbar->Write();
  h_dphis_ss->Write();
  h_dphis_sbarbar->Write();
  p_delta_s->Write();


  outFile->Close();
}