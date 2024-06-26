#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TColor.h"
#include "TLegend.h"


void Draw() {
  int ci[4];
  TColor* color[4];
  ci[0] = TColor::GetFreeColorIndex();
  color[0] = new TColor(ci[0],   0/255.,  24/255., 113/255.);//dark blue
  ci[1] = TColor::GetFreeColorIndex();
  color[1] = new TColor(ci[1], 255/255.,  88/255.,  93/255.);//red
  ci[2] = TColor::GetFreeColorIndex();
  color[2] = new TColor(ci[2], 255/255., 181/255.,  73/255.);//yellow
  ci[3] = TColor::GetFreeColorIndex();
  color[3] = new TColor(ci[3], 65/255.,  182/255., 230/255.);//light blue

  TFile *f = new TFile("obvs_coalHadrons-FinalMerged.root");

  TH1D* h_mult = (TH1D*)f->Get("h_mult");
  TH1D* h_pt = (TH1D*)f->Get("h_pt");
  TH1D* h_eta = (TH1D*)f->Get("h_eta");
  TH1D* h_phi = (TH1D*)f->Get("h_phi");
  //Observable delta gamma

  TH1D* h_delta_lam_pro[4];
  TH1D* h_gamma_lam_pro[4];

  TProfile* p_delta_lam_pro[4];
  TProfile* p_gamma_lam_pro[4];
  TProfile* p_v2_pt_lam[2];
  TProfile* p_v2_pt_pro[2];

  for (int i = 0; i < 4; i++) {
    p_delta_lam_pro[i] = (TProfile*)f->Get(Form("p_delta_lam_pro_%d", i));
    h_delta_lam_pro[i] = new TH1D(Form("h_delta_lam_pro_%d", i), "" , 4, 0., 4.);
    h_delta_lam_pro[i]->SetBinContent(i+1, p_delta_lam_pro[i]->GetBinContent(1));
    h_delta_lam_pro[i]->SetBinError(i+1, p_delta_lam_pro[i]->GetBinError(1));
    h_delta_lam_pro[i]->SetMarkerStyle(20);
    h_delta_lam_pro[i]->SetMarkerSize(1.3);
    h_delta_lam_pro[i]->SetMarkerColor(ci[i]);
    h_delta_lam_pro[i]->SetLineColor(ci[i]);

    p_gamma_lam_pro[i] = (TProfile*)f->Get(Form("p_gamma_lam_pro_%d", i));
    h_gamma_lam_pro[i] = new TH1D(Form("h_gamma_lam_pro_%d", i), "" , 4, 0., 4.);
    h_gamma_lam_pro[i]->SetBinContent(i+1, p_gamma_lam_pro[i]->GetBinContent(1));
    h_gamma_lam_pro[i]->SetBinError(i+1, p_gamma_lam_pro[i]->GetBinError(1));
    h_gamma_lam_pro[i]->SetMarkerStyle(20);
    h_gamma_lam_pro[i]->SetMarkerSize(1.3);
    h_gamma_lam_pro[i]->SetMarkerColor(ci[i]);
    h_gamma_lam_pro[i]->SetLineColor(ci[i]);
  }

  p_v2_pt_lam[0] = (TProfile*)f->Get("p_v2_pt_lam");
  p_v2_pt_pro[0] = (TProfile*)f->Get("p_v2_pt_pro");
  p_v2_pt_lam[1] = (TProfile*)f->Get("p_v2_pt_antilam");
  p_v2_pt_pro[1] = (TProfile*)f->Get("p_v2_pt_antipro");

  for (int i = 0; i < 2; i++) {
    p_v2_pt_lam[i]->SetMarkerColor(kRed);
    p_v2_pt_pro[i]->SetLineColor(kBlue);
    p_v2_pt_pro[i]->SetMarkerColor(kBlue);
    p_v2_pt_lam[i]->SetLineColor(kRed);
  }

  p_v2_pt_lam[0]->SetMarkerStyle(kFullCircle);
  p_v2_pt_pro[0]->SetMarkerStyle(kFullCircle);
  p_v2_pt_lam[1]->SetMarkerStyle(kOpenSquare);
  p_v2_pt_pro[1]->SetMarkerStyle(kOpenSquare);

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->SetNColumns(2);
  leg->AddEntry(p_v2_pt_lam[0], "#Lambda", "pl");
  leg->AddEntry(p_v2_pt_pro[0], "p", "pl");
  leg->AddEntry(p_v2_pt_lam[1], "#bar{#Lambda}", "pl");
  leg->AddEntry(p_v2_pt_pro[1], "#bar{p}", "pl");

  TCanvas *c1 = new TCanvas("cflow","cflow",400,400);
  c1->cd()->DrawFrame(0., 0., 3, 0.2, ";p_{T} [GeV/c];v_{2}");
  c1->cd()->SetGrid();
  for (int i = 0; i < 2; i++) {
    p_v2_pt_lam[i]->Draw("same P");
    p_v2_pt_pro[i]->Draw("same P");
  }
  leg->Draw();

  TCanvas *c2 = new TCanvas("cObv","cObv",800,400);
  c2->Divide(2,1);
  c2->cd(1);
  c2->cd(1)->SetGrid();
  TH2D* dummyDelta = new TH2D("dummyDelta", ";;#delta", 4, 0., 4.,1,-5e-4,5e-4);
  dummyDelta->SetStats(0);
  dummyDelta->GetXaxis()->SetBinLabel(1, "#Lambda-p");
  dummyDelta->GetXaxis()->SetBinLabel(2, "#Lambda-#bar{p}");
  dummyDelta->GetXaxis()->SetBinLabel(3, "#bar{#Lambda}-p");
  dummyDelta->GetXaxis()->SetBinLabel(4, "#bar{#Lambda}-#bar{p}");
  dummyDelta->Draw();

  c2->cd(2)->SetGrid();
  TH2D* dummyGamma = new TH2D("dummyGamma", ";;#gamma", 4, 0., 4.,1,-5e-4,5e-4);
  dummyGamma->SetStats(0);
  dummyGamma->GetXaxis()->SetBinLabel(1, "#Lambda-p");
  dummyGamma->GetXaxis()->SetBinLabel(2, "#Lambda-#bar{p}");
  dummyGamma->GetXaxis()->SetBinLabel(3, "#bar{#Lambda}-p");
  dummyGamma->GetXaxis()->SetBinLabel(4, "#bar{#Lambda}-#bar{p}");
  dummyGamma->Draw();
  for (int i = 0; i < 4; i++) {
    c2->cd(1);
    h_delta_lam_pro[i]->Draw("same");
    c2->cd(2);
    h_gamma_lam_pro[i]->Draw("same");
  }

}