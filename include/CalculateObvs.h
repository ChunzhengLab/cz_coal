#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "Particle.h"

class CalculateObvs {
  private:
  std::unique_ptr<TFile> file;
  //QA
  TH1D* h_mult;
  TH1D* h_pt;
  TH1D* h_eta;
  TH1D* h_phi;
  // spectra
  TH1D* h_pt_piLike[3];
  TH1D* h_pt_lam[2];
  TH1D* h_pt_pro[2];
  // flow
  TProfile* p_v2_pt_piLike[3];
  TProfile* p_v2_pt_lam[2];
  TProfile* p_v2_pt_pro[2];
  //Observable delta gamma
  TProfile* p_delta_lam_pro[4];
  TProfile* p_gamma_lam_pro[4];
  //Observable C(Δφ), C(sumφ)
  TH1D* h_dphi_lam_pro[4];
  TH1D* h_sphi_lam_pro[4];

  public:
  CalculateObvs(std::string obvsFile) {
    file = std::unique_ptr<TFile>(new TFile(obvsFile.c_str(), "RECREATE"));

    //QA
    h_mult = new TH1D("h_mult", "Multiplicity", 500, 0, 10000);
    h_pt = new TH1D("h_pt", "p_{T}", 100, 0, 10);
    h_eta = new TH1D("h_eta", "#eta", 100, -2.5, 2.5);
    h_phi = new TH1D("h_phi", "#phi", 100, -M_PI, M_PI);

    // spectra
    h_pt_piLike[0] = new TH1D("h_pt_piLike_0", "p_{T} #pi+ like", 200, 0, 10);
    h_pt_piLike[1] = new TH1D("h_pt_piLike_1", "p_{T} #pi0 like", 200, 0, 10);
    h_pt_piLike[2] = new TH1D("h_pt_piLike_2", "p_{T} #pi- like", 200, 0, 10);
    h_pt_lam[0] = new TH1D("h_pt_lam_0", "p_{T} #Lambda", 200, 0, 10);
    h_pt_lam[1] = new TH1D("h_pt_lam_1", "p_{T} #bar{#Lambda}", 200, 0, 10);
    h_pt_pro[0] = new TH1D("h_pt_pro_0", "p_{T} p", 200, 0, 10);
    h_pt_pro[1] = new TH1D("h_pt_pro_1", "p_{T} #bar{p}", 200, 0, 10);

    // elliptic flow
    p_v2_pt_piLike[0] = new TProfile("p_v2_pt_piLike_0", "v_{2} #pi+ like", 50, 0., 10.);
    p_v2_pt_piLike[1] = new TProfile("p_v2_pt_piLike_1", "v_{2} #pi0 like", 50, 0., 10.);
    p_v2_pt_piLike[2] = new TProfile("p_v2_pt_piLike_2", "v_{2} #pi- like", 50, 0., 10.);
    p_v2_pt_lam[0] = new TProfile("p_v2_pt_lam_0", "v_{2} #Lambda", 50, 0., 10.);
    p_v2_pt_lam[1] = new TProfile("p_v2_pt_lam_1", "v_{2} #bar{#Lambda}", 50, 0., 10.);
    p_v2_pt_pro[0] = new TProfile("p_v2_pt_pro_0", "v_{2} p", 50, 0., 10.);
    p_v2_pt_pro[1] = new TProfile("p_v2_pt_pro_1", "v_{2} #bar{p}", 50, 0., 10.);

    //Observable delta gamma
    // 0 -> lambda - proton
    // 1 -> lambda - anti-proton
    // 2 -> anti-lambda - proton
    // 3 -> anti-lambda - anti-proton
    std::vector<std::string> obvsName = {"#Lambda-p", "#Lambda-#bar{p}", "#bar{#Lambda}-p", "#bar{#Lambda}-#bar{p}"};
    for (int i = 0; i < 4; i++) {
      p_delta_lam_pro[i] = new TProfile(Form("p_delta_lam_pro_%d", i), Form("#delta %s", obvsName[i].c_str()), 1, 0., 1.);
      p_gamma_lam_pro[i] = new TProfile(Form("p_gamma_lam_pro_%d", i), Form("#gamma %s", obvsName[i].c_str()), 1, 0., 1.);
    }
    //Observable C(Δφ), C(sumφ)
    for (int i = 0; i < 4; i++) {
      h_dphi_lam_pro[i] = new TH1D(Form("h_dphi_lam_pro_%d", i), Form("C(#Delta#phi) %s", obvsName[i].c_str()), 100, -0.5*M_PI, 1.5*M_PI);
      h_sphi_lam_pro[i] = new TH1D(Form("h_sphi_lam_pro_%d", i), Form("C(s#phi) %s", obvsName[i].c_str()), 100, -0.5*M_PI, 1.5*M_PI);
    }
  }

  inline float RangeDPhi(float dphi) {
    while (dphi < -0.5 * M_PI) dphi += 2 * M_PI;
    while (dphi > 1.5 * M_PI) dphi -= 2 * M_PI;
    return dphi;
  }

  ~CalculateObvs() {
    file->cd();
    //QA
    h_mult->Write();
    h_pt->Write();
    h_eta->Write();
    h_phi->Write();
    // spectra
    for (int i = 0; i < 3; i++) h_pt_piLike[i]->Write();
    for (int i = 0; i < 2; i++) h_pt_lam[i]->Write();
    for (int i = 0; i < 2; i++) h_pt_pro[i]->Write();
    // elliptic flow
    for (int i = 0; i < 3; i++) p_v2_pt_piLike[i]->Write();
    for (int i = 0; i < 2; i++) p_v2_pt_lam[i]->Write();
    for (int i = 0; i < 2; i++) p_v2_pt_pro[i]->Write();
    //Observable delta gamma
    for (int i = 0; i < 4; i++) {
      p_delta_lam_pro[i]->Write();
      p_gamma_lam_pro[i]->Write();
    }
    //Observable C(Δφ), C(sumφ)
    for (int i = 0; i < 4; i++) {
      h_dphi_lam_pro[i]->Write();
      h_sphi_lam_pro[i]->Write();
    }
    file->Close();
    // NOTE：这里需要delete吗？
    // h_mult->Delete();
    // h_pt->Delete();
    // h_eta->Delete();
    // h_phi->Delete();
    // for (int i = 0; i < 4; i++) {
    //   p_delta_lam_pro[i]->Delete();
    //   p_gamma_lam_pro[i]->Delete();
    // }
  }

  void Process(std::vector<Hadron>& hadrons);
  void Print() {
    std::cout <<"--------------------------" << std::endl;
    std::cout << "CalculateObvs:" << std::endl;
    std::cout << "result saved in " << file->GetName() << std::endl;
    std::cout << "-----------------" << std::endl;
  }

  CalculateObvs(const CalculateObvs&) = delete;
  CalculateObvs& operator=(const CalculateObvs&) = delete;
};