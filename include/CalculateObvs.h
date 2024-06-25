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
  //Observable delta gamma
  TProfile* p_delta_lam_pro[4];
  TProfile* p_gamma_lam_pro[4];
  TProfile* p_v2_pt_lam[2];
  TProfile* p_v2_pt_pro[2];

  public:
  CalculateObvs(std::string obvsFile) {
    file = std::unique_ptr<TFile>(new TFile(obvsFile.c_str(), "RECREATE"));
    h_mult = new TH1D("h_mult", "Multiplicity", 500, 0, 10000);
    h_pt = new TH1D("h_pt", "p_{T}", 100, 0, 10);
    h_eta = new TH1D("h_eta", "#eta", 100, -2.5, 2.5);
    h_phi = new TH1D("h_phi", "#phi", 100, -M_PI, M_PI);
    // FIXME: phi range
    // 0 -> lambda - proton
    // 1 -> lambda - anti-proton
    // 2 -> anti-lambda - proton
    // 3 -> anti-lambda - anti-proton
    std::vector<std::string> obvsName = {"#Lambda-p", "#Lambda-#bar{p}", "#bar{#Lambda}-p", "#bar{#Lambda}-#bar{p}"};
    for (int i = 0; i < 4; i++) {
      //FIXME
      p_delta_lam_pro[i] = new TProfile(Form("p_delta_lam_pro_%d", i), Form("#delta %s", obvsName[i].c_str()), 1, 0., 1.);
      p_gamma_lam_pro[i] = new TProfile(Form("p_gamma_lam_pro_%d", i), Form("#gamma %s", obvsName[i].c_str()), 1, 0., 1.);
    }
    p_v2_pt_lam[0] = new TProfile("p_v2_pt_lam", "v_{2} #Lambda", 25, 0.2, 5.2);
    p_v2_pt_lam[1] = new TProfile("p_v2_pt_pro", "v_{2} p", 25, 0.2, 5.2);
    p_v2_pt_pro[0] = new TProfile("p_v2_pt_antilam", "v_{2} #bar{#Lambda}", 25, 0.2, 5.2);
    p_v2_pt_pro[1] = new TProfile("p_v2_pt_antipro", "v_{2} #bar{p}", 25, 0.2, 5.2);
  }

  ~CalculateObvs() {
    file->cd();
    h_mult->Write();
    h_pt->Write();
    h_eta->Write();
    h_phi->Write();
    for (int i = 0; i < 4; i++) {
      p_delta_lam_pro[i]->Write();
      p_gamma_lam_pro[i]->Write();
    }
    for (int i = 0; i < 2; i++) {
      p_v2_pt_lam[i]->Write();
      p_v2_pt_pro[i]->Write();
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