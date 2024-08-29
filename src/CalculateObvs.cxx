#include "CalculateObvs.h"
#include "Event.h"
#include "TBits.h"
#include "Par.h"
#include <TLorentzVector.h>


namespace par {
  std::map<int, float> pdgToBinCenter = {
    //pi
    {211, 0.5},
    {-211, 1.5},
    {111, 2.5},
    //K
    {321, 3.5},
    {-321, 4.5},
    {311, 5.5},
    {-311, 6.5},
    //rho
    {213, 7.5},
    {-213, 8.5},
    {113, 9.5},
    //eta
    {221, 10.5},
    //omega
    {223, 11.5},
    //phi
    {333, 12.5},
    //p
    {2212, 13.5},
    {-2212, 14.5},
    //n
    {2112, 15.5},
    {-2112, 16.5},
    //Delta
    {2224, 17.5},
    {-2224, 18.5},
    {2214, 19.5},
    {-2214, 20.5},
    {2114, 21.5},
    {-2114, 22.5},
    {1114, 23.5},
    {-1114, 24.5},
    //Lambda
    {3122, 25.5},
    {-3122, 26.5},
    //Sigma
    {3222, 27.5},
    {-3222, 28.5},
    {3212, 29.5},
    {-3212, 30.5},
    {3112, 31.5},
    {-3112, 32.5},
    //Xi
    {3322, 33.5},
    {-3322, 34.5},
    {3312, 35.5},
    {-3312, 36.5},
    //Omega
    {3334, 37.5},
    {-3334, 38.5}
  };
}

void CalculateObvs::Process(std::vector<Hadron>& hadrons) {
  h_mult->Fill(hadrons.size());
  for (int i = 0; i < hadrons.size(); i++) {
    //只要lambda 和 proton
    float pt = hadrons[i].Pt();
    float eta = hadrons[i].Eta();
    float phi = hadrons[i].Phi();
    int pdg = hadrons[i].PDG();
    h_pt->Fill(pt);
    h_eta->Fill(eta);
    h_phi->Fill(phi);
    if (par::pdgToBinCenter.find(pdg) != par::pdgToBinCenter.end()) {
      h_pid->Fill(par::pdgToBinCenter[pdg]);
    }

    if(pt < 0.2 || pt > 10.0) continue;
    if(abs(eta) > 0.8) continue;

    if (pdg == 211) {  // pi+
      h_pt_piLike[0]->Fill(pt);
      p_v2_pt_piLike[0]->Fill(pt, cos(2. * phi));
    } else if (pdg == -211) { // pi-
      h_pt_piLike[1]->Fill(pt);
      p_v2_pt_piLike[1]->Fill(pt, cos(2. * phi));
    } else if (pdg == 111) { // pi0
      h_pt_piLike[2]->Fill(pt);
      p_v2_pt_piLike[2]->Fill(pt, cos(2. * phi));
    } else if (pdg == 3122) { // lambda
      h_pt_la[0]->Fill(pt);
      p_v2_pt_la[0]->Fill(pt, cos(2. * phi));
    } else if (pdg == -3122) { // anti-lambda
      h_pt_la[1]->Fill(pt);
      p_v2_pt_la[1]->Fill(pt, cos(2. * phi));
    } else if (pdg == 2112) { // proton
      h_pt_pr[0]->Fill(pt);
      p_v2_pt_pr[0]->Fill(pt, cos(2. * phi));
    } else if (pdg == -2112) { // anti-proton
      h_pt_pr[1]->Fill(pt);
      p_v2_pt_pr[1]->Fill(pt, cos(2. * phi));
    }

    if (abs(pdg) != 3122) continue; // 第一个粒子只选lambda(anti-lambda)
    //快度Lambda |y| < 0.5
    TLorentzVector lv;
    lv.SetPtEtaPhiM(pt, eta, phi, 1.115683);
    float y = lv.Rapidity();
    if (abs(y) > 0.5) continue;

    for (int j = 0; j < hadrons.size(); j++) {
      if (i == j) continue;
      if (abs(hadrons[j].PDG()) != 211 && abs(hadrons[j].PDG()) != 2112 && abs(hadrons[j].PDG()) != 3122) continue;
      float pt_j = hadrons[j].Pt();
      float eta_j = hadrons[j].Eta();
      float phi_j = hadrons[j].Phi();
      float pdg_j = hadrons[j].PDG();
      // 第二个粒子只选pion(anti-pion)和proton(anti-proton)和lambda(anti-lambda)
      if(pt_j < 0.2 || pt_j > 10.0) continue;
      if(abs(eta_j) > 0.8) continue;

      //delta = <cos(phi_0 - phi_1)>
      //gamma = <cos(phi_0 + phi_1)>
      float delta = cos(phi - phi_j);
      float gamma = cos(phi + phi_j);

      TBits bits_la_pi(4);
      TBits bits_la_pr(4);
      TBits bits_la_la(4);
      // lambda - pion
      bits_la_pi.SetBitNumber(0, pdg == 3122 && pdg_j == 211);
      bits_la_pi.SetBitNumber(1, pdg == 3122 && pdg_j == -211);
      bits_la_pi.SetBitNumber(2, pdg == -3122 && pdg_j == 211);
      bits_la_pi.SetBitNumber(3, pdg == -3122 && pdg_j == -211);
      // lambda - proton
      bits_la_pr.SetBitNumber(0, pdg == 3122 && pdg_j == 2112);
      bits_la_pr.SetBitNumber(1, pdg == 3122 && pdg_j == -2112);
      bits_la_pr.SetBitNumber(2, pdg == -3122 && pdg_j == 2112);
      bits_la_pr.SetBitNumber(3, pdg == -3122 && pdg_j == -2112);
      // lambda - lambda
      bits_la_la.SetBitNumber(0, pdg == 3122 && pdg_j == 3122);
      bits_la_la.SetBitNumber(1, pdg == 3122 && pdg_j == -3122);
      bits_la_la.SetBitNumber(2, pdg == -3122 && pdg_j == 3122);
      bits_la_la.SetBitNumber(3, pdg == -3122 && pdg_j == -3122);

      for (int iBit = 0; iBit < 4; iBit++) {
        if (bits_la_pi.TestBitNumber(iBit)) {
          p_delta_la_pi[iBit]->Fill(0.5, delta);
          p_gamma_la_pi[iBit]->Fill(0.5, gamma);
          h_dphi_la_pi[iBit]->Fill(RangeDPhi(phi - phi_j));
          h_sphi_la_pi[iBit]->Fill(RangeDPhi(phi + phi_j));
        }
        if (bits_la_pr.TestBitNumber(iBit)) {
          p_delta_la_pr[iBit]->Fill(0.5, delta);
          p_gamma_la_pr[iBit]->Fill(0.5, gamma);
          h_dphi_la_pr[iBit]->Fill(RangeDPhi(phi - phi_j));
          h_sphi_la_pr[iBit]->Fill(RangeDPhi(phi + phi_j));
        }
        if (bits_la_la.TestBitNumber(iBit)) {
          p_delta_la_la[iBit]->Fill(0.5, delta);
          p_gamma_la_la[iBit]->Fill(0.5, gamma);
          h_dphi_la_la[iBit]->Fill(RangeDPhi(phi - phi_j));
          h_sphi_la_la[iBit]->Fill(RangeDPhi(phi + phi_j));
        }
      }
    }
  }
}