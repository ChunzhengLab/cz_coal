#include "CalculateObvs.h"
#include "Event.h"
#include "TBits.h"

void CalculateObvs::Process(std::vector<Hadron>& hadrons) {
  h_mult->Fill(hadrons.size());
  for (int i = 0; i < hadrons.size(); i++) {
    //只要lambda 和 proton
    float pt = hadrons[i].Pt();
    float eta = hadrons[i].Eta();
    float phi = hadrons[i].Phi();
    h_pt->Fill(pt);
    h_eta->Fill(eta);
    h_phi->Fill(phi);

    if(pt < 0.2 || pt > 10.0) continue;
    if(abs(eta) > 0.8) continue;

    if (hadrons[i].PDG() == 211) {  // pi+
      h_pt_piLike[0]->Fill(pt);
      p_v2_pt_piLike[0]->Fill(pt, cos(2. * phi));
    } else if (hadrons[i].PDG() == -211) { // pi-
      h_pt_piLike[1]->Fill(pt);
      p_v2_pt_piLike[1]->Fill(pt, cos(2. * phi));
    } else if (hadrons[i].PDG() == 111) { // pi0
      h_pt_piLike[2]->Fill(pt);
      p_v2_pt_piLike[2]->Fill(pt, cos(2. * phi));
    } else if (hadrons[i].PDG() == 3122) { // lambda
      h_pt_lam[0]->Fill(pt);
      p_v2_pt_lam[0]->Fill(pt, cos(2. * phi));
    } else if (hadrons[i].PDG() == -3122) { // anti-lambda
      h_pt_lam[1]->Fill(pt);
      p_v2_pt_lam[1]->Fill(pt, cos(2. * phi));
    } else if (hadrons[i].PDG() == 2212) { // proton
      h_pt_pro[0]->Fill(pt);
      p_v2_pt_pro[0]->Fill(pt, cos(2. * phi));
    } else if (hadrons[i].PDG() == -2212) { // anti-proton
      h_pt_pro[1]->Fill(pt);
      p_v2_pt_pro[1]->Fill(pt, cos(2. * phi));
    }

    if (abs(hadrons[i].PDG()) != 3122) continue; // 第一个粒子只选lambda(anti-lambda)
    for (int j = i + 1; j < hadrons.size(); j++) {
      float pt_j = hadrons[j].Pt();
      float eta_j = hadrons[j].Eta();
      float phi_j = hadrons[j].Phi();

      if (abs(hadrons[j].PDG()) != 2212) continue; // 第二个粒子只选proton(anti-proton)
      if(pt_j < 0.2 || pt_j > 10.0) continue;
      if(abs(eta_j) > 0.8) continue;

      //delta = <cos(phi_0 - phi_1)>
      //gamma = <cos(phi_0 + phi_1)>
      float delta = cos(phi - phi_j);
      float gamma = cos(hadrons[i].Phi() + hadrons[j].Phi());

      TBits bits(4);
      bits.SetBitNumber(0, hadrons[i].PDG() == 3122 && hadrons[j].PDG() == 2212);
      bits.SetBitNumber(1, hadrons[i].PDG() == 3122 && hadrons[j].PDG() == -2212);
      bits.SetBitNumber(2, hadrons[i].PDG() == -3122 && hadrons[j].PDG() == 2212);
      bits.SetBitNumber(3, hadrons[i].PDG() == -3122 && hadrons[j].PDG() == -2212);

      for (int iBit = 0; iBit < 4; iBit++) {
        if (bits.TestBitNumber(iBit)) {
          p_delta_lam_pro[iBit]->Fill(0.5, delta);
          p_gamma_lam_pro[iBit]->Fill(0.5, gamma);
          h_dphi_lam_pro[iBit]->Fill(RangeDPhi(hadrons[i].Phi() - hadrons[j].Phi()));
          h_sphi_lam_pro[iBit]->Fill(RangeDPhi(hadrons[i].Phi() + hadrons[j].Phi()));
        }
      }
    }
  }
}