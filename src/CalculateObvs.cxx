#include "CalculateObvs.h"
#include "Event.h"
#include "TBits.h"

void CalculateObvs::Process(std::vector<Hadron>& hadrons) {
  h_mult->Fill(hadrons.size());
  for (int i = 0; i < hadrons.size(); i++) {
    //只要lambda 和 proton
    h_pt->Fill(hadrons[i].Pt());
    h_eta->Fill(hadrons[i].Eta());
    h_phi->Fill(hadrons[i].Phi());

    if (hadrons[i].PDG() == 3122) {
      p_v2_pt_lam[0]->Fill(hadrons[i].Pt(), cos(2. * hadrons[i].Phi()));
    } else if (hadrons[i].PDG() == -3122) {
      p_v2_pt_lam[1]->Fill(hadrons[i].Pt(), cos(2. * hadrons[i].Phi()));
    } else if (hadrons[i].PDG() == 2212) {
      p_v2_pt_pro[0]->Fill(hadrons[i].Pt(), cos(2. * hadrons[i].Phi()));
    } else if (hadrons[i].PDG() == -2212) {
      p_v2_pt_pro[1]->Fill(hadrons[i].Pt(), cos(2. * hadrons[i].Phi()));
    }

    if (abs(hadrons[i].PDG()) != 3122) continue;

    for (int j = i + 1; j < hadrons.size(); j++) {
      if (abs(hadrons[j].PDG()) != 2212) continue;
      //delta = <cos(phi_0 - phi_1)>
      //gamma = <cos(phi_0 + phi_1)>
      float delta = cos(hadrons[i].Phi() - hadrons[j].Phi());
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
        }
      }
    }
  }
}