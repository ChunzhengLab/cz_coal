#include "EventsWriter.h"
#include <TTree.h>
#include "Par.h"

void EventsWriter::Event2HadronEventStruct(const Event<Hadron>& event) {
  //清空之前事件的信息
  hadronEventStruct.Clear();
  //填充新的事件信息
  hadronEventStruct.nSeries = event.GetSerial();
  hadronEventStruct.nTracks = event.GetNTrks();
  std::vector<Hadron> particles = event.GetParticles();
  for (auto particle : particles) {
    hadronEventStruct.PDG.emplace_back(particle.PDG());
    hadronEventStruct.Pt.emplace_back(particle.Pt());
    hadronEventStruct.Eta.emplace_back(particle.Eta());
    hadronEventStruct.Phi.emplace_back(particle.Phi());
    hadronEventStruct.x.emplace_back(particle.X());
    hadronEventStruct.y.emplace_back(particle.Y());
    hadronEventStruct.z.emplace_back(particle.Z());
    hadronEventStruct.dis.emplace_back(particle.Distance());
    hadronEventStruct.quark0.emplace_back(particle.GetPartonSerial0());
    hadronEventStruct.quark1.emplace_back(particle.GetPartonSerial1());
    hadronEventStruct.quark2.emplace_back(particle.GetPartonSerial2());
    if (par::isDebug) {
      float x, y, z;
      particle.GetParton0Position(x, y, z);
      hadronEventStruct.quark0_x.emplace_back(x);
      hadronEventStruct.quark0_y.emplace_back(y);
      hadronEventStruct.quark0_z.emplace_back(z);
      particle.GetParton1Position(x, y, z);
      hadronEventStruct.quark1_x.emplace_back(x);
      hadronEventStruct.quark1_y.emplace_back(y);
      hadronEventStruct.quark1_z.emplace_back(z);
      particle.GetParton2Position(x, y, z);
      hadronEventStruct.quark2_x.emplace_back(x);
      hadronEventStruct.quark2_y.emplace_back(y);
      hadronEventStruct.quark2_z.emplace_back(z);
    }
  }
}

void EventsWriter::WriteEvent(const Event<Hadron>& event) {
  Event2HadronEventStruct(event);
  tree->Fill();
}

