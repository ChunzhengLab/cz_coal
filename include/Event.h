#ifndef EVENT_H
#define EVENT_H
#include "Particle.h"

// AMPT Event
struct PartonEventStruct {
  int nevent;
  int nparton;
  int   ID[30000];   //[nparton]
  float Px[30000], Py[30000], Pz[30000];   //[nparton]
  float X[30000], Y[30000], Z[30000];   //[nparton]
  float Time[30000];   //[nparton]
};

struct HadronEventStruct {
  int nSeries;
  int nTracks;
  std::vector<int> PDG;
  std::vector<float> Pt, Eta, Phi;
  std::vector<float> x, y, z;
  std::vector<float> dis;
  std::vector<float> time;
  std::vector<int> quark0, quark1, quark2;
  std::vector<float> quark0_x, quark0_y, quark0_z;
  std::vector<float> quark1_x, quark1_y, quark1_z;
  std::vector<float> quark2_x, quark2_y, quark2_z;

  void Clear() {
    nSeries = 0;
    nTracks = 0;
    PDG.clear();
    Pt.clear(), Eta.clear(), Phi.clear();
    x.clear(), y.clear(), z.clear();
    dis.clear();
    time.clear();
    quark0.clear(), quark1.clear(), quark2.clear();
    quark0_x.clear(), quark0_y.clear(), quark0_z.clear();
    quark1_x.clear(), quark1_y.clear(), quark1_z.clear();
    quark2_x.clear(), quark2_y.clear(), quark2_z.clear();
  }
};

template <typename T>
class Event {
private:
  //构造函数
  int nSerial; // unique serial number for the event, 从1开始
  int nTrks;
  std::vector<T> Particles;
public:
  Event(): nSerial(0), nTrks(0) {}
  Event(int nSerial, int nTrks, std::vector<Hadron>&& particles): nSerial(nSerial), nTrks(nTrks), Particles(std::move(particles)) {}
  Event(const PartonEventStruct& partonEventStruct) {
    nSerial = partonEventStruct.nevent;
    nTrks = partonEventStruct.nparton;
    Particles.reserve(partonEventStruct.nparton);
    for (int i = 0; i < partonEventStruct.nparton; i++) {
        Particles.emplace_back(i, partonEventStruct.ID[i], 
                   partonEventStruct.X[i], partonEventStruct.Y[i], partonEventStruct.Z[i],
                   partonEventStruct.Px[i], partonEventStruct.Py[i], partonEventStruct.Pz[i],
                   partonEventStruct.Time[i]);
    }
  }
  ~Event() {}
  const std::vector<T>& GetParticles() const { return Particles; }
  void SetSerial(int nSerial) { this->nSerial = nSerial; }
  int GetSerial() const { return nSerial; }
  void SetNTrks(int nTrks) { this->nTrks = nTrks; }
  int GetNTrks() const { return nTrks; }

  void Print() {
    std::cout << "This is event " << nSerial << " with " << nTrks << " particles" << std::endl;
  }
};

#endif // EVENT_H


