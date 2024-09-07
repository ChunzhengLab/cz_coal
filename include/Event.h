#ifndef EVENT_H
#define EVENT_H
#include "Par.h"
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
  std::vector<T> Particles;
public:
  Event(): nSerial(0) {}
  Event(int nSerial, std::vector<Hadron>&& particles): nSerial(nSerial), Particles(std::move(particles)) {}
  Event(const PartonEventStruct& partonEventStruct) {
    nSerial = partonEventStruct.nevent;
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
  void wash() {
    //如果存在夸克的PDG是0（代表没有PDG信息），那么直接移除并报错
    Particles.erase(std::remove_if(Particles.begin(), Particles.end(), [](const T& p) { if(p.PDG() == 0) std::cerr << p.GetSerial() << "th parton has PDG code 0, which is invalid." << std::endl; return p.PDG() == 0; }), Particles.end());

    // 如果是重夸克c = 4, b = 5, t = 6，直接移除这个quark，用Lambda表达式，std::remove_if
    if(par::isRemoveHFQuarks) {
      Particles.erase(std::remove_if(Particles.begin(), Particles.end(), [](const T& p) { return abs(p.PDG()) > 3; }), Particles.end());
    }

    // 如果不启用quark move on，那么将所有parton的时间设置为0, 也就是所有parton都在同一时刻
    if(!par::isEnableQuarkMoveOn) {
      for (int i = 0; i < Particles.size(); i++) Particles[i].SetTime(0);
    }

    // 如果是要求忘记Z坐标或者LocalDraw，那么将所有parton的Z坐标设置为0
    if(par::isForgetZ || par::isCreateAnimation) {
      for (int i = 0; i < Particles.size(); i++) Particles[i].SetPosition(Particles[i].X(), Particles[i].Y(), 0);
    }

    //如果要求平衡夸克数，那么将正负夸克数平衡, 这里借用了LabalAsUsed(),希望以后可以改进！
    if (par::isBalanceQuarkNumber) {
      // 如果启用了平衡夸克数，那么将正负夸克数平衡
      // 1. 保存多余的正夸克或者负夸克的index
      std::vector<int> indexPositiveQuark, indexNegativeQuark;
      for (int i = 0; i < Particles.size(); i++) {
        if (Particles[i].PDG() > 0) indexPositiveQuark.push_back(i);
        else if (Particles[i].PDG() < 0) indexNegativeQuark.push_back(i);
        else std::cerr << "Warning: Quark PDG code is 0, please check the input." << std::endl;
      }

      // 2. 统计正负夸克数
      int nPositiveQuark = indexPositiveQuark.size();
      int nNegativeQuark = indexNegativeQuark.size();
      int nQuarkToBeRemoved = std::abs(nPositiveQuark - nNegativeQuark);

      // 3. 随机删除数量多的正夸克或者负夸克
      if (nPositiveQuark > nNegativeQuark) {
        std::shuffle(indexPositiveQuark.begin(), indexPositiveQuark.end(), par::gen);
        // 只标记前 nQuarkToBeRemoved 个多余的正夸克
        for (int i = 0; i < nQuarkToBeRemoved; ++i) {
            Particles[indexPositiveQuark[i]].LabelAsUsed();
        }
      } else if (nPositiveQuark < nNegativeQuark) {
        std::shuffle(indexNegativeQuark.begin(), indexNegativeQuark.end(), par::gen);
        // 只标记前 nQuarkToBeRemoved 个多余的负夸克
        for (int i = 0; i < nQuarkToBeRemoved; ++i) {
            Particles[indexNegativeQuark[i]].LabelAsUsed();
        }
      }
      Particles.erase(std::remove_if(Particles.begin(), Particles.end(), [](const T& p) { return p.IsUsed(); }), Particles.end());
    }
  }


  void Print() {
    std::cout << "This is event " << nSerial << " with " << Particles.size() << " particles" << std::endl;
  }
};

#endif // EVENT_H


