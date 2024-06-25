#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <cmath>

class Particle
{
private:
  int nSerial; // unique serial number for the particle in this event
  int pdg; // PDG code
  float x, y, z; // position
  float px, py, pz; // momentum
public:
  //构造函数
  Particle(int nSerial = 0, int pdg = 0, float x = 0, float y = 0, float z = 0, float px = 0, float py = 0, float pz = 0)
      : nSerial(nSerial), pdg(pdg), x(x), y(y), z(z), px(px), py(py), pz(pz) {}
  virtual ~Particle() {}
  void SetSerial(int nSerial) { this->nSerial = nSerial; }
  int GetSerial() const { return nSerial; }
  void SetMomentum(float px, float py, float pz) { this->px = px; this->py = py; this->pz = pz; }
  void SetPosition(float x, float y, float z) { this->x = x; this->y = y; this->z = z; }
  void GetMomentum(float &px, float &py, float &pz) const { px = this->px; py = this->py; pz = this->pz; }
  void GetPosition(float &x, float &y, float &z) const { x = this->x; y = this->y; z = this->z; }
  void SetPDG(int pdg) { this->pdg = pdg; }
  int PDG() const { return pdg; }
  float X() const { return x; }
  float Y() const { return y; }
  float Z() const { return z; }
  float Pt() const { return sqrt(px * px + py * py); }
  float Eta() const { 
    float p = sqrt(px * px + py * py + pz * pz);
    return 0.5 * log((p + pz) / (p - pz));
  }
  float Phi() const { return atan2(py, px); }
  float PhiS() const { return atan2(y, x); }

  virtual void Print() const {
    std::cout << "This is particle " << nSerial << " with PDG code " << pdg << std::endl;
    std::cout << "Position: (" << x << ", " << y << ", " << z << ")" << std::endl;
    std::cout << "Momentum: (" << px << ", " << py << ", " << pz << ")" << std::endl;
  }

};

// Parton和Hadron的实现
class Parton : public Particle {
  public:
  //构造函数
  Parton(int nSerial = 0, int pdg = 0, float x = 0, float y = 0, float z = 0, float px = 0, float py = 0, float pz = 0)
      : Particle(nSerial, pdg, x, y, z, px, py, pz) {}
};

class Hadron : public Particle {
  private:
  float distance; //距离
  int nSerial0, nSerial1, nSerial2; //三个parton的序号
  float x0, y0, z0; // parton0的位置
  float x1, y1, z1; // parton1的位置
  float x2, y2, z2; // parton2的位置

  public:
  //构造函数
  Hadron(int nSerial = 0, int pdg = 0, float x = 0, float y = 0, float z = 0, float px = 0, float py = 0, float pz = 0, float distance = 0, int nSerial0 = 0, int nSerial1 = 0, int nSerial2 = 0)
      : Particle(nSerial, pdg, x, y, z, px, py, pz), distance(distance), nSerial0(nSerial0), nSerial1(nSerial1), nSerial2(nSerial2) {}
  float Distance() const { return distance; }
  void GetPartonSerials(int &nSerial0, int &nSerial1, int &nSerial2) const { nSerial0 = this->nSerial0; nSerial1 = this->nSerial1; nSerial2 = this->nSerial2; }
  void GetPartonSerials(int &nSerial0, int &nSerial1) const { nSerial0 = this->nSerial0; nSerial1 = this->nSerial1; }
  int GetPartonSerial0() const { return nSerial0; }
  int GetPartonSerial1() const { return nSerial1; }
  int GetPartonSerial2() const { return nSerial2; }
  void SetParton0Position(float x0, float y0, float z0) { this->x0 = x0; this->y0 = y0; this->z0 = z0; }
  void SetParton1Position(float x1, float y1, float z1) { this->x1 = x1; this->y1 = y1; this->z1 = z1; }
  void SetParton2Position(float x2, float y2, float z2) { this->x2 = x2; this->y2 = y2; this->z2 = z2; }
  void GetParton0Position(float &x0, float &y0, float &z0) const { x0 = this->x0; y0 = this->y0; z0 = this->z0; }
  void GetParton1Position(float &x1, float &y1, float &z1) const { x1 = this->x1; y1 = this->y1; z1 = this->z1; }
  void GetParton2Position(float &x2, float &y2, float &z2) const { x2 = this->x2; y2 = this->y2; z2 = this->z2; }

  //距离比较函数
  bool operator>(const Hadron &rhs) const {
    return distance > rhs.distance;
  }

  void Print() const {
    std::cout << "This is hadron " << GetSerial() << " with PDG code " << PDG() << std::endl;
    std::cout << "Distance: " << distance << std::endl;
    std::cout << "Formed by parton " << nSerial0 << ", parton " << nSerial1 << " and parton " << nSerial2 << std::endl;
  }
};

#endif // PARTICLE_H