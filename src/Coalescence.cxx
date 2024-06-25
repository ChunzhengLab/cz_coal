#include <queue>
#include <set>
#include <iostream>
#include "Coalescence.h"
#include "DistanceFun.h"
#include "Par.h"
#include <random>
#include "TStopwatch.h"

bool Coalescence::initialized = false;
std::map<BaryonCombination, int> Coalescence::baryonLookupTable;
std::map<MesonCombination, int> Coalescence::mesonLookupTable;

void Coalescence::Process(std::vector<Parton> const &partons, std::vector<Hadron> &hadrons) {
  switch (coalescenceAlgorithm) {
    case CoalescenceAlgorithm::kClassic:
      ProcessClassic(partons, hadrons);
      break;
    case CoalescenceAlgorithm::kFromParton:
      ProcessFromParton(partons, hadrons);
      break;
    default:
      std::cerr << "Unknown coalescence algorithm" << std::endl;
      break;
  }
}

void Coalescence::ProcessClassic(std::vector<Parton> const &partons, std::vector<Hadron> &hadrons) {
  TStopwatch timer;
  timer.Start();

  if(par::isDebug) std::cout<<"How many partons in this event for coalescence: "<<partons.size()<<std::endl;
  
  // 将Hadron对象放入优先队列，使用优先队列来排序距离
  std::priority_queue<Hadron, std::vector<Hadron>, std::greater<>> queHadronCadidates;

  //计算所有parton之间的空间距离 -> 介子
  int nSerialHadron = 0;
  for (int i = 0; i < partons.size(); i++) {
    int pdg_quark_0 = partons[i].PDG();
    int nSerial_0 = partons[i].GetSerial();
    float x0, y0, z0;
    float px0, py0, pz0;
    partons[i].GetPosition(x0, y0, z0);
    partons[i].GetMomentum(px0, py0, pz0);

    // 从i+1开始，避免重复计算
    for (int j = i + 1; j < partons.size(); j++) {
      int pdg_quark_1 = partons[j].PDG();
      //如果前两个夸克同号，则已经不可能是介子，可以直接跳过
      //采用异或运算符判断两个数是否异号
      if (pdg_quark_0 * pdg_quark_1 > 0) {
        continue;
      }
      int nSerial_1 = partons[j].GetSerial();

      int pdg_meson = 0;
      pdg_meson = LookupSpecies(partons[i].PDG(), partons[j].PDG(), 0);
      if (pdg_meson == 0) continue;

      float x1, y1, z1;
      float px1, py1, pz1;
      partons[j].GetPosition(x1, y1, z1);
      partons[j].GetMomentum(px1, py1, pz1);

      float d_meson = distance3D(x0, y0, z0, x1, y1, z1);
      // 介子不乘系数
      float d = d_meson;

      nSerialHadron++;
      Hadron hadronCadi(nSerialHadron, pdg_meson, (x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2, (px0 + px1) / 2, (py0 + py1) / 2, (pz0 + pz1) / 2, d, nSerial_0, nSerial_1, 0);
      queHadronCadidates.push(hadronCadi);
      // std::cout << "Distance between parton " << nSerial_0 << " and parton " << nSerial_1 << " is " << d << std::endl;
    }
  }
  std::cout << "queHadronCadidates size: " << queHadronCadidates.size() << std::endl;

  timer.Stop();
  if(par::isDebug) std::cout<<"Time for calculating meson distance: "<<timer.RealTime()<<" s"<<std::endl;
  timer.Start();

  // 计算三个parton的费马点和费马距离 -> 重子
  for (int i = 0; i < partons.size(); i++) {
    int pdg_quark_0 = partons[i].PDG();
    int nSerial_0 = partons[i].GetSerial();
    float x0, y0, z0;
    float px0, py0, pz0;
    partons[i].GetPosition(x0, y0, z0);
    partons[i].GetMomentum(px0, py0, pz0);
    
    // 从i+1开始，避免重复计算
    for (int j = i + 1; j < partons.size(); j++) {
      int pdg_quark_1 = partons[j].PDG();
      //如果前两个夸克异号，则已经不可能是重子，可以直接跳过
      if (pdg_quark_0 * pdg_quark_1 < 0) {
        continue;
      }
      int nSerial_1 = partons[j].GetSerial();
      float x1, y1, z1;
      float px1, py1, pz1;
      partons[j].GetPosition(x1, y1, z1);
      partons[j].GetMomentum(px1, py1, pz1);

      // 从j+1开始，避免重复计算
      for (int k = j + 1; k < partons.size(); k++) {
        int pdg_quark_2 = partons[k].PDG();
        //如果第三个夸克和前两个不同号，则已经不可能是重子，可以直接跳过
        if (pdg_quark_0 * pdg_quark_1 * pdg_quark_2 < 0) {
          continue;
        }
        int nSerial_2 = partons[k].GetSerial();
        float x2, y2, z2;
        float px2, py2, pz2;
        partons[k].GetPosition(x2, y2, z2);
        partons[k].GetMomentum(px2, py2, pz2);
        int pdg_baryon = 0;
        pdg_baryon = LookupSpecies(partons[i].PDG(), partons[j].PDG(), partons[k].PDG());
        if (pdg_baryon == 0) continue;

        float x, y, z;
        fermatPoint(x0, y0, z0, x1, y1, z1, x2, y2, z2, x, y, z);
        float d_baryon = distance3D(x0, y0, z0, x, y, z);
        //重子乘以系数
        float d = d_baryon * r_bm;

        nSerialHadron++;
        Hadron hadronCadi(nSerialHadron, pdg_baryon, x, y, z, (px0 + px1 + px2) / 3, (py0 + py1 + py2) / 3, (pz0 + pz1 + pz2) / 3, d, nSerial_0, nSerial_1, nSerial_2);
        queHadronCadidates.push(hadronCadi);
        // std::cout << "Fermat point of parton " << nSerial_0 << ", parton " << nSerial_1 << " and parton " << nSerial_2 << " is (" << x << ", " << y << ", " << z << ") with distance " << d << std::endl;
      }
    }
  }
  timer.Stop();
  if(par::isDebug) std::cout<<"Time for calculating baryon distance: "<<timer.RealTime()<<" s"<<std::endl;
  timer.Start();

  // 存储已被选中的 parton 的序列号
  std::set<int> selectedPartons;
  bool isTimeToBreak = false;
  // FIXME: 这里的条件可能有问题
  // 我们来想一下，理论上说，上面的proiority queue里面，应该已经写入了所有的hadron了，
  // 这也意味着，所有可能的parton组合都已经被写入了
  // 这时候会存在某些parton没有被选中，的可能吗？

  while (!queHadronCadidates.empty() && selectedPartons.size() < partons.size()) {
    // 利用top()函数获取队列中的第一个元素
    Hadron hc = queHadronCadidates.top();
    // 利用pop()函数将其从队列中删除
    queHadronCadidates.pop();

    // 确保所有的parton都只被选中一次
    int nSerial0, nSerial1, nSerial2;
    // 但是如果nSerial2 = 0 说明是介子，不需要检查nSerial2
    hc.GetPartonSerials(nSerial0, nSerial1, nSerial2);
    if (nSerial2 == 0) {
      // nSerial2 == 0 介子，只需要检查nSerial0和nSerial1
      if (selectedPartons.find(nSerial0) != selectedPartons.end() || selectedPartons.find(nSerial1) != selectedPartons.end()) {
        continue;
      }
    } else {
      // nSerial2 != 0 重子，需要检查nSerial0, nSerial1和nSerial2
      if (selectedPartons.find(nSerial0) != selectedPartons.end() || selectedPartons.find(nSerial1) != selectedPartons.end() || selectedPartons.find(nSerial2) != selectedPartons.end()) {
        continue;
      } 
    }
    // 如果没有被选中，那么将这个hadron加入到hadrons中
    selectedPartons.insert(nSerial0);
    selectedPartons.insert(nSerial1);
    if (nSerial2 != 0) selectedPartons.insert(nSerial2);

    hadrons.emplace_back(hc);
  }

  timer.Stop();
  if(par::isDebug) std::cout<<"Time for selecting hadrons: "<<timer.RealTime()<<" s"<<std::endl;

  std::cout << "Hadron size: " << hadrons.size() << std::endl;
}


int Coalescence::LookupSpecies(int pdg_quark_0, int pdg_quark_1, int pdg_quark_2) {
  if (pdg_quark_2 == 0) {
    return LookupMesonSpecies(pdg_quark_0, pdg_quark_1);
  } else {
    return LookupBaryonSpecies(pdg_quark_0, pdg_quark_1, pdg_quark_2);
  }
}

void Coalescence::InitMesonLookupTable() {
  // π介子
  // π^+ (u, anti-d) -> 211
  mesonLookupTable[std::make_tuple(-1, 2)] = 211;
  // π^- (d, anti-u) -> -211
  mesonLookupTable[std::make_tuple(-2, 1)] = -211;
  // π^0 (u, anti-u) or (d, anti-d) -> 111
  mesonLookupTable[std::make_tuple(-2, 2)] = 111;
  mesonLookupTable[std::make_tuple(-1, 1)] = 111;

  // K介子
  // K^+ (u, anti-s) -> 321
  mesonLookupTable[std::make_tuple(-3, 2)] = 321;
  // K^- (s, anti-u) -> -321
  mesonLookupTable[std::make_tuple(-2, 3)] = -321;
  // K^0 (d, anti-s) -> 311
  mesonLookupTable[std::make_tuple(-3, 1)] = 311;
  // anti-K^0 (s, anti-d) -> -311
  mesonLookupTable[std::make_tuple(-1, 3)] = -311;

  // η介子 (暂时注释掉)
  // η (u, anti-u) or (d, anti-d) or (s, anti-s) -> 221
  // mesonLookupTable[std::make_tuple(-2, 2)] = 221;
  // mesonLookupTable[std::make_tuple(-1, 1)] = 221;
  // mesonLookupTable[std::make_tuple(-3, 3)] = 221;

  // ρ介子 (暂时注释掉)
  // ρ^+ (u, anti-d) -> 213
  // mesonLookupTable[std::make_tuple(-1, 2)] = 213;
  // ρ^- (d, anti-u) -> -213
  // mesonLookupTable[std::make_tuple(-2, 1)] = -213;
  // ρ^0 (u, anti-u) or (d, anti-d) -> 113
  // mesonLookupTable[std::make_tuple(-2, 2)] = 113;
  // mesonLookupTable[std::make_tuple(-1, 1)] = 113;

  // φ介子
  // φ (s, anti-s) -> 333
  mesonLookupTable[std::make_tuple(-3, 3)] = 333;

  // J/ψ介子
  // J/ψ (c, anti-c) -> 443
  mesonLookupTable[std::make_tuple(-4, 4)] = 443;

  // D介子
  // D^0 (c, anti-u) -> 421
  mesonLookupTable[std::make_tuple(-2, 4)] = 421;
  // anti-D^0 (u, anti-c) -> -421
  mesonLookupTable[std::make_tuple(-4, 2)] = -421;
  // D^+ (c, anti-d) -> 411
  mesonLookupTable[std::make_tuple(-1, 4)] = 411;
  // D^- (d, anti-c) -> -411
  mesonLookupTable[std::make_tuple(-4, 1)] = -411;

  // B介子
  // B^0 (d, anti-b) -> 511
  mesonLookupTable[std::make_tuple(-5, 1)] = 511;
  // anti-B^0 (b, anti-d) -> -511
  mesonLookupTable[std::make_tuple(-1, 5)] = -511;
  // B^+ (u, anti-b) -> 521
  mesonLookupTable[std::make_tuple(-5, 2)] = 521;
  // B^- (b, anti-u) -> -521
  mesonLookupTable[std::make_tuple(-2, 5)] = -521;
}

void Coalescence::InitBaryonLookupTable() {
  // 质子 (uud) -> 2212
  baryonLookupTable[std::make_tuple(1, 2, 2)] = 2212;
  // 反质子 (anti-u, anti-u, anti-d) -> -2212
  baryonLookupTable[std::make_tuple(-2, -2, -1)] = -2212;

  // 中子 (udd) -> 2112
  baryonLookupTable[std::make_tuple(1, 1, 2)] = 2112;
  // 反中子 (anti-d, anti-d, anti-u) -> -2112
  baryonLookupTable[std::make_tuple(-2, -1, -1)] = -2112;

  // Λ重子 (uds) -> 3122
  baryonLookupTable[std::make_tuple(1, 2, 3)] = 3122;
  // 反Λ重子 (anti-u, anti-d, anti-s) -> -3122
  baryonLookupTable[std::make_tuple(-3, -2, -1)] = -3122;

  // Σ^+重子 (uus) -> 3222
  baryonLookupTable[std::make_tuple(2, 2, 3)] = 3222;
  // 反Σ^+重子 (anti-u, anti-u, anti-s) -> -3222
  baryonLookupTable[std::make_tuple(-3, -2, -2)] = -3222;

  // // Σ^0重子 (uds) -> 3212 (暂时注释掉)
  // baryonLookupTable[std::make_tuple(1, 2, 3)] = 3212;
  // // 反Σ^0重子 (anti-u, anti-d, anti-s) -> -3212
  // baryonLookupTable[std::make_tuple(-3, -2, -1)] = -3212;

  // Σ^-重子 (dds) -> 3112
  baryonLookupTable[std::make_tuple(1, 1, 3)] = 3112;
  // 反Σ^-重子 (anti-d, anti-d, anti-s) -> -3112
  baryonLookupTable[std::make_tuple(-3, -1, -1)] = -3112;

  // Ξ^0重子 (uss) -> 3322
  baryonLookupTable[std::make_tuple(2, 3, 3)] = 3322;
  // 反Ξ^0重子 (anti-u, anti-s, anti-s) -> -3322
  baryonLookupTable[std::make_tuple(-3, -3, -2)] = -3322;

  // Ξ^-重子 (dss) -> 3312
  baryonLookupTable[std::make_tuple(1, 3, 3)] = 3312;
  // 反Ξ^-重子 (anti-d, anti-s, anti-s) -> -3312
  baryonLookupTable[std::make_tuple(-3, -3, -1)] = -3312;

  // Ω^-重子 (sss) -> 3334
  baryonLookupTable[std::make_tuple(3, 3, 3)] = 3334;
  // 反Ω^-重子 (anti-s, anti-s, anti-s) -> -3334
  baryonLookupTable[std::make_tuple(-3, -3, -3)] = -3334;

  // 其他重子的组合和对应的PDG代码
}

//pdg_code for a quark
// u = 2, d = 1, s = 3, c = 4, b = 5, t = 6
// ubar = -2, dbar = -1, sbar = -3, cbar = -4, bbar = -5, tbar = -6

int Coalescence::LookupMesonSpecies(int pdg_quark_0, int pdg_quark_1) {
  // std::cout << "r_bm = " << r_bm << std::endl;
  // std::cout << "Coalescence to hadron" << std::endl;
  if (pdg_quark_0 == 0 || pdg_quark_1 == 0) {
    std::cout << "There is a quark without species, return 0" << std::endl;
    return 0;
  }
  // There are two quarks with the same sign, return 0
  if (pdg_quark_0 * pdg_quark_1 > 0) return 0;

  //如果第三个夸克不存在，那么只能是介子，直接查找介子信息
  int quarks[] = {pdg_quark_0, pdg_quark_1};
  //从小到大排序
  std::sort(quarks, quarks + 2);
  MesonCombination key = std::make_tuple(quarks[0], quarks[1]);
  if (mesonLookupTable.find(key) != mesonLookupTable.end()) {
    return mesonLookupTable[key];
  } else {
    return 0;
  }
  return 0;
}

int Coalescence::LookupBaryonSpecies(int pdg_quark_0, int pdg_quark_1, int pdg_quark_2) {
  // std::cout << "r_bm = " << r_bm << std::endl;
  // std::cout << "Coalescence to hadron" << std::endl;
  if (pdg_quark_0 == 0 || pdg_quark_1 == 0 || pdg_quark_2 == 0) {
    std::cout << "There is a quark without species, return 0" << std::endl;
    return 0;
  }
  // 三个夸克一定全是同号，否则返回0
  if (pdg_quark_0 * pdg_quark_1 < 0) return 0;
  if (pdg_quark_0 * pdg_quark_2 < 0) return 0;
  if (pdg_quark_1 * pdg_quark_2 < 0) return 0;

  //如果第三个夸克存在，那么一定是重子，直接查找重子信息
  int quarks[] = {pdg_quark_0, pdg_quark_1, pdg_quark_2};
  std::sort(quarks, quarks + 3);
  BaryonCombination key = std::make_tuple(quarks[0], quarks[1], quarks[2]);
  if (baryonLookupTable.find(key) != baryonLookupTable.end()) {
    return baryonLookupTable[key];
  } else {
    return 0;
  }
  return 0;
}

void Coalescence::ProcessFromParton(std::vector<Parton> const &partons0, std::vector<Hadron> &hadrons) {
  // 这种算法，每次循环必出现一个hadron
  // 随机化parton vector

  auto& partons = const_cast<std::vector<Parton>&>(partons0);
  std::shuffle(partons.begin(), partons.end(), par::gen);

  bool isThisPartonUsed[partons.size()];
  int nPartons = partons.size();
  float x[nPartons], y[nPartons], z[nPartons];
  float px[nPartons], py[nPartons], pz[nPartons];
  int pdg_quark[nPartons];

  for (int i = 0; i < nPartons; i++) {
    partons[i].GetPosition(x[i], y[i], z[i]);
    partons[i].GetMomentum(px[i], py[i], pz[i]);
    pdg_quark[i] = partons[i].PDG();
    isThisPartonUsed[i] = false;
    // 如果是重夸克，直接标记为已使用
    if (par::isRemoveHFQuarks && (abs(pdg_quark[i]) == 4 || abs(pdg_quark[i]) == 5 || abs(pdg_quark[i]) == 6)) {
      isThisPartonUsed[i] = true;
    }
  }

  int nHadronSerial = 0;
  for (int iParton = 0; iParton < nPartons; iParton++) {
    if (isThisPartonUsed[iParton]) continue;
    int pdg_quark_0 = pdg_quark[iParton];
    if (pdg_quark_0 == 0) continue;

    float dis = -1;
    float x0 = x[iParton];
    float y0 = y[iParton];
    float z0 = z[iParton];

    float d_meson_min = 1e5;
    float d_diquark_min = 1e5;
    float d_baryon_min = 1e5;
    
    int meson_quark_label[2] = {iParton, -1};
    int diquark_quark_label[2] = {iParton, -1};
    int baryon_quark_label[3] = {iParton, -1, -1};

    std::vector<bool> isPartonUsedAsDiquarkThisLoop(nPartons, false);
    int nSerialLastDiquark = -1;
    for (int jParton = iParton + 1; jParton < nPartons; jParton++) {
      if (isThisPartonUsed[jParton]) continue;
      int pdg_quark_1 = pdg_quark[jParton];
      if (pdg_quark_1 == 0) continue;
      float x1 = x[jParton];
      float y1 = y[jParton];
      float z1 = z[jParton];

      float d = distance3D(x0, y0, z0, x1, y1, z1);
      float d_meson = d; // 这里是为了下面的pi0的特殊处理
      if (pdg_quark_0 * pdg_quark_1 < 0) {
        int pdg_meson = LookupMesonSpecies(pdg_quark_0, pdg_quark_1);
        // TODO：如果说是 u-ubar 或者 d-dbar，那么查表都会得到一个pi0
        // 我们设置50%的概率为p0，50%的概率这次不生成，即让距离变得无限大
        if (pdg_meson == 111 && static_cast<bool>(par::zero_or_one(par::gen))) {
          d_meson = 1e6;
        }
        if (pdg_meson != 0 && d_meson < d_meson_min) {
          d_meson_min = d_meson;
          meson_quark_label[1] = jParton;
        }
      } else if (pdg_quark_0 * pdg_quark_1 > 0) {
        // 只可能是di-quark
        if(d < d_diquark_min) {
          d_diquark_min = d;
          diquark_quark_label[1] = jParton;
          if (nSerialLastDiquark != -1) isPartonUsedAsDiquarkThisLoop[nSerialLastDiquark] = false;
          isPartonUsedAsDiquarkThisLoop[jParton] = true;
          nSerialLastDiquark = jParton;
        }
      } else continue;
    }
    if (diquark_quark_label[1] == -1 || meson_quark_label[1] == -1) continue;

    // std::cout<<"now we get the diquark candidate with quark label: "<<diquark_quark_label[0]<<", "<<diquark_quark_label[1]<<std::endl;
    // std::cout<<"now we get the meson candidate with quark label: "<<meson_quark_label[0]<<", "<<meson_quark_label[1]<<std::endl;

    // 对于di-quark，位置是
    float x_di = (x0 + x[diquark_quark_label[1]]) / 2;
    float y_di = (y0 + y[diquark_quark_label[1]]) / 2;
    float z_di = (z0 + z[diquark_quark_label[1]]) / 2;

    //这时候di-quark已经找到了
    for (int kParton = iParton + 1; kParton < nPartons; kParton++) {
      if (isThisPartonUsed[kParton]) continue;
      if (isPartonUsedAsDiquarkThisLoop[kParton]) continue;
      int pdg_quark_2 = pdg_quark[kParton];
      if (pdg_quark_0 * pdg_quark_2 < 0) continue;
      float x2 = x[kParton];
      float y2 = y[kParton];
      float z2 = z[kParton];

      float d = distance3D(x_di, y_di, z_di, x2, y2, z2);
      int pdg_baryon = LookupBaryonSpecies(pdg_quark_0, pdg_quark[diquark_quark_label[1]], pdg_quark_2);

      if (pdg_baryon !=0 && d < d_baryon_min) {
        d_baryon_min = d;
        baryon_quark_label[1] = diquark_quark_label[1];
        baryon_quark_label[2] = kParton;
      }
    }

    // 在这里，我们已经找到了一个hadron，我们需要将这个hadron加入到hadrons中，根据b_meson 和 r_bm * b_baryon的 大小关系

    int pdg_hadron = 0;
    float x_hadron = 0, y_hadron = 0, z_hadron = 0;
    float px_hadron = 0, py_hadron = 0, pz_hadron = 0;

    if (d_meson_min < r_bm * d_baryon_min) {
      dis = d_meson_min;
      pdg_hadron = LookupMesonSpecies(pdg_quark_0, pdg_quark[meson_quark_label[1]]);
      x_hadron = (x0 + x[meson_quark_label[1]]) / 2;
      y_hadron = (y0 + y[meson_quark_label[1]]) / 2;
      z_hadron = (z0 + z[meson_quark_label[1]]) / 2;
      px_hadron = (px[meson_quark_label[0]] + px[meson_quark_label[1]]) / 2;
      py_hadron = (py[meson_quark_label[0]] + py[meson_quark_label[1]]) / 2;
      pz_hadron = (pz[meson_quark_label[0]] + pz[meson_quark_label[1]]) / 2;
      //hadrons.emplace_back(iSerial++, pdg_hadron, x_hadron, y_hadron, z_hadron, px_hadron, py_hadron, pz_hadron, dis, meson_quark_label[0], meson_quark_label[1], -1);
      Hadron hc(nHadronSerial++, pdg_hadron, x_hadron, y_hadron, z_hadron, px_hadron, py_hadron, pz_hadron, dis, partons[meson_quark_label[0]].GetSerial(), partons[meson_quark_label[1]].GetSerial(), -1);
      if(par::isDebug) {
        hc.SetParton0Position(x0, y0, z0);
        hc.SetParton1Position(x[meson_quark_label[1]], y[meson_quark_label[1]], z[meson_quark_label[1]]);
      } 
      hadrons.emplace_back(hc);
      isThisPartonUsed[meson_quark_label[0]] = true;
      isThisPartonUsed[meson_quark_label[1]] = true;
    } else {
      dis = r_bm * d_baryon_min;
      pdg_hadron = LookupBaryonSpecies(pdg_quark_0, pdg_quark[diquark_quark_label[1]], pdg_quark[baryon_quark_label[2]]);
      x_hadron = (x[baryon_quark_label[0]] + x[baryon_quark_label[1]] + x[baryon_quark_label[2]]) / 3;
      y_hadron = (y[baryon_quark_label[0]] + y[baryon_quark_label[1]] + y[baryon_quark_label[2]]) / 3;
      z_hadron = (x[baryon_quark_label[0]] + z[baryon_quark_label[1]] + z[baryon_quark_label[2]]) / 3;
      px_hadron = (px[baryon_quark_label[0]] + px[baryon_quark_label[1]] + px[baryon_quark_label[2]]) / 3;
      py_hadron = (py[baryon_quark_label[0]] + py[baryon_quark_label[1]] + py[baryon_quark_label[2]]) / 3;
      pz_hadron = (pz[baryon_quark_label[0]] + pz[baryon_quark_label[1]] + pz[baryon_quark_label[2]]) / 3;
      //hadrons.emplace_back(iSerial++, pdg_hadron, x_hadron, y_hadron, z_hadron, px_hadron, py_hadron, pz_hadron, dis, baryon_quark_label[0], baryon_quark_label[1], baryon_quark_label[2]);
      Hadron hc(nHadronSerial++, pdg_hadron, x_hadron, y_hadron, z_hadron, px_hadron, py_hadron, pz_hadron, dis, partons[baryon_quark_label[0]].GetSerial(), partons[baryon_quark_label[1]].GetSerial(), partons[baryon_quark_label[2]].GetSerial());
      if(par::isDebug) {
        hc.SetParton0Position(x0, y0, z0);
        hc.SetParton1Position(x[diquark_quark_label[1]], y[diquark_quark_label[1]], z[diquark_quark_label[1]]);
        hc.SetParton2Position(x[baryon_quark_label[2]], y[baryon_quark_label[2]], z[baryon_quark_label[2]]);
      }
      hadrons.emplace_back(hc);
      isThisPartonUsed[baryon_quark_label[0]] = true;
      isThisPartonUsed[baryon_quark_label[1]] = true;
      isThisPartonUsed[baryon_quark_label[2]] = true;
    }
  }
}

