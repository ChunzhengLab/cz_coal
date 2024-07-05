#include <queue>
#include <set>
#include <random>
#include <iostream>
#include <algorithm>
#include "Coalescence.h"
#include "DistanceFun.h"
#include "TStopwatch.h"
#include "Par.h"

bool Coalescence::initialized = false;
std::map<BaryonCombination, int> Coalescence::baryonLookupTable;
std::map<MesonCombination, int> Coalescence::mesonLookupTable;

void Coalescence::Process(std::vector<Parton> const &partons, std::vector<Hadron> &hadrons) {
  nPartonsThisEvent = partons.size();
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
  std::cout<<"ProcessClassic: "<<std::endl;
  std::cout<<"It will take a really long time, suggest parton number < 10"<<std::endl;
  TStopwatch timer;
  timer.Start();

  if(par::isDebug) std::cout<<"How many partons in this event for coalescence: "<<partons.size()<<std::endl;
  
  // 将Hadron对象放入优先队列，使用优先队列来排序距离
  std::priority_queue<Hadron, std::vector<Hadron>, std::greater<Hadron>> queHadronCadidates;

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

void Coalescence::ProcessFromParton(std::vector<Parton> const &partons0, std::vector<Hadron> &hadrons, int nLastHadronSerial) {
  // 这种算法，每次循环必出现一个hadron
  // 随机化parton vector

  auto& partons = const_cast<std::vector<Parton>&>(partons0);
  std::shuffle(partons.begin(), partons.end(), par::gen);

  bool isThisPartonUsed[partons.size()];
  int nPartons = partons.size();

  for (int i = 0; i < nPartons; i++) {
    isThisPartonUsed[i] = false;
    // 如果是重夸克c = 4, b = 5, t = 6，直接标记为已使用
    if (par::isRemoveHFQuarks && (abs(partons[i].PDG()) > 3)) {
      isThisPartonUsed[i] = true;
    }
  }

  int nHadronSerial = nLastHadronSerial;
  for (int iParton = 0; iParton < nPartons; iParton++) {
    if (isThisPartonUsed[iParton]) continue;
    // 第0个parton
    float x0 = 0, y0 = 0, z0 = 0;
    float px0 = 0, py0 = 0, pz0 = 0;
    partons[iParton].GetPosition(x0, y0, z0);
    partons[iParton].GetMomentum(px0, py0, pz0);
    float t0 = partons[iParton].Time();
    float pdg0 = partons[iParton].PDG();
    if (pdg0 == 0) continue;

    // meson
    float x_me = 0, y_me = 0, z_me = 0;
    float px_me = 0, py_me = 0, pz_me = 0;
    float t_me = 0;
    int pdg_me = 0;

    // baryon
    float x_ba = 0, y_ba = 0, z_ba = 0;
    float px_ba = 0, py_ba = 0, pz_ba = 0;
    float t_ba = 0;
    int pdg_ba = 0;

    float d_meson_min = 1e6;
    float d_diquark_min = 1e6;
    float d_baryon_min = 1e6;
    
    int meson_quark_label[2] = {iParton, -1};
    int diquark_quark_label[2] = {iParton, -1};
    int baryon_quark_label[3] = {iParton, -1, -1};
    std::vector<bool> isPartonUsedAsDiquarkThisLoop(nPartons, false);

    int nSerialLastDiquark = -1; // 上一个di-quark的序列号

    for (int jParton = iParton + 1; jParton < nPartons; jParton++) {
      float x1 = 0, y1 = 0, z1 = 0;
      float px1 = 0, py1 = 0, pz1 = 0;
      float t1 = 0;
      int pdg1 = 0;
      if (isThisPartonUsed[jParton]) continue;
      partons[jParton].GetPosition(x1, y1, z1);
      partons[jParton].GetMomentum(px1, py1, pz1);
      t1 = partons[jParton].Time();
      pdg1 = partons[jParton].PDG();
      if (pdg1 == 0) continue;

      // std::cout<<"x0 "<<x0<<", y0 "<<y0<<", z0 "<<z0<<", x1 "<<x1<<", y1 "<<y1<<", z1 "<<z1<<std::endl;

      float x0_tmp = x0, y0_tmp = y0, z0_tmp = z0;
      float x1_tmp = x1, y1_tmp = y1, z1_tmp = z1;

      float d = distance3DMoveOn(x0_tmp, y0_tmp, z0_tmp, x1_tmp, y1_tmp, z1_tmp, px0, py0, pz0, px1, py1, pz1, t0, t1);
      float d_meson = d; // 这里是为了下面的pi0的特殊处理
      float d_diquark = d;

      if (pdg0 * pdg1 < 0) {
        // 如果说是 u-ubar 或者 d-dbar，那么查表都会得到一个pi0
        if (pdg0 == -2 && pdg1 == 2 || pdg0 == -1 && pdg1 == 1 || pdg0 == 2 && pdg1 == -2 || pdg0 == 1 && pdg1 == -1) {
          // 我们设置50%的概率为pi0，50%的概率这次不生成，即让距离变得无限大
          if(static_cast<bool>(par::zero_or_one(par::gen))) d_meson = std::numeric_limits<float>::max();
        }
        if (d_meson < d_meson_min) {
          d_meson_min = d_meson;
          meson_quark_label[1] = jParton;
          pdg_me = LookupMesonSpecies(pdg0, pdg1);
          x_me = (x0_tmp + x1_tmp) / 2, y_me = (y0_tmp + y1_tmp) / 2, z_me = (z0_tmp + z1_tmp) / 2;
          px_me = (px0 + px1), py_me = (py0 + py1), pz_me = (pz0 + pz1);
          t_me = t0 > t1 ? t0 : t1; // 取最大的时间
        }
      } else if (pdg0 * pdg1 > 0) {
        // 只可能是di-quark
        if(d_diquark < d_diquark_min) {
          d_diquark_min = d_diquark;
          diquark_quark_label[1] = jParton;
          if (nSerialLastDiquark != -1) isPartonUsedAsDiquarkThisLoop[nSerialLastDiquark] = false;
          isPartonUsedAsDiquarkThisLoop[jParton] = true;
          // 好像把label记下来就行了，不需要更新x_di, y_di, z_di, t_di
          nSerialLastDiquark = jParton;
        }
      } else continue;
    }
    if (diquark_quark_label[1] == -1 || meson_quark_label[1] == -1) continue;

    //这时候di-quark已经找到了,我们需要找到一个quark来组成重子
    for (int kParton = iParton + 1; kParton < nPartons; kParton++) {
      if (isThisPartonUsed[kParton]) continue;
      if (isPartonUsedAsDiquarkThisLoop[kParton]) continue;
      // 第2个parton
      int pdg2 = 0;
      float x2 = 0, y2 = 0, z2 = 0;
      float px2 = 0, py2 = 0, pz2 = 0;
      float t2 = 0;
      pdg2 = partons[kParton].PDG();
      if (pdg0 * pdg2 < 0) continue; // 重子的三个夸克必须同号,pdg0和pdg1已经同号
      partons[kParton].GetPosition(x2, y2, z2);
      partons[kParton].GetMomentum(px2, py2, pz2);
      t2 = partons[kParton].Time();

      // 第0个parton
      float x0_tmp = x0, y0_tmp = y0, z0_tmp = z0;

      // 第1个parton
      int pdg1 = 0;
      float x1_tmp = 0, y1_tmp = 0, z1_tmp = 0;
      float px1 = 0, py1 = 0, pz1 = 0;
      float t1 = 0;
      // 根据记录的di-quark的label，找到第1个quark的位置
      pdg1 = partons[diquark_quark_label[1]].PDG();
      partons[diquark_quark_label[1]].GetPosition(x1_tmp, y1_tmp, z1_tmp);
      partons[diquark_quark_label[1]].GetMomentum(px1, py1, pz1);
      t1 = partons[diquark_quark_label[1]].Time();

      float d_baryon = perimeterMoveOn(x0_tmp, y0_tmp, z0_tmp, x1_tmp, y1_tmp, z1_tmp, x2, y2, z2, px0, py0, pz0, px1, py1, pz1, px2, py2, pz2, t0, t1, t2);

      d_baryon = d_baryon / 3.; // 周长的距离除以3，得到平均距离
      // std::cout<<"pdg0: "<<pdg0<<", pdg1: "<<pdg1<<", pdg2: "<<pdg2<<", d_baryon: "<<d_baryon<<", d_TStruct: "<<d_TStruct<<std::endl;

      if (d_baryon < d_baryon_min) {
        d_baryon_min = d_baryon;
        baryon_quark_label[1] = diquark_quark_label[1];
        baryon_quark_label[2] = kParton;
        pdg_ba = LookupBaryonSpecies(pdg0, pdg1, pdg2);
        x_ba = (x0_tmp + x1_tmp + x2) / 3, y_ba = (y0_tmp + y1_tmp + y2) / 3, z_ba = (z0_tmp + z1_tmp + z2) / 3;
        px_ba = px0 + px1 + px2, py_ba = py0 + py1 + py2, pz_ba = pz0 + pz1 + pz2;
        t_ba = (t0 > t1) ? ((t0 > t2) ? t0 : t2) : ((t1 > t2) ? t1 : t2); // 取最大的时间
      }
    }

    // 在这里，我们已经找到了一个hadron，我们需要将这个hadron加入到hadrons中，根据b_meson 和 r_bm * b_baryon的 大小关系

    if (d_meson_min < r_bm * d_baryon_min) {
      hadrons.emplace_back(nHadronSerial++, pdg_me, x_me, y_me, z_me, px_me, py_me, pz_me, t_me, d_meson_min, partons[meson_quark_label[0]].GetSerial(), partons[meson_quark_label[1]].GetSerial(), -9999);
      if(par::isDebug) {
        hadrons.back().SetParton0Position(x0, y0, z0);
        hadrons.back().SetParton1Position(partons[meson_quark_label[1]].X(), partons[meson_quark_label[1]].Y(), partons[meson_quark_label[1]].Z());
        hadrons.back().SetParton2Position(-9999, -9999, -9999);
      } 
      isThisPartonUsed[meson_quark_label[0]] = true;
      isThisPartonUsed[meson_quark_label[1]] = true;
    } else {
      hadrons.emplace_back(nHadronSerial++, pdg_ba, x_ba, y_ba, z_ba, px_ba, py_ba, pz_ba, t_ba, d_baryon_min, partons[baryon_quark_label[0]].GetSerial(), partons[baryon_quark_label[1]].GetSerial(), partons[baryon_quark_label[2]].GetSerial());
      if(par::isDebug) {
        hadrons.back().SetParton0Position(x0, y0, z0);
        hadrons.back().SetParton1Position(partons[baryon_quark_label[1]].X(), partons[baryon_quark_label[1]].Y(), partons[baryon_quark_label[1]].Z());
        hadrons.back().SetParton2Position(partons[baryon_quark_label[2]].X(), partons[baryon_quark_label[2]].Y(), partons[baryon_quark_label[2]].Z());
      }
      isThisPartonUsed[baryon_quark_label[0]] = true;
      isThisPartonUsed[baryon_quark_label[1]] = true;
      isThisPartonUsed[baryon_quark_label[2]] = true;
    }
  }

  // 如果还有没有被使用的parton，这些parton的数量理论上应该很少
  // 可以读取isThisPartonUsed数组，找到没有被使用的parton，然后将这些parton打包成一个vector<Parton>
  // 进行递归
  std::vector<Parton> partonsUnused;
  for (int i = 0; i < nPartons; i++) {
    if (!isThisPartonUsed[i]) {
      partonsUnused.emplace_back(partons[i]);
    }
  }
  std::cout<<"Partons left with PDG code: ";
  for (auto& p : partonsUnused) {
    std::cout<<p.PDG()<<" ";
  }
  std::cout<<std::endl;

  // bool isNeedRecursion = true;
  // if (nRecursionThisEvent > 5) {
  //   // 如果递归超过5次，那么直接结束递归
  //   isNeedRecursion = false;
  //   std::cout<<"Recursion times exceed 5, recursion ends."<<std::endl;
  // } else {
  //   if (partonsUnused.size() == 0) {
  //     // 如果没有没有被使用的parton，那么递归结束
  //     isNeedRecursion = false;
  //   } else if (partonsUnused.size() == 1) {
  //     // 如果只有一个没有被使用的parton，那么递归结束
  //     if(par::isDebug) {
  //       std::cout<<"One parton with PDG code: "<<partonsUnused[0].PDG()<<" is left, recursion ends."<<std::endl;
  //     }
  //     isNeedRecursion = false;
  //   } else if (partonsUnused.size() == 2) {
  //     // 递归必然结束
  //     isNeedRecursion = false;
  //     int pdg_lookup = LookupMesonSpecies(partonsUnused[0].PDG(), partonsUnused[1].PDG());
  //     // 如果刚好能够形成一个介子，打包成一个介子
  //     if(pdg_lookup != 0) {
  //       if(par::isDebug) std::cout<<"Two partons (No."<<partonsUnused[0].GetSerial()<<", No."<<partonsUnused[1].GetSerial()<<") left with PDG code: "<<partonsUnused[0].PDG()<<", "<<partonsUnused[1].PDG()<<std::endl;
  //       // 如果有两个没有被使用的parton，而且这两个parton的PDG code同号，那么递归结束
  //       float x, y, z, px, py, pz, t;
  //       float x0, y0, z0, px0, py0, pz0, t0 = partonsUnused[0].Time();
  //       float x1, y1, z1, px1, py1, pz1, t1 = partonsUnused[1].Time();
  //       partonsUnused[0].GetPosition(x0, y0, z0);
  //       partonsUnused[0].GetMomentum(px0, py0, pz0);
  //       partonsUnused[1].GetPosition(x1, y1, z1);
  //       partonsUnused[1].GetMomentum(px1, py1, pz1);

  //       if (t0 > t1) moveOn(x1, y1, z1, px1, py1, pz1, t0 - t1);
  //       else if (t0 < t1) moveOn(x0, y0, z0, px0, py0, pz0, t1 - t0);

  //       px = px0 + px1, py = py0 + py1, pz = pz0 + pz1;
  //       float d = distance3D(x0, y0, z0, x1, y1, z1);

  //       hadrons.emplace_back(nHadronSerial++, pdg_lookup, x, y, z, px, py, pz, d, partonsUnused[0].GetSerial(), partonsUnused[1].GetSerial(), -1);
  //       if(par::isDebug) {
  //         hadrons.back().SetParton0Position(x0, y0, z0);
  //         hadrons.back().SetParton1Position(x1, y1, z1);
  //       }
  //       if(par::isDebug) std::cout<<"------These two partons can form a meson: "<<pdg_lookup<<std::endl;
  //       partonsUnused.clear();
  //       if(par::isDebug) std::cout<<"------Used partons vector cleared."<<std::endl;
  //     } else {
  //       // 无法形成介子，递归结束
  //       std::cout<<"Two partons left with pdg code: "<<partonsUnused[0].PDG()<<", "<<partonsUnused[1].PDG()<<std::endl;
  //       std::cout<<"------These two partons cannot form a meson."<<std::endl;
  //     }
  //   } else if (partonsUnused.size() == 3) {
  //     int pdg_lookup = LookupBaryonSpecies(partonsUnused[0].PDG(), partonsUnused[1].PDG(), partonsUnused[2].PDG());
  //     if (pdg_lookup != 0) {
  //       if(par::isDebug) std::cout<<"Three partons (No."<<partonsUnused[0].GetSerial()<<", No."<<partonsUnused[1].GetSerial()<<", No."<<partonsUnused[2].GetSerial()<<") left with PDG code: "<<partonsUnused[0].PDG()<<", "<<partonsUnused[1].PDG()<<", "<<partonsUnused[2].PDG()<<std::endl;
  //       // 可以形成重子，直接打包成一个hadron
  //       float x, y, z, px, py, pz, t;
  //       float x0, y0, z0, px0, py0, pz0, t0 = partonsUnused[0].Time();
  //       float x1, y1, z1, px1, py1, pz1, t1 = partonsUnused[1].Time();
  //       float x2, y2, z2, px2, py2, pz2, t2 = partonsUnused[2].Time();
  //       partonsUnused[0].GetPosition(x0, y0, z0);
  //       partonsUnused[0].GetMomentum(px0, py0, pz0);
  //       partonsUnused[1].GetPosition(x1, y1, z1);
  //       partonsUnused[1].GetMomentum(px1, py1, pz1);
  //       partonsUnused[2].GetPosition(x2, y2, z2);
  //       partonsUnused[2].GetMomentum(px2, py2, pz2);

  //       float d = perimeterMoveOn(x0, y0, z0, px0, py0, pz0, t0, x1, y1, z1, px1, py1, pz1, t1, x2, y2, z2, px2, py2, pz2, t2);
  //       x = (x0 + x1 + x2) / 3, y = (y0 + y1 + y2) / 3, z = (z0 + z1 + z2) / 3;
  //       px = px0 + px1 + px2, py = py0 + py1 + py2, pz = pz0 + pz1 + pz2;
  //       t = t0 > t1 ? t0 : t1;
  //       t = t > t2 ? t : t2;

  //       hadrons.emplace_back(nHadronSerial++, pdg_lookup, x, y, z, px, py, pz, d, partonsUnused[0].GetSerial(), partonsUnused[1].GetSerial(), partonsUnused[2].GetSerial());

  //       if(par::isDebug) {
  //         hadrons.back().SetParton0Position(partonsUnused[0].X(), partonsUnused[0].Y(), partonsUnused[0].Z());
  //         hadrons.back().SetParton1Position(partonsUnused[1].X(), partonsUnused[1].Y(), partonsUnused[1].Z());
  //         hadrons.back().SetParton2Position(partonsUnused[2].X(), partonsUnused[2].Y(), partonsUnused[2].Z());
  //       }
  //       if(par::isDebug) std::cout<<"------These three partons can form a baryon: "<<pdg_lookup<<std::endl;
  //       partonsUnused.clear();
  //       if(par::isDebug) std::cout<<"------Used partons vector cleared."<<std::endl;
  //       isNeedRecursion = false;
  //     }
  //   } else {
  //     // 其他情况，需要递归
  //     isNeedRecursion = true;
  //   }
  // }

  // // 递归
  // if (isNeedRecursion) {
  //   nRecursionThisEvent++;
  //   if(par::isDebug) {
  //     std::cout<<"Recursion "<<nRecursionThisEvent<<" starts."<<std::endl;
  //   }
  //   ProcessFromParton(partonsUnused, hadrons, nHadronSerial);
  // } else {
  //   if(par::isDebug) {
  //     std::cout<<"----------------------------------------"<<std::endl;
  //     std::cout<<"Recursion ends."<<std::endl;
  //     std::cout<<"Number of recursion in this event: "<<nRecursionThisEvent<<std::endl;
  //     std::cout<<"----------------------------------------"<<std::endl;
  //   }
  //   // 如果最终partonsUnused大于0，说明递归结束后还有一些parton没有被使用
  //   if (partonsUnused.size() > 0) {
  //     if (par::isDebug) {
  //       if (partonsUnused.size() > par::flavourBreakTolerance * partons.size()) {
  //         std::cerr<<"Error: "<<partonsUnused.size()<<" partons are left unused after recursion, which is more than "<<par::flavourBreakTolerance / 100. <<"% of the total partons."<<std::endl;
  //         std::cerr<<"This event will be saved as empty."<<std::endl;
  //         std::cerr<<"Clearing hadrons array."<<std::endl;
  //         std::vector<Hadron>().swap(hadrons);
  //       } else {
  //         std::cout<<"Warning: "<<partonsUnused.size()<<" partons are left unused after recursion."<<std::endl;
  //       }
  //     }
  //   }
  //   // 递归结束, Reset
  //   ResetRecursionForNextEvent();
  // }
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
  // N重子
  // 质子 (uud) -> 2212
  baryonLookupTable[std::make_tuple(1, 2, 2)] = 2212;
  // 反质子 (anti-u, anti-u, anti-d) -> -2212
  baryonLookupTable[std::make_tuple(-2, -2, -1)] = -2212;
  // 中子 (udd) -> 2112
  baryonLookupTable[std::make_tuple(1, 1, 2)] = 2112;
  // 反中子 (anti-d, anti-d, anti-u) -> -2112
  baryonLookupTable[std::make_tuple(-2, -1, -1)] = -2112;
  
  // Δ重子
  // Δ^++ (uuu) -> 2224
  baryonLookupTable[std::make_tuple(2, 2, 2)] = 2224;
  // 反Δ^++ (anti-u, anti-u, anti-u) -> -2224
  baryonLookupTable[std::make_tuple(-2, -2, -2)] = -2224;
  // // Δ^+ (uud) -> 2214 // 组分与质子相同，暂时注释掉
  // baryonLookupTable[std::make_tuple(1, 2, 2)] = 2214;
  // // Δ^0 (udd) -> 2114 // 组分与中子相同，暂时注释掉
  // baryonLookupTable[std::make_tuple(1, 1, 2)] = 2114;
  // Δ^- (ddd) -> 1114
  baryonLookupTable[std::make_tuple(1, 1, 1)] = 1114;
  // 反Δ^- (anti-d, anti-d, anti-d) -> -1114
  baryonLookupTable[std::make_tuple(-1, -1, -1)] = -1114;

  // Λ重子
  // Λ0 (uds) -> 3122
  baryonLookupTable[std::make_tuple(1, 2, 3)] = 3122;
  // 反Λ0 (anti-u, anti-d, anti-s) -> -3122
  baryonLookupTable[std::make_tuple(-3, -2, -1)] = -3122;

  // Σ 重子
  // Σ^+重子 (uus) -> 3222
  baryonLookupTable[std::make_tuple(2, 2, 3)] = 3222;
  // 反Σ^+重子 (anti-u, anti-u, anti-s) -> -3222
  baryonLookupTable[std::make_tuple(-3, -2, -2)] = -3222;
  // // Σ^0重子 (uds) -> 3212 // 组分与Λ相同，暂时注释掉
  // baryonLookupTable[std::make_tuple(1, 2, 3)] = 3212;
  // // 反Σ^0重子 (anti-u, anti-d, anti-s) -> -3212 // 组分与反Λ相同，暂时注释掉
  // baryonLookupTable[std::make_tuple(-3, -2, -1)] = -3212;
  // Σ^-重子 (dds) -> 3112
  baryonLookupTable[std::make_tuple(1, 1, 3)] = 3112;
  // 反Σ^-重子 (anti-d, anti-d, anti-s) -> -3112
  baryonLookupTable[std::make_tuple(-3, -1, -1)] = -3112;
   
  // Ξ重子
  // Ξ^0重子 (uss) -> 3322
  baryonLookupTable[std::make_tuple(2, 3, 3)] = 3322;
  // 反Ξ^0重子 (anti-u, anti-s, anti-s) -> -3322
  baryonLookupTable[std::make_tuple(-3, -3, -2)] = -3322;
  // Ξ^-重子 (dss) -> 3312
  baryonLookupTable[std::make_tuple(1, 3, 3)] = 3312;
  // 反Ξ^-重子 (anti-d, anti-s, anti-s) -> -3312
  baryonLookupTable[std::make_tuple(-3, -3, -1)] = -3312;
  
  // Ω重子
  // Ω^-重子 (sss) -> 3334
  baryonLookupTable[std::make_tuple(3, 3, 3)] = 3334;
  // 反Ω^-重子 (anti-s, anti-s, anti-s) -> -3334
  baryonLookupTable[std::make_tuple(-3, -3, -3)] = -3334;

  // 其他重子的组合和对应的PDG代码
}


