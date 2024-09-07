#include <queue>
#include <set>
#include <random>
#include <iostream>
#include <algorithm>
#include "Coalescence.h"
#include "DistanceFun.h"
#include "TStopwatch.h"
#include "Par.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TCanvas.h"

bool Coalescence::initialized = false;
std::map<BaryonCombination, int> Coalescence::baryonLookupTable;
std::map<MesonCombination, int> Coalescence::mesonLookupTable;

void Coalescence::Process(std::vector<Parton> const &partons, std::vector<Hadron> &hadrons) {
  nPartonsThisEvent = partons.size();
  switch (coalescenceAlgorithm) {
    case CoalescenceAlgorithm::kClassic:
      // ProcessClassic(partons, hadrons);
      break;
    case CoalescenceAlgorithm::kFromParton:
      ProcessFromParton(partons, hadrons);
      break;
    default:
      std::cerr << "Unknown coalescence algorithm" << std::endl;
      break;
  }

  if (par::isWriteCoalQA) {
    h_stats->Fill(0.5);
    if (isRejectByFlavourTolerance) h_stats->Fill(1.5);
    if (nRecursionThisEvent == 1) h_stats->Fill(2.5);
    if (nRecursionThisEvent == 2) h_stats->Fill(3.5);
    if (nRecursionThisEvent == 3) h_stats->Fill(4.5);
    if (nRecursionThisEvent > 3) h_stats->Fill(5.5);
    h_track_mass_reject->Fill(nTrackRejectByMass);
  }
  std::cout<<nTrackRejectByMass<<std::endl;

  ResetRecursionForNextEvent();
}

// void Coalescence::ProcessClassic(std::vector<Parton> const &partons, std::vector<Hadron> &hadrons) {
//   std::cout<<"ProcessClassic: "<<std::endl;
//   std::cout<<"It will take a really long time, suggest parton number < 10"<<std::endl;
//   TStopwatch timer;
//   timer.Start();

//   if(par::isDebug) std::cout<<"How many partons in this event for coalescence: "<<partons.size()<<std::endl;
  
//   // 将Hadron对象放入优先队列，使用优先队列来排序距离
//   std::priority_queue<Hadron, std::vector<Hadron>, std::greater<Hadron>> queHadronCadidates;

//   //计算所有parton之间的空间距离 -> 介子
//   int nSerialHadron = 0;
//   for (int i = 0; i < partons.size(); i++) {
//     int pdg_quark_0 = partons[i].PDG();
//     int nSerial_0 = partons[i].GetSerial();
//     float x0, y0, z0;
//     float px0, py0, pz0;
//     partons[i].GetPosition(x0, y0, z0);
//     partons[i].GetMomentum(px0, py0, pz0);

//     // 从i+1开始，避免重复计算
//     for (int j = i + 1; j < partons.size(); j++) {
//       int pdg_quark_1 = partons[j].PDG();
//       //如果前两个夸克同号，则已经不可能是介子，可以直接跳过
//       //采用异或运算符判断两个数是否异号
//       if (pdg_quark_0 * pdg_quark_1 > 0) {
//         continue;
//       }
//       int nSerial_1 = partons[j].GetSerial();

//       int pdg_meson = 0;
//       pdg_meson = LookupSpecies(partons[i].PDG(), partons[j].PDG(), 0);
//       if (pdg_meson == 0) continue;

//       float x1, y1, z1;
//       float px1, py1, pz1;
//       partons[j].GetPosition(x1, y1, z1);
//       partons[j].GetMomentum(px1, py1, pz1);

//       float d_meson = distance3D(x0, y0, z0, x1, y1, z1);
//       // 介子不乘系数
//       float d = d_meson;

//       nSerialHadron++;
//       Hadron hadronCadi(nSerialHadron, pdg_meson, (x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2, (px0 + px1) / 2, (py0 + py1) / 2, (pz0 + pz1) / 2, d, nSerial_0, nSerial_1, 0);
//       queHadronCadidates.push(hadronCadi);
//       // std::cout << "Distance between parton " << nSerial_0 << " and parton " << nSerial_1 << " is " << d << std::endl;
//     }
//   }

//   std::cout << "queHadronCadidates size: " << queHadronCadidates.size() << std::endl;
//   timer.Stop();
//   if(par::isDebug) std::cout<<"Time for calculating meson distance: "<<timer.RealTime()<<" s"<<std::endl;
//   timer.Start();

//   // 计算三个parton的费马点和费马距离 -> 重子
//   for (int i = 0; i < partons.size(); i++) {
//     int pdg_quark_0 = partons[i].PDG();
//     int nSerial_0 = partons[i].GetSerial();
//     float x0, y0, z0;
//     float px0, py0, pz0;
//     partons[i].GetPosition(x0, y0, z0);
//     partons[i].GetMomentum(px0, py0, pz0);
    
//     // 从i+1开始，避免重复计算
//     for (int j = i + 1; j < partons.size(); j++) {
//       int pdg_quark_1 = partons[j].PDG();
//       //如果前两个夸克异号，则已经不可能是重子，可以直接跳过
//       if (pdg_quark_0 * pdg_quark_1 < 0) {
//         continue;
//       }
//       int nSerial_1 = partons[j].GetSerial();
//       float x1, y1, z1;
//       float px1, py1, pz1;
//       partons[j].GetPosition(x1, y1, z1);
//       partons[j].GetMomentum(px1, py1, pz1);

//       // 从j+1开始，避免重复计算
//       for (int k = j + 1; k < partons.size(); k++) {
//         int pdg_quark_2 = partons[k].PDG();
//         //如果第三个夸克和前两个不同号，则已经不可能是重子，可以直接跳过
//         if (pdg_quark_0 * pdg_quark_1 * pdg_quark_2 < 0) {
//           continue;
//         }
//         int nSerial_2 = partons[k].GetSerial();
//         float x2, y2, z2;
//         float px2, py2, pz2;
//         partons[k].GetPosition(x2, y2, z2);
//         partons[k].GetMomentum(px2, py2, pz2);
//         int pdg_baryon = 0;
//         pdg_baryon = LookupSpecies(partons[i].PDG(), partons[j].PDG(), partons[k].PDG());
//         if (pdg_baryon == 0) continue;

//         float x, y, z;
//         fermatPoint(x0, y0, z0, x1, y1, z1, x2, y2, z2, x, y, z);
//         float d_baryon = distance3D(x0, y0, z0, x, y, z);
//         //重子乘以系数
//         float d = d_baryon * r_bm;

//         nSerialHadron++;
//         Hadron hadronCadi(nSerialHadron, pdg_baryon, x, y, z, (px0 + px1 + px2) / 3, (py0 + py1 + py2) / 3, (pz0 + pz1 + pz2) / 3, d, nSerial_0, nSerial_1, nSerial_2);
//         queHadronCadidates.push(hadronCadi);
//         // std::cout << "Fermat point of parton " << nSerial_0 << ", parton " << nSerial_1 << " and parton " << nSerial_2 << " is (" << x << ", " << y << ", " << z << ") with distance " << d << std::endl;
//       }
//     }
//   }
//   timer.Stop();
//   if(par::isDebug) std::cout<<"Time for calculating baryon distance: "<<timer.RealTime()<<" s"<<std::endl;
//   timer.Start();

//   // 存储已被选中的 parton 的序列号
//   std::set<int> selectedPartons;
//   bool isTimeToBreak = false;
//   // FIXME: 这里的条件可能有问题
//   // 我们来想一下，理论上说，上面的proiority queue里面，应该已经写入了所有的hadron了，
//   // 这也意味着，所有可能的parton组合都已经被写入了
//   // 这时候会存在某些parton没有被选中，的可能吗？

//   while (!queHadronCadidates.empty() && selectedPartons.size() < partons.size()) {
//     // 利用top()函数获取队列中的第一个元素
//     Hadron hc = queHadronCadidates.top();
//     // 利用pop()函数将其从队列中删除
//     queHadronCadidates.pop();

//     // 确保所有的parton都只被选中一次
//     int nSerial0, nSerial1, nSerial2;
//     // 但是如果nSerial2 = 0 说明是介子，不需要检查nSerial2
//     hc.GetPartonSerials(nSerial0, nSerial1, nSerial2);
//     if (nSerial2 == 0) {
//       // nSerial2 == 0 介子，只需要检查nSerial0和nSerial1
//       if (selectedPartons.find(nSerial0) != selectedPartons.end() || selectedPartons.find(nSerial1) != selectedPartons.end()) {
//         continue;
//       }
//     } else {
//       // nSerial2 != 0 重子，需要检查nSerial0, nSerial1和nSerial2
//       if (selectedPartons.find(nSerial0) != selectedPartons.end() || selectedPartons.find(nSerial1) != selectedPartons.end() || selectedPartons.find(nSerial2) != selectedPartons.end()) {
//         continue;
//       } 
//     }
//     // 如果没有被选中，那么将这个hadron加入到hadrons中
//     selectedPartons.insert(nSerial0);
//     selectedPartons.insert(nSerial1);
//     if (nSerial2 != 0) selectedPartons.insert(nSerial2);

//     hadrons.emplace_back(hc);
//   }

//   timer.Stop();
//   if(par::isDebug) std::cout<<"Time for selecting hadrons: "<<timer.RealTime()<<" s"<<std::endl;

//   std::cout << "Hadron size: " << hadrons.size() << std::endl;
// }


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
    std::cout << "Looking up meson species, there is a quark without species, return 0" << std::endl;
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
    std::cout << "Looking up baryon species, there is a quark don't know it's species, return 0" << std::endl;
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
  // 这种算法，每次循环一般必出现一个hadron
  auto partons = const_cast<std::vector<Parton>&>(partons0);
  // ================================================================================
  // For CreateAnimation
  // ================================================================================
  std::unique_ptr<TFile> file;
  std::unique_ptr<TCanvas> canvas;
  // 一直保留的元素
  std::unique_ptr<TGraph> g_quark_all;
  std::unique_ptr<TGraph> g_anti_quark_all;
  // 每次循环都要重新画的元素
  std::unique_ptr<TGraph> g_quark;
  std::unique_ptr<TGraph> g_anti_quark;
  std::unique_ptr<TGraph> g_quark_start;
  std::unique_ptr<TGraph> g_meson;
  std::unique_ptr<TGraph> g_baryon;
  std::unique_ptr<TGraph> g_anti_baryon;
  std::vector<std::unique_ptr<TGraph>> g_meson_shape;
  std::vector<std::unique_ptr<TGraph>> g_baryon_shape;
  int frame_number = 1;
  std::unique_ptr<TH2> dummy = std::unique_ptr<TH2>(new TH2F("", ";x;y", 1, -10, 10, 1, -10, 10));
  if (par::isCreateAnimation) {
    file = std::unique_ptr<TFile>(new TFile("particle_positon.root", "RECREATE"));
    canvas = std::unique_ptr<TCanvas>(new TCanvas("frame_0", "frame_0", 800, 800));
    g_quark_all = std::unique_ptr<TGraph>(new TGraph());
    g_quark_all->SetMarkerStyle(kFullCircle);
    g_quark_all->SetMarkerColor(kRed);
    g_quark_all->SetMarkerSize(0.5);

    g_anti_quark_all = std::unique_ptr<TGraph>(new TGraph());
    g_anti_quark_all->SetMarkerStyle(kFullCircle);
    g_anti_quark_all->SetMarkerColor(kBlue);
    g_anti_quark_all->SetMarkerSize(0.5);

    g_quark = std::unique_ptr<TGraph>(new TGraph());
    g_quark->SetMarkerStyle(kOpenSquare);
    g_quark->SetMarkerColor(kRed);
    g_quark->SetMarkerSize(0.8);

    g_anti_quark = std::unique_ptr<TGraph>(new TGraph());
    g_anti_quark->SetMarkerStyle(kOpenSquare);
    g_anti_quark->SetMarkerColor(kBlue);
    g_anti_quark->SetMarkerSize(0.8);

    g_quark_start = std::unique_ptr<TGraph>(new TGraph());
    g_quark_start->SetMarkerStyle(kStar);
    g_quark_start->SetMarkerColor(kBlack);
    g_quark_start->SetMarkerSize(1.5);

    g_meson = std::unique_ptr<TGraph>(new TGraph());
    g_meson->SetMarkerStyle(kFullCircle);
    g_meson->SetMarkerColor(kGreen);
    g_meson->SetMarkerSize(1);

    g_baryon = std::unique_ptr<TGraph>(new TGraph());
    g_baryon->SetMarkerStyle(kFullCircle);
    g_baryon->SetMarkerColor(kRed + 2);
    g_baryon->SetMarkerSize(1.2);

    g_anti_baryon = std::unique_ptr<TGraph>(new TGraph());
    g_anti_baryon->SetMarkerStyle(kFullCircle);
    g_anti_baryon->SetMarkerColor(kBlue + 2);
    g_anti_baryon->SetMarkerSize(1.2);

    // 如果需要制作动画，那么限制parton数量为60
    partons.resize(60);
  }
  // ================================================================================

  // 随机化parton vector，以保证夸克开始聚合的顺序是随机的
  if(!par::isCreateAnimation && !par::isDebug) {
    std::shuffle(partons.begin(), partons.end(), par::gen);
  }
  int nPartons = partons.size();
  
  // 为了递归之后的hadron序列号连续
  int nHadronSerial = nLastHadronSerial;


  // 支持夸克的xy坐标重采样
  if(par::isResampleQuarkXY) {
    std::unique_ptr<TH1F> h_x = std::unique_ptr<TH1F>(new TH1F("h_x", "h_x", 100, -10, 10));
    std::unique_ptr<TH1F> h_y = std::unique_ptr<TH1F>(new TH1F("h_y", "h_y", 100, -10, 10));
    for (int i = 0; i < nPartons; i++) {
      h_x->Fill(partons[i].X());
      h_y->Fill(partons[i].Y());
    }
    for (int i = 0; i < nPartons; i++) {
      partons[i].SetXY(h_x->GetRandom(), h_y->GetRandom());
    }
  }
  // 支持夸克的xy坐标随机旋转
  if(par::isRandomRotateQuarkXY) {
    for (int i = 0; i < nPartons; i++) {
      float x, y;
      partons[i].GetXY(x, y);
      float dphi = std::uniform_real_distribution<float>(0, 2 * M_PI)(par::gen);
      float x_new = x * cos(dphi) - y * sin(dphi);
      float y_new = x * sin(dphi) + y * cos(dphi);
      partons[i].SetXY(x_new, y_new);
    }
  }
  // 支持夸克的xy坐标饼状抽样
  if(par::isPiesampleXY) {
    std::unique_ptr<TH1F> h_r = std::unique_ptr<TH1F>(new TH1F("h_r", "h_r", 100, 0, 10));
    for (int i = 0; i < nPartons; i++) {
      float x, y;
      partons[i].GetXY(x, y);
      float r = sqrt(x * x + y * y);
      h_r->Fill(r);
    }
    float r_mean = h_r->GetMean();
    float r_sigma = h_r->GetRMS();
    for (int i = 0; i < nPartons; i++) {
      float r = std::normal_distribution<float>(r_mean, r_sigma)(par::gen);
      float phi = std::uniform_real_distribution<float>(0, 2 * M_PI)(par::gen);
      float x_new = sqrt(r) * cos(phi);
      float y_new = sqrt(r) * sin(phi);
      partons[i].SetXY(x_new, y_new);
    }
  }

  // 支持夸克的px, py重采样
  if(par::isResampleQuarkPxPy) {
    std::unique_ptr<TH1F> h_px = std::unique_ptr<TH1F>(new TH1F("h_px", "h_px", 100, -10, 10));
    std::unique_ptr<TH1F> h_py = std::unique_ptr<TH1F>(new TH1F("h_py", "h_py", 100, -10, 10));
    for (int i = 0; i < nPartons; i++) {
      h_px->Fill(partons[i].Px());
      h_py->Fill(partons[i].Py());
    }
    for (int i = 0; i < nPartons; i++) {
      partons[i].SetPxPy(h_px->GetRandom(), h_py->GetRandom());
    }
  }
  // 支持夸克的px, py随机旋转
  if(par::isRandomRotateQuarkPxPy) {
    for (int i = 0; i < nPartons; i++) {
      float px, py;
      partons[i].GetPxPy(px, py);
      float dphi = std::uniform_real_distribution<float>(0, 2 * M_PI)(par::gen);
      float px_new = px * cos(dphi) - py * sin(dphi);
      float py_new = px * sin(dphi) + py * cos(dphi);
      partons[i].SetPxPy(px_new, py_new);
    }
  }
  // 支持夸克的px, py饼状抽样
  if(par::isPiesamplePxPy) {
    std::unique_ptr<TH1F> h_pT = std::unique_ptr<TH1F>(new TH1F("h_pT", "h_pT", 100, 0, 10));
    for (int i = 0; i < nPartons; i++) {
      float pT = sqrt(partons[i].Px() * partons[i].Px() + partons[i].Py() * partons[i].Py());
      h_pT->Fill(pT);
    }
    float pT_mean = h_pT->GetMean();
    float pT_sigma = h_pT->GetRMS();
    for (int i = 0; i < nPartons; i++) {
      float r = std::normal_distribution<float>(pT_mean, pT_sigma)(par::gen);
      float phi = std::uniform_real_distribution<float>(0, 2 * M_PI)(par::gen);
      float px_new = r * cos(phi);
      float py_new = r * sin(phi);
      partons[i].SetPxPy(px_new, py_new);
    }
  }

  // ================================================================================
  // For CreateAnimation
  if (par::isCreateAnimation) {
    int nQuark = 0;
    int nAntiQuark = 0;
    for (int i = 0; i < nPartons; i++) {
      if (partons[i].PDG() > 0) {
        // g_quark_all->AddPoint(partons[i].X(), partons[i].Y());
        // 改成setpoint, 以支持旧版本的ROOT
        g_quark_all->SetPoint(nQuark, partons[i].X(), partons[i].Y());
        nQuark++;
      } else {
        // g_anti_quark_all->AddPoint(partons[i].X(), partons[i].Y());
        // 改成setpoint
        g_anti_quark_all->SetPoint(nAntiQuark, partons[i].X(), partons[i].Y());
        nAntiQuark++;
      }
    }
    dummy->Draw();
    g_quark_all->Draw("same P");
    g_anti_quark_all->Draw("same P");
    canvas->Update();
    canvas->Write();
  }
  // ================================================================================


  // *=*=*=*=*=*=*=*=*=*=*=*=*=*=
  // Start Coalescence
  // *=*=*=*=*=*=*=*=*=*=*=*=*=*=
  for (int iParton = 0; iParton < nPartons; iParton++) {
    if (partons[iParton].IsUsed()) continue;
    // 第0个parton
    float pdg0 = partons[iParton].PDG();
    float x0 = partons[iParton].X(), y0 = partons[iParton].Y(), z0 = partons[iParton].Z();
    float px0 = partons[iParton].Px(), py0 = partons[iParton].Py(), pz0 = partons[iParton].Pz();
    float t0 = partons[iParton].Time();

    if (par::isCreateAnimation) {
      g_quark_start->SetPoint(g_quark_start->GetN(), x0, y0);
    }

    // meson
    int pdg_me = 0;
    float x_me = 0, y_me = 0, z_me = 0;
    float px_me = 0, py_me = 0, pz_me = 0;
    float t_me = 0;
    float d_meson_min = std::numeric_limits<float>::max();

    // diquark
    float d_diquark_min = std::numeric_limits<float>::max();

    // baryon
    int pdg_ba = 0;
    float x_ba = 0, y_ba = 0, z_ba = 0;
    float px_ba = 0, py_ba = 0, pz_ba = 0;
    float t_ba = 0;
    float d_baryon_min = std::numeric_limits<float>::max();

    // 存储用于生成hadron的quark的脚标
    int meson_quark_label[2] = {iParton, -1};
    int diquark_quark_label[2] = {iParton, -1};
    int baryon_quark_label[3] = {iParton, -1, -1};

    int nSerialLastDiquark = -1; // 上一个di-quark的序列号

    for (int jParton = iParton + 1; jParton < nPartons; jParton++) {
      if (partons[jParton].IsUsed()) continue;

      int pdg1 = partons[jParton].PDG();
      float x1 = partons[jParton].X(), y1 = partons[jParton].Y(), z1 = partons[jParton].Z();
      float px1 = partons[jParton].Px(), py1 = partons[jParton].Py(), pz1 = partons[jParton].Pz();
      float t1 = partons[jParton].Time();

      //临时的x0, y0, z0, x1, y1, z1，用来做move on，move on会改变这些值
      float x0_tmp = x0, y0_tmp = y0, z0_tmp = z0;
      float x1_tmp = x1, y1_tmp = y1, z1_tmp = z1;

      float d = distance3DMoveOn(x0_tmp, y0_tmp, z0_tmp, x1_tmp, y1_tmp, z1_tmp, px0, py0, pz0, px1, py1, pz1, t0, t1);

      float d_meson = d; // 这里是为了下面的pi0的特殊处理
      float d_diquark = d;

      if (pdg0 * pdg1 < 0) {
        // 如果说是 u-ubar 或者 d-dbar，那么查表可以得到一个介子的pdg
        int pdg_me_tmp = LookupMesonSpecies(pdg0, pdg1);
        
        // 临时的px, py, pz，用来做mass varify，mass varify会改变这些值
        float px_me_tmp = 0, py_me_tmp = 0, pz_me_tmp = 0;
        bool isMassValid = DeriveHadronPxPyPz(pdg_me_tmp, pdg0, pdg1, px_me_tmp, py_me_tmp, pz_me_tmp, px0, py0, pz0, px1, py1, pz1);
        if (!isMassValid && par::isWriteCoalQA) nTrackRejectByMass++;

        // 设置50%的概率为pi0，50%的概率这次不生成，即让距离变得无限大
        if (pdg_me_tmp == 111) {
          if(static_cast<bool>(par::zero_or_one(par::gen))) d_meson = std::numeric_limits<float>::max();
        }

        if (isMassValid && d_meson < d_meson_min) {
          d_meson_min = d_meson; // 更新最小距离
          meson_quark_label[1] = jParton;

          // 本次生成的介子：
          pdg_me = pdg_me_tmp;
          x_me = (x0_tmp + x1_tmp) / 2, y_me = (y0_tmp + y1_tmp) / 2, z_me = (z0_tmp + z1_tmp) / 2;
          px_me = px_me_tmp, py_me = py_me_tmp, pz_me = pz_me_tmp;
          t_me = t0 > t1 ? t0 : t1; // 取最大的时间
        }

      } else if (pdg0 * pdg1 > 0) {

        // 只可能是di-quark
        if(d_diquark < d_diquark_min) {
          d_diquark_min = d_diquark; // 更新最小距离
          diquark_quark_label[1] = jParton;

          // 找到了新的di-quark，将jParton标记为已使用
          partons[jParton].LabelAsUsedByDiQuark();
          if (nSerialLastDiquark != -1) {
            // 将上一个备选di-quark中的标记清除
            partons[nSerialLastDiquark].ClearLabelAsUsedByDiQuark();
          }
          nSerialLastDiquark = jParton;
        }
      } else {
        continue; // 不会出现这种情况
      }
    }

    // 查看是否找到了meson
    bool isThereMeson = meson_quark_label[1] != -1;
    // 查看是否找到了di-quark
    bool isThereDiquark = diquark_quark_label[1] != -1;
    bool isThereBaryon = false;

    int pdg1 = 0;
    float x1 = 0, y1 = 0, z1 = 0;
    float px1 = 0, py1 = 0, pz1 = 0;
    float t1 = 0;

    // 找到di-quark或者meson的quark
    if (isThereDiquark) {
      //一个简单的检查,确保iParton和diquark_quark_label[0]是同一个quark
      if (partons[iParton].GetSerial() != partons[diquark_quark_label[0]].GetSerial()) {
        std::cerr<<"Error: Di-quark quark0 is not the same as partons[diquark_quark_label[0]]"<<std::endl;
        std::abort();
      }
      // 第1个parton(用于生成di-quark的parton,用diquark_quark_label[1]标记了)
      pdg1 = partons[diquark_quark_label[1]].PDG();
      x1 = partons[diquark_quark_label[1]].X(), y1 = partons[diquark_quark_label[1]].Y(), z1 = partons[diquark_quark_label[1]].Z();
      px1 = partons[diquark_quark_label[1]].Px(), py1 = partons[diquark_quark_label[1]].Py(), pz1 = partons[diquark_quark_label[1]].Pz();
      t1 = partons[diquark_quark_label[1]].Time();
    }

    if (isThereMeson) {
      if (partons[iParton].GetSerial() != partons[meson_quark_label[0]].GetSerial()) {
        std::cerr<<"Error: quark0 is not the same as partons[meson_quark_label[0]]"<<std::endl;
        std::abort();
      }
      // 对于meson事实上已经不需要读取第1个quark的信息了
    }


    if (isThereDiquark) {
      //如果找到了di-quark，那么还需要找到第2个quark以组成一个baryon
      for (int kParton = iParton + 1; kParton < nPartons; kParton++) {
        if (partons[kParton].IsUsed()) continue;
        if (partons[kParton].IsUsedAsDiQuark()) continue;

        // 第0个parton(在iParton循环中的parton)
        float x0_tmp = x0, y0_tmp = y0, z0_tmp = z0;
        // 第1个parton(被diquark_quark_label[1]标记的parton)
        float x1_tmp = x1, y1_tmp = y1, z1_tmp = z1;
  
        // 第2个parton(这个循环中的parton)
        int pdg2 = partons[kParton].PDG();
        if (pdg2 * pdg1 < 0) continue; // 保证第2个quark和第1个quark是同号的
        float x2 = partons[kParton].X(), y2 = partons[kParton].Y(), z2 = partons[kParton].Z();
        float px2 = partons[kParton].Px(), py2 = partons[kParton].Py(), pz2 = partons[kParton].Pz();
        float t2 = partons[kParton].Time();

        int pdg_ba_tmp = LookupBaryonSpecies(pdg0, pdg1, pdg2);
        float px_ba_tmp = 0, py_ba_tmp = 0, pz_ba_tmp = 0;
        bool isMassValid = DeriveHadronPxPyPz(pdg_ba_tmp, pdg0, pdg1, pdg2, px_ba_tmp, py_ba_tmp, pz_ba_tmp, px0, py0, pz0, px1, py1, pz1, px2, py2, pz2);
        if (!isMassValid && par::isWriteCoalQA) nTrackRejectByMass++;
  
        float d_baryon = perimeterMoveOn(x0_tmp, y0_tmp, z0_tmp, x1_tmp, y1_tmp, z1_tmp, x2, y2, z2, px0, py0, pz0, px1, py1, pz1, px2, py2, pz2, t0, t1, t2);
        d_baryon = d_baryon / 3.; // 周长的距离除以3，得到平均距离

        if (isMassValid && d_baryon < d_baryon_min) {
          d_baryon_min = d_baryon;
          baryon_quark_label[1] = diquark_quark_label[1];
          baryon_quark_label[2] = kParton;

          pdg_ba = pdg_ba_tmp;
          px_ba = px_ba_tmp, py_ba = py_ba_tmp, pz_ba = pz_ba_tmp;
          x_ba = (x0_tmp + x1_tmp + x2) / 3, y_ba = (y0_tmp + y1_tmp + y2) / 3, z_ba = (z0_tmp + z1_tmp + z2) / 3;
          t_ba = (t0 > t1) ? ((t0 > t2) ? t0 : t2) : ((t1 > t2) ? t1 : t2); // 取最大的时间
          isThereBaryon = true;
        }
      }
      partons[diquark_quark_label[1]].ClearLabelAsUsedByDiQuark();
    } else {
      // 如果没有找到di-quark，那么也不可能找到baryon
      isThereBaryon = false;
    }

    if(!isThereMeson && !isThereBaryon) {
      // 如果没有找到meson和diquark，那么这个parton就是一个孤立的parton，无法形成hadron
      if(par::isDebug) std::cout<<"No meson or baryon found: this parton is isolated, and should left by a unsuccessful pi0 coalescence"<<std::endl;
    }

    // if(par::isDebug) {
    //   std::cout<<"d_meson_min = " << d_meson_min << std::endl;
    //   std::cout<<"r_bm * d_baryon_min = " << r_bm * d_baryon_min << std::endl;
    // }

    // 在这里，我们已经找到了一个meson或者一个baryon
    // 我们需要根据b_meson 和 r_bm * b_baryon的 大小关系，选择一个距离最小的，然后生成hadron将这个hadron加入到hadrons中
    if (d_meson_min < r_bm * d_baryon_min) {
      if (isThereMeson) {
        hadrons.emplace_back(nHadronSerial++, pdg_me, x_me, y_me, z_me, px_me, py_me, pz_me, t_me, d_meson_min, partons[meson_quark_label[0]].GetSerial(), partons[meson_quark_label[1]].GetSerial(), -9999);
        if(par::isDebug) {
          hadrons.back().SetParton0Position(x0, y0, z0);
          hadrons.back().SetParton1Position(partons[meson_quark_label[1]].X(), partons[meson_quark_label[1]].Y(), partons[meson_quark_label[1]].Z());
          hadrons.back().SetParton2Position(-9999, -9999, -9999);
        }
        partons[meson_quark_label[0]].LabelAsUsed();
        partons[meson_quark_label[1]].LabelAsUsed();

        if(par::isCreateAnimation) {
          if (partons[meson_quark_label[0]].PDG() < 0) {
            g_anti_quark->SetPoint(g_anti_quark->GetN(), partons[meson_quark_label[0]].X(), partons[meson_quark_label[0]].Y());
          } else {
            g_quark->SetPoint(g_quark->GetN(), partons[meson_quark_label[0]].X(), partons[meson_quark_label[0]].Y());
          }
          if (partons[meson_quark_label[1]].PDG() < 0) {
            g_anti_quark->SetPoint(g_anti_quark->GetN(), partons[meson_quark_label[1]].X(), partons[meson_quark_label[1]].Y());
          } else {
            g_quark->SetPoint(g_quark->GetN(), partons[meson_quark_label[1]].X(), partons[meson_quark_label[1]].Y());
          }
          g_meson->SetPoint(g_meson->GetN(), x_me, y_me);
          std::unique_ptr<TGraph> g_meson_shape_tmp = std::unique_ptr<TGraph>(new TGraph());
          g_meson_shape_tmp->SetPoint(0, partons[meson_quark_label[0]].X(), partons[meson_quark_label[0]].Y());
          g_meson_shape_tmp->SetPoint(1, partons[meson_quark_label[1]].X(), partons[meson_quark_label[1]].Y());
          g_meson_shape.emplace_back(std::move(g_meson_shape_tmp));
        }
        if (par::isWriteCoalQA) h_coal_dis_meson->Fill(d_meson_min);
      }
    } else if (d_meson_min > r_bm * d_baryon_min) {
      if (isThereBaryon) {
        hadrons.emplace_back(nHadronSerial++, pdg_ba, x_ba, y_ba, z_ba, px_ba, py_ba, pz_ba, t_ba, d_baryon_min, partons[baryon_quark_label[0]].GetSerial(), partons[baryon_quark_label[1]].GetSerial(), partons[baryon_quark_label[2]].GetSerial());
        if(par::isDebug) {
          hadrons.back().SetParton0Position(x0, y0, z0);
          hadrons.back().SetParton1Position(partons[baryon_quark_label[1]].X(), partons[baryon_quark_label[1]].Y(), partons[baryon_quark_label[1]].Z());
          hadrons.back().SetParton2Position(partons[baryon_quark_label[2]].X(), partons[baryon_quark_label[2]].Y(), partons[baryon_quark_label[2]].Z());
        }
        partons[baryon_quark_label[0]].LabelAsUsed();
        partons[baryon_quark_label[1]].LabelAsUsed();
        partons[baryon_quark_label[2]].LabelAsUsed();

        if(par::isCreateAnimation) {
          if (partons[baryon_quark_label[0]].PDG() < 0) {
            g_anti_quark->SetPoint(g_anti_quark->GetN(), partons[baryon_quark_label[0]].X(), partons[baryon_quark_label[0]].Y());
          } else {
            g_quark->SetPoint(g_quark->GetN(), partons[baryon_quark_label[0]].X(), partons[baryon_quark_label[0]].Y());
          }
          if (partons[baryon_quark_label[1]].PDG() < 0) {
            g_anti_quark->SetPoint(g_anti_quark->GetN(), partons[baryon_quark_label[1]].X(), partons[baryon_quark_label[1]].Y());
          } else {
            g_quark->SetPoint(g_quark->GetN(), partons[baryon_quark_label[1]].X(), partons[baryon_quark_label[1]].Y());
          }
          if (partons[baryon_quark_label[2]].PDG() < 0) {
            g_anti_quark->SetPoint(g_anti_quark->GetN(), partons[baryon_quark_label[2]].X(), partons[baryon_quark_label[2]].Y());
          } else {
            g_quark->SetPoint(g_quark->GetN(), partons[baryon_quark_label[2]].X(), partons[baryon_quark_label[2]].Y());
          }
          if (pdg_ba > 0) {
            g_baryon->SetPoint(g_baryon->GetN(), x_ba, y_ba);
          } else {
            g_anti_baryon->SetPoint(g_anti_baryon->GetN(), x_ba, y_ba);
          }
          std::unique_ptr<TGraph> g_baryon_shape_tmp = std::unique_ptr<TGraph>(new TGraph());
          g_baryon_shape_tmp->SetPoint(0, partons[baryon_quark_label[0]].X(), partons[baryon_quark_label[0]].Y());
          g_baryon_shape_tmp->SetPoint(1, partons[baryon_quark_label[1]].X(), partons[baryon_quark_label[1]].Y());
          g_baryon_shape_tmp->SetPoint(2, partons[baryon_quark_label[2]].X(), partons[baryon_quark_label[2]].Y());
          g_baryon_shape_tmp->SetPoint(3, partons[baryon_quark_label[0]].X(), partons[baryon_quark_label[0]].Y());
          g_baryon_shape.emplace_back(std::move(g_baryon_shape_tmp));
        }
        if (par::isWriteCoalQA) h_coal_dis_baryon->Fill(d_baryon_min);
      }
    }

    if(par::isCreateAnimation) {
      canvas->SetName(Form("frame_%d", frame_number));
      dummy->Draw();
      //Draw all partons
      // 只有当GetN() > 0时才画，否则会报错
      if (g_quark_all->GetN() > 0) g_quark_all->Draw("same P");
      if (g_anti_quark_all->GetN() > 0) g_anti_quark_all->Draw("same P");
      if (g_quark_start->GetN() > 0) g_quark_start->Draw("same P");
      if (g_quark->GetN() > 0) g_quark->Draw("same P");
      if (g_anti_quark->GetN() > 0) g_anti_quark->Draw("same P");
      if (g_meson->GetN() > 0) g_meson->Draw("same P");
      if (g_baryon->GetN() > 0) g_baryon->Draw("same P");
      if (g_anti_baryon->GetN() > 0) g_anti_baryon->Draw("same P");
      for (int i = 0; i < g_meson_shape.size(); i++) {
        g_meson_shape[i]->Draw("same L");
      }
      for (int i = 0; i < g_baryon_shape.size(); i++) {
        g_baryon_shape[i]->Draw("same L");
      }
      frame_number++;
      canvas->Update();
      canvas->Write();
    }

    if(par::isDebug) {
      if (!partons[iParton].IsUsed()) {
        std::cout<<"Parton "<<iParton<<" is not used, its PDG is "<<partons[iParton].PDG()<<std::endl;
        std::cout<<"now d_meson_min: "<<d_meson_min<<" and d_baryon_min: "<<d_baryon_min<<std::endl;
      }
    }

  }

  if(par::isCreateAnimation) {
    canvas->Update();
    canvas->Write();
    file->Close();
  }

  // 如果还有没有被使用的parton，这些parton的数量理论上应该很少
  // 可以读取isThisPartonUsed数组，找到没有被使用的parton，然后将这些parton打包成一个vector<Parton>
  // 进行递归
  std::vector<Parton> partonsUnused;
  for (auto parton : partons) {
    if (!parton.IsUsed()) {
      partonsUnused.push_back(parton);
    }
  }
  if(par::isDebug) {
    std::cout<<"Number of partons left: "<<partonsUnused.size()<<std::endl;
    for (int i = 0; i < partonsUnused.size(); i++) {
      std::cout<<"Parton "<<partonsUnused[i].GetSerial()<<" with PDG code: "<<partonsUnused[i].PDG()<<" is left."<<std::endl;
    }
  }
  // 如果是CreateAnimation模式，那么递归结束
  if(par::isCreateAnimation) {
    std::cout<<"CreateAnimation mode is on, recursion ends."<<std::endl;
    return;
  }
  // 是否通过输出参数直接关闭递归
  if(par::isTurnOffRecursion && partonsUnused.size() > 0) {
    std::cout<<"Recursion is turned off, recursion ends."<<std::endl;
    return;
  }

  if (partonsUnused.size() == 2) {
    // 如果只有两个个没有被使用的parton，可以尝试直接生成一个介子
    float x, y, z, px, py, pz, t, d;
    float x0, y0, z0, px0, py0, pz0, t0 = partonsUnused[0].Time();
    float x1, y1, z1, px1, py1, pz1, t1 = partonsUnused[1].Time();
    partonsUnused[0].GetPosition(x0, y0, z0);
    partonsUnused[0].GetMomentum(px0, py0, pz0);
    partonsUnused[1].GetPosition(x1, y1, z1);
    partonsUnused[1].GetMomentum(px1, py1, pz1);
    int pdg_lookup = LookupMesonSpecies(partonsUnused[0].PDG(), partonsUnused[1].PDG());
    // 如果刚好能够形成一个介子, 而且可以通过DeriveHadronPxPyPz函数验证质量（已经包含了pdg组合验证和mass验证）, 那么打包成一个介子,递归结束
    bool isMassValid = DeriveHadronPxPyPz(pdg_lookup, partonsUnused[0].PDG(), partonsUnused[1].PDG(), px, py, pz, px0, py0, pz0, px1, py1, pz1);
    if (!isMassValid && par::isWriteCoalQA) nTrackRejectByMass++;
    if(isMassValid) {
      if (par::isWriteCoalQA) h_coal_dis_meson->Fill(d);
      if(par::isDebug) std::cout<<"Two partons (No."<<partonsUnused[0].GetSerial()<<", No."<<partonsUnused[1].GetSerial()<<") left with PDG code: "<<partonsUnused[0].PDG()<<", "<<partonsUnused[1].PDG()<<std::endl;
      d = distance3DMoveOn(x0, y0, z0, x1, y1, z1, px0, py0, pz0, px1, py1, pz1, t0, t1);
      x = (x0 + x1) / 2, y = (y0 + y1) / 2, z = (z0 + z1) / 2;
      t = t0 > t1 ? t0 : t1;
      hadrons.emplace_back(nHadronSerial++, pdg_lookup, x, y, z, px, py, pz, t, d, partonsUnused[0].GetSerial(), partonsUnused[1].GetSerial(), 0);
      if(par::isDebug) {
        hadrons.back().SetParton0Position(x0, y0, z0);
        hadrons.back().SetParton1Position(x1, y1, z1);
      }
      if(par::isDebug) std::cout<<"------These two partons can form a meson: "<<pdg_lookup<<std::endl;
      partonsUnused.clear();
      if(par::isDebug) std::cout<<"------Used partons vector cleared."<<std::endl;
    } else {
      // 无法形成介子
      if(par::isDebug) {
        std::cout<<"Two partons left with pdg code: "<<partonsUnused[0].PDG()<<", "<<partonsUnused[1].PDG()<<std::endl;
        std::cout<<"------These two partons cannot form a meson."<<std::endl;
      }
    }
  } else if (partonsUnused.size() == 3) {
    // 如果只有三个个没有被使用的parton
    float x, y, z, px, py, pz, t, d;
    float x0, y0, z0, px0, py0, pz0, t0 = partonsUnused[0].Time();
    float x1, y1, z1, px1, py1, pz1, t1 = partonsUnused[1].Time();
    float x2, y2, z2, px2, py2, pz2, t2 = partonsUnused[2].Time();
    partonsUnused[0].GetPosition(x0, y0, z0);
    partonsUnused[0].GetMomentum(px0, py0, pz0);
    partonsUnused[1].GetPosition(x1, y1, z1);
    partonsUnused[1].GetMomentum(px1, py1, pz1);
    partonsUnused[2].GetPosition(x2, y2, z2);
    partonsUnused[2].GetMomentum(px2, py2, pz2);
    // 如果刚好能够形成一个重子, 而且可以通过MassVarify函数验证质量（MassVarify已经包含了堆pdg的验证）, 那么打包成一个重子,递归结束
    int pdg_lookup = LookupBaryonSpecies(partonsUnused[0].PDG(), partonsUnused[1].PDG(), partonsUnused[2].PDG());
    bool isMassValid = DeriveHadronPxPyPz(pdg_lookup, partonsUnused[0].PDG(), partonsUnused[1].PDG(), partonsUnused[2].PDG(), px, py, pz, px0, py0, pz0, px1, py1, pz1, px2, py2, pz2);
    if (!isMassValid && par::isWriteCoalQA) nTrackRejectByMass++;
    if (isMassValid) {
      if (par::isWriteCoalQA) h_coal_dis_baryon->Fill(d);
      if(par::isDebug) std::cout<<"Three partons (No."<<partonsUnused[0].GetSerial()<<", No."<<partonsUnused[1].GetSerial()<<", No."<<partonsUnused[2].GetSerial()<<") left with PDG code: "<<partonsUnused[0].PDG()<<", "<<partonsUnused[1].PDG()<<", "<<partonsUnused[2].PDG()<<std::endl;
      d = perimeterMoveOn(x0, y0, z0, px0, py0, pz0, t0, x1, y1, z1, px1, py1, pz1, t1, x2, y2, z2, px2, py2, pz2, t2);
      x = (x0 + x1 + x2) / 3, y = (y0 + y1 + y2) / 3, z = (z0 + z1 + z2) / 3;
      t = (t0 > t1) ? ((t0 > t2) ? t0 : t2) : ((t1 > t2) ? t1 : t2);
      hadrons.emplace_back(nHadronSerial++, pdg_lookup, x, y, z, px, py, pz, t, d, partonsUnused[0].GetSerial(), partonsUnused[1].GetSerial(), partonsUnused[2].GetSerial());
      if(par::isDebug) {
        hadrons.back().SetParton0Position(partonsUnused[0].X(), partonsUnused[0].Y(), partonsUnused[0].Z());
        hadrons.back().SetParton1Position(partonsUnused[1].X(), partonsUnused[1].Y(), partonsUnused[1].Z());
        hadrons.back().SetParton2Position(partonsUnused[2].X(), partonsUnused[2].Y(), partonsUnused[2].Z());
      }
      if(par::isDebug) std::cout<<"------These three partons can form a baryon: "<<pdg_lookup<<std::endl;
      partonsUnused.clear();
      if(par::isDebug) std::cout<<"------Used partons vector cleared."<<std::endl;
    }
  }

  bool isNeedRecursion = false;
  partonsUnused.size() > 0 ? isNeedRecursion = true : isNeedRecursion = false;
  if (nRecursionThisEvent > 3) {
    // 如果递归超过3次，那么直接结束递归
    isNeedRecursion = false;
    if(par::isDebug) std::cout<<"Recursion times exceed 3, recursion ends."<<std::endl;
  }

  // 递归
  if (isNeedRecursion) {
    nRecursionThisEvent++;
    if(par::isDebug) {
      std::cout<<"Recursion "<<nRecursionThisEvent<<" starts."<<std::endl;
    }
    ProcessFromParton(partonsUnused, hadrons, nHadronSerial);
  } else {
    if(par::isDebug) {
      std::cout<<"----------------------------------------"<<std::endl;
      std::cout<<"No need for recursion again."<<std::endl;
      std::cout<<"Number of recursion in this event: "<<nRecursionThisEvent<<std::endl;
      std::cout<<"----------------------------------------"<<std::endl;
    }
    if (partonsUnused.size() > par::flavourBreakTolerance * nPartonsThisEvent) {
      std::cout<<"Error: "<<partonsUnused.size()<<" partons are left unused after recursion, which is more than "<<par::flavourBreakTolerance * 100. <<"% of the total partons."<<std::endl;
      std::cout<<"This event will be saved as empty."<<std::endl;
      std::cout<<"Clearing hadrons array."<<std::endl;
      isRejectByFlavourTolerance = true;
      std::vector<Hadron>().swap(hadrons);
    }
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

  // ω介子 (暂时注释掉)
  // ω (u, anti-u) or (d, anti-d)
  // mesonLookupTable[std::make_tuple(-2, 2)] = 223;
  // mesonLookupTable[std::make_tuple(-1, 1)] = 223;

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
  // // 反Δ^+ (anti-u, anti-u, anti-d) -> -2214 // 组分与反质子相同，暂时注释掉
  // baryonLookupTable[std::make_tuple(-2, -2, -1)] = -2214;
  // // Δ^0 (udd) -> 2114 // 组分与中子相同，暂时注释掉
  // baryonLookupTable[std::make_tuple(1, 1, 2)] = 2114;
  // // 反Δ^0 (anti-d, anti-d, anti-u) -> -2114 // 组分与反中子相同，暂时注释掉
  // baryonLookupTable[std::make_tuple(-2, -1, -1)] = -2114;

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


bool Coalescence::DeriveHadronPxPyPz(const int genPdg, const int pdg0, const int pdg1, float& genpx, float& genpy, float& genpz, const float px0, const float py0, const float pz0, const float px1, const float py1, const float pz1) {
  // 计算两个粒子的能量
  if (par::mass.find(genPdg) == par::mass.end() || par::mass.find(pdg0) == par::mass.end() || par::mass.find(pdg1) == par::mass.end()) return false;
  // 不做质量验证
  if (!par::isEnableMassVarify) {
    genpx = px0 + px1;
    genpy = py0 + py1;
    genpz = pz0 + pz1;
    return true;
  }
  float E0 = sqrt(px0 * px0 + py0 * py0 + pz0 * pz0 + par::mass[pdg0] * par::mass[pdg0]);
  float E1 = sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + par::mass[pdg1] * par::mass[pdg1]);
  float E = E0 + E1;
  // 检查总能量是否小于目标粒子的静质量
  if (E < par::mass[genPdg]) return false;
  // 计算方向
  float px = px0 + px1;
  float py = py0 + py1;
  float pz = pz0 + pz1;
  float p = sqrt(px * px + py * py + pz * pz);
  // 计算总动量的大小
  // 确保动量方向保持一致并且重新归一化
  if (p > 1.e-6) {
    genpx = sqrt((E * E) - par::mass[genPdg] * par::mass[genPdg]) * px / p;
    genpy = sqrt((E * E) - par::mass[genPdg] * par::mass[genPdg]) * py / p;
    genpz = sqrt((E * E) - par::mass[genPdg] * par::mass[genPdg]) * pz / p;
  } else return false;
  return true;
}

bool Coalescence::DeriveHadronPxPyPz(const int genPdg, const int pdg0, const int pdg1, const int pdg2, float& genpx, float& genpy, float& genpz, const float px0, const float py0, const float pz0, const float px1, const float py1, const float pz1, const float px2, const float py2, const float pz2) {
  if (par::mass.find(genPdg) == par::mass.end() || par::mass.find(pdg0) == par::mass.end() || par::mass.find(pdg1) == par::mass.end() || par::mass.find(pdg2) == par::mass.end()) return false;
  //不做质量验证
  if (!par::isEnableMassVarify) {
    genpx = px0 + px1 + px2;
    genpy = py0 + py1 + py2;
    genpz = pz0 + pz1 + pz2;
    return true;
  }
  // 计算两个粒子的能量
  float E0 = sqrt(px0 * px0 + py0 * py0 + pz0 * pz0 + par::mass[pdg0] * par::mass[pdg0]);
  float E1 = sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + par::mass[pdg1] * par::mass[pdg1]);
  float E2 = sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + par::mass[pdg2] * par::mass[pdg2]);
  float E = E0 + E1 + E2;
  // 检查总能量是否小于目标粒子的静质量
  if (E < par::mass[genPdg]) return false;
  // 计算方向
  float px = px0 + px1 + px2;
  float py = py0 + py1 + py2;
  float pz = pz0 + pz1 + pz2;
  // 计算总动量的大小
  float p = sqrt(px * px + py * py + pz * pz);
  // 确保动量方向保持一致并且重新归一化
  if (p > 1.e-6) {
    genpx = sqrt((E * E) - par::mass[genPdg] * par::mass[genPdg]) * px / p;
    genpy = sqrt((E * E) - par::mass[genPdg] * par::mass[genPdg]) * py / p;
    genpz = sqrt((E * E) - par::mass[genPdg] * par::mass[genPdg]) * pz / p;
  } else return false;
  return true;
}