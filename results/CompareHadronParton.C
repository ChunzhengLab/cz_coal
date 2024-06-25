#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"

using namespace std;

#include <iostream>
#include <vector>
#include <unordered_set>

#include <iostream>
#include <vector>
#include <unordered_set>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <set>

void printUniqueElements(const std::vector<int>& v1 , const std::vector<int>& v2) {
  std::set<int> set1(v1.begin() , v1.end());
  std::set<int> set2(v2.begin() , v2.end());
  std::vector<int> diff;

  // Find elements in set1 not in set2
  std::set_difference(set1.begin() , set1.end() , set2.begin() , set2.end() , std::back_inserter(diff));
  // Find elements in set2 not in set1
  std::set_difference(set2.begin() , set2.end() , set1.begin() , set1.end() , std::back_inserter(diff));

  // Output the unique elements
  std::cout << "Non-used labels: ";
  for (const auto& elem : diff) {
    std::cout << elem << " ";
  }
  std::cout << std::endl;
}

std::vector<int> mergeVectors(const std::vector<int>& v1 , const std::vector<int>& v2 , const std::vector<int>& v3) {
  std::unordered_set<int> set;

  // 插入所有向量的元素到set中，重复的自动被忽略
  set.insert(v1.begin() , v1.end());
  set.insert(v2.begin() , v2.end());
  set.insert(v3.begin() , v3.end());

  // 将set的元素转回vector
  std::vector<int> merged(set.begin() , set.end());
  return merged;
}

void printCommonElements(const std::vector<int>& v1 , const std::vector<int>& v2) {
  std::unordered_set<int> seen; // 用于存储v1中的元素
  std::unordered_set<int> commonElements; // 用于存储找到的共同元素，避免重复输出

  // 将v1中的元素加入到哈希集中
  for (int num : v1) {
    seen.insert(num);
  }

  // 检查v2中的元素是否在哈希集中
  for (int num : v2) {
    if (seen.find(num) != seen.end()) {
      commonElements.insert(num);
    }
  }

  // 输出所有找到的共同元素
  std::cout << " ";
  for (int num : commonElements) {
    if (num != -1) { // 排除-1，因为它已经被特殊处理
      std::cout << " " << num;
    }
  }
  std::cout << std::endl;
}

struct PartonEvent {
  int nevent;
  int nparton;
  int   ID[30000];   //[nparton]
  float Px[30000];   //[nparton]
  float Py[30000];   //[nparton]
  float Pz[30000];   //[nparton]
  float X[30000];   //[nparton]
  float Y[30000];   //[nparton]
  float Z[30000];   //[nparton]
};

struct HadronEvent {
  int nSeries;
  int nTracks;
  vector<int>* PDG = nullptr;
  vector<float>* x = nullptr;
  vector<float>* y = nullptr;
  vector<float>* z = nullptr;
  vector<int>* quark0 = nullptr;
  vector<int>* quark1 = nullptr;
  vector<int>* quark2 = nullptr;
  vector<float>* quark0_x = nullptr;
  vector<float>* quark0_y = nullptr;
  vector<float>* quark0_z = nullptr;
  vector<float>* quark1_x = nullptr;
  vector<float>* quark1_y = nullptr;
  vector<float>* quark1_z = nullptr;
  vector<float>* quark2_x = nullptr;
  vector<float>* quark2_y = nullptr;
  vector<float>* quark2_z = nullptr;
};


void CompareHadronParton() {
  gStyle->SetOptStat(0);

  TFile* fileParton = new TFile("../rand_partons_ampt_sampling_np30.root");
  TFile* fileHadron = new TFile("../coal_hadrons_rndparton_np30.root");

  TTree* treeParton = (TTree*)fileParton->Get("AMPT");
  TTree* treeHadron = (TTree*)fileHadron->Get("CoalHadrons");

  HadronEvent hadron;
  treeHadron->SetBranchAddress("nSeries" , &hadron.nSeries);
  treeHadron->SetBranchAddress("nTracks" , &hadron.nTracks);
  treeHadron->SetBranchAddress("PDG" , &(hadron.PDG));
  treeHadron->SetBranchAddress("x" , &(hadron.x));
  treeHadron->SetBranchAddress("y" , &(hadron.y));
  treeHadron->SetBranchAddress("z" , &(hadron.z));
  treeHadron->SetBranchAddress("quark0" , &(hadron.quark0));
  treeHadron->SetBranchAddress("quark1" , &(hadron.quark1));
  treeHadron->SetBranchAddress("quark2" , &(hadron.quark2));
  treeHadron->SetBranchAddress("quark0_x" , &(hadron.quark0_x));
  treeHadron->SetBranchAddress("quark0_y" , &(hadron.quark0_y));
  treeHadron->SetBranchAddress("quark0_z" , &(hadron.quark0_z));
  treeHadron->SetBranchAddress("quark1_x" , &(hadron.quark1_x));
  treeHadron->SetBranchAddress("quark1_y" , &(hadron.quark1_y));
  treeHadron->SetBranchAddress("quark1_z" , &(hadron.quark1_z));
  treeHadron->SetBranchAddress("quark2_x" , &(hadron.quark2_x));
  treeHadron->SetBranchAddress("quark2_y" , &(hadron.quark2_y));
  treeHadron->SetBranchAddress("quark2_z" , &(hadron.quark2_z));

  PartonEvent parton;
  treeParton->SetBranchAddress("Event" , &parton);
  treeParton->SetBranchAddress("ID" , parton.ID);
  treeParton->SetBranchAddress("Px" , parton.Px);
  treeParton->SetBranchAddress("Py" , parton.Py);
  treeParton->SetBranchAddress("Pz" , parton.Pz);
  treeParton->SetBranchAddress("X" , parton.X);
  treeParton->SetBranchAddress("Y" , parton.Y);
  treeParton->SetBranchAddress("Z" , parton.Z);

  TGraph* xyHadron = new TGraph();
  TGraph* xyParton = new TGraph();

  TGraph* xyMeson = new TGraph();
  TGraph* xyBaryon = new TGraph();
  TGraph* xyQuark = new TGraph();
  TGraph* xyAntiQuark = new TGraph();
  TGraph* xyQuarkFromHadron = new TGraph();

  vector<int> seriesHadronParton0;
  vector<int> seriesHadronParton1;
  vector<int> seriesHadronParton2;
  vector<int> seriesParton;

  std::vector<TGraph2D*> graphHadronForms;
  std::vector<TGraph*> xyHadronShapes; //这里面放小三角形或者小线段

  //std::vector<TGraph2D*> graphPartonForms;
  int nHadronEvents = treeHadron->GetEntries();
  int nPartonEvents = treeParton->GetEntries();
  if (nHadronEvents != nPartonEvents) {
    std::cout << "Number of events in Hadron and Parton files do not match!" << std::endl;
    return;
  }

  for (int iEvent = 3; iEvent < 4; iEvent++) {
    treeHadron->GetEntry(iEvent);
    treeParton->GetEntry(iEvent);
    int nHadrons = hadron.nTracks;
    int nPartons = parton.nparton;

    for (int iHadron = 0; iHadron < nHadrons; iHadron++) {

      cout << "HadronPDG: " << hadron.PDG->at(iHadron) << endl;
      cout << "HadronXYZ: " << hadron.x->at(iHadron) << " " << hadron.y->at(iHadron) << " " << hadron.z->at(iHadron) << endl;
      cout << "Quark0XYZ: " << hadron.quark0_x->at(iHadron) << " " << hadron.quark0_y->at(iHadron) << " " << hadron.quark0_z->at(iHadron) << endl;
      cout << "Quark1XYZ: " << hadron.quark1_x->at(iHadron) << " " << hadron.quark1_y->at(iHadron) << " " << hadron.quark1_z->at(iHadron) << endl;
      // if(hadron.quark2->at(iHadron)!=-1){
      //   cout<<"Quark2XYZ: "<<hadron.quark2_x->at(iHadron)<<" "<<hadron.quark2_y->at(iHadron)<<" "<<hadron.quark2_z->at(iHadron)<<endl;
      // }
      // //查看在 hadron quark0 真实对应的 parton 的位置
      // cout<<"---------------------"<<endl;
      // cout<<"Real Parton0 XYZ : "<<parton.X[hadron.quark0->at(iHadron)]<<" "<<parton.Y[hadron.quark0->at(iHadron)]<<" "<<parton.Z[hadron.quark0->at(iHadron)]<<endl;
      // cout<<"Real Parton1 XYZ : "<<parton.X[hadron.quark1->at(iHadron)]<<" "<<parton.Y[hadron.quark1->at(iHadron)]<<" "<<parton.Z[hadron.quark1->at(iHadron)]<<endl;
      // // if(hadron.quark2->at(iHadron)!=-1){
      // //   cout<<"Real Parton2 XYZ : "<<parton.X[hadron.quark2->at(iHadron)]<<" "<<parton.Y[hadron.quark2->at(iHadron)]<<" "<<parton.Z[hadron.quark2->at(iHadron)]<<endl;
      // // }

      cout << "=====================" << endl;
      //z-x-y
      xyHadron->AddPoint(hadron.x->at(iHadron) , hadron.y->at(iHadron));
      if (abs(hadron.PDG->at(iHadron)) < 1000) { // meson
        xyMeson->AddPoint(hadron.x->at(iHadron) , hadron.y->at(iHadron));
      }
      else { // baryon
        xyBaryon->AddPoint(hadron.x->at(iHadron) , hadron.y->at(iHadron));
      }

      seriesHadronParton0.push_back(hadron.quark0->at(iHadron));
      seriesHadronParton1.push_back(hadron.quark1->at(iHadron));
      seriesHadronParton2.push_back(hadron.quark2->at(iHadron));

      // 这里画出Hadron的形状
      TGraph* xyHadronShape = new TGraph();
      xyHadronShape->AddPoint(parton.X[hadron.quark0->at(iHadron)] , parton.Y[hadron.quark0->at(iHadron)]);
      xyHadronShape->AddPoint(parton.X[hadron.quark1->at(iHadron)] , parton.Y[hadron.quark1->at(iHadron)]);

      xyQuarkFromHadron->AddPoint(parton.X[hadron.quark1->at(iHadron)] , parton.Y[hadron.quark1->at(iHadron)]);
      xyQuarkFromHadron->AddPoint(parton.X[hadron.quark0->at(iHadron)] , parton.Y[hadron.quark0->at(iHadron)]);
      if (hadron.quark2->at(iHadron) == -1) {
        xyHadronShape->AddPoint(parton.X[hadron.quark0->at(iHadron)] , parton.Y[hadron.quark0->at(iHadron)]);
      }
      else {
        xyHadronShape->AddPoint(parton.X[hadron.quark2->at(iHadron)] , parton.Y[hadron.quark2->at(iHadron)]);
        xyHadronShape->AddPoint(parton.X[hadron.quark0->at(iHadron)] , parton.Y[hadron.quark0->at(iHadron)]);

        xyQuarkFromHadron->AddPoint(parton.X[hadron.quark2->at(iHadron)] , parton.Y[hadron.quark2->at(iHadron)]);
      }
      xyHadronShapes.push_back(xyHadronShape);
    }

    for (int iParton = 0; iParton < nPartons; iParton++) {
      xyParton->AddPoint(parton.X[iParton] , parton.Y[iParton]);
      seriesParton.push_back(iParton);
      if (parton.ID[iParton] < 0) {
        xyQuark->AddPoint(parton.X[iParton] , parton.Y[iParton]);
      }
      else if (parton.ID[iParton] > 0) {
        xyAntiQuark->AddPoint(parton.X[iParton] , parton.Y[iParton]);
      }
      else {
        cout << "Quark ID: " << parton.ID[iParton] << endl;
      }
    }

    cout << "======================================" << endl;
    cout << "Same label used between parton0 and parton1" << endl;
    printCommonElements(seriesHadronParton0 , seriesHadronParton1);
    cout << "Same label used between parton1 and parton2" << endl;
    printCommonElements(seriesHadronParton1 , seriesHadronParton2);
    cout << "Same label used between parton0 and parton2" << endl;
    printCommonElements(seriesHadronParton0 , seriesHadronParton2);
    vector<int> merged = mergeVectors(seriesHadronParton0 , seriesHadronParton1 , seriesHadronParton2);
    cout << "Different elements between all partons and partons forme to a hadron: " << endl;
    printUniqueElements(seriesParton , merged);

    //输出这些vector的内
    cout << "seriesHadronParton0: ";
    for (int i = 0; i < seriesHadronParton0.size(); i++) {
      cout << seriesHadronParton0[i] << " ";
    }
    cout << endl;
    cout << "seriesHadronParton1: ";
    for (int i = 0; i < seriesHadronParton1.size(); i++) {
      cout << seriesHadronParton1[i] << " ";
    }
    cout << endl;
    cout << "seriesHadronParton2: ";
    for (int i = 0; i < seriesHadronParton2.size(); i++) {
      cout << seriesHadronParton2[i] << " ";
    }
    cout << endl;
    cout << "seriesParton: ";
    for (int i = 0; i < seriesParton.size(); i++) {
      cout << seriesParton[i] << " ";
    }
    cout << endl;

    seriesHadronParton0.clear();
    seriesHadronParton1.clear();
    seriesHadronParton2.clear();
    seriesParton.clear();

    cout << "======================================" << endl;
  }

  TCanvas* c3 = new TCanvas("c3" , "c3" , 800 , 800);
  TH2D* dummy = new TH2D("dummy" , "Partile Positions;x[fm];y[fm]" , 1 , -10, 10 , 1 , -10 , 10);
  c3->cd();
  dummy->Draw();

  for (auto xyHadronShape : xyHadronShapes) {
    xyHadronShape->SetLineColor(kBlack);
    xyHadronShape->Draw("PL SAME");
  }

  xyMeson->SetMarkerColor(kRed);
  xyMeson->SetMarkerStyle(20);
  xyMeson->Draw("P SAME");

  xyBaryon->SetMarkerColor(kBlue);
  xyBaryon->SetMarkerStyle(20);
  xyBaryon->Draw("P SAME");

  xyQuark->SetMarkerColor(kBlack);
  xyQuark->SetMarkerStyle(20);
  xyQuark->Draw("P SAME");

  xyAntiQuark->SetMarkerColor(kBlack);
  xyAntiQuark->SetMarkerStyle(kOpenCircle);
  xyAntiQuark->Draw("P SAME");

  xyQuarkFromHadron->SetMarkerColor(kOrange);
  xyQuarkFromHadron->SetMarkerStyle(kFullCrossX);
  xyQuarkFromHadron->SetMarkerSize(0.7);
  xyQuarkFromHadron->Draw("P SAME");

  TLegend* legend = new TLegend(0.1 , 0.7 , 0.3 , 0.9);
  legend->AddEntry(xyQuark , "Quark" , "P");
  legend->AddEntry(xyAntiQuark , "Anti-Quark" , "P");
  legend->AddEntry(xyMeson , "Meson" , "P");
  legend->AddEntry(xyBaryon , "Baryon" , "P");
  legend->Draw();
}

