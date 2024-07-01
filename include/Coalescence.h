#ifndef COALESCENCE_H
#define COALESCENCE_H

#include <vector>
#include <map>
#include <tuple>
#include "Par.h" // Include the Par header
#include "Particle.h" // Include the Particle header

// 使用 std::tuple 表示夸克的组合
typedef std::tuple<int, int, int> BaryonCombination;
typedef std::tuple<int, int> MesonCombination;

class Coalescence {
  private:
  // 有一个参数对所有的coalescence对象共享 r_bm
  static bool initialized;
  CoalescenceAlgorithm coalescenceAlgorithm;
  float r_bm;
  //递归次数
  int nRecursionThisEvent = 0;
  int nPartonsThisEvent = 0;
  // 重子信息的查找表
  static std::map<BaryonCombination, int> baryonLookupTable;
  // 介子信息的查找表
  static std::map<MesonCombination, int> mesonLookupTable;
  public:
  explicit Coalescence(float r_bm, CoalescenceAlgorithm coalescenceAlgorithm) : r_bm(r_bm) , coalescenceAlgorithm(coalescenceAlgorithm) {
    if (!initialized) {
      if (baryonLookupTable.empty()) InitBaryonLookupTable();
      if (mesonLookupTable.empty()) InitMesonLookupTable();
      initialized = true;
    }
  }
  bool IsInitialized() const { return initialized; }
  void Process(std::vector<Parton> const &partons, std::vector<Hadron> &hadrons);
  void ProcessFromParton(std::vector<Parton> const &partons0, std::vector<Hadron> &hadrons, int nLastHadronSerial = 0);
  void ProcessClassic(std::vector<Parton> const &partons, std::vector<Hadron> &hadrons);
  int LookupSpecies(int pdg_quark_0, int pdg_quark_1, int pdg_quark_2);
  int LookupMesonSpecies(int pdg_quark_0, int pdg_quark_1);
  int LookupBaryonSpecies(int pdg_quark_0, int pdg_quark_1, int pdg_quark_2);
  void InitBaryonLookupTable();
  void InitMesonLookupTable();
  void ResetRecursionForThisEvent() { nRecursionThisEvent = 0; }
  void Print() const {
    std::cout << "--------------------------" << std::endl;
    std::cout << "Coalescence:" << std::endl;
    std::cout << "Coalenscence process with r_bm = " << r_bm << std::endl;
    std::cout << "--------------------------" << std::endl;
  }
  ~Coalescence() {
    // Do nothing
  }

  //TODO:如果效率不够，可能需要考虑从parton开始loop，而不是从距离开始loop
  //效率太差了，AMPT将会读如～10k个parton，每个parton都要和其他parton计算距离，还要计算和第三个parton的距离，这样的复杂度是O(n^3)
  //如果从parton开始loop，可以将复杂度降低到O(n^2)，方法简单一些，但是有没有办法进一步降低复杂度？
  //可以考虑使用KD树，KD树是一种二叉树，每个节点是一个k维空间中的一个点，每个节点的左子树和右子树分别是根据一个维度划分的，这样可以快速找到最近邻的点
  //KD树的构建复杂度是O(nlogn)，查找复杂度是O(logn)，所以总的复杂度是O(nlogn)
  //如果要使用KD树，可能要调用第三方库，比如FLANN、nanoflann等，哪一个更好呢？我的答案是nanoflann，因为它更小，更快，更简单
};

#endif // COALESCENCE_H
