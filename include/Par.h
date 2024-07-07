#pragma once
#include <random>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>

enum class CoalescenceAlgorithm {
  kClassic,
  kFromParton
};

enum class EventType {
  kAMPT,
  kRandom
};

void parseConfig(const std::string& line);
void printConfig();

namespace par {
  extern bool isDebug;
  extern bool isWriteEvents;
  extern bool isCalculateObvs;
  extern bool isRemoveHFQuarks;
  extern EventType eventType;
  //r_bm的默认值
  extern float r_bm;
  extern float flavourBreakTolerance;
  //聚合算法
  extern CoalescenceAlgorithm coalescenceAlgorithm;
  extern std::string inputFile;
  extern std::string outputFile;
  extern std::string obvsFile;

  //维护一个全局的随机数生成器
  extern std::mt19937 gen;
  extern std::uniform_int_distribution<> zero_or_one;

  extern std::map<int, float> mass;
}
