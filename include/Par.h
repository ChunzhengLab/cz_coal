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

namespace par {
  extern bool isDebug;
  extern bool isCreateAnimation;
  extern bool isWriteCoalQA;
  extern bool isWriteEvents;
  extern bool isCalculateObvs;
  extern bool isRemoveHFQuarks;
  extern bool isEnableMassVarify;
  extern bool isEnableQuarkMoveOn;
  extern bool isBalanceQuarkNumber;
  extern bool isTurnOffRecursion;

  extern bool isRandomRotateQuarkXY;
  extern bool isResampleQuarkXY;
  extern bool isPiesampleXY;
  extern bool isForgetZ;

  extern bool isRandomRotateQuarkPxPy;
  extern bool isResampleQuarkPxPy;
  extern bool isPiesamplePxPy;

  extern EventType eventType;
  extern float r_bm;
  extern float flavourBreakTolerance;
  extern CoalescenceAlgorithm coalescenceAlgorithm;
  extern std::string inputFile;
  extern std::string outputFile;
  extern std::string obvsFile;

  extern std::mt19937 gen;
  extern std::uniform_int_distribution<> zero_or_one;

  extern std::map<int, float> mass;

  void printConfig();
  void initConfig(int argc, char** argv); // 新增的配置初始化函数
}